import os
import argparse
import glob
import re
import subprocess
import sys
import time
import torch
import esm
import warnings

# Suppress DeepSpeed's deprecated PyTorch AMP warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="deepspeed.*")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Batch ESMFold + Foldseek pipeline with telemetry, VRAM safety, and dynamic ETA estimation.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--input", type=str, required=True, help="Path to input FASTA file, directory, or glob pattern")
    parser.add_argument("--db", type=str, default="../reference_dbs/pdb_db", help="Path to local Foldseek database")
    parser.add_argument("--pdb_dir", type=str, default="../foldinator_outputs/pdbs", help="Directory where PDBs are stored")
    parser.add_argument("--output_dir", type=str, default="../foldinator_outputs/foldseek_results", help="Directory to save Foldseek TSV results")
    parser.add_argument("--threshold", type=int, default=500, help="Amino acid length threshold triggering warning (default: 500)")
    parser.add_argument("--min_length", type=int, default=50, help="Automatically skip proteins below this amino acid length (default: 50)")
    parser.add_argument("--limit", type=int, default=None, help="Maximum number of newly folded targets to process in this run")
    parser.add_argument("--max_length", type=int, default=None, help="Automatically skip proteins exceeding this length")
    parser.add_argument("--max_hours", type=float, default=None, help="Maximum runtime hours before graceful exit")
    parser.add_argument("--threads", type=int, default=12, help="Set PyTorch intra-op thread count")
    return parser.parse_args()

def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', s)]

def format_time(seconds):
    seconds = int(seconds)
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    if h > 0:
        return f"{h}h {m}m {s}s"
    elif m > 0:
        return f"{m}m {s}s"
    return f"{s}s"

def estimate_remaining_time(completed_folds, remaining_lengths):
    if not completed_folds or len(completed_folds) < 3 or not remaining_lengths:
        return None

    total_time = sum(t for _, t in completed_folds)
    total_aa = sum(n for n, _ in completed_folds)

    if total_aa == 0:
        return None

    seconds_per_aa = total_time / total_aa
    raw_eta = seconds_per_aa * sum(remaining_lengths)

    if raw_eta > 172800:  # 48 hour cap guard
        return None

    return raw_eta

def translate_dna(dna_seq):
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    dna_seq = dna_seq.upper()
    protein = []
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3]
        protein.append(table.get(codon, 'X'))
    return "".join(protein).rstrip('_')

def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if name: yield name, "".join(seq)
                name = line[1:]
                seq = []
            else:
                seq.append(line)
        if name: yield name, "".join(seq)

def get_clean_target_name(header, filename):
    target_part = header.split("||")[-1] if "||" in header else header
    if "_start=" in target_part:
        clean_name = target_part.split("_start=")[0]
    else:
        clean_name = target_part.split("_#_")[0].replace(" ", "_")

    if not clean_name or len(clean_name) > 30:
        clean_name = os.path.splitext(os.path.basename(filename))[0].split("_start=")[0]

    return clean_name

def main():
    args = parse_args()
    torch.set_num_threads(args.threads)

    pipeline_start = time.time()
    max_seconds = args.max_hours * 3600.0 if args.max_hours else None

    if os.path.isdir(args.input):
        input_files = sorted(
            glob.glob(os.path.join(args.input, "*.fna")) +
            glob.glob(os.path.join(args.input, "*.faa")) +
            glob.glob(os.path.join(args.input, "*.fasta")),
            key=natural_sort_key
        )
    else:
        input_files = sorted(glob.glob(args.input), key=natural_sort_key)

    if not input_files:
        print(f"[ERROR] No valid sequence files (.fna, .faa, .fasta) found matching: {args.input}")
        sys.exit(1)

    os.makedirs(args.pdb_dir, exist_ok=True)
    os.makedirs(args.output_dir, exist_ok=True)

    all_candidates = []
    no_hit_targets = []
    prev_fin_candidates = []
    not_folded_candidates = []

    for file_path in input_files:
        try:
            for header, seq in read_fasta(file_path):
                clean_name = get_clean_target_name(header, file_path)
                all_candidates.append((clean_name, seq, file_path, header))
        except Exception as e:
            print(f"  [ERROR] Failed to read FASTA file {file_path}: {e}")

    total_candidates = len(all_candidates)
    if total_candidates == 0:
        print(f"[ERROR] No valid sequence records found in input files.")
        sys.exit(1)


    print("=" * 70)
    print(f"  [BATCH PIPELINE INITIALIZED] Total Candidate Files : {len(input_files)}")
    print(f"  [BATCH PIPELINE INITIALIZED] Total Sequence Records : {total_candidates}")
    print(f"  [COMPUTE] PyTorch Intra-op Threads : {args.threads}")
    if args.limit is not None:
        print(f"  [QUOTA] Target Limit (Fresh Folds): {args.limit}")
    if args.max_hours is not None:
        print(f"  [TIME-BUDGET] Max Runtime Window: {args.max_hours} hours")
    print("=" * 70)

    model = None
    model_loaded = False
    model_load_duration = 0.0
    performance_log = []
    completed_fold_records = []
    completed_target_names = []

    last_target_expected = None
    batch_tsvs = []
    folded_count = 0

# TODO: #1 in-place updating instead of QC REJECTED spam for short peptides
    for idx, (clean_name, seq, file_path, header) in enumerate(all_candidates, 1):
        if max_seconds is not None and (time.time() - pipeline_start) >= max_seconds:
            print(f"\n  [TIME-BUDGET EXPIRED] Reached time limit of {args.max_hours} hours. Ending batch gracefully...")
            break

        if args.limit is not None and folded_count >= args.limit:
            print(f"\n  [QUOTA REACHED] Successfully completed requested limit of {args.limit} fresh folds. Halting batch.")
            break

        target_start = time.time()
        print(f"\n--- Scanning Target {idx}/{total_candidates}: {clean_name} (from {os.path.basename(file_path)}) ---")

        output_pdb = os.path.join(args.pdb_dir, f"{clean_name}.pdb")
        results_tsv = os.path.join(args.output_dir, f"{clean_name}.tsv")

        is_dna = set(seq.upper()).issubset({"A", "C", "G", "T", "N"})
        if is_dna:
            seq = translate_dna(seq)

        aa_len = len(seq)

        if aa_len == 0:
            print(f"  [ERROR] Target '{clean_name}' resulted in a 0-length sequence. Skipping.")
            performance_log.append((clean_name, 0, "FAILED (Empty Seq)", "Skipped", 0.0))
            continue

        if aa_len < args.min_length:
            print(f" -> [QC REJECTED] Length ({aa_len} aa) is below minimum threshold of {args.min_length} aa. Skipped.")
            performance_log.append((clean_name, aa_len, "Skipped (Too Short)", "Skipped", 0.0))
            continue

        print(f"  -> Clean Name  : {clean_name}")
        print(f"  -> Length      : {aa_len} amino acids")

        pdb_exists = os.path.exists(output_pdb) and os.path.getsize(output_pdb) > 0
        tsv_exists = os.path.exists(results_tsv)

        if pdb_exists and tsv_exists:
            print(f"  -> [CACHE HIT] PDB and Foldseek results already exist. Skipping computation.")
            batch_tsvs.append(results_tsv)
            completed_target_names.append(clean_name)
            performance_log.append((clean_name, aa_len, "Cached (0s)", "Cached (0s)", 0.0))
            continue

        fold_time = 0.0
        if pdb_exists:
            print(f"  -> [CACHE HIT] PDB structure already exists. Skipping folding.")
            fold_status = "Cached"
        else:
            if args.max_length is not None and aa_len > args.max_length:
                print(f"  -> [LENGTH FILTER] Target exceeds max length ({aa_len}aa > {args.max_length}aa). Skipped.")
                performance_log.append((clean_name, aa_len, "Skipped (Length)", "Skipped", 0.0))
                continue

            if aa_len > args.threshold:
                print("\a\a\a", end="", flush=True)
                print("!" * 66)
                print(f"  [ATTENTION REQUIRED] Target '{clean_name}' is {aa_len} aa")
                print(f"  exceeds your threshold of {args.threshold} aa.")
                print("!" * 66)

                choice = input(f"  Proceed with folding '{clean_name}'? [y/n]: ").strip().lower()
                if choice not in ['y', 'yes']:
                    print("  -> Skipped by user.")
                    performance_log.append((clean_name, aa_len, "Skipped", "Skipped", 0.0))
                    continue

            folded_count += 1
            print(f"  -> [COMPUTE SLOT {folded_count}/{args.limit if args.limit else 'All'}] Running fresh ESMFold prediction.")

            next_target_expected_secs = None
            if completed_fold_records:
                remaining_unpredicted_lens = []
                true_future_candidates = all_candidates[idx:]

                slots_remaining_in_batch = args.limit - folded_count if args.limit is not None else float('inf')

                for f_clean, f_seq, _, _ in true_future_candidates:
                    if len(remaining_unpredicted_lens) >= slots_remaining_in_batch:
                        break
                    f_is_dna = set(f_seq.upper()).issubset({"A", "C", "G", "T", "N"})
                    if f_is_dna:
                        f_seq = translate_dna(f_seq)
                    f_len = len(f_seq)
                    f_pdb = os.path.join(args.pdb_dir, f"{f_clean}.pdb")

                    if not (os.path.exists(f_pdb) and os.path.getsize(f_pdb) > 0):
                        remaining_unpredicted_lens.append(f_len)

                if remaining_unpredicted_lens:
                    recent_records = completed_fold_records[-5:] if len(completed_fold_records) > 1 else completed_fold_records
                    eta_secs = estimate_remaining_time(recent_records, remaining_unpredicted_lens)
                    if eta_secs is not None and eta_secs > 0:
                        print(f"     [ETA] Projected time for remaining unpredicted targets: ~{format_time(eta_secs)}")
                    next_target_expected_secs = estimate_remaining_time(recent_records, [remaining_unpredicted_lens[0]])
                else:
                    print(f"     [ETA] Final target in queue. Wrapping up computation...")

            device = torch.device("cuda" if torch.cuda.is_available() else ("xpu" if hasattr(torch, "xpu") and torch.xpu.is_available() else "cpu"))

            if not model_loaded:
                print("  -> Loading ESMFold model into memory...")
                load_start = time.time()
                try:
                    model = esm.pretrained.esmfold_v1()
                    model = model.eval().to(device)
                except Exception as e:
                    print(f"  [CRITICAL ERROR] Failed to load ESMFold model: {e}")
                    sys.exit(1)
                model_load_duration = time.time() - load_start
                print(f"  -> Model loaded successfully in {format_time(model_load_duration)}")
                model_loaded = True

            t_fold_start = time.time()
            print(f"  -> Folding via ESMFold...")

            try:
                with torch.no_grad():
                    output = model.infer(seq)

                pdb_lines = model.output_to_pdb(output)
                pdb_string = "\n".join(pdb_lines) if isinstance(pdb_lines, list) else pdb_lines

                with open(output_pdb, "w") as f:
                    f.write(pdb_string)

                fold_time = time.time() - t_fold_start
                fold_status = "Computed"
                completed_fold_records.append((aa_len, fold_time))
                completed_target_names.append(clean_name)

                if fold_time > 0:
                    rate_sec = aa_len / fold_time
                    rate_str = f"{rate_sec:.2f} aa/s" if rate_sec >= 1.0 else f"{rate_sec * 60.0:.2f} aa/min"
                    print(f"  -> Saved PDB: {output_pdb}")
                    print(f"     [Metrics] Length: {aa_len} aa | Time: {format_time(fold_time)} | Rate: {rate_str}")

                    if last_target_expected is not None:
                        variance = fold_time - last_target_expected
                        variance_pct = (abs(variance) / fold_time) * 100 if fold_time > 0 else 0
                        direction = "slower than expected" if variance > 0 else "faster than expected"

                        print(f"     [Accuracy] Predicted: {format_time(last_target_expected)} | Actual: {format_time(fold_time)}")
                        if abs(variance) < 2:
                            print(f"     [Accuracy] Prediction was spot on! (Within 2 seconds)")
                        else:
                            print(f"     [Accuracy] Model was {format_time(abs(variance))} {direction} ({variance_pct:.1f}% error)")
                else:
                    print(f"  -> Saved PDB: {output_pdb} (Completed in < 1s)")

                last_target_expected = next_target_expected_secs

            except (torch.cuda.OutOfMemoryError, RuntimeError) as oom_err:
                print("\a", end="", flush=True)
                print("!" * 70)
                print(f"  [VRAM EXHAUSTION / OOM ERROR] Target '{clean_name}' ({aa_len}aa) failed!")
                print(f"  Details: {oom_err}")
                print("!" * 70)
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
                performance_log.append((clean_name, aa_len, "FAILED (OOM)", "Skipped", time.time() - target_start))
                last_target_expected = None
                continue
            except Exception as gen_err:
                print(f"  [ERROR] Unexpected folding failure for {clean_name}: {gen_err}")
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
                performance_log.append((clean_name, aa_len, "FAILED (Error)", "Skipped", time.time() - target_start))
                last_target_expected = None
                continue

        # Search Phase
        search_time = 0.0
        if tsv_exists:
            print(f"  -> [CACHE HIT] Foldseek results already exist. Skipping alignment.")
            search_status = "Cached"
        else:
            print(f"  -> Running Foldseek search...")
            t_search_start = time.time()
            tmp_dir = f"tmp_foldseek_{clean_name}"
            foldseek_cmd = ["foldseek", "easy-search", output_pdb, args.db, results_tsv, tmp_dir, "--format-output", "query,target,pident,evalue,bits,theader,alntmscore"]

            try:
                f_res = subprocess.run(foldseek_cmd, capture_output=True, text=True, check=False)
                search_status = "Computed" if f_res.returncode == 0 else "Failed"
            except FileNotFoundError:
                print(f"  [CRITICAL ERROR] 'foldseek' binary not found in system PATH!")
                search_status = "Failed (Missing Binary)"
            except Exception as e:
                print(f"  [ERROR] Unexpected error during Foldseek: {e}")
                search_status = "Failed (Exception)"

            if os.path.exists(tmp_dir):
                try:
                    subprocess.run(["rm", "-rf", tmp_dir], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                except Exception:
                    pass

            search_time = time.time() - t_search_start

            if search_status == "Computed" and os.path.exists(results_tsv) and os.path.getsize(results_tsv) == 0:
                os.remove(results_tsv)
                no_hit_targets.append(clean_name)
                print(f"  -> Foldseek search complete: 0 structural matches found. Discarded empty TSV.")

            else:
                print(f"  -> Foldseek search complete. Saved to {results_tsv} (Completed in {format_time(search_time)})")

        batch_tsvs.append(results_tsv)
        total_target_time = time.time() - target_start
        performance_log.append((clean_name, aa_len, f"{fold_status} ({format_time(fold_time)})", f"{search_status} ({format_time(search_time)})", total_target_time))

    if args.limit is not None and folded_count < args.limit:
        print(f"\n  [NOTICE] Reached end of candidate pool. Processed {folded_count} of {args.limit} requested fresh folds.")

    # Write out batch TSVs manifest
    manifest_path = os.path.join(os.path.dirname(args.output_dir), "latest_batch_tsvs.txt")
    with open(manifest_path, "w") as mf:
        for tsv in batch_tsvs:
            mf.write(f"{tsv}\n")

    # Write out completion tracking log for parent orchestrator
    completed_log_path = os.path.join(os.path.dirname(args.output_dir), "recently_completed.txt")
    with open(completed_log_path, "w") as clf:
        for c_name in completed_target_names:
            clf.write(f"{c_name}\n")

    total_duration = time.time() - pipeline_start
    print("\n" + "=" * 80)
    print(f"  [PIPELINE PERFORMANCE TELEMETRY SUMMARY]")
    print("=" * 80)
    if model_loaded:
        print(f"  * ESMFold Model Cold-Load Time : {format_time(model_load_duration)}")
    print(f"  * Total Batch Elapsed Time     : {format_time(total_duration)}")
    print("-" * 80)
    print(f"  {'Target Name':<20} | {'Length':<8} | {'Folding Status/Time':<22} | {'Search Status/Time':<22}")
    print("-" * 80)
    for name, length, f_stat, s_stat, _ in performance_log:
        if f_stat.startswith("Computed"):
            print(f"  {name:<20} | {length:<8} | {f_stat:<22} | {s_stat:<22}")
    print("=" * 80)

if __name__ == "__main__":
    main()
