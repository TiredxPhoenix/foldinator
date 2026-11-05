import os

# Force CPU fallback for OpenFold C++ extensions before any library imports
os.environ["OPENFOLD_USE_CPP_OPS"] = "0"

import sys
import subprocess
import time
import datetime
import platform
import glob
import argparse
import shutil

# Import our dedicated QC module
from qc_cop import audit_and_clean_tsv

def get_tool_versions():
    """Dynamically queries underlying tool and library versions."""
    versions = {
        "Python": platform.python_version(),
        "PyTorch": "Unknown",
        "PyMOL": "Unknown",
        "Foldseek": "Unknown"
    }

    try:
        import torch
        versions["PyTorch"] = torch.__version__
    except ImportError:
        pass

    try:
        import pymol
        versions["PyMOL"] = getattr(pymol, "__version__", "Open-Source")
    except ImportError:
        pass

    try:
        res = subprocess.run(["foldseek", "version"], capture_output=True, text=True, timeout=2)
        if res.returncode == 0:
            versions["Foldseek"] = res.stdout.strip()
    except Exception:
        pass

    return versions

def get_bin_counts(lengths, num_bins=10):
    """Computes raw counts per bin for delta tracking."""
    if not lengths:
        return []
    min_val = min(lengths)
    max_val = max(lengths)
    if min_val == max_val:
        return [len(lengths)]
    bin_width = (max_val - min_val) / num_bins
    bins = [0] * num_bins
    for length in lengths:
        bin_idx = int((length - min_val) / bin_width)
        if bin_idx >= num_bins:
            bin_idx = num_bins - 1
        bins[bin_idx] += 1
    return bins

def print_terminal_histogram(lengths, initial_bins=None, title="DATASET PROTEIN LENGTH DISTRIBUTION", num_bins=10, max_bar_width=40):
    """Renders a text-based ASCII histogram of protein lengths directly to the terminal with optional progress deltas."""
    if not lengths:
        print(f"  [Histogram] No sequences available to plot.")
        return

    min_val = min(lengths)
    max_val = max(lengths)

    if min_val == max_val:
        print(f"  [Histogram] All proteins share the exact same length: {min_val} aa")
        return

    bin_width = (max_val - min_val) / num_bins
    bins = [0] * num_bins

    for length in lengths:
        bin_idx = int((length - min_val) / bin_width)
        if bin_idx >= num_bins:
            bin_idx = num_bins - 1
        bins[bin_idx] += 1

    max_count = max(bins) if bins else 1

    print("\n" + "  " + "=" * 60)
    print(f"  [{title}] (Total Active Remaining: {len(lengths)})")
    print("  " + "=" * 60)

    for i in range(num_bins):
        bin_start = int(min_val + (i * bin_width))
        bin_end = int(min_val + ((i + 1) * bin_width))
        count = bins[i]

        bar_len = int((count / max_count) * max_bar_width) if max_count > 0 else 0
        bar = "█" * bar_len

        if initial_bins and i < len(initial_bins):
            diff = count - initial_bins[i]
            if diff != 0:
                count_str = f"{initial_bins[i]} {diff}"
            else:
                count_str = f"{count}"
        else:
            count_str = f"{count}"

        print(f"  {bin_start:5d} - {bin_end:5d} aa | {bar} ({count_str})")

    print("  " + "=" * 60)

def display_banner():
    v = get_tool_versions()
    print("=" * 80)
    print("  [THE FOLDINATOR] Automated Structural Pipeline Orchestrator")
    print("=" * 80)
    print("  Environment & Tool Versions:")
    print(f"    - Python   : {v['Python']}")
    print(f"    - PyTorch  : {v['PyTorch']}")
    print(f"    - PyMOL    : {v['PyMOL']}")
    print(f"    - Foldseek : {v['Foldseek']}")
    print("-" * 80)
    print("  Pipeline Architecture:")
    print("    [Raw FASTA Inputs / OrfanNapper]")
    print("        ├── Prodigal Metagenomic Gene Prediction (CDS)")
    print("        ├── Rational Lineage & Header Standardization")
    print("        └── MMseqs2 Redundancy Sieve & Representative Filtering")
    print("        │")
    print("        ▼")
    print("    [Tool A: fold_and_seek.py]")
    print("        ├── VRAM-cached ESMFold Structure Prediction")
    print("        └── Foldseek Homology Search (with TM-scores)")
    print("        │")
    print("        ▼")
    print("    [Tool B: qc_cop.py] (Statistical Gate & PDB Quarantine)")
    print("        │")
    print("        ▼")
    print("    [Tool C: verifold.py]")
    print("        └── Executive Summary Audit Report")
    print("=" * 80)

def sniff_input_type(input_path):
    target_file = input_path
    if "*" in input_path:
        matches = glob.glob(input_path)
        if not matches:
            print(f"[ERROR] No files match pattern: {input_path}")
            sys.exit(1)
        target_file = matches[0]

    if not os.path.exists(target_file):
        print(f"[ERROR] Path does not exist on disk: {target_file}")
        sys.exit(1)

    ext = os.path.splitext(target_file)[1].lower()

    if ext in [".fna", ".fasta", ".fa"]:
        with open(target_file, "r") as f:
            for line in f:
                if not line.startswith(">"):
                    seq = line.strip().upper()
                    aa_chars = set("EQFLIPM")
                    if any(char in seq for char in aa_chars) and len(set(seq).intersection(set("ATCGN"))) < 4:
                        return "protein_unclustered"
                    break
        return "nucleotide_genomes"
    elif ext in [".faa"]:
        return "protein_unclustered"
    else:
        return "nucleotide_genomes"

def get_clean_target_name(header):
    target_part = header.split("||")[-1] if "||" in header else header
    if "_start=" in target_part:
        clean_name = target_part.split("_start=")[0]
    else:
        clean_name = target_part.split("_#_")[0].replace(" ", "_")
    return clean_name

def process_and_filter_manifest(fasta_path, pdb_dir, max_length_filter=None):
    sequences = []
    current_header = None
    current_seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_header:
                    sequences.append((current_header, "".join(current_seq)))
                current_header = line.strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_header:
            sequences.append((current_header, "".join(current_seq)))

    total_raw_count = len(sequences)

    existing_pdbs = set()
    if os.path.exists(pdb_dir):
        for fname in os.listdir(pdb_dir):
            if fname.endswith(".pdb"):
                existing_pdbs.add(fname[:-4])

    filtered_sequences = []
    qc_rejected_count = 0
    cached_skipped_count = 0

    for header, seq in sequences:
        protein_id = get_clean_target_name(header.strip(">"))

        if protein_id in existing_pdbs:
            cached_skipped_count += 1
            continue

        if len(seq) < 30:
            qc_rejected_count += 1
            continue

        if max_length_filter is not None and len(seq) > max_length_filter:
            continue

        filtered_sequences.append((header, seq))

    filtered_sequences.sort(key=lambda x: len(x[1]))
    active_lengths = [len(seq) for _, seq in filtered_sequences]
    histogram_title = "SUB-500AA LENGTH DISTRIBUTION (CAPPED)" if max_length_filter == 500 else "FULL DATASET LENGTH DISTRIBUTION"

    print_terminal_histogram(active_lengths, initial_bins=None, title=histogram_title, num_bins=10, max_bar_width=35)

    active_manifest_path = fasta_path.replace(".faa", "_active_sorted.faa")
    with open(active_manifest_path, "w") as f:
        for header, seq in filtered_sequences:
            f.write(f"{header}\n")
            for i in range(0, len(seq), 60):
                f.write(f"{seq[i:i+60]}\n")

    return active_manifest_path, total_raw_count, len(filtered_sequences), qc_rejected_count, cached_skipped_count, active_lengths

def main():
    parser = argparse.ArgumentParser(description="Foldinator Bioinformatics Pipeline")
    parser.add_argument("--input", type=str, default=None, help="Input path or glob pattern")
    parser.add_argument("--output_dir", type=str, default=None, help="Base output directory")
    parser.add_argument("--max_targets", type=int, default=None, help="Max targets to process")
    parser.add_argument("--max_hours", type=float, default=None, help="Max runtime hours before graceful exit")
    parser.add_argument("--threads", type=int, default=12, help="Set PyTorch thread count")
    parser.add_argument("--large-proteins", type=str, choices=["prompt", "skip", "fold"], default="prompt", help="How to handle proteins >500aa")

    args = parser.parse_args()

    display_banner()

    MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_ROOT = os.path.dirname(os.path.dirname(MODULE_DIR))

    FOLD_SEEK_PATH = os.path.join(MODULE_DIR, "fold_and_seek.py")
    VERIFOLD_PATH = os.path.join(MODULE_DIR, "verifold.py")

    default_input = os.path.join(PROJECT_ROOT, "DATA", "extracted_fastas", "full_*.fna")
    if args.input:
        input_target = args.input
    else:
        input_prompt = input(f"\n  -> Enter input path / glob pattern [default: DATA/extracted_fastas/full_*.fna]: ").strip()
        input_target = input_prompt if input_prompt else default_input

    # Updated default output folder name to foldinator_outputs
    default_output = os.path.join(PROJECT_ROOT, "foldinator_outputs")
    if args.output_dir:
        base_output = args.output_dir
    else:
        output_prompt = input(f"  -> Enter base output directory [default: foldinator_outputs]: ").strip()
        base_output = output_prompt if output_prompt else default_output

    pdb_dir = os.path.join(base_output, "pdbs")
    foldseek_results_dir = os.path.join(base_output, "foldseek_results")
    verifold_manifest = os.path.join(base_output, "latest_batch_tsvs.txt")

    print("\n" + "=" * 80)
    print("  [PHASE 1/3] Inspecting Input & Executing Upstream Prep")
    print("=" * 80)

    input_provenance = sniff_input_type(input_target)
    print(f"  -> Detected input provenance: {input_provenance}")

    prep_output_dir = os.path.join(base_output, "upstream_prep")

    if input_provenance in ["nucleotide_genomes", "protein_unclustered"]:
        if input_provenance == "nucleotide_genomes":
            print("  -> Routing through full OrfanNapper pipeline (Prodigal + MMseqs2)...")
        else:
            print("  -> Raw proteins detected. Routing directly to MMseqs2 sieve...")

        from orfan_napper import run_orfan_napper
        raw_manifest = run_orfan_napper(
            input_fasta=input_target,
            output_dir=prep_output_dir,
            seq_identity=0.85
        )
    else:
        raw_manifest = input_target

    if args.max_hours is not None:
        max_hrs_val = args.max_hours
    else:
        hours_input = input(f"\n  -> Maximum runtime in hours (e.g. 0.25, 0.5, 2)? [Press Enter for unconstrained]: ").strip()
        max_hrs_val = float(hours_input) if hours_input else None

    if max_hrs_val is not None:
        max_length_filter = 500
        print(f" -> [AUTO-SKIP] Time limit chosen: ({max_hrs_val} hrs). Auto-skipping proteins >500aa.")
    else:
        while True:
            skip_input = input("  -> Skip proteins over 500 amino acids? (y/n): ").strip().lower()
            if skip_input in ['y', 'yes']:
                max_length_filter = 500
                break
            elif skip_input in ['n', 'no']:
                max_length_filter = None
                break
            else:
                print("  [ERROR] Please enter 'y' or 'n'.")

    active_manifest, total_raw, active_count, qc_rejected, cached_skipped, initial_lengths = process_and_filter_manifest(
        raw_manifest, pdb_dir, max_length_filter=max_length_filter
    )
    initial_bins = get_bin_counts(initial_lengths)

    if args.max_targets is not None:
        limit_val = args.max_targets
    elif max_hrs_val is not None:
        limit_val = None
    else:
        while True:
            limit_input = input(f"\n  -> How many proteins to process? (Enter a number up to {active_count}, or 'all'): ").strip().lower()
            if limit_input == 'all':
                limit_val = None
                break
            elif limit_input.isdigit():
                limit_val = int(limit_input)
                break
            else:
                print("  [ERROR] Please enter a valid integer or 'all'.")

    print("\n" + "-" * 80)
    print(f"  Configuration Confirmed:")
    print(f"    - Input Target      : {input_target}")
    print(f"    - Active Queue Size : {active_count} proteins ready")
    print(f"    - Output Folder     : {base_output}")
    print(f"    - Target Limit      : {limit_val if limit_val else 'All available'}")
    print(f"    - Time-Budget       : {f'{max_hrs_val} hours' if max_hrs_val else 'None'}")
    print(f"    - PyTorch Threads   : {args.threads}")
    print("-" * 80 + "\n")

    pipeline_start_time = time.time()

    try:
        fold_cmd = [
            sys.executable,
            FOLD_SEEK_PATH,
            "--input", active_manifest,
            "--pdb_dir", pdb_dir,
            "--output_dir", foldseek_results_dir
        ]

        if limit_val is not None:
            fold_cmd.extend(["--limit", str(limit_val)])
        if max_length_filter is not None:
            fold_cmd.extend(["--max_length", str(max_length_filter)])
        if max_hrs_val is not None:
            fold_cmd.extend(["--max_hours", str(max_hrs_val)])
        if args.threads is not None:
            fold_cmd.extend(["--threads", str(args.threads)])

        print("[PHASE 2/4] Launching Fold & Search Core (fold_and_seek.py)...")
        phase1_start = time.time()
        result1 = subprocess.run(fold_cmd)
        phase1_duration = time.time() - phase1_start

        if result1.returncode != 0:
            print("[FATAL ERROR] Fold and Seek phase failed. Aborting pipeline.")
            sys.exit(1)

        # PHASE 3: QC Cop Audit, Noise Purge & PDB Quarantine
        print(f"\n[PHASE 3/4] Launching QC Cop Audit & Quarantine (qc_cop.py)...")
        phase3_start = time.time()
        audited_count = 0

        if os.path.exists(verifold_manifest):
            with open(verifold_manifest, "r") as mf:
                for line in mf:
                    tsv_path = line.strip()
                    if tsv_path:
                        audit_and_clean_tsv(tsv_path, pdb_dir=pdb_dir)
                        audited_count += 1

        phase3_duration = time.time() - phase3_start
        print(f"  -> QC Cop audited {audited_count} result files in {phase3_duration:.2f} seconds.")

        # PHASE 4: Executive Structural Audit Reporting
        verifold_cmd = [
            sys.executable,
            VERIFOLD_PATH,
            "--batch_file", verifold_manifest,
            "--output_dir", base_output
        ]

        print(f"\n[PHASE 4/4] Launching Structural Audit Core (verifold.py)...")
        phase4_start = time.time()
        result2 = subprocess.run(verifold_cmd)
        phase4_duration = time.time() - phase4_start

        if result2.returncode != 0:
            print("[WARNING] Verifold audit completed with non-zero exit code.")
            sys.exit(result2.returncode)

        # PHASE 5: Compile Interactive HTML Atlas (foldinator_atlas.py)
        ATLAS_PATH = os.path.join(MODULE_DIR, "foldinator_atlas.py")
        output_html_path = os.path.join(base_output, "foldinator_atlas.html")

        atlas_cmd = [
            sys.executable,
            ATLAS_PATH,
            "--output_dir", base_output,
            "--out_html", output_html_path
        ]

        print(f"\n[PHASE 5/5] Compiling Interactive Atlas (foldinator_atlas.py)...")
        phase5_start = time.time()
        result_atlas = subprocess.run(atlas_cmd)
        phase5_duration = time.time() - phase5_start

        if result_atlas.returncode != 0:
            print("[WARNING] Atlas compilation completed with non-zero exit code.")
        else:
            print(f"  -> Atlas compiled successfully in {phase5_duration:.2f} seconds.")


        total_duration = time.time() - pipeline_start_time

        print("\n" + "=" * 80)
        print(f"  [FOLDINATOR TELEMETRY & BATCH SUMMARY]")
        print("=" * 80)
        print(f"  * Total Source Proteins Discovered : {total_raw}")
        print(f"  * QC Skipped (<30aa micro-peptides)  : {qc_rejected}")
        print(f"  * Previously Completed (Cached PDBs) : {cached_skipped}")
        print(f"  * Actively Processed This Run      : {min(limit_val, active_count) if limit_val else active_count}")
        print("-" * 80)
        print(f"  * Phase 2 (Fold & Search) Duration : {phase1_duration / 60:.2f} minutes")
        print(f"  * Phase 3 (QC Cop Audit)           : {phase3_duration:.2f} seconds")
        print(f"  * Phase 4 (Verifold Report)        : {phase4_duration:.2f} seconds")
        print(f"  * Total Pipeline Runtime           : {total_duration / 60:.2f} minutes")
        print("-" * 80)
        print(f"  * Audit Report Generated           : {os.path.join(base_output, 'latest_audit_report.md')}")
        print(f"  -> Explore your interactive structural atlas!")
        print(f"     Open file://{os.path.abspath(output_html_path)} in your browser.")
        print("-" * 80)
        print(f"  * Audit Report Generated           : {os.path.join(base_output, 'latest_audit_report.md')}")
        print(f"  * Completion Timestamp             : {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print("=" * 80)
        print("\a\a\a")

    except KeyboardInterrupt:
        print("\n\n[INTERRUPTED] Ctrl+C detected. Shutting down pipeline gracefully...")
        for tmp_folder in glob.glob("tmp_foldseek_*"):
            try:
                shutil.rmtree(tmp_folder)
            except Exception:
                pass
        print("  -> Cleanup complete. Exiting safely.")
        sys.exit(0)

if __name__ == "__main__":
    main()
