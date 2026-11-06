import os
import subprocess
import sys
import shutil
import time

def run_prodigal_with_file_monitor(cmd, output_faa):
    print("  -> Launching Prodigal gene prediction pipeline...")

    # Let Prodigal write directly to files without trapping pipes
    process = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    last_count = 0
    while process.poll() is None:
        if os.path.exists(output_faa):
            try:
                with open(output_faa, "r") as f:
                    current_count = sum(1 for line in f if line.startswith(">"))
                if current_count != last_count:
                    last_count = current_count
                    sys.stdout.write(f"\r  -> Progress: {last_count} protein sequences predicted... \033[K")
                    sys.stdout.flush()
            except Exception:
                pass
        time.sleep(0.1)

    process.wait()

    final_count = 0
    if os.path.exists(output_faa):
        with open(output_faa, "r") as f:
            final_count = sum(1 for line in f if line.startswith(">"))

    sys.stdout.write(f"\r  -> Progress: {final_count} protein sequences predicted. Complete!          \033[K\n")
    sys.stdout.flush()

    if process.returncode != 0:
        print("[FATAL ERROR] Prodigal execution failed.")
        sys.exit(1)

    print(f"  -> Prodigal complete. Total CDS predicted: {final_count}")

def run_orfan_napper(input_fasta, output_dir, seq_identity=0.85):
    """
    OrfanNapper: Slices raw vOTU genomes into predicted proteins via Prodigal,
    enforces rational lineage naming, and slashes redundancy via MMseqs2.
    """
    os.makedirs(output_dir, exist_ok=True)

    raw_proteins = os.path.join(output_dir, "prodigal_raw_proteins.faa")
    rational_proteins = os.path.join(output_dir, "prodigal_rational_proteins.faa")
    cluster_prefix = os.path.join(output_dir, "mmseqs_cluster")
    rep_proteins = f"{cluster_prefix}_rep_seq.fasta"

    # Verify external dependencies are available in PATH
    for binary in ["prodigal", "mmseqs"]:
        if subprocess.run(["which", binary], capture_output=True).returncode != 0:
            print(f"[FATAL ERROR] Required binary '{binary}' not found in system PATH.")
            print(f"Please ensure {binary} is installed and active in your environment.")
            sys.exit(1)

    print(f"\n[ORFAN-NAPPER: PHASE 1A] Kidnapping ORFs via Prodigal on: {os.path.basename(input_fasta)}")

    # Clean, standard metagenomic Prodigal command
    prodigal_cmd = [
        "prodigal",
        "-i", input_fasta,
        "-a", raw_proteins,
        "-p", "meta"
    ]

    run_prodigal_with_file_monitor(prodigal_cmd, raw_proteins)

    print(f"[ORFAN-NAPPER: PHASE 1B] Enforcing rational lineage naming conventions & building translation map...")
    parsed_count = 0
    mapping_records = []

    with open(raw_proteins, "r") as inf, open(rational_proteins, "w") as outf:
        gene_counter = 0
        for line in inf:
            if line.startswith(">"):
                original_header = line[1:].strip()
                header = original_header
                # Extract base vOTU contig name prior to Prodigal coordinate markers
                votu_base = header.split(" # ")[0].rsplit("_#_", 1)[0].strip()
                if not votu_base or votu_base == header:
                    votu_base = header.split()[0]

                gene_counter += 1
                clean_name = f"{votu_base}_gene_{gene_counter}"
                new_header = f">{clean_name}\n"
                outf.write(new_header)

                mapping_records.append((clean_name, original_header))
                parsed_count += 1
            else:
                outf.write(line)

    # Write out header translation lookup map
    mapping_path = os.path.join(output_dir, "header_translation_map.tsv")
    with open(mapping_path, "w") as map_f:
        map_f.write("Clean_Name\tOriginal_Header\n")
        for c_name, o_head in mapping_records:
            map_f.write(f"{c_name}\t{o_head}\n")

    print(f"  -> Successfully processed and re-tagged {parsed_count} protein coding sequences.")
    print(f"  -> Header translation map saved to: {mapping_path}")

    print(f"[ORFAN-NAPPER: PHASE 1C] Running MMseqs2 redundancy sieve (Clustering at {seq_identity} amino acid identity to remove near-duplicates)...")
    mmseqs_cmd = [
        "mmseqs", "easy-cluster", rational_proteins, cluster_prefix,
        os.path.join(output_dir, "tmp_mmseqs"),
        "--min-seq-id", str(seq_identity),
        "--cov-mode", "1",
        "-v", "0"
    ]
    res = subprocess.run(mmseqs_cmd)
    if res.returncode != 0:
        print("[ERROR] MMseqs2 clustering failed.")
        sys.exit(1)

    # Clean up temporary MMseqs scratch directory to save space
    tmp_mmseqs_dir = os.path.join(output_dir, "tmp_mmseqs")
    if os.path.exists(tmp_mmseqs_dir):
        try:
            shutil.rmtree(tmp_mmseqs_dir)
        except Exception:
            pass

    if os.path.exists(rep_proteins):
        rep_count = sum(1 for line in open(rep_proteins) if line.startswith(">"))
        print(f"  -> All orfans kidnapped and in the bag!")
        print(f"  -> Compressed redundant protein dataset down to {rep_count} unique representative folds.")
        print(f"  -> Clean manifest handed off to: {rep_proteins}")
        return rep_proteins
    else:
        print("[ERROR] OrfanNapper failed to generate the expected MMseqs representative file.")
        sys.exit(1)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="OrfanNapper: Upstream vOTU gene prediction and clustering.")
    parser.add_argument("--input", required=True, help="Path to input vOTU genomic fasta (viruses.fna)")
    parser.add_argument("--output_dir", required=True, help="Directory to store intermediate and final outputs")
    parser.add_argument("--identity", type=float, default=0.85, help="MMseqs2 sequence identity threshold (default: 0.85)")

    args = parser.parse_args()
    run_orfan_napper(args.input, args.output_dir, args.identity)
