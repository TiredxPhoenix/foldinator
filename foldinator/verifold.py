import os
import argparse
import datetime
import sys

def main():
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_ROOT = os.path.dirname(BASE_DIR)

    default_manifest = os.path.join(PROJECT_ROOT, "OUTPUTS", "latest_batch_tsvs.txt")
    default_output_dir = os.path.join(PROJECT_ROOT, "OUTPUTS")

    parser = argparse.ArgumentParser(description="Audit Foldseek results and generate structural dark matter report.")
    parser.add_argument("--batch_file", type=str, default=default_manifest, help="Path to batch manifest of TSVs")
    parser.add_argument("--output_dir", type=str, default=default_output_dir, help="Directory to save audit reports")
    parser.add_argument("--evalue", type=float, default=1e-3, help="E-value threshold for significant structural hits")
    args = parser.parse_args()

    # Guardrail check using the correct argument name
    if not os.path.exists(args.batch_file):
        print(f"[ERROR] Batch manifest not found at: {args.batch_file}")
        print("Run your fold-and-seek pipeline first to generate alignment targets.")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = os.path.join(args.output_dir, f"verifold_report_{timestamp}.tsv")

    # Load TSV paths from manifest
    tsv_files = []
    with open(args.batch_file, "r") as f:
        for line in f:
            path = line.strip()
            if path:
                tsv_files.append(path)

    if not tsv_files:
        print(f"[WARNING] Manifest file is empty: {args.batch_file}")
        sys.exit(0)

    print(f"=" * 65)
    print(f"  [VERIFOLD AUDIT] Processing {len(tsv_files)} target results...")
    print(f"=" * 65)

    rows = []
    for tsv_path in tsv_files:
        target_name = os.path.basename(tsv_path).replace("foldseek_results_", "").replace(".tsv", "")

        # Check if file exists and contains data (handles structural dark matter / 0-byte files)
        if not os.path.exists(tsv_path) or os.path.getsize(tsv_path) == 0:
            rows.append({
                "Target": target_name,
                "Hit_ID": "NO_HIT",
                "E_value": "N/A",
                "Bitscore": "N/A",
                "Description": "Structural dark matter (no alignment detected)"
            })
            continue

        # Parse top hit from Foldseek easy-search output
        best_hit = None
        min_evalue = float('inf')

        try:
            with open(tsv_path, "r") as tf:
                for line in tf:
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        hit_id = parts[1]

                        try:
                            evalue = float(parts[3])
                        except ValueError:
                            evalue = 1.0

                        bitscore = parts[4] if len(parts) > 4 else "0.0"
                        description = parts[5] if len(parts) > 5 else "No description"

                        if evalue < min_evalue:
                            min_evalue = evalue
                            best_hit = {
                                "Target": target_name,
                                "Hit_ID": hit_id,
                                "E_value": f"{evalue:.2e}" if evalue != 1.0 else "1.0",
                                "Bitscore": bitscore,
                                "Description": description
                            }
        except Exception as e:
            print(f"  [WARNING] Error parsing {tsv_path}: {e}")
            continue

        if best_hit and min_evalue <= args.evalue:
            rows.append(best_hit)
        else:
            rows.append({
                "Target": target_name,
                "Hit_ID": "BELOW_THRESHOLD" if best_hit else "NO_ALIGNMENT",
                "E_value": f"{min_evalue:.2e}" if best_hit else "N/A",
                "Bitscore": "N/A",
                "Description": f"Best hit exceeded threshold ({args.evalue})" if best_hit else "No significant structural match"
            })

    # Write clean TSV report with threshold in column header
    with open(report_path, "w", encoding="utf-8") as rf:
        rf.write(f"Target\tTop_Hit_ID\tE-value ({args.evalue:g})\tBitscore\tProtein_Description\n")

        for r in rows:
            rf.write(f"{r['Target']}\t{r['Hit_ID']}\t{r['E_value']}\t{r['Bitscore']}\t{r['Description']}\n")

    print(f"  -> Audit complete.")
    print(f"  -> Report saved to: {report_path}")
    print(f"=" * 65)

if __name__ == "__main__":
    main()
