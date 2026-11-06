#!/usr/bin/env python3
"""
Foldinator Structural Atlas Generator
Generates a self-contained, high-performance interactive HTML dashboard
for viral metagenome contigs, structural predictions, and dark matter loci.
"""

import argparse
import base64
import glob
import json
import os
import sys


def generate_dataset_from_pipeline(output_dir):
    pdbs_dir = os.path.join(output_dir, "pdbs")
    tsv_dir = os.path.join(output_dir, "foldseek_results")

    data = []
    if not os.path.exists(pdbs_dir):
        print("Could not load pdb directory (line 486 foldinator_atlas.py)")
        return None

    # Recursively find all PDBs (including quarantined ones if needed)
    pdb_files = []
    for root, _, files in os.walk(pdbs_dir):
        for f in files:
            if f.endswith(".pdb"):
                pdb_files.append(os.path.join(root, f))

    if not pdb_files:
        return None

    print(f"[INFO] Compiling atlas dataset from {len(pdb_files)} PDB structures...")

    for pdb_path in pdb_files:
        target_id = os.path.splitext(os.path.basename(pdb_path))[0]

        plddt_vals = []
        aa_count = 0
        pdb_content = ""
        try:
            with open(pdb_path, "r", encoding="utf-8") as f:
                pdb_content = f.read()

            for line in pdb_content.splitlines():
                if line.startswith("ATOM") and (len(line) < 14 or line[12:16].strip() == "CA"):
                    aa_count += 1
                    try:
                        bfactor = float(line[60:66].strip())
                        plddt_vals.append(bfactor)
                    except ValueError:
                        pass
        except Exception:
            continue

        mean_plddt = f"{sum(plddt_vals)/len(plddt_vals):.1f}" if plddt_vals else "N/A"

        tsv_path = os.path.join(tsv_dir, f"{target_id}.tsv")
        homology_status = "Structural Dark Matter"
        best_evalue = "None"

        if os.path.exists(tsv_path) and os.path.getsize(tsv_path) > 0:
            try:
                with open(tsv_path, "r", encoding="utf-8") as tf:
                    first_line = tf.readline().strip()
                    if first_line:
                        parts = first_line.split("\t")
                        if len(parts) >= 4:
                            best_evalue = parts[3]
                            desc = parts[5] if len(parts) > 5 else parts[1]
                            homology_status = f"Hit: {desc[:30]}..."
            except Exception:
                pass

        pdb_b64 = base64.b64encode(pdb_content.encode('utf-8')).decode('utf-8')

        data.append({
            "id": target_id,
            "length_aa": aa_count if aa_count > 0 else 100,
            "mean_plddt": mean_plddt,
            "checkv_quality": homology_status,
            "dark_matter_loci": best_evalue,
            "description": f"Protein target {target_id} modeled via ESMFold & Foldseek.",
            "has_pdb": True,
            "pdb_b64": pdb_b64
        })

    return data if data else None

def main():
    parser = argparse.ArgumentParser(description="Generate Foldinator Interactive Structural Atlas HTML.")
    parser.add_argument("--output_dir", help="Directory to output the HTML file.", default="../foldinator_outputs")
    parser.add_argument("--out_html", help="Path to output HTML file (pipeline flag).", default=None)
    args = parser.parse_args()

    output_path = "foldinator_structural_atlas.html"
    if args.out_html:
        output_path = args.out_html
    elif args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)
        output_path = os.path.join(args.output_dir, "foldinator_atlas.html")

    target_dir = args.output_dir if args.output_dir else os.path.dirname(output_path)
    data = dataset = generate_dataset_from_pipeline(target_dir)

    if not data:
        print("[WARNING] No PDBs found to compile into the atlas.")
        data = []

    # Read the local ngl.js file and embed it directly into the template
    ngl_js_path = os.path.join(os.path.dirname(__file__), "assets", "ngl.js")
    if not os.path.exists(ngl_js_path):
        # Fallback search relative to working directory if script path differs
        ngl_js_path = "assets/ngl.js"

    try:
        with open(ngl_js_path, "r", encoding="utf-8") as f:
            ngl_script_content = f.read()
    except Exception as e:
        print(f"[ERROR] Could not read local ngl.js from {ngl_js_path}: {e}")
        ngl_script_content = "// ERROR: ngl.js could not be loaded during compilation."

    # Path to your new external HTML template file
    template_path = os.path.join(os.path.dirname(__file__), "atlas_template.html")
    if not os.path.exists(template_path):
        template_path = "atlas_template.html" # Fallback relative check

    try:
        with open(template_path, "r", encoding="utf-8") as tf:
            html_template = tf.read()
    except Exception as e:
        print(f"[ERROR] Could not read atlas_template.html from {template_path}: {e}")
        sys.exit(1)

    data_json_str = json.dumps(data)
    final_html = html_template.replace("__EMBEDDED_NGL_JS__", ngl_script_content)
    final_html = final_html.replace("__ATLAS_DATA_JSON__", data_json_str)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(final_html)

    print(f"Successfully generated structural atlas: {os.path.abspath(output_path)}")

if __name__ == "__main__":
    main()
