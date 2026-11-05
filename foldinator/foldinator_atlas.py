#!/usr/bin/env python3
"""
Foldinator Structural Atlas Generator
Generates a self-contained, high-performance interactive HTML dashboard
for deep-sea Arctic viral contigs, structural predictions, and dark matter loci.
"""

import argparse
import base64
import glob
import json
import os
import sys

HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Foldinator Structural Atlas - Arctic Viral Metagenomics</title>
    <!-- NGL Viewer Library -->
    <script src="https://cdn.jsdelivr.net/npm/ngl@2.0.0/dist/ngl.js"></script>
    <style>
        :root {
            --bg-dark: #0f172a;
            --bg-panel: #1e293b;
            --bg-card: #334155;
            --accent: #38bdf8;
            --accent-hover: #0ea5e9;
            --text-main: #f8fafc;
            --text-muted: #94a3b8;
            --border: #475569;
            --success: #22c55e;
            --warning: #f59e0b;
        }

        * {
            box-sizing: border-box;
            margin: 0;
            padding: 0;
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
        }

        body {
            background-color: var(--bg-dark);
            color: var(--text-main);
            display: flex;
            height: 100vh;
            overflow: hidden;
        }

        aside {
            width: 380px;
            background-color: var(--bg-panel);
            border-right: 1px solid var(--border);
            display: flex;
            flex-direction: column;
            height: 100%;
        }

        .sidebar-header {
            padding: 20px;
            border-bottom: 1px solid var(--border);
        }

        .sidebar-header h1 {
            font-size: 1.25rem;
            color: var(--accent);
            font-weight: 700;
            letter-spacing: -0.025em;
        }

        .sidebar-header p {
            font-size: 0.85rem;
            color: var(--text-muted);
            margin-top: 4px;
        }

        .search-box {
            padding: 12px 20px;
            border-bottom: 1px solid var(--border);
        }

        .search-box input {
            width: 100%;
            padding: 8px 12px;
            background-color: var(--bg-dark);
            border: 1px solid var(--border);
            border-radius: 6px;
            color: var(--text-main);
            font-size: 0.9rem;
        }

        .search-box input:focus {
            outline: none;
            border-color: var(--accent);
        }

        .contig-list {
            flex: 1;
            overflow-y: auto;
            padding: 10px;
        }

        .contig-item {
            padding: 12px 14px;
            border-radius: 8px;
            cursor: pointer;
            transition: background 0.2s, border-color 0.2s;
            margin-bottom: 8px;
            border: 1px solid transparent;
            background-color: var(--bg-card);
        }

        .contig-item:hover {
            background-color: #3b82f622;
            border-color: var(--accent);
        }

        .contig-item.active {
            background-color: #3b82f633;
            border-color: var(--accent);
        }

        .contig-title {
            font-weight: 600;
            font-size: 0.95rem;
            color: var(--text-main);
            display: flex;
            justify-content: space-between;
        }

        .contig-meta {
            font-size: 0.8rem;
            color: var(--text-muted);
            margin-top: 4px;
            display: flex;
            gap: 12px;
        }

        main {
            flex: 1;
            display: flex;
            flex-direction: column;
            height: 100%;
            position: relative;
        }

        #viewport-container {
            flex: 1;
            width: 100%;
            height: 100%;
            position: relative;
            background: radial-gradient(circle at center, #1e293b 0%, #0f172a 100%);
        }

        .viewer-overlay {
            position: absolute;
            top: 20px;
            left: 20px;
            z-index: 10;
            background: rgba(30, 41, 59, 0.85);
            backdrop-filter: blur(8px);
            border: 1px solid var(--border);
            padding: 16px;
            border-radius: 10px;
            max-width: 450px;
            box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.3);
        }

        .viewer-overlay h2 {
            font-size: 1.1rem;
            margin-bottom: 8px;
            color: var(--accent);
        }

        .stats-grid {
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 8px;
            margin-top: 10px;
            font-size: 0.85rem;
        }

        .stat-card {
            background: rgba(15, 23, 42, 0.6);
            padding: 8px;
            border-radius: 6px;
            border: 1px solid var(--border);
        }

        .stat-label {
            color: var(--text-muted);
            font-size: 0.75rem;
        }

        .stat-value {
            font-weight: 600;
            color: var(--text-main);
            margin-top: 2px;
        }

        .controls-bar {
            position: absolute;
            bottom: 20px;
            left: 50%;
            transform: translateX(-50%);
            z-index: 10;
            background: rgba(30, 41, 59, 0.9);
            backdrop-filter: blur(8px);
            border: 1px solid var(--border);
            padding: 10px 20px;
            border-radius: 30px;
            display: flex;
            gap: 12px;
            box-shadow: 0 10px 25px -5px rgba(0, 0, 0, 0.3);
        }

        button {
            background-color: var(--bg-card);
            color: var(--text-main);
            border: 1px solid var(--border);
            padding: 6px 14px;
            border-radius: 6px;
            cursor: pointer;
            font-size: 0.85rem;
            transition: all 0.2s;
        }

        button:hover {
            background-color: var(--accent);
            color: var(--bg-dark);
            border-color: var(--accent);
        }

        .no-structure {
            position: absolute;
            inset: 0;
            display: flex;
            align-items: center;
            justify-content: center;
            background: var(--bg-dark);
            color: var(--text-muted);
            font-size: 1.1rem;
            z-index: 5;
            display: none;
        }

        .loading-spinner {
            position: absolute;
            inset: 0;
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            background: rgba(15, 23, 42, 0.75);
            backdrop-filter: blur(4px);
            z-index: 20;
            color: var(--accent);
            font-size: 0.9rem;
            gap: 12px;
        }

        .loading-spinner.hidden {
            display: none;
        }

        .spinner {
            width: 40px;
            height: 40px;
            border: 3px solid rgba(56, 189, 248, 0.2);
            border-top-color: var(--accent);
            border-radius: 50%;
            animation: spin 0.8s linear infinite;
        }

        @keyframes spin {
            to { transform: rotate(360deg); }
        }
    </style>
</head>
<body>

    <aside>
        <div class="sidebar-header">
            <h1>Foldinator Atlas</h1>
            <p>Arctic Deep-Sea Viral Metagenome Pipeline</p>
        </div>
        <div class="search-box">
            <input type="text" id="searchInput" placeholder="Filter targets (e.g., k141)..." oninput="filterContigs()">
        </div>
        <div class="contig-list" id="contigList"></div>
    </aside>

    <main>
        <div id="viewport-container">
            <div id="loadingSpinner" class="loading-spinner hidden">
                <div class="spinner"></div>
                <span>Decoding PDB & Rendering Structure...</span>
            </div>
        </div>

        <div id="noStructureMsg" class="no-structure">No structural model (PDB) available for this target.</div>

        <div class="viewer-overlay" id="infoPanel">
            <h2 id="activeTitle">Select a target</h2>
            <p id="activeDesc" style="font-size: 0.85rem; color: var(--text-muted);">Choose a protein from the sidebar to inspect structural predictions and homology matches.</p>
            <div class="stats-grid" id="statsGrid" style="display: none;">
                <div class="stat-card">
                    <div class="stat-label">Length</div>
                    <div class="stat-value" id="statLen">-</div>
                </div>
                <div class="stat-card">
                    <div class="stat-label">Mean pLDDT</div>
                    <div class="stat-value" id="statPlddt">-</div>
                </div>
                <div class="stat-card">
                    <div class="stat-label">Homology Status</div>
                    <div class="stat-value" id="statCheckv">-</div>
                </div>
                <div class="stat-card">
                    <div class="stat-label">Best E-value</div>
                    <div class="stat-value" id="statDark">-</div>
                </div>
            </div>
        </div>

        <div class="controls-bar" id="controlsBar" style="display: none;">
            <button onclick="setRepresentation('cartoon')">Cartoon</button>
            <button onclick="setRepresentation('surface')">Surface</button>
            <button onclick="setRepresentation('ball+stick')">Ball & Stick</button>
            <button onclick="resetView()">Reset View</button>
        </div>
    </main>

    <script>
        const ATLAS_DATA = __ATLAS_DATA_JSON__;

        let nglStage = null;
        let currentComponent = null;

        document.addEventListener("DOMContentLoaded", () => {
            nglStage = new NGL.Stage("viewport-container", { backgroundColor: "#0f172a" });

            window.addEventListener("resize", () => {
                nglStage.handleResize();
            });

            renderContigList(ATLAS_DATA);
            if (ATLAS_DATA.length > 0) {
                selectContig(0);
            }
        });

        function renderContigList(data) {
            const listEl = document.getElementById("contigList");
            listEl.innerHTML = "";
            data.forEach((item, index) => {
                const div = document.createElement("div");
                div.className = `contig-item ${index === 0 ? 'active' : ''}`;
                div.id = `contig-item-${index}`;
                div.onclick = () => selectContig(index);
                div.innerHTML = `
                    <div class="contig-title">
                        <span>${item.id}</span>
                        <span style="font-size: 0.75rem; color: var(--accent);">${item.length_aa} aa</span>
                    </div>
                    <div class="contig-meta">
                        <span>pLDDT: ${item.mean_plddt || 'N/A'}</span>
                        <span>${item.checkv_quality || 'Analyzed'}</span>
                    </div>
                `;
                listEl.appendChild(div);
            });
        }

        function filterContigs() {
            const query = document.getElementById("searchInput").value.toLowerCase();
            const filtered = ATLAS_DATA.filter(item => item.id.toLowerCase().includes(query));
            const listEl = document.getElementById("contigList");
            listEl.innerHTML = "";

            filtered.forEach((item) => {
                const originalIndex = ATLAS_DATA.findIndex(d => d.id === item.id);
                const div = document.createElement("div");
                div.className = `contig-item`;
                div.id = `contig-item-${originalIndex}`;
                div.onclick = () => selectContig(originalIndex);
                div.innerHTML = `
                    <div class="contig-title">
                        <span>${item.id}</span>
                        <span style="font-size: 0.75rem; color: var(--accent);">${item.length_aa} aa</span>
                    </div>
                    <div class="contig-meta">
                        <span>pLDDT: ${item.mean_plddt || 'N/A'}</span>
                        <span>${item.checkv_quality || 'Analyzed'}</span>
                    </div>
                `;
                listEl.appendChild(div);
            });
        }

        function selectContig(index) {
            document.querySelectorAll('.contig-item').forEach(el => el.classList.remove('active'));
            const activeEl = document.getElementById(`contig-item-${index}`);
            if (activeEl) activeEl.classList.add('active');

            const item = ATLAS_DATA[index];
            document.getElementById("activeTitle").innerText = item.id;
            document.getElementById("activeDesc").innerText = item.description || "Deep-sea viral protein target analyzed via Foldinator pipeline.";

            document.getElementById("statsGrid").style.display = "grid";
            document.getElementById("controlsBar").style.display = "flex";
            document.getElementById("statLen").innerText = item.length_aa + " aa";
            document.getElementById("statPlddt").innerText = item.mean_plddt || "N/A";
            document.getElementById("statCheckv").innerText = item.checkv_quality || "Unknown";
            document.getElementById("statDark").innerText = item.dark_matter_loci || "N/A";

            loadStructure(item);
        }

        function loadStructure(item) {
            const spinner = document.getElementById("loadingSpinner");
            const noStructMsg = document.getElementById("noStructureMsg");

            nglStage.removeAllComponents();

            if (!item.has_pdb || !item.pdb_b64) {
                spinner.classList.add("hidden");
                noStructMsg.style.display = "flex";
                return;
            }

            noStructMsg.style.display = "none";
            spinner.classList.remove("hidden");

            setTimeout(() => {
                try {
                    const decodedPdb = atob(item.pdb_b64);
                    const blob = new Blob([decodedPdb], { type: 'text/plain' });

                    nglStage.loadFile(blob, { ext: "pdb" }).then(function(component) {
                        currentComponent = component;
                        component.addRepresentation("cartoon", { colorScheme: "plddt" });
                        component.autoView();
                        spinner.classList.add("hidden");
                    }).catch(err => {
                        console.error("Error parsing PDB structure in NGL viewer:", err);
                        spinner.classList.add("hidden");
                        noStructMsg.style.display = "flex";
                        noStructMsg.innerText = "Error parsing structure coordinate file.";
                    });
                } catch (e) {
                    console.error("Base64 decoding failed:", e);
                    spinner.classList.add("hidden");
                    noStructMsg.style.display = "flex";
                    noStructMsg.innerText = "Failed to decode PDB payload.";
                }
            }, 30);
        }

        function setRepresentation(type) {
            if (!currentComponent) return;
            currentComponent.removeAllRepresentations();
            currentComponent.addRepresentation(type, { colorScheme: "plddt" });
        }

        function resetView() {
            if (nglStage) {
                nglStage.autoView();
            }
        }
    </script>
</body>
</html>
"""

def generate_dataset_from_pipeline(output_dir):
    pdbs_dir = os.path.join(output_dir, "pdbs")
    tsv_dir = os.path.join(output_dir, "foldseek_results")

    data = []
    if not os.path.exists(pdbs_dir):
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
            "description": f"Deep-sea Arctic viral protein target {target_id} modeled via ESMFold & Foldseek.",
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
    data = generate_dataset_from_pipeline(target_dir)

    if not data:
        print("[WARNING] No PDBs found to compile into the atlas.")
        data = []

    data_json_str = json.dumps(data)
    final_html = HTML_TEMPLATE.replace("__ATLAS_DATA_JSON__", data_json_str)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(final_html)

    print(f"Successfully generated structural atlas: {os.path.abspath(output_path)}")

if __name__ == "__main__":
    main()
