#!/usr/bin/env python3
"""
foldinator_atlas.py - Phase 5 Compiler for Foldinator Atlas
Ingests processed pipeline outputs (Foldseek TSVs, ESMFold PDBs, and Quarantine folders)
to compile a lightning-fast, interactive HTML atlas with embedded 3Dmol.js WebGL models.
"""

import os
import sys
import glob
import json
import argparse
import base64
import pandas as pd
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Compile Foldinator Atlas Dashboard")
    parser.add_argument("--output_dir", default="foldinator_outputs", help="Root output directory of the pipeline")
    parser.add_argument("--out_html", default="foldinator_atlas.html", help="Path to write final HTML dashboard")
    return parser.parse_args()

def normalize_id(identifier):
    """Clean pipeline identifiers by removing common extension noise and tags."""
    if not identifier:
        return ""
    clean = str(identifier).strip()
    for suffix in [".pdb", ".tsv", "_unfolded", "_folded", "_dark_matter"]:
        if clean.endswith(suffix):
            clean = clean[:-len(suffix)]
    return clean

def extract_gene_tag(raw_str):
    """Extract standard gene tag (e.g. full_1_gene_6054) from complex headers or filenames."""
    cleaned = normalize_id(raw_str)
    if "||" in cleaned:
        cleaned = cleaned.split("||")[1]
    if "_" in cleaned and "full_" in cleaned:
        parts = cleaned.split("_")
        for i in range(len(parts) - 3):
            if parts[i] == "full" and parts[i+2] == "gene":
                return f"{parts[i]}_{parts[i+1]}_{parts[i+2]}_{parts[i+3]}"
    return cleaned

def scan_pipeline_data(output_dir):
    print(f"[*] Scanning processed pipeline outputs in: {output_dir}")

    proteins_data = {}
    similarity_edges = []
    processed_ids = set()
    pdb_map = {}
    gene_to_pdb = {}

    # 1. Map PDB files from novel_dark_matter
    pdb_dir = os.path.join(output_dir, "pdbs", "novel_dark_matter")
    if not os.path.exists(pdb_dir):
        pdb_dir = os.path.join(output_dir, "pdbs")

    pdb_files = glob.glob(os.path.join(pdb_dir, "**", "*.pdb"), recursive=True)

    for pdb_path in pdb_files:
        if "quarantine_corrupt" in pdb_path:
            continue

        raw_name = os.path.basename(pdb_path).replace(".pdb", "")
        clean_id = normalize_id(raw_name)
        gtag = extract_gene_tag(raw_name)

        processed_ids.add(clean_id)
        pdb_map[clean_id] = pdb_path
        pdb_map[raw_name] = pdb_path
        gene_to_pdb[gtag] = pdb_path

    print(f"[*] Discovered {len(pdb_map)} valid PDB files.")

    # 2. Map Foldseek results
    foldseek_dir = os.path.join(output_dir, "foldseek_results")
    tsv_files = glob.glob(os.path.join(foldseek_dir, "**", "*.tsv"), recursive=True)
    if not tsv_files:
        tsv_files = glob.glob(os.path.join(output_dir, "*.tsv"))

    foldseek_results = {}

    for tsv in tsv_files:
        try:
            df = pd.read_csv(tsv, sep="\t", header=None)
            if df.shape[1] >= 7:
                for _, row in df.iterrows():
                    q_raw = str(row[0])
                    q = normalize_id(q_raw)
                    q_gtag = extract_gene_tag(q_raw)
                    target_id = normalize_id(str(row[1]))

                    foldseek_results[q] = {
                        "hit": str(row[5]),
                        "pident": float(row[2]),
                        "evalue": float(row[3]),
                        "bits": float(row[4]),
                        "tm_score": float(row[6])
                    }
                    foldseek_results[q_gtag] = foldseek_results[q]

                    if float(row[3]) <= 1e-3:
                        similarity_edges.append({
                            "source": q,
                            "target": target_id,
                            "value": float(row[4])
                        })
        except Exception as e:
            print(f"[!] Warning: Error parsing TSV {tsv}: {e}")

    # 3. Load sequences from FAA into a quick lookup dictionary first
    seq_lookup = {}
    faa_files = [
        os.path.join(output_dir, "upstream_prep", "prodigal_raw_proteins.faa"),
        os.path.join(output_dir, "upstream_prep", "prodigal_rational_proteins.faa")
    ]
    faa_files = [f for f in faa_files if os.path.exists(f)]
    if not faa_files:
        faa_files = glob.glob(os.path.join(output_dir, "**", "*.faa"), recursive=True)

    for faa in faa_files:
        try:
            for record in SeqIO.parse(faa, "fasta"):
                seq_lookup[record.id] = record
                seq_lookup[normalize_id(record.id)] = record
                seq_lookup[extract_gene_tag(record.id)] = record
        except Exception as e:
            print(f"[!] Warning parsing FAA {faa}: {e}")

    # 4. Build protein records STRICTLY for genes that have a valid PDB file on disk
    for pid, pdb_path in pdb_map.items():
        gtag = extract_gene_tag(pid)

        rec = seq_lookup.get(pid) or seq_lookup.get(normalize_id(pid)) or seq_lookup.get(gtag)

        if rec:
            rec_id = rec.id
            seq_str = str(rec.seq)
            parts = rec.id.split(" # ")
            header_raw = parts[0]
            if "||" in header_raw:
                contig_id, gene_tag = header_raw.split("||", 1)
            else:
                contig_id, gene_tag = header_raw, rec.id

            start_pos = int(parts[1]) if len(parts) > 1 else 0
            end_pos = int(parts[2]) if len(parts) > 2 else 0
            strand = parts[3] if len(parts) > 3 else "1"
        else:
            rec_id = pid
            seq_str = "N/A"
            contig_id = pid.split("||")[0] if "||" in pid else "Unassigned_Contig"
            gene_tag = pid
            start_pos, end_pos, strand = 0, 0, "1"

        rec_id_clean = normalize_id(rec_id)
        fs = foldseek_results.get(rec_id_clean, foldseek_results.get(gtag, {
            "hit": "Structural Dark Matter (Orphan Fold)",
            "pident": 0.0,
            "evalue": 1.0,
            "bits": 0.0,
            "tm_score": 0.0
        }))

        b_factors = []
        pdb_base64 = ""
        try:
            with open(pdb_path, "r") as pf:
                pdb_content = pf.read()
                pdb_base64 = base64.b64encode(pdb_content.encode('utf-8')).decode('utf-8')
                for line in pdb_content.splitlines():
                    if line.startswith("ATOM"):
                        try:
                            b_factors.append(float(line[60:66].strip()))
                        except ValueError:
                            pass
            plddt = round(sum(b_factors) / len(b_factors), 1) if b_factors else 0.0
        except Exception:
            plddt = 0.0

        proteins_data[rec_id] = {
            "id": rec_id,
            "contig_id": contig_id,
            "gene_tag": gene_tag,
            "start": start_pos,
            "end": end_pos,
            "strand": strand,
            "seq": seq_str,
            "length": len(seq_str) if seq_str != "N/A" else 0,
            "foldseek_hit": fs["hit"],
            "evalue": fs["evalue"],
            "pident": fs["pident"],
            "bits": fs["bits"],
            "tm_score": fs["tm_score"],
            "plddt": plddt,
            "has_pdb": True,
            "pdb_b64": pdb_base64,
            "pdb_path": pdb_path
        }

    contigs_map = {}
    for pid, p in proteins_data.items():
        cid = p["contig_id"]
        if cid not in contigs_map:
            contigs_map[cid] = {
                "contig_id": cid,
                "genes": [],
                "taxonomy": "Processed Metagenomic Bin"
            }
        if pid not in contigs_map[cid]["genes"]:
            contigs_map[cid]["genes"].append(pid)

    print(f"[*] Successfully compiled {len(proteins_data)} folded proteins across {len(contigs_map)} contigs.")
    return list(contigs_map.values()), list(proteins_data.values()), similarity_edges

def generate_html_dashboard(contigs, proteins, edges, output_html):
    print(f"[*] Compiling cumulative interactive dashboard to {output_html}...")

    payload = {
        "contigs": contigs,
        "proteins": proteins,
        "edges": edges[:500]
    }

    html_template = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Foldinator Atlas — Structural Genomics Explorer</title>
    <script src="https://cdn.jsdelivr.net/npm/@tailwindcss/browser@4"></script>
    <script type="text/javascript" src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js"></script>
    <style>
        body {{ background-color: #0f172a; color: #f8fafc; font-family: ui-sans-serif, system-ui, sans-serif; }}
        .custom-scrollbar::-webkit-scrollbar {{ width: 6px; }}
        .custom-scrollbar::-webkit-scrollbar-track {{ background: #1e293b; }}
        .custom-scrollbar::-webkit-scrollbar-thumb {{ background: #475569; border-radius: 3px; }}
    </style>
</head>
<body class="h-screen flex flex-col overflow-hidden">

    <header class="bg-slate-900 border-b border-slate-800 px-6 py-4 flex justify-between items-center shrink-0">
        <div class="flex items-center space-x-3">
            <div class="w-3 h-3 bg-cyan-500 rounded-full animate-pulse"></div>
            <h1 class="text-xl font-bold tracking-wider text-slate-100">FOLDINATOR ATLAS</h1>
            <span class="text-xs bg-slate-800 text-cyan-400 px-2.5 py-1 rounded-full border border-slate-700">Cumulative Landscape</span>
        </div>
        <div class="flex space-x-6 text-sm">
            <div class="bg-slate-800/80 px-4 py-1.5 rounded-lg border border-slate-700">
                <span class="text-slate-400 text-xs block">Active Contigs</span>
                <span id="kpi-contigs" class="font-bold text-slate-200">0</span>
            </div>
            <div class="bg-slate-800/80 px-4 py-1.5 rounded-lg border border-slate-700">
                <span class="text-slate-400 text-xs block">Processed Proteins</span>
                <span id="kpi-proteins" class="font-bold text-slate-200">0</span>
            </div>
            <div class="bg-slate-800/80 px-4 py-1.5 rounded-lg border border-slate-700">
                <span class="text-slate-400 text-xs block">Mean pLDDT Confidence</span>
                <span id="kpi-plddt" class="font-bold text-cyan-400">0.0</span>
            </div>
        </div>
    </header>

    <div class="flex-1 flex overflow-hidden">
        <div class="w-2/5 border-r border-slate-800 flex flex-col bg-slate-900/50">
            <div class="flex border-b border-slate-800 bg-slate-900">
                <button onclick="switchTab('tree')" id="btn-tree" class="flex-1 py-3 text-sm font-semibold text-cyan-400 border-b-2 border-cyan-500 cursor-pointer">Processed Contig Tree</button>
                <button onclick="switchTab('network')" id="btn-network" class="flex-1 py-3 text-sm font-semibold text-slate-400 hover:text-slate-200 cursor-pointer">Similarity Network</button>
            </div>

            <div id="view-tree" class="flex-1 overflow-y-auto custom-scrollbar p-4 space-y-3">
                <div class="relative">
                    <input type="text" id="search-box" placeholder="Search processed contigs or ORFs..." oninput="filterTree()" class="w-full bg-slate-800 border border-slate-700 rounded-lg px-4 py-2 text-sm text-slate-200 focus:outline-none focus:border-cyan-500">
                </div>
                <div id="contig-list-container" class="space-y-2"></div>
            </div>

            <div id="view-network" class="flex-1 hidden relative bg-slate-950">
                <div id="network-canvas" class="absolute inset-0"></div>
                <div class="absolute bottom-3 left-3 bg-slate-900/80 backdrop-blur px-3 py-1.5 rounded border border-slate-800 text-xs text-slate-400">
                    <span class="inline-block w-2 h-2 rounded-full bg-cyan-400 mr-1.5"></span>Known Homolog / <span class="inline-block w-2 h-2 rounded-full bg-amber-400 mr-1.5 ml-2"></span>Structural Dark Matter
                </div>
            </div>
        </div>

        <div class="w-3/5 flex flex-col bg-slate-950 overflow-y-auto custom-scrollbar p-6">
            <div id="inspector-welcome" class="h-full flex flex-col items-center justify-center text-center text-slate-500 space-y-3">
                <svg class="w-16 h-16 stroke-slate-700" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                    <path stroke-linecap="round" stroke-linejoin="round" stroke-width="1.5" d="M19 11H5m14 0a2 2 0 012 2v6a2 2 0 01-2 2H5a2 2 0 01-2-2v-6a2 2 0 012-2m14 0V9a2 2 0 00-2-2M5 11V9a2 2 0 012-2m0 0V5a2 2 0 012-2h6a2 2 0 012 2v2M7 7h10"/>
                </svg>
                <p class="text-base font-medium">Select a processed protein from the tree or network to inspect its structure.</p>
            </div>

            <div id="inspector-content" class="hidden space-y-6">
                <div class="flex justify-between items-start border-b border-slate-800 pb-4">
                    <div>
                        <span id="ins-badge" class="text-xs bg-cyan-950 text-cyan-400 border border-cyan-800 px-2.5 py-1 rounded-full uppercase tracking-wider">ORF Inspector</span>
                        <h2 id="ins-title" class="text-2xl font-bold text-slate-100 mt-2">--</h2>
                        <p id="ins-subtitle" class="text-sm text-slate-400 mt-0.5">--</p>
                    </div>
                    <div id="ins-plddt-badge" class="text-right">
                        <span class="text-xs text-slate-400 block">pLDDT Confidence</span>
                        <span id="ins-plddt-val" class="text-lg font-bold text-cyan-400">--</span>
                    </div>
                </div>

                <div class="grid grid-cols-1 lg:grid-cols-2 gap-6">
                    <div class="bg-slate-900 border border-slate-800 rounded-xl p-4 flex flex-col">
                        <div class="flex justify-between items-center mb-3">
                            <span class="text-xs font-semibold uppercase tracking-wider text-slate-400">ESMFold 3D Structure</span>
                            <span class="text-xs text-slate-500">WebGL View (pLDDT Spectrum)</span>
                        </div>
                        <div id="viewer-container" class="relative w-full h-72 bg-slate-950 rounded-lg border border-slate-800 overflow-hidden flex items-center justify-center">
                            <span id="no-pdb-msg" class="text-xs text-slate-500 hidden">No PDB structure loaded</span>
                        </div>
                    </div>

                    <div class="bg-slate-900 border border-slate-800 rounded-xl p-4 flex flex-col justify-between">
                        <div>
                            <span class="text-xs font-semibold uppercase tracking-wider text-slate-400 block mb-3">Foldseek Structural Match</span>
                            <div class="space-y-3 text-sm">
                                <div>
                                    <span class="text-xs text-slate-500 block">Top Database Hit</span>
                                    <span id="ins-hit" class="font-medium text-slate-200">--</span>
                                </div>
                                <div class="grid grid-cols-4 gap-2 pt-2 border-t border-slate-800">
                                    <div>
                                        <span class="text-xs text-slate-500 block">TM-Score</span>
                                        <span id="ins-tm" class="font-mono text-cyan-400">--</span>
                                    </div>
                                    <div>
                                        <span class="text-xs text-slate-500 block">E-Value</span>
                                        <span id="ins-evalue" class="font-mono text-cyan-400">--</span>
                                    </div>
                                    <div>
                                        <span class="text-xs text-slate-500 block">Identity %</span>
                                        <span id="ins-pident" class="font-mono text-cyan-400">--</span>
                                    </div>
                                    <div>
                                        <span class="text-xs text-slate-500 block">Bit Score</span>
                                        <span id="ins-bits" class="font-mono text-cyan-400">--</span>
                                    </div>
                                </div>
                            </div>
                        </div>
                        <div class="pt-4 border-t border-slate-800 mt-4">
                            <span class="text-xs text-slate-500 block">Genomic Coordinates</span>
                            <span id="ins-coords" class="font-mono text-xs text-slate-300">--</span>
                        </div>
                    </div>
                </div>

                <div class="border border-slate-800 rounded-xl bg-slate-900 overflow-hidden">
                    <button onclick="toggleAccordion()" class="w-full px-4 py-3 flex justify-between items-center text-left text-sm font-medium text-slate-300 hover:bg-slate-800/50 cursor-pointer">
                        <span>[+] Raw Telemetry & Sequence Data</span>
                        <svg id="accordion-icon" class="w-4 h-4 transform transition-transform" fill="none" viewBox="0 0 24 24" stroke="currentColor">
                            <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M19 9l-7 7-7-7"/>
                        </svg>
                    </button>
                    <div id="accordion-body" class="hidden p-4 border-t border-slate-800 bg-slate-950 font-mono text-xs text-slate-400 space-y-2">
                        <div>
                            <span class="text-slate-500 block">Full Identifier:</span>
                            <span id="raw-id" class="text-slate-300">--</span>
                        </div>
                        <div>
                            <span class="text-slate-500 block">Amino Acid Sequence:</span>
                            <div id="raw-seq" class="bg-slate-900 p-3 rounded border border-slate-800 overflow-x-auto max-h-32 text-slate-300">--</div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </div>

    <script>
        const ATLAS_DATA = {json.dumps(payload)};
        let viewer3D = null;

        document.addEventListener("DOMContentLoaded", () => {{
            initDashboard();
        }});

        function initDashboard() {{
            document.getElementById("kpi-contigs").innerText = ATLAS_DATA.contigs.length;
            document.getElementById("kpi-proteins").innerText = ATLAS_DATA.proteins.length;

            const plddts = ATLAS_DATA.proteins.filter(p => p.plddt > 0).map(p => p.plddt);
            const meanPlddt = plddts.length ? (plddts.reduce((a, b) => a + b, 0) / plddts.length).toFixed(1) : "0.0";
            document.getElementById("kpi-plddt").innerText = meanPlddt;

            renderTree(ATLAS_Data.contigs);
            initNetwork();

            let element = document.getElementById("viewer-container");
            viewer3D = $3Dmol.createViewer(element, {{ backgroundColor: "#020617" }});

            // Explicitly inject scene lighting so WebGL shaders have a light source to render
            viewer3D.addLight();
        }}

        function renderTree(contigs) {{
            const container = document.getElementById("contig-list-container");
            container.innerHTML = "";

            contigs.forEach(c => {{
                const div = document.createElement("div");
                div.className = "bg-slate-900 border border-slate-800 rounded-lg p-3 hover:border-slate-700 transition";

                let genesHtml = c.genes.map(gid => {{
                    const prot = ATLAS_DATA.proteins.find(p => p.id === gid);
                    const hasHit = prot && prot.evalue <= 1e-3;
                    const badgeColor = hasHit ? "bg-cyan-950 text-cyan-400 border-cyan-800" : "bg-amber-950 text-amber-400 border-amber-800";
                    const label = hasHit ? "Known Homolog" : "Structural Dark Matter";

                    return `<div onclick="inspectProtein('${{gid}}')" class="flex justify-between items-center py-1.5 px-2 rounded hover:bg-slate-800 cursor-pointer text-xs">
                        <span class="font-mono text-slate-300 truncate w-3/5">${{gid}}</span>
                        <span class="px-2 py-0.5 rounded border text-[10px] ${{badgeColor}}">${{label}}</span>
                    </div>`;
                }}).join("");

                div.innerHTML = `
                    <div class="flex justify-between items-center cursor-pointer" onclick="toggleContig('${{c.contig_id}}')">
                        <div>
                            <span class="font-mono text-xs font-bold text-slate-200">${{c.contig_id}}</span>
                            <span class="text-[10px] text-slate-400 block">${{c.genes.length}} Processed ORFs</span>
                        </div>
                        <span class="text-xs text-cyan-400">Expand ▾</span>
                    </div>
                    <div id="contig-genes-${{c.contig_id}}" class="hidden mt-2 pt-2 border-t border-slate-800 space-y-1 max-h-48 overflow-y-auto custom-scrollbar">
                        ${{genesHtml}}
                    </div>
                `;
                container.appendChild(div);
            }});
        }}

        function toggleContig(cid) {{
            const el = document.getElementById(`contig-genes-${{cid}}`);
            if (el) el.classList.toggle("hidden");
        }}

        function switchTab(tab) {{
            document.getElementById("view-tree").classList.toggle("hidden", tab !== "tree");
            document.getElementById("view-network").classList.toggle("hidden", tab !== "network");
            document.getElementById("btn-tree").className = tab === "tree" ? "flex-1 py-3 text-sm font-semibold text-cyan-400 border-b-2 border-cyan-500 cursor-pointer" : "flex-1 py-3 text-sm font-semibold text-slate-400 hover:text-slate-200 cursor-pointer";
            document.getElementById("btn-network").className = tab === "network" ? "flex-1 py-3 text-sm font-semibold text-cyan-400 border-b-2 border-cyan-500 cursor-pointer" : "flex-1 py-3 text-sm font-semibold text-slate-400 hover:text-slate-200 cursor-pointer";

            if (tab === "network" && window.networkInstance) {{
                window.networkInstance.redraw();
            }}
        }}

        function initNetwork() {{
            const nodes = [];
            const edges = [];

            ATLAS_DATA.proteins.forEach(p => {{
                const hasHit = p.evalue <= 1e-3;
                nodes.push({{
                    id: p.id,
                    label: p.id.split("_").pop(),
                    color: hasHit ? "#06b6d4" : "#f59e0b",
                    size: hasHit ? 12 : 15
                }});
            }});

            ATLAS_DATA.edges.forEach(e => {{
                edges.push({{ from: e.source, to: e.target, color: "#334155" }});
            }});

            const container = document.getElementById("network-canvas");
            const data = {{ nodes: new vis.DataSet(nodes), edges: new vis.DataSet(edges) }};
            const options = {{
                nodes: {{ shape: "dot", font: {{ color: "#cbd5e1", size: 10 }}}},
                physics: {{ barnesHut: {{ gravitationalConstant: -2000, centralGravity: 0.3, springLength: 95 }}}}
            }};
            window.networkInstance = new vis.Network(container, data, options);

            window.networkInstance.on("selectNode", params => {{
                if (params.nodes.length > 0) inspectProtein(params.nodes[0]);
            }});
        }}

        function inspectProtein(pid) {{
            const prot = ATLAS_DATA.proteins.find(p => p.id === pid);
            if (!prot) return;

            document.getElementById("inspector-welcome").classList.add("hidden");
            document.getElementById("inspector-content").classList.remove("hidden");

            document.getElementById("ins-title").innerText = prot.id;
            document.getElementById("ins-subtitle").innerText = "Contig: " + prot.contig_id + " • Length: " + prot.length + " aa";
            document.getElementById("ins-plddt-val").innerText = prot.plddt > 0 ? prot.plddt : "N/A";

            document.getElementById("ins-hit").innerText = prot.foldseek_hit;
            document.getElementById("ins-tm").innerText = prot.tm_score > 0 ? prot.tm_score.toFixed(3) : "N/A";
            document.getElementById("ins-evalue").innerText = prot.evalue;
            document.getElementById("ins-pident").innerText = prot.pident + "%";
            document.getElementById("ins-bits").innerText = prot.bits;
            document.getElementById("ins-coords").innerText = "Start: " + prot.start + " | End: " + prot.end + " | Strand: " + prot.strand;

            document.getElementById("raw-id").innerText = prot.id;
            document.getElementById("raw-seq").innerText = prot.seq && prot.seq.length > 0 ? prot.seq : "No sequence data available";

            if (viewer3D) {{
                viewer3D.clear();
                const msgEl = document.getElementById("no-pdb-msg");
                if (prot.has_pdb && prot.pdb_b64) {{
                    msgEl.classList.remove("hidden");

                    const rawPdb = atob(prot.pdb_b64);
                    const cleanPdb = rawPdb
                        .split('\n')
                        .filter(line => !line.startsWith('PARENT'))
                        .join('\n');

                    viewer3D.addModel(cleanPdb, "pdb");

                    // Force 3Dmol to calculate backbone geometry and secondary structure for ESMFold outputs
                    viewer3D.calculateSecondaryStructure();

                    viewer3D.setStyle({{}}, {{
                        stick: {{
                            colorscheme: 'element'
                        }}
                    }});

                    viewer3D.zoomTo();
                    viewer3D.resize();
                    viewer3D.render();
                }} else {{
                    msgEl.classList.remove("hidden");
                }}
            }}
        }}

        function toggleAccordion() {{
            const body = document.getElementById("accordion-body");
            const icon = document.getElementById("accordion-icon");
            body.classList.toggle("hidden");
            icon.classList.toggle("rotate-180");
        }}

        function filterTree() {{
            const query = document.getElementById("search-box").value.toLowerCase();
            const filteredContigs = ATLAS_DATA.contigs.filter(c =>
                c.contig_id.toLowerCase().includes(query) ||
                c.genes.some(g => g.toLowerCase().includes(query))
            );
            renderTree(filteredContigs);
        }}
    </script>
</body>
</html>
    """

    with open(output_html, "w") as f:
        f.write(html_template)
    print(f"[+] Successfully compiled cumulative atlas to {output_html}!")

def main():
    args = parse_args()
    contigs, proteins, edges = scan_pipeline_data(args.output_dir)
    generate_html_dashboard(contigs, proteins, edges, args.out_html)

if __name__ == "__main__":
    main()
