import os
import sqlite3
import shutil

def init_db(db_path="foldinator.db"):
    """Initializes the SQLite pipeline ledger."""
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS pipeline_runs (
            target_name TEXT PRIMARY KEY,
            search_status TEXT,
            evalue_threshold REAL,
            high_confidence INTEGER,
            moderate_fold INTEGER,
            low_confidence INTEGER,
            discarded_count INTEGER,
            file_status TEXT,
            last_updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    """)
    conn.commit()
    conn.close()

def audit_and_clean_tsv(tsv_path, pdb_dir, db_path="foldinator.db", evalue_threshold=1e-3):
    """
    Acts as the QC cop: validates 7-column Foldseek TSV, bins valid hits by TM-score,
    purges noise files, manages PDB quarantine, and logs telemetry directly to SQLite.
    """
    target_name = os.path.basename(tsv_path).replace(".tsv", "")
    init_db(db_path)

    bins = {
        "high_confidence": 0, # TM-score >= 0.5
        "moderate_fold": 0,   # 0.3 <= TM-score < 0.5
        "low_confidence": 0   # TM-score < 0.3
    }
    discarded_count = 0
    file_status = "valid"
    has_any_hits = False

    # Check for missing or 0-byte structural dark matter
    if not os.path.exists(tsv_path) or os.path.getsize(tsv_path) == 0:
        if os.path.exists(tsv_path):
            os.remove(tsv_path)  # Purge useless empty file
        file_status = "structural_dark_matter"
        _log_to_db(db_path, target_name, "No Alignment", evalue_threshold, bins, discarded_count, file_status)
        manage_pdb_quarantine(target_name, pdb_dir, file_status, has_hits=False)
        return {"status": file_status, "bins": bins, "discarded": discarded_count}

    valid_hits_found = False

    with open(tsv_path, "r", encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            has_any_hits = True

            try:
                # Column layout: [0]query, [1]target, [2]pident, [3]evalue, [4]bits, [5]theader, [6]alntmscore
                evalue = float(parts[3])
                tm_score = float(parts[6]) if len(parts) > 6 else 0.0
            except (ValueError, IndexError):
                discarded_count += 1
                continue

            # Gate 1: E-value threshold check
            if evalue > evalue_threshold:
                discarded_count += 1
                continue

            valid_hits_found = True

            # Gate 2: TM-Score Categorization
            if tm_score >= 0.5:
                bins["high_confidence"] += 1
            elif tm_score >= 0.3:
                bins["moderate_fold"] += 1
            else:
                bins["low_confidence"] += 1

    # If file had bytes but zero hits passed statistical thresholds, purge it
    if not valid_hits_found:
        os.remove(tsv_path)
        file_status = "purged_noise_file"

    _log_to_db(db_path, target_name, "Computed", evalue_threshold, bins, discarded_count, file_status)
    manage_pdb_quarantine(target_name, pdb_dir, file_status, has_hits=has_any_hits)

    return {"status": file_status, "bins": bins, "discarded": discarded_count}

def manage_pdb_quarantine(target_name, pdb_dir, file_status, has_hits):
    pdb_path = os.path.join(pdb_dir, f"{target_name}.pdb")
    if not os.path.exists(pdb_path):
        return

    if file_status in ["structural_dark_matter", "purged_noise_file"]:
        if not has_hits:
            # Valid PDB, but completely unknown fold (True Structural Dark Matter)
            target_sub = os.path.join(pdb_dir, "novel_dark_matter")
        else:
            # Had hits, but all failed statistical thresholds (Noise / Garbage)
            target_sub = os.path.join(pdb_dir, "quarantine_corrupt")

        os.makedirs(target_sub, exist_ok=True)
        shutil.move(pdb_path, os.path.join(target_sub, f"{target_name}.pdb"))

def _log_to_db(db_path, target_name, status, evalue_threshold, bins, discarded_count, file_status):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute("""
        INSERT OR REPLACE INTO pipeline_runs
        (target_name, search_status, evalue_threshold, high_confidence, moderate_fold, low_confidence, discarded_count, file_status, last_updated)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
    """, (
        target_name, status, evalue_threshold,
        bins["high_confidence"], bins["moderate_fold"], bins["low_confidence"],
        discarded_count, file_status
    ))
    conn.commit()
    conn.close()
