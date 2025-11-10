#!/usr/bin/env python3
import re
import time
import gzip
import requests
import pandas as pd
from io import BytesIO
from typing import Any, Dict, List

# --------------- CONFIG ---------------
GSE_LIST = [
    "GSE131685", "GSE107585", "GSE211785", "GSE131882",
    "GSE114530"
]
OUTPUT_XLSX = "kidney_metadata.xlsx"
SLEEP_SEC = 0.5   # be nicer to NCBI
# --------------------------------------


def series_prefix(gse: str) -> str:
    gse = gse.strip()
    base = gse[:-3]
    return f"{base}nnn"


def download_soft(gse: str) -> List[str]:
    prefix = series_prefix(gse)
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{gse}/soft/{gse}_family.soft.gz"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    with gzip.GzipFile(fileobj=BytesIO(r.content)) as gz:
        text = gz.read().decode("utf-8", errors="replace")
    return text.splitlines()


def parse_samples_from_soft(lines: List[str]) -> List[Dict[str, Any]]:
    """
    Pull GSM, title, source_name_ch1, characteristics_ch1, relations.
    Organism is handled via HTML scraping.
    """
    samples = []
    current = None
    for line in lines:
        if line.startswith("^SAMPLE"):
            if current:
                samples.append(current)
            acc = line.split("=", 1)[1].strip()
            current = {
                "gsm": acc,
                "title": None,
                "source_name_ch1": [],
                "characteristics_ch1": [],
                "relations": [],
            }
        elif current is not None and line.startswith("!Sample_"):
            key, _, val = line.partition("=")
            val = val.strip()
            if key == "!Sample_title":
                current["title"] = val
            elif key == "!Sample_source_name_ch1":
                current["source_name_ch1"].append(val)
            elif key == "!Sample_characteristics_ch1":
                current["characteristics_ch1"].append(val)
            elif key == "!Sample_relation":
                current["relations"].append(val)
    if current:
        samples.append(current)
    return samples


# ----------------- HTML SCRAPER -----------------
def fetch_gsm_page(gsm: str) -> str:
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.text


def clean_html_text(s: str) -> str:
    s = re.sub(r"<br\s*/?>", "; ", s, flags=re.IGNORECASE)
    s = re.sub(r"<.*?>", " ", s)
    s = re.sub(r"\s+", " ", s)
    return s.strip()


def scrape_gsm_meta_from_html(html: str) -> Dict[str, str]:
    """
    Parse GSM HTML for:
      - Title
      - Organism
      - Source name
      - Characteristics (possibly multiple)
      - Extraction protocol
      - Instrument model
      - Library strategy
      - Library source
    """
    title = ""
    organism = ""
    source_name = ""
    characteristics: List[str] = []
    extraction_protocol = ""
    instrument_model = ""
    library_strategy = ""
    library_source = ""

    pattern_th = re.compile(
        r"<tr[^>]*>\s*<th[^>]*>([^<]+)</th>\s*<td[^>]*>(.*?)</td>\s*</tr>",
        re.IGNORECASE | re.DOTALL,
    )
    for m in pattern_th.finditer(html):
        label = m.group(1).strip()
        val = clean_html_text(m.group(2))
        ll = label.lower()
        if ll.startswith("title"):
            title = val
        elif ll.startswith("organism"):
            organism = val
        elif ll.startswith("source name"):
            source_name = val
        elif ll.startswith("characteristics"):
            characteristics.append(val)
        elif ll.startswith("extraction protocol"):
            extraction_protocol = val
        elif ll.startswith("instrument model"):
            instrument_model = val
        elif ll.startswith("library strategy"):
            library_strategy = val
        elif ll.startswith("library source"):
            library_source = val

    pattern_td = re.compile(
        r"<tr[^>]*>\s*<td[^>]*>([^<]+)</td>\s*<td[^>]*>(.*?)</td>\s*</tr>",
        re.IGNORECASE | re.DOTALL,
    )
    for m in pattern_td.finditer(html):
        label = m.group(1).strip()
        val = clean_html_text(m.group(2))
        ll = label.lower()
        if ll.startswith("title") and not title:
            title = val
        elif ll.startswith("organism") and not organism:
            organism = val
        elif ll.startswith("source name") and not source_name:
            source_name = val
        elif ll.startswith("characteristics"):
            characteristics.append(val)
        elif ll.startswith("extraction protocol") and not extraction_protocol:
            extraction_protocol = val
        elif ll.startswith("instrument model") and not instrument_model:
            instrument_model = val
        elif ll.startswith("library strategy") and not library_strategy:
            library_strategy = val
        elif ll.startswith("library source") and not library_source:
            library_source = val

    return {
        "title": title,
        "organism": organism,
        "source_name": source_name,
        "characteristics": "; ".join([c for c in characteristics if c]),
        "extraction_protocol": extraction_protocol,
        "instrument_model": instrument_model,
        "library_strategy": library_strategy,
        "library_source": library_source,
    }


# ----------------- SRA -----------------
def fetch_runinfo_csv(acc: str) -> pd.DataFrame:
    url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra/?sp=runinfo&acc={acc}"
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=120)
            r.raise_for_status()
            text = r.text.strip()
            if not text:
                return pd.DataFrame()
            return pd.read_csv(BytesIO(text.encode("utf-8")))
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
            print(f"[retry] runinfo timeout for {acc}, attempt {attempt+1}/3")
            time.sleep(5 * (attempt + 1))
        except Exception as e:
            print(f"[warn] fetch_runinfo_csv failed for {acc}: {e}")
            return pd.DataFrame()
    print(f"[fail] runinfo for {acc} after retries")
    return pd.DataFrame()


def runs_from_relations(relations: List[str]) -> List[str]:
    for rel in relations:
        m = re.search(r"SRA:\s*([A-Z0-9_\.]+)", rel)
        if not m:
            continue
        acc = m.group(1).strip()
        df = fetch_runinfo_csv(acc)
        if not df.empty and "Run" in df.columns:
            runs = df["Run"].astype(str).tolist()
            return list(dict.fromkeys(runs))
    return []


def get_runs_for_gsm(sample: Dict[str, Any]) -> List[str]:
    runs = runs_from_relations(sample.get("relations", []))
    if runs:
        return runs
    df = fetch_runinfo_csv(sample["gsm"])
    if not df.empty and "Run" in df.columns:
        return list(dict.fromkeys(df["Run"].astype(str).tolist()))
    return []


# ----------------- HELPERS FOR CLEANUP -----------------
def strip_trailing_semicolon(s: str) -> str:
    return re.sub(r";\s*$", "", s.strip()) if s else s


def extract_disease_only(char_string: str) -> str:
    if not char_string:
        return ""

    parts = [p.strip() for p in char_string.split(";") if p.strip()]
    if not parts:
        return ""

    # 1) 'disease...'
    for p in parts:
        low = p.lower()
        if "disease" in low:
            if ":" in p:
                val = p.split(":", 1)[1].strip()
            else:
                m = re.search(r"disease[^\w]*(.+)", p, flags=re.IGNORECASE)
                val = m.group(1).strip(" :=") if m else p
            return strip_trailing_semicolon(val)

    # 2) 'diagnosis' / 'state' / 'stage'
    for p in parts:
        low = p.lower()
        if any(k in low for k in ["diagnosis", "state", "stage"]):
            if ":" in p:
                val = p.split(":", 1)[1].strip()
            else:
                m = re.search(r"(diagnosis|state|stage)[^\w]*(.+)", p, flags=re.IGNORECASE)
                val = m.group(2).strip(" :=") if m else p
            return strip_trailing_semicolon(val)

    return ""


# ---------- CLASSIFIER: sc / sn / spatial / bulk ----------

CATEGORY_KEYWORDS = {
    "spatial": [
        "spatial transcriptomics",
        "visium",
        "slide-seq",
        "slideseq",
        "merfish",
        "seqfish",
        "spatial gene expression",
        "spatial profiling",
        "st data",
        "ffpe",  # FFPE often indicates spatial
    ],
    "sn": [
        "single-nucleus",
        "single nucleus",
        "snrna",
        "sn-rna",
        "snrna-seq",
        "snatac",
        "sn-atac",
        "snatacseq",
        "single nucleus rna",
        "single nucleus atac",
        "single-nuc",
    ],
    "sc": [
        "single-cell",
        "single cell",
        "scrna",
        "sc-rna",
        "scrna-seq",
        "single cell rna",
        "single-cell rna",
        "scatac",
        "sc-atac",
        "scatac-seq",
        "single-cell atac",
        "single cell atac",
        "10x genomics chromium single cell",
        "10x chromium",
    ],
    "bulk": [
        "bulk rna",
        "bulk-rna",
        "bulk atac",
        "bulk-atac",
        "whole tissue",
        "bulk transcriptome",
        "whole kidney",
        "dissociated kidney",
    ],
}

SECONDARY_BULK_HINTS = [
    "rna-seq",
    "rna seq",
    "mrna-seq",
    "mrna seq",
    "transcriptomic",
    "atac-seq",
    "atac seq",
]


def classify_from_title(title: str) -> str:
    """
    Use *only* the GSM title to infer:
      - patterns with 'ST' or '_ST' or 'FFPE'  -> 'spatial'
      - patterns with 'sn' prefix    -> 'sn'
      - patterns with 'sc' prefix    -> 'sc'
      - patterns with 'whole' -> 'bulk'
    If nothing matches, return '' and let caller decide.
    """
    if not title:
        return ""
    t = title.lower()
    
    # Check for spatial: _ST or just ST as word boundary, or FFPE
    if re.search(r"[_\-]st\b|^st\b|\bst[_\-]|\bst$|\bffpe\b", t):
        return "spatial"
    
    # Check for single-nucleus: snRNA, sn-RNA, snRNAseq, snATAC, etc.
    # Look for 'sn' followed by 'rna' or 'atac' (with optional separators)
    if re.search(r"sn[_\-]?rna|sn[_\-]?atac|single[_\-]?nuc", t):
        return "sn"
    
    # Check for single-cell: scRNA, sc-RNA, scRNAseq, scATAC, etc.
    # Look for 'sc' followed by 'rna' or 'atac' (with optional separators)
    if re.search(r"sc[_\-]?rna|sc[_\-]?atac|single[_\-]?cell", t):
        return "sc"
    
    # Check for bulk indicators: "whole kidney", "whole tissue", "dissociated"
    if re.search(r"\bwhole\s+\w+|\bdissociated\b|\bbulk\b", t):
        return "bulk"

    return ""


def classify_sc_sn_spatial_bulk(
    title: str,
    library_strategy: str,
    library_source: str,
    extraction_protocol: str,
    html_text: str,
) -> str:
    """
    Logic:

      1. Build a big text blob from all metadata.
      2. See which categories (spatial/sn/sc/bulk) have *any* keyword hit.
      3. If 0 categories -> if only generic RNA-seq hints, call 'bulk', else ''.
      4. If exactly 1 category -> return that category.
      5. If >1 categories -> use TITLE ONLY as tie-breaker.
         If title doesn't help, return '' (ambiguous).
    """
    # 1) Combine all relevant text
    wide_text = " ".join(
        [
            title or "",
            library_strategy or "",
            library_source or "",
            extraction_protocol or "",
            html_text or "",
        ]
    ).lower()

    # 2) Which categories have any hit?
    categories_present: List[str] = []
    for cat, kws in CATEGORY_KEYWORDS.items():
        if any(k in wide_text for k in kws):
            categories_present.append(cat)

    # 3) No categories → maybe generic bulk, else unknown
    if not categories_present:
        if any(h in wide_text for h in SECONDARY_BULK_HINTS):
            return "bulk"
        return ""

    # 4) Exactly one category → just use it
    if len(categories_present) == 1:
        return categories_present[0]

    # 5) Multiple categories → use TITLE ONLY as tie-breaker
    label_from_title = classify_from_title(title or "")
    if label_from_title:
        return label_from_title

    # Title is useless & multiple categories matched → don't guess
    print(f"[WARN] Multiple categories {categories_present} found for title '{title}', cannot resolve")
    return ""


# --------- TECHNOLOGY DETECTOR (3′ / 5′) ----------
def detect_library_prep(
    extraction_protocol: str,
    title: str,
    library_strategy: str,
) -> str:
    pieces = [extraction_protocol or "", title or "", library_strategy or ""]
    t = " ".join(pieces).lower()
    if not t:
        return ""
    if "5'" in t or "5 prime" in t or "5′" in t:
        return "5′"
    if "3'" in t or "3 prime" in t or "3′" in t:
        return "3′"
    return ""


def fetch_gse_page(gse: str) -> str:
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}"
    r = requests.get(url, timeout=60)
    r.raise_for_status()
    return r.text


# --------- RAW COUNT DATA CHECK ----------
def infer_raw_count_flag(series_html: str) -> str:
    if not series_html:
        return ""

    t = series_html.lower()

    raw_like_keywords = [
        "matrix.mtx",
        "matrix.mtx.gz",
        "barcodes.tsv",
        "genes.tsv",
        "features.tsv",
        "count matrix",
        "counts.mtx",
        "umi counts",
        "umi_count",
        "umi matrix",
        "dge.txt",
        "dge.mtx",
        "filtered_gene_bc_matrices",
        "_raw",
        "raw_counts",
        "rawcount",
        "raw count matrix",
    ]

    processed_patterns = [
        ("reads and fpkm", "reads and fpkm"),
        ("fpkm", "fpkm"),
        ("tpm", "tpm"),
        ("rpkm", "rpkm"),
        ("normalized expression", "normalized expression"),
        ("log-normalized", "log-normalized"),
        ("lognormalized", "log-normalized"),
        ("vst counts", "vst counts"),
        ("cpm", "cpm"),
        ("processed data", "processed data"),
    ]

    has_raw_like = any(k in t for k in raw_like_keywords)

    matched_processed_labels = []
    for pat, label in processed_patterns:
        if pat in t:
            matched_processed_labels.append(label)

    matched_processed_labels = sorted(set(matched_processed_labels))

    if has_raw_like:
        return "yes"
    if matched_processed_labels:
        return "no, " + ", ".join(sorted(matched_processed_labels))

    return ""


# ----------------- MAIN PIPE -----------------
def process_gse(gse: str) -> List[Dict[str, Any]]:
    # inspect GSE page once
    try:
        series_html = fetch_gse_page(gse)
        raw_flag = infer_raw_count_flag(series_html)
    except Exception as e:
        print(f"[warn] could not inspect GSE page for {gse}: {e}")
        raw_flag = ""

    lines = download_soft(gse)
    samples = parse_samples_from_soft(lines)
    rows: List[Dict[str, Any]] = []

    for s in samples:
        gsm = s["gsm"]

        try:
            html = fetch_gsm_page(gsm)
            meta = scrape_gsm_meta_from_html(html)
            html_clean = clean_html_text(html)
        except Exception as e:
            print(f"[warn] scrape failed for {gsm}: {e}")
            meta = {
                "title": "",
                "organism": "",
                "source_name": "",
                "characteristics": "",
                "extraction_protocol": "",
                "instrument_model": "",
                "library_strategy": "",
                "library_source": "",
            }
            html_clean = ""

        organism = meta.get("organism", "").strip()
        if "homo sapiens" not in organism.lower():
            continue

        raw_source = "; ".join(s.get("source_name_ch1", [])) or meta.get("source_name", "") or ""
        source_name = strip_trailing_semicolon(raw_source)

        raw_chars = "; ".join(s.get("characteristics_ch1", [])) or meta.get("characteristics", "") or ""
        disease_diag = extract_disease_only(raw_chars)
        if not disease_diag.strip():
            disease_diag = "Normal"

        title = s.get("title") or meta.get("title", "") or ""
        lib_strategy = meta.get("library_strategy", "")
        lib_source = meta.get("library_source", "")
        extraction_protocol = meta.get("extraction_protocol", "")
        sequencer = meta.get("instrument_model", "")

        sc_sn = classify_sc_sn_spatial_bulk(
            title=title,
            library_strategy=lib_strategy,
            library_source=lib_source,
            extraction_protocol=extraction_protocol,
            html_text=html_clean,
        )

        technology = detect_library_prep(
            extraction_protocol=extraction_protocol,
            title=title,
            library_strategy=lib_strategy,
        )

        time.sleep(SLEEP_SEC)
        try:
            runs = get_runs_for_gsm(s)
        except Exception as e:
            print(f"[warn] run fetch failed for {gsm}: {e}")
            runs = []

        base_row = {
            "Dataset": gse,
            "Sample": gsm,
            "source_name": source_name,
            "Disease_diagnosis": disease_diag,
            "sc_sn": sc_sn,
            "Technology": technology,
            "Library_strategy": lib_strategy,
            "Sequencer": sequencer,
            "raw_count_mtx_available": raw_flag,
        }

        if runs:
            for run_id in runs:
                r = dict(base_row)
                r["run_id"] = run_id
                rows.append(r)
        else:
            r = dict(base_row)
            r["run_id"] = ""
            rows.append(r)

    return rows


def main():
    all_rows: List[Dict[str, Any]] = []
    for gse in GSE_LIST:
        print(f"[+] processing {gse}")
        rows = process_gse(gse)
        all_rows.extend(rows)

    if not all_rows:
        print("No human rows found. Then these GSEs genuinely have no GSM with Organism = Homo sapiens.")
        return

    df = pd.DataFrame(all_rows)
    df = df[
        [
            "Dataset",
            "Sample",
            "source_name",
            "Disease_diagnosis",
            "sc_sn",
            "Technology",
            "Library_strategy",
            "Sequencer",
            "run_id",
            "raw_count_mtx_available",
        ]
    ]
    df.to_excel(OUTPUT_XLSX, index=False)
    print(f"Wrote {len(df)} rows to {OUTPUT_XLSX}")


if __name__ == "__main__":
    main()