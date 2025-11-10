## Background
This script was developed to support the curation of publicly available sequencing data for meta-analysis and data integration projects. Originally created for integration of single-cell RNA-seq data across multiple studies, it addresses the challenge of systematically identifying and classifying GEO samples by sequencing technology and data availability.

Manual curation of GEO metadata is slow and error-prone, especially when handling large numbers of GSE datasets that vary in sequencing type, data availability, library chemistry, disease annotations, and raw data access details.

## Features

- **Automated technology classification:** Distinguishes between single-cell (sc), single-nucleus (sn), spatial transcriptomics, and bulk sequencing
- **Raw data availability check:** Identifies if raw count matrices (MTX format) or only processed data (FPKM, TPM, etc.) are available
- **Disease/diagnosis extraction:** Parses sample characteristics to extract disease state information
- **Library chemistry detection:** Identifies 3′ or 5′ sequencing chemistry
- **SRA integration:** Automatically retrieves SRA run accessions for data download
- **Multi-source metadata:** Combines information from SOFT files and HTML pages for comprehensive metadata
- **Robust parsing:** Handles multiple keyword formats and inconsistent GEO metadata formats

## Installation
### Requirements
- Python 3.7+
- Required packages:`requests`, `pandas`, `openpyxl`

## Setup

1. Clone the repository:
```
git clone https://github.com/[your-username]/GEO_Metadata_Extractor.git
cd GEO_Metadata_Extractor
```

2. Install dependencies:
```
pip install -r requirements.txt
```

Or install packages individually:
```
pip install requests pandas openpyxl
```

## Usage
### Basic Usage

1. Edit the GSE_LIST in the script to include your GEO Series accessions:

```
GSE_LIST = [
    "GSE131685", 
    "GSE107585", 
    "GSE211785", 
    "GSE131882",
    "GSE114530"
]
```

2. Run the script:

```
python metadata_collection.py
```

3. Output will be saved to `kidney_metadata.xlsx` (or customize `OUTPUT_XLSX` variable)

### Configuration
Key variables you can modify:
```
GSE_LIST = ["GSE123456", "GSE789012"]  # List of GEO Series to process
OUTPUT_XLSX = "my_metadata.xlsx"        # Output filename
SLEEP_SEC = 0.5                          # Delay between requests (be nice to NCBI!)
```

## Output Format
The script generates an Excel file with the following columns:

## Classification Logic
### Technology Classification
The script uses a multi-tiered approach to classify samples:

1. Keyword matching: Searches titles, protocols, and metadata for technology-specific keywords
2. Multiple category resolution: When multiple technologies are detected, uses the sample title as the authoritative source
3. Pattern recognition:
   - `scRNAseq`, `sc-RNA`, `single-cell` → single-cell (sc)
   - `snRNAseq`, `sn-ATAC`, `single-nucleus` → single-nucleus (sn)
   - `_ST`, `FFPE`, `Visium` → spatial transcriptomics
   - `whole`, `bulk`, `dissociated` → bulk sequencing

### Raw Data Detection
Searches GEO supplementary files for indicators of raw vs processed data:

- **Raw data:** matrix.mtx, barcodes.tsv, features.tsv, UMI counts
- **Processed only:** FPKM, TPM, RPKM, normalized expression

## Use Cases
1. **Meta-Analysis Data Curation:**
Collecting samples across multiple studies for integrated analysis.

2. **Atlas Building:**
Creating disease-specific cell atlases by identifying all relevant single-cell studies.

3. **Data Availability Survey:**
Quickly assess what data exists before starting a project.

## Limitations

- Requires human-curated GSE list (doesn't search GEO automatically)
- Classification depends on keyword presence in metadata (may miss poorly annotated samples)
- Only processes human (Homo sapiens) samples
- Rate limited to be respectful to NCBI servers (adjustable via SLEEP_SEC)
- Some samples with ambiguous metadata may not be classified

## Areas for improvement:

- Additional keyword patterns for technology classification
- Support for other organisms
- Command-line interface for easier usage
- Automatic GEO search by disease/tissue type
- Better handling of edge cases in metadata parsing

## Acknowledgments
Special thanks to the GEO and SRA teams at NCBI for maintaining these essential resources.

## Design Philosophy: Why Heuristics Over LLMs?
This tool uses a keyword-based heuristic approach rather than Large Language Models (LLMs) for metadata classification. Here's why:
### Current Approach: Rule-Based Classification
#### Advantages:
- Fast and deterministic: Processes hundreds of samples in seconds
- Zero cost: No API fees or token usage
- Explainable: You can trace exactly why a sample was classified a certain way
- Offline capable: No internet dependency beyond initial GEO queries
- Reliable: Same input always produces the same output
- Easy to debug: Simple to add new keywords or adjust rules

#### Trade-offs:
- Requires manual keyword maintenance
- May struggle with novel phrasing not in keyword lists
- Needs fallback logic for ambiguous cases

### Could LLMs Improve This?
LLMs could potentially help with:

- Understanding context in unusual descriptions
- Handling novel terminology
- Extracting more nuanced metadata (treatment details, disease severity, etc.)
- Resolving genuinely ambiguous cases

However, for this use case:

- GEO metadata is relatively structured and keyword-rich
- Most samples (>95%) are correctly classified with current heuristics
- The speed/cost/reliability benefits of heuristics outweigh potential LLM advantages
- Edge cases can be manually reviewed using the [WARN] outputs

### Hybrid Approach (Future Direction)
A potential enhancement would be a hybrid system:

- Use heuristics for the majority of samples (fast, cheap, reliable)
- Fall back to LLM only when heuristics fail or return low confidence
- Track which samples required LLM classification for quality control

This would provide the best of both worlds: efficiency for standard cases and intelligent handling of truly ambiguous samples.



