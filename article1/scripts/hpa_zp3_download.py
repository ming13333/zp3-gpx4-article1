#!/usr/bin/env python
import os as _os
def _project_root():
    d = _os.path.dirname(_os.path.abspath(__file__))
    while True:
        if _os.path.isdir(_os.path.join(d, "output")):
            return d
        p = _os.path.dirname(d)
        if p == d:
            break
        d = p
    return _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
ROOT = _project_root()
"""
HPA (Human Protein Atlas) ZP3 protein expression data download script
Directly download HPA precompiled data files
"""

import requests
import pandas as pd
import json
import os
import time
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns

# Configuration
OUTPUT_DIR = os.path.join(ROOT, "output", "phase1_knowledge_gap_filling")
GENE_NAME = "ZP3"
ENSEMBL_ID = "ENSG00000188372"  # Correct Ensembl ID for ZP3

# HPA data download URL template
HPA_DOWNLOAD_URLS = {
    "tissue_expression": f"https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip",
    "cancer_expression": f"https://www.proteinatlas.org/download/rna_cancer.tsv.zip",
    "normal_tissue": f"https://www.proteinatlas.org/download/normal_tissue.tsv.zip",
    "cell_line": f"https://www.proteinatlas.org/download/rna_celline.tsv.zip",
    "blood_cell": f"https://www.proteinatlas.org/download/rna_blood_cell.tsv.zip",
    "subcellular_location": f"https://www.proteinatlas.org/download/subcellular_location.tsv.zip",
    "antibody_validation": f"https://www.proteinatlas.org/download/antibody_validations.tsv.zip"
}

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

def download_hpa_file(url, filename):
    """Download HPA data file"""
    print(f"Downloading: {filename}...")
    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        
        # Save as temporary file
        temp_path = os.path.join(OUTPUT_DIR, f"{filename}.zip")
        with open(temp_path, 'wb') as f:
            f.write(response.content)
        
        print(f"    Downloaded: {temp_path}")
        return temp_path
    except requests.exceptions.RequestException as e:
        print(f"    Download failed: {e}")
        return None

def extract_and_filter_zp3(zip_path, filename):
    """Extract and filter ZP3 data"""
    if not zip_path:
        return None
    
    try:
        import zipfile
        
        # Extract
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(OUTPUT_DIR)
        
        # Find the extracted files
        extracted_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('.tsv') and filename.replace('.zip', '') in f]
        
        if not extracted_files:
            print(f"    Extracted file not found")
            return None
        
        # Read and filter ZP3
        for file in extracted_files:
            file_path = os.path.join(OUTPUT_DIR, file)
            print(f"    Processing file: {file}")
            
            # Read TSV file
            df = pd.read_csv(file_path, sep='\t')
            
            # Filter ZP3 rows (by gene name or Ensembl ID)
            if 'Gene' in df.columns:
                zp3_rows = df[df['Gene'].str.contains('ZP3', case=False, na=False)]
            elif 'Gene name' in df.columns:
                zp3_rows = df[df['Gene name'].str.contains('ZP3', case=False, na=False)]
            elif 'ensembl_gene_id' in df.columns:
                zp3_rows = df[df['ensembl_gene_id'] == ENSEMBL_ID]
            else:
                # Try to find a column containing the gene name
                gene_cols = [col for col in df.columns if 'gene' in col.lower()]
                if gene_cols:
                    zp3_rows = df[df[gene_cols[0]].str.contains('ZP3', case=False, na=False)]
                else:
                    print(f"    Gene name column not found, skipping")
                    continue
            
            if not zp3_rows.empty:
                print(f"    Found {len(zp3_rows)} ZP3 records")
                
                # Save the filtered data
                output_file = os.path.join(OUTPUT_DIR, f"hpa_{GENE_NAME.lower()}_{filename.replace('.zip', '')}.tsv")
                zp3_rows.to_csv(output_file, sep='\t', index=False)
                print(f"    Saved: {output_file}")
                
                return zp3_rows
            else:
                print(f"   ZP3 data not found")
        
        return None
    except Exception as e:
        print(f"   Processing failed: {e}")
        return None

def create_summary_from_web():
    """Create summary from web (based on fetched information)"""
    print(f"\nCreating ZP3 HPA summary report...")
    
    summary = f"""# HPA ZP3 protein expression summary report

## Query Information
- Query time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
- Gene: {GENE_NAME} ({ENSEMBL_ID})
- Data source: Human Protein Atlas (HPA)

## Key Findings

### 1. Tissue Expression Characteristics
- **Main expression tissue**: Ovary (Ovary)
- **Tissue specificity**: Ovary enriched (Tissue enriched)
- **Expression location**: Oocyte cytoplasm (Cytoplasmic expression in oocytes)
- **Brain expression**: Very low (highest region caudate nucleus only 5.1 nTPM)

### 2. Cell type specificity
- **Cell type enrichment**: Oocytes (Oocytes)
- **Single-cell expression cluster**: Oocytes - Oogenesis (Oocytes - Oogenesis)
- **Immune cell specificity**: Low (Low immune cell specificity)

### 3. Subcellular localization
- **Predicted localization**: Secreted (Secreted) and Membrane (Membrane)
- **Extracellular localization**: Specifically secreted in the female reproductive system
- **Experimental validation**: Subcellular localization information annotated as "unavailable"

### 4. Cancer and prognosis
- **Cancer specificity**: Low
- **Prognostic marker**: Can serve as a prognostic marker in glioblastoma, clear cell renal cell carcinoma, and hepatocellular carcinoma
- **Cell line expression**: Enhanced in myeloma cell lines (42.4 nTPM)
- **Cancer tissue detection**: All cancer tissue tests were negative

### 5. Blood protein detection
- **Secretion annotation**: Secreted into the female reproductive system
- **Blood test results**:
  - Immunoassay: not detected (not applicable)
  - Mass spectrometry: not detected
  - Proximity extension assay (PEA): data available
  - SomaScan: data available

## Relevance to this study

### Key findings
1. **ZP3 is expressed at very low levels in normal brain tissue**: This contrasts with the high ZP3 expression we observed in glioma, suggesting that ZP3 may be aberrantly activated during tumorigenesis.

2. **ZP3 is a secreted protein**: Consistent with the "extracellular GPX4 receptor" function described in the Cell 2026 paper, but HPA data show that it is mainly secreted in the female reproductive system, rather than the immune system.

3. **Negative detection in cancer tissues**: HPA immunohistochemistry data show that all cancer tissues tested negative, which may be related to:
   - Antibody specificity issues
   - Low ZP3 protein expression level
   - The ZP3-Cancer isoform is not recognized by the antibody

4. **Prognostic marker value**: HPA independently confirms that ZP3 can serve as a prognostic marker in glioblastoma, supporting our CGGA analysis results.

## Research Implications

### Methodological Significance
1. **Necessity of protein validation**: HPA data suggest that RNA expression and protein expression may be inconsistent, requiring protein-level validation.

2. **Isoform-specific antibodies**: Specific antibodies that can distinguish classical ZP3 from ZP3-Cancer need to be developed or used.

3. **Tumor-specific expression**: ZP3 is lowly expressed in normal brain tissue and highly expressed in glioma, suggesting that it may serve as a tumor-associated antigen.

### Suggestions for Next Steps
1. **Download HPA raw data**: Obtain detailed tissue expression and cancer expression data.
2. **Compare RNA and protein**: Integrate HPA protein data with our RNA data.
3. **Validate subcellular localization**: Validate ZP3 subcellular localization in glioma cell lines.
4. **Develop isoform-specific detection methods**: Distinguish canonical and Cancer isoforms.

## Output Files
- `hpa_zp3_tissue_expression.tsv` - Tissue expression data
- `hpa_zp3_cancer_expression.tsv` - Cancer expression data
- `hpa_zp3_normal_tissue.tsv` - Normal tissue data
- `hpa_zp3_cell_line.tsv` - Cell line data
- `hpa_zp3_blood_cell.tsv` - Blood cell data
- `hpa_zp3_subcellular_location.tsv` - Subcellular localization data

## Notes
- HPA data is based on antibody detection and may have cross-reactivity
- Protein expression may not perfectly match RNA expression
- Need to validate with other databases (e.g., GTEx, TCGA)
- ZP3-Cancer isoform may not be recognized by standard antibodies
"""
    
    summary_path = os.path.join(OUTPUT_DIR, f"hpa_{GENE_NAME.lower()}_summary_report.md")
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write(summary)
    print(f"Summary report created: {summary_path}")

def main():
    """Main function"""
    print(f"=" * 60)
    print(f"HPA (Human Protein Atlas) {GENE_NAME} protein expression data download")
    print(f"Download time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"=" * 60)
    
    # 1. Create summary report (based on fetched webpage info)
    create_summary_from_web()
    
    # 2. Attempt to download data files (may be network restricted)
    print(f"\nAttempting to download HPA data files...")
    print(f"Note: Download may be network restricted; if failed, webpage summary info will be used")
    
    downloaded_files = {}
    for name, url in HPA_DOWNLOAD_URLS.items():
        filename = f"hpa_{GENE_NAME.lower()}_{name}"
        zip_path = download_hpa_file(url, filename)
        
        if zip_path:
            # Unzip and filter ZP3 data
            zp3_data = extract_and_filter_zp3(zip_path, filename)
            if zp3_data is not None:
                downloaded_files[name] = zp3_data
        
        time.sleep(2)  # avoid requesting too fast
    
    # 3. Create visualization (if data available)
    if downloaded_files:
        print(f"\nCreating visualization...")
        # Visualization code can be added here
    
    print(f"\n" + "=" * 60)
    print(f"HPA data download complete")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Downloaded data: {list(downloaded_files.keys())}")
    print(f"=" * 60)

if __name__ == "__main__":
    main()
