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
HPA (Human Protein Atlas) ZP3 protein expression query script
Query ZP3 protein expression data in normal tissues and tumor tissues
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
HPA_API_BASE = "https://www.proteinatlas.org/api"
OUTPUT_DIR = os.path.join(ROOT, "article1", "results", "hpa")
GENE_NAME = "ZP3"
ENSEMBL_ID = "ENSG00000188010"  # Ensembl ID of ZP3

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

def hpa_api_request(endpoint, params=None):
    """Send HPA API request"""
    url = f"{HPA_API_BASE}/{endpoint}"
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
    }
    
    try:
        response = requests.get(url, params=params, headers=headers, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"API request failed: {e}")
        return None

def get_gene_info():
    """Get basic gene information"""
    print(f"1. Fetching {GENE_NAME} gene information...")
    endpoint = f"gene/{ENSEMBL_ID}"
    data = hpa_api_request(endpoint)
    
    if data:
        print(f"   Gene name: {data.get('gene_name', 'N/A')}")
        print(f"   Ensembl ID: {data.get('ensembl_id', 'N/A')}")
        print(f"   Protein name: {data.get('protein_name', 'N/A')}")
        print(f"   Gene description: {data.get('gene_description', 'N/A')[:100]}...")
        return data
    return None

def get_tissue_expression():
    """Get tissue expression data"""
    print(f"\n2. Fetching {GENE_NAME} tissue expression data...")
    
    # HPA tissue expression endpoint
    endpoint = f"gene/{ENSEMBL_ID}/tissue_expression"
    data = hpa_api_request(endpoint)
    
    if data:
        print(f"   Retrieved {len(data)} tissue expression records")
        return data
    else:
        print("   Trying fallback endpoint...")
        # Try other possible endpoints
        endpoints_to_try = [
            f"expression/{ENSEMBL_ID}/tissue",
            f"tissue/{ENSEMBL_ID}",
            f"gene/{ENSEMBL_ID}/expression"
        ]
        
        for ep in endpoints_to_try:
            data = hpa_api_request(ep)
            if data:
                print(f"   Using endpoint: {ep}")
                return data
        
        print("   Unable to retrieve tissue expression data")
        return None

def get_cancer_expression():
    """Get cancer tissue expression data"""
    print(f"\n3. Fetching {GENE_NAME} cancer expression data...")
    
    endpoints_to_try = [
        f"gene/{ENSEMBL_ID}/cancer_expression",
        f"cancer/{ENSEMBL_ID}",
        f"gene/{ENSEMBL_ID}/expression/cancer"
    ]
    
    for ep in endpoints_to_try:
        data = hpa_api_request(ep)
        if data:
            print(f"   Using endpoint: {ep}")
            return data
    
    print("   Unable to retrieve cancer expression data")
    return None

def get_subcellular_localization():
    """Get subcellular localization data"""
    print(f"\n4. Fetching {GENE_NAME} subcellular localization data...")
    
    endpoints_to_try = [
        f"gene/{ENSEMBL_ID}/subcellular_location",
        f"subcellular/{ENSEMBL_ID}",
        f"gene/{ENSEMBL_ID}/location"
    ]
    
    for ep in endpoints_to_try:
        data = hpa_api_request(ep)
        if data:
            print(f"   Using endpoint: {ep}")
            return data
    
    print("   Unable to retrieve subcellular localization data")
    return None

def get_antibody_data():
    """Get antibody data (protein validation)"""
    print(f"\n5. Fetching {GENE_NAME} antibody data...")
    
    endpoints_to_try = [
        f"gene/{ENSEMBL_ID}/antibody",
        f"antibody/{ENSEMBL_ID}",
        f"gene/{ENSEMBL_ID}/antibodies"
    ]
    
    for ep in endpoints_to_try:
        data = hpa_api_request(ep)
        if data:
            print(f"   Using endpoint: {ep}")
            return data
    
    print("   Unable to retrieve antibody data")
    return None

def save_data(data, filename):
    """Save data as JSON and CSV"""
    if not data:
        return
    
    # Save JSON
    json_path = os.path.join(OUTPUT_DIR, f"{filename}.json")
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"   Saved: {json_path}")
    
    # Try converting to CSV (if data structure allows)
    try:
        if isinstance(data, list) and len(data) > 0:
            df = pd.json_normalize(data)
            csv_path = os.path.join(OUTPUT_DIR, f"{filename}.csv")
            df.to_csv(csv_path, index=False, encoding='utf-8-sig')
            print(f"   Saved: {csv_path}")
            return df
        elif isinstance(data, dict):
            # Try to flatten the dictionary
            df = pd.json_normalize(data)
            csv_path = os.path.join(OUTPUT_DIR, f"{filename}.csv")
            df.to_csv(csv_path, index=False, encoding='utf-8-sig')
            print(f"   Saved: {csv_path}")
            return df
    except Exception as e:
        print(f"   CSV conversion failed: {e}")
    
    return None

def visualize_tissue_expression(tissue_data):
    """Visualize tissue expression data"""
    if not tissue_data:
        return
    
    print(f"\n6. Generating visualization...")
    
    try:
        # Try to create DataFrame
        if isinstance(tissue_data, list):
            df = pd.json_normalize(tissue_data)
        elif isinstance(tissue_data, dict):
            df = pd.json_normalize(tissue_data)
        else:
            print("   Unable to parse data structure")
            return
        
        # Look for possible expression value columns
        value_cols = [col for col in df.columns if 'expression' in col.lower() or 'value' in col.lower() or 'level' in col.lower()]
        tissue_cols = [col for col in df.columns if 'tissue' in col.lower() or 'organ' in col.lower()]
        
        if value_cols and tissue_cols:
            # Create bar chart
            plt.figure(figsize=(12, 6))
            sns.barplot(data=df, x=tissue_cols[0], y=value_cols[0])
            plt.xticks(rotation=45, ha='right')
            plt.title(f'{GENE_NAME} Protein Expression Across Tissues (HPA)')
            plt.ylabel('Expression Level')
            plt.tight_layout()
            
            fig_path = os.path.join(OUTPUT_DIR, f"fig_{GENE_NAME.lower()}_tissue_expression.png")
            plt.savefig(fig_path, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"   Saved visualization: {fig_path}")
        else:
            print(f"   No suitable columns found for visualization. Available columns: {list(df.columns)[:10]}...")
    except Exception as e:
        print(f"   Visualization failed: {e}")

def main():
    """Main function"""
    print(f"=" * 60)
    print(f"HPA (Human Protein Atlas) {GENE_NAME} protein expression query")
    print(f"Query time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"=" * 60)
    
    # 1. Get gene information
    gene_info = get_gene_info()
    if gene_info:
        save_data(gene_info, f"hpa_{GENE_NAME.lower()}_gene_info")
    
    time.sleep(1)  # Avoid requesting too fast
    
    # 2. Get tissue expression data
    tissue_data = get_tissue_expression()
    if tissue_data:
        df_tissue = save_data(tissue_data, f"hpa_{GENE_NAME.lower()}_tissue_expression")
        visualize_tissue_expression(tissue_data)
    
    time.sleep(1)
    
    # 3. Get cancer expression data
    cancer_data = get_cancer_expression()
    if cancer_data:
        save_data(cancer_data, f"hpa_{GENE_NAME.lower()}_cancer_expression")
    
    time.sleep(1)
    
    # 4. Get subcellular localization data
    location_data = get_subcellular_localization()
    if location_data:
        save_data(location_data, f"hpa_{GENE_NAME.lower()}_subcellular_location")
    
    time.sleep(1)
    
    # 5. Get antibody data
    antibody_data = get_antibody_data()
    if antibody_data:
        save_data(antibody_data, f"hpa_{GENE_NAME.lower()}_antibody")
    
    print(f"\n" + "=" * 60)
    print(f"HPA query completed")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"=" * 60)
    
    # Create summary report
    create_summary_report()

def create_summary_report():
    """Create summary report"""
    report = f"""# HPA {GENE_NAME} Protein Expression Query Report

## Query Information
- Query time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
- Gene: {GENE_NAME} ({ENSEMBL_ID})
- Data source: Human Protein Atlas (HPA)

## Query Contents
1. Basic gene information
2. Normal tissue expression data
3. Tumor tissue expression data
4. Subcellular localization data
5. Antibody validation data

## Output files
- `hpa_{GENE_NAME.lower()}_gene_info.json` - Gene information
- `hpa_{GENE_NAME.lower()}_tissue_expression.json/csv` - Tissue expression data
- `hpa_{GENE_NAME.lower()}_cancer_expression.json/csv` - Cancer expression data
- `hpa_{GENE_NAME.lower()}_subcellular_location.json/csv` - Subcellular localization
- `hpa_{GENE_NAME.lower()}_antibody.json/csv` - Antibody data
- `fig_{GENE_NAME.lower()}_tissue_expression.png` - Tissue expression visualization

## Analysis points
1. **Normal tissue expression**: In which normal tissues is ZP3 expressed? (Focus on brain tissue)
2. **Tumor tissue expression**: Compared to normal, how does ZP3 expression change in tumors?
3. **Subcellular localization**: Membrane localization vs cytoplasmic localization (distinguishing classic vs Cancer isoforms)
4. **Protein validation**: Does HPA immunohistochemistry data support the RNA expression results?

## Notes
- HPA data is based on antibody detection and may have cross-reactivity
- Protein expression may not be fully consistent with RNA expression
- Needs to be validated in combination with other databases (e.g., GTEx, TCGA)
"""
    
    report_path = os.path.join(OUTPUT_DIR, f"hpa_{GENE_NAME.lower()}_query_report.md")
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"Summary report created: {report_path}")

if __name__ == "__main__":
    main()
