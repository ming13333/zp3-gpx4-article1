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
HPA (Human Protein Atlas) ZP3 蛋白表达查询脚本
查询 ZP3 蛋白在正常组织和肿瘤组织中的表达数据
"""

import requests
import pandas as pd
import json
import os
import time
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns

# 配置
HPA_API_BASE = "https://www.proteinatlas.org/api"
OUTPUT_DIR = os.path.join(ROOT, "article1", "results", "hpa")
GENE_NAME = "ZP3"
ENSEMBL_ID = "ENSG00000188010"  # ZP3 的 Ensembl ID

# 创建输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)

def hpa_api_request(endpoint, params=None):
    """发送 HPA API 请求"""
    url = f"{HPA_API_BASE}/{endpoint}"
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
    }
    
    try:
        response = requests.get(url, params=params, headers=headers, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"API 请求失败: {e}")
        return None

def get_gene_info():
    """获取基因基本信息"""
    print(f"1. 获取 {GENE_NAME} 基因信息...")
    endpoint = f"gene/{ENSEMBL_ID}"
    data = hpa_api_request(endpoint)
    
    if data:
        print(f"   基因名称: {data.get('gene_name', 'N/A')}")
        print(f"   Ensembl ID: {data.get('ensembl_id', 'N/A')}")
        print(f"   蛋白名称: {data.get('protein_name', 'N/A')}")
        print(f"   基因描述: {data.get('gene_description', 'N/A')[:100]}...")
        return data
    return None

def get_tissue_expression():
    """获取组织表达数据"""
    print(f"\n2. 获取 {GENE_NAME} 组织表达数据...")
    
    # HPA 的 tissue expression endpoint
    endpoint = f"gene/{ENSEMBL_ID}/tissue_expression"
    data = hpa_api_request(endpoint)
    
    if data:
        print(f"   获取到 {len(data)} 个组织表达记录")
        return data
    else:
        print("   尝试备用 endpoint...")
        # 尝试其他可能的 endpoint
        endpoints_to_try = [
            f"expression/{ENSEMBL_ID}/tissue",
            f"tissue/{ENSEMBL_ID}",
            f"gene/{ENSEMBL_ID}/expression"
        ]
        
        for ep in endpoints_to_try:
            data = hpa_api_request(ep)
            if data:
                print(f"   使用 endpoint: {ep}")
                return data
        
        print("   无法获取组织表达数据")
        return None

def get_cancer_expression():
    """获取肿瘤组织表达数据"""
    print(f"\n3. 获取 {GENE_NAME} 肿瘤表达数据...")
    
    endpoints_to_try = [
        f"gene/{ENSEMBL_ID}/cancer_expression",
        f"cancer/{ENSEMBL_ID}",
        f"gene/{ENSEMBL_ID}/expression/cancer"
    ]
    
    for ep in endpoints_to_try:
        data = hpa_api_request(ep)
        if data:
            print(f"   使用 endpoint: {ep}")
            return data
    
    print("   无法获取肿瘤表达数据")
    return None

def get_subcellular_localization():
    """获取亚细胞定位数据"""
    print(f"\n4. 获取 {GENE_NAME} 亚细胞定位数据...")
    
    endpoints_to_try = [
        f"gene/{ENSEMBL_ID}/subcellular_location",
        f"subcellular/{ENSEMBL_ID}",
        f"gene/{ENSEMBL_ID}/location"
    ]
    
    for ep in endpoints_to_try:
        data = hpa_api_request(ep)
        if data:
            print(f"   使用 endpoint: {ep}")
            return data
    
    print("   无法获取亚细胞定位数据")
    return None

def get_antibody_data():
    """获取抗体数据（蛋白验证）"""
    print(f"\n5. 获取 {GENE_NAME} 抗体数据...")
    
    endpoints_to_try = [
        f"gene/{ENSEMBL_ID}/antibody",
        f"antibody/{ENSEMBL_ID}",
        f"gene/{ENSEMBL_ID}/antibodies"
    ]
    
    for ep in endpoints_to_try:
        data = hpa_api_request(ep)
        if data:
            print(f"   使用 endpoint: {ep}")
            return data
    
    print("   无法获取抗体数据")
    return None

def save_data(data, filename):
    """保存数据为 JSON 和 CSV"""
    if not data:
        return
    
    # 保存 JSON
    json_path = os.path.join(OUTPUT_DIR, f"{filename}.json")
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"   已保存: {json_path}")
    
    # 尝试转换为 CSV（如果数据结构允许）
    try:
        if isinstance(data, list) and len(data) > 0:
            df = pd.json_normalize(data)
            csv_path = os.path.join(OUTPUT_DIR, f"{filename}.csv")
            df.to_csv(csv_path, index=False, encoding='utf-8-sig')
            print(f"   已保存: {csv_path}")
            return df
        elif isinstance(data, dict):
            # 尝试展平字典
            df = pd.json_normalize(data)
            csv_path = os.path.join(OUTPUT_DIR, f"{filename}.csv")
            df.to_csv(csv_path, index=False, encoding='utf-8-sig')
            print(f"   已保存: {csv_path}")
            return df
    except Exception as e:
        print(f"   转换 CSV 失败: {e}")
    
    return None

def visualize_tissue_expression(tissue_data):
    """可视化组织表达数据"""
    if not tissue_data:
        return
    
    print(f"\n6. 生成可视化...")
    
    try:
        # 尝试创建 DataFrame
        if isinstance(tissue_data, list):
            df = pd.json_normalize(tissue_data)
        elif isinstance(tissue_data, dict):
            df = pd.json_normalize(tissue_data)
        else:
            print("   无法解析数据结构")
            return
        
        # 查找可能的表达值列
        value_cols = [col for col in df.columns if 'expression' in col.lower() or 'value' in col.lower() or 'level' in col.lower()]
        tissue_cols = [col for col in df.columns if 'tissue' in col.lower() or 'organ' in col.lower()]
        
        if value_cols and tissue_cols:
            # 创建柱状图
            plt.figure(figsize=(12, 6))
            sns.barplot(data=df, x=tissue_cols[0], y=value_cols[0])
            plt.xticks(rotation=45, ha='right')
            plt.title(f'{GENE_NAME} Protein Expression Across Tissues (HPA)')
            plt.ylabel('Expression Level')
            plt.tight_layout()
            
            fig_path = os.path.join(OUTPUT_DIR, f"fig_{GENE_NAME.lower()}_tissue_expression.png")
            plt.savefig(fig_path, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"   已保存可视化: {fig_path}")
        else:
            print(f"   未找到合适的列进行可视化。可用列: {list(df.columns)[:10]}...")
    except Exception as e:
        print(f"   可视化失败: {e}")

def main():
    """主函数"""
    print(f"=" * 60)
    print(f"HPA (Human Protein Atlas) {GENE_NAME} 蛋白表达查询")
    print(f"查询时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"=" * 60)
    
    # 1. 获取基因信息
    gene_info = get_gene_info()
    if gene_info:
        save_data(gene_info, f"hpa_{GENE_NAME.lower()}_gene_info")
    
    time.sleep(1)  # 避免请求过快
    
    # 2. 获取组织表达数据
    tissue_data = get_tissue_expression()
    if tissue_data:
        df_tissue = save_data(tissue_data, f"hpa_{GENE_NAME.lower()}_tissue_expression")
        visualize_tissue_expression(tissue_data)
    
    time.sleep(1)
    
    # 3. 获取肿瘤表达数据
    cancer_data = get_cancer_expression()
    if cancer_data:
        save_data(cancer_data, f"hpa_{GENE_NAME.lower()}_cancer_expression")
    
    time.sleep(1)
    
    # 4. 获取亚细胞定位数据
    location_data = get_subcellular_localization()
    if location_data:
        save_data(location_data, f"hpa_{GENE_NAME.lower()}_subcellular_location")
    
    time.sleep(1)
    
    # 5. 获取抗体数据
    antibody_data = get_antibody_data()
    if antibody_data:
        save_data(antibody_data, f"hpa_{GENE_NAME.lower()}_antibody")
    
    print(f"\n" + "=" * 60)
    print(f"HPA 查询完成")
    print(f"输出目录: {OUTPUT_DIR}")
    print(f"=" * 60)
    
    # 创建摘要报告
    create_summary_report()

def create_summary_report():
    """创建摘要报告"""
    report = f"""# HPA {GENE_NAME} 蛋白表达查询报告

## 查询信息
- 查询时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
- 基因: {GENE_NAME} ({ENSEMBL_ID})
- 数据源: Human Protein Atlas (HPA)

## 查询内容
1. 基因基本信息
2. 正常组织表达数据
3. 肿瘤组织表达数据
4. 亚细胞定位数据
5. 抗体验证数据

## 输出文件
- `hpa_{GENE_NAME.lower()}_gene_info.json` - 基因信息
- `hpa_{GENE_NAME.lower()}_tissue_expression.json/csv` - 组织表达数据
- `hpa_{GENE_NAME.lower()}_cancer_expression.json/csv` - 肿瘤表达数据
- `hpa_{GENE_NAME.lower()}_subcellular_location.json/csv` - 亚细胞定位
- `hpa_{GENE_NAME.lower()}_antibody.json/csv` - 抗体数据
- `fig_{GENE_NAME.lower()}_tissue_expression.png` - 组织表达可视化

## 分析要点
1. **正常组织表达**: ZP3 在哪些正常组织中表达？（重点关注脑组织）
2. **肿瘤组织表达**: 与正常相比，肿瘤中 ZP3 表达如何变化？
3. **亚细胞定位**: 膜定位 vs 胞质定位（区分经典 vs Cancer 异构体）
4. **蛋白验证**: HPA 免疫组化数据是否支持 RNA 表达结果？

## 注意事项
- HPA 数据基于抗体检测，可能存在交叉反应
- 蛋白表达与 RNA 表达可能不完全一致
- 需要结合其他数据库（如GTEx、TCGA）进行验证
"""
    
    report_path = os.path.join(OUTPUT_DIR, f"hpa_{GENE_NAME.lower()}_query_report.md")
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(report)
    print(f"已创建摘要报告: {report_path}")

if __name__ == "__main__":
    main()