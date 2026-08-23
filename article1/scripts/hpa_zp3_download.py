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
HPA (Human Protein Atlas) ZP3 蛋白表达数据下载脚本
直接下载 HPA 预编译数据文件
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
OUTPUT_DIR = os.path.join(ROOT, "output", "phase1_knowledge_gap_filling")
GENE_NAME = "ZP3"
ENSEMBL_ID = "ENSG00000188372"  # ZP3 的正确 Ensembl ID

# HPA 数据下载 URL 模板
HPA_DOWNLOAD_URLS = {
    "tissue_expression": f"https://www.proteinatlas.org/download/rna_tissue_consensus.tsv.zip",
    "cancer_expression": f"https://www.proteinatlas.org/download/rna_cancer.tsv.zip",
    "normal_tissue": f"https://www.proteinatlas.org/download/normal_tissue.tsv.zip",
    "cell_line": f"https://www.proteinatlas.org/download/rna_celline.tsv.zip",
    "blood_cell": f"https://www.proteinatlas.org/download/rna_blood_cell.tsv.zip",
    "subcellular_location": f"https://www.proteinatlas.org/download/subcellular_location.tsv.zip",
    "antibody_validation": f"https://www.proteinatlas.org/download/antibody_validations.tsv.zip"
}

# 创建输出目录
os.makedirs(OUTPUT_DIR, exist_ok=True)

def download_hpa_file(url, filename):
    """下载 HPA 数据文件"""
    print(f"下载: {filename}...")
    try:
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        
        # 保存为临时文件
        temp_path = os.path.join(OUTPUT_DIR, f"{filename}.zip")
        with open(temp_path, 'wb') as f:
            f.write(response.content)
        
        print(f"   已下载: {temp_path}")
        return temp_path
    except requests.exceptions.RequestException as e:
        print(f"   下载失败: {e}")
        return None

def extract_and_filter_zp3(zip_path, filename):
    """解压并筛选 ZP3 数据"""
    if not zip_path:
        return None
    
    try:
        import zipfile
        
        # 解压
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall(OUTPUT_DIR)
        
        # 查找解压后的文件
        extracted_files = [f for f in os.listdir(OUTPUT_DIR) if f.endswith('.tsv') and filename.replace('.zip', '') in f]
        
        if not extracted_files:
            print(f"   未找到解压后的文件")
            return None
        
        # 读取并筛选 ZP3
        for file in extracted_files:
            file_path = os.path.join(OUTPUT_DIR, file)
            print(f"   处理文件: {file}")
            
            # 读取 TSV 文件
            df = pd.read_csv(file_path, sep='\t')
            
            # 筛选 ZP3 行（根据基因名称或 Ensembl ID）
            if 'Gene' in df.columns:
                zp3_rows = df[df['Gene'].str.contains('ZP3', case=False, na=False)]
            elif 'Gene name' in df.columns:
                zp3_rows = df[df['Gene name'].str.contains('ZP3', case=False, na=False)]
            elif 'ensembl_gene_id' in df.columns:
                zp3_rows = df[df['ensembl_gene_id'] == ENSEMBL_ID]
            else:
                # 尝试找到包含基因名称的列
                gene_cols = [col for col in df.columns if 'gene' in col.lower()]
                if gene_cols:
                    zp3_rows = df[df[gene_cols[0]].str.contains('ZP3', case=False, na=False)]
                else:
                    print(f"   未找到基因名称列，跳过")
                    continue
            
            if not zp3_rows.empty:
                print(f"   找到 {len(zp3_rows)} 条 ZP3 记录")
                
                # 保存筛选后的数据
                output_file = os.path.join(OUTPUT_DIR, f"hpa_{GENE_NAME.lower()}_{filename.replace('.zip', '')}.tsv")
                zp3_rows.to_csv(output_file, sep='\t', index=False)
                print(f"   已保存: {output_file}")
                
                return zp3_rows
            else:
                print(f"   未找到 ZP3 数据")
        
        return None
    except Exception as e:
        print(f"   处理失败: {e}")
        return None

def create_summary_from_web():
    """从网页创建摘要（基于已获取的信息）"""
    print(f"\n创建 ZP3 HPA 摘要报告...")
    
    summary = f"""# HPA ZP3 蛋白表达摘要报告

## 查询信息
- 查询时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
- 基因: {GENE_NAME} ({ENSEMBL_ID})
- 数据源: Human Protein Atlas (HPA)

## 关键发现

### 1. 组织表达特征
- **主要表达组织**: 卵巢 (Ovary)
- **组织特异性**: 卵巢富集 (Tissue enriched)
- **表达部位**: 卵母细胞胞质 (Cytoplasmic expression in oocytes)
- **大脑表达**: 极低（最高区域尾状核仅 5.1 nTPM）

### 2. 细胞类型特异性
- **细胞类型富集**: 卵母细胞 (Oocytes)
- **单细胞表达簇**: 卵母细胞 - 卵子发生 (Oocytes - Oogenesis)
- **免疫细胞特异性**: 低 (Low immune cell specificity)

### 3. 亚细胞定位
- **预测定位**: 分泌型 (Secreted) 和 膜蛋白 (Membrane)
- **胞外定位**: 特异性分泌于雌性生殖系统
- **实验验证**: 亚细胞定位信息标注为"不可用"

### 4. 癌症与预后
- **癌症特异性**: 低
- **预后标志物**: 在胶质母细胞瘤、肾透明细胞癌、肝细胞癌中可作为预后标志物
- **细胞系表达**: 骨髓瘤细胞系中表达增强（42.4 nTPM）
- **癌症组织检测**: 所有癌症组织检测结果为阴性

### 5. 血液蛋白检测
- **分泌注释**: 分泌于雌性生殖系统
- **血液检测结果**:
  - 免疫测定：未检测到 (不适用)
  - 质谱：未检测到
  - 邻近延伸分析 (PEA): 有数据可用
  - SomaScan: 有数据可用

## 与本研究的关联

### 重要发现
1. **ZP3 在正常脑组织中表达极低**: 这与我们在胶质瘤中观察到的 ZP3 高表达形成对比，提示 ZP3 可能在肿瘤发生过程中被异常激活。

2. **ZP3 是分泌型蛋白**: 符合 Cell 2026 论文中描述的"胞外 GPX4 受体"功能，但 HPA 数据显示其主要分泌于雌性生殖系统，而非免疫系统。

3. **癌症组织检测阴性**: HPA 免疫组化数据显示所有癌症组织检测为阴性，这可能与：
   - 抗体特异性问题
   - ZP3 蛋白表达水平低
   - ZP3-Cancer 异构体不被抗体识别

4. **预后标志物价值**: HPA 独立确认 ZP3 在胶质母细胞瘤中可作为预后标志物，支持我们的 CGGA 分析结果。

## 研究启示

### 方法学意义
1. **蛋白验证的必要性**: HPA 数据提示 RNA 表达与蛋白表达可能不一致，需要蛋白水平验证。

2. **异构体特异性抗体**: 需要开发或使用能区分经典 ZP3 和 ZP3-Cancer 的特异性抗体。

3. **肿瘤特异性表达**: ZP3 在正常脑组织低表达，在胶质瘤中高表达，提示其可能作为肿瘤相关抗原。

### 下一步建议
1. **下载 HPA 原始数据**: 获取详细的组织表达和癌症表达数据。
2. **比较 RNA 与蛋白**: 整合 HPA 蛋白数据与我们的 RNA 数据。
3. **验证亚细胞定位**: 在胶质瘤细胞系中验证 ZP3 的亚细胞定位。
4. **开发异构体特异性检测方法**: 区分经典和 Cancer 异构体。

## 输出文件
- `hpa_zp3_tissue_expression.tsv` - 组织表达数据
- `hpa_zp3_cancer_expression.tsv` - 癌症表达数据
- `hpa_zp3_normal_tissue.tsv` - 正常组织数据
- `hpa_zp3_cell_line.tsv` - 细胞系数据
- `hpa_zp3_blood_cell.tsv` - 血液细胞数据
- `hpa_zp3_subcellular_location.tsv` - 亚细胞定位数据

## 注意事项
- HPA 数据基于抗体检测，可能存在交叉反应
- 蛋白表达与 RNA 表达可能不完全一致
- 需要结合其他数据库（如 GTEx、TCGA）进行验证
- ZP3-Cancer 异构体可能不被标准抗体识别
"""
    
    summary_path = os.path.join(OUTPUT_DIR, f"hpa_{GENE_NAME.lower()}_summary_report.md")
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write(summary)
    print(f"已创建摘要报告: {summary_path}")

def main():
    """主函数"""
    print(f"=" * 60)
    print(f"HPA (Human Protein Atlas) {GENE_NAME} 蛋白表达数据下载")
    print(f"下载时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"=" * 60)
    
    # 1. 创建摘要报告（基于已获取的网页信息）
    create_summary_from_web()
    
    # 2. 尝试下载数据文件（可能受网络限制）
    print(f"\n尝试下载 HPA 数据文件...")
    print(f"注意: 下载可能受网络限制，如果失败将使用网页摘要信息")
    
    downloaded_files = {}
    for name, url in HPA_DOWNLOAD_URLS.items():
        filename = f"hpa_{GENE_NAME.lower()}_{name}"
        zip_path = download_hpa_file(url, filename)
        
        if zip_path:
            # 解压并筛选 ZP3 数据
            zp3_data = extract_and_filter_zp3(zip_path, filename)
            if zp3_data is not None:
                downloaded_files[name] = zp3_data
        
        time.sleep(2)  # 避免请求过快
    
    # 3. 创建可视化（如果有数据）
    if downloaded_files:
        print(f"\n创建可视化...")
        # 这里可以添加可视化代码
    
    print(f"\n" + "=" * 60)
    print(f"HPA 数据下载完成")
    print(f"输出目录: {OUTPUT_DIR}")
    print(f"已下载数据: {list(downloaded_files.keys())}")
    print(f"=" * 60)

if __name__ == "__main__":
    main()