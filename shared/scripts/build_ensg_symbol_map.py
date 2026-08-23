# -*- coding: utf-8 -*-
"""
构建 ensg -> symbol 全基因组映射（缓存到本地）
用于 GSEA 排序列表的基因名转换。
数据源：Ensembl REST /lookup/id 批量查询（分批）。
缓存：ensg_symbol_map.json
"""
import os
import json
import time
import urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(BASE, "ensg_symbol_map.json")


def build_map(ensg_list, batch=1000, delay=0.2):
    ensg_list = sorted(set(ensg_list))
    out = {}
    url = "https://rest.ensembl.org/lookup/id"
    for i in range(0, len(ensg_list), batch):
        chunk = ensg_list[i:i + batch]
        req = urllib.request.Request(
            url,
            data=json.dumps({"ids": chunk}).encode(),
            headers={"Content-Type": "application/json",
                     "Accept": "application/json"})
        ok = False
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req, timeout=60) as r:
                    d = json.load(r)
                for k, v in d.items():
                    if v:
                        out[k] = v.get("display_name") or k
                ok = True
                break
            except Exception as e:
                print(f"  批次 {i//batch} 重试 {attempt+1}: {str(e)[:50]}")
                time.sleep(2)
        if not ok:
            print(f"  批次 {i//batch} 失败")
        print(f"  进度: {min(i+batch, len(ensg_list))}/{len(ensg_list)}")
        time.sleep(delay)
    return out


def main():
    # 从 ranklist 读取 ensg
    rnk_path = os.path.join(BASE, "zp3_gsea_results", "fl_vs_ri_ranklist.csv")
    import pandas as pd
    rnk = pd.read_csv(rnk_path)
    ensgs = [g for g in rnk["symbol"] if g.startswith("ENSG")]
    print(f"需要映射的 ensg 数: {len(ensgs)}")
    m = build_map(ensgs)
    print(f"映射成功: {len(m)}/{len(ensgs)}")
    with open(CACHE, "w") as f:
        json.dump(m, f)
    print(f"缓存已存: {CACHE}")
    # 抽查
    for g in ["ENSG00000188372", "ENSG00000141510"]:
        print(f"  {g} -> {m.get(g)}")


if __name__ == "__main__":
    main()
