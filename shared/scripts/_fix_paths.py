# -*- coding: utf-8 -*-
"""一次性工具：把脚本中的硬编码项目根路径替换为基于 __file__ 计算的 ROOT。
仅做字符串级替换，不改动逻辑。干跑模式打印改动，确认后去掉 DRY_RUN=False 应用。"""
import os, re, sys

ROOT_LITERAL = "C:/D/workbuddy/科研/细胞外GPX4免疫抑制"
SUBDIRS = ["h2_bulk", "phase1_knowledge_gap_filling", "h1_pilot",
           "cgga_validation", "gse91061_validation", "immunotherapy_validation",
           "tcga_pancan", "comprehensive_lit_search", "pubmed_probe",
           "pubmed_probe2", "output"]

DRY_RUN = (len(sys.argv) > 1 and sys.argv[1] == "apply")

def root_def():
    # 向上寻找含 'output' 子目录的项目根
    return (
        'import os as _os\n'
        'def _project_root():\n'
        '    d = _os.path.dirname(_os.path.abspath(__file__))\n'
        '    while True:\n'
        '        if _os.path.isdir(_os.path.join(d, "output")):\n'
        '            return d\n'
        '        p = _os.path.dirname(d)\n'
        '        if p == d:\n'
        '            break\n'
        '        d = p\n'
        '    return _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))\n'
        'ROOT = _project_root()\n'
    )

def transform(text):
    changed = False
    # 1) 完整子目录路径字面量 -> os.path.join(ROOT, "output", d)
    for d in SUBDIRS:
        full = ROOT_LITERAL + "/output/" + d
        for q in ('"', "'"):
            needle = q + full + q
            repl = 'os.path.join(ROOT, "output", "' + d + '")'
            if needle in text:
                text = text.replace(needle, repl)
                changed = True
    # 2) 其余裸根字面量 -> ROOT
    for q in ('"', "'"):
        needle = q + ROOT_LITERAL + q
        if needle in text:
            text = text.replace(needle, 'ROOT')
            changed = True
    # 3) 形如 S + "/output/h2_bulk/..." 已在上一步变 ROOT + "..."，合法。
    return text, changed

def main():
    base = "C:/D/workbuddy/科研/细胞外GPX4免疫抑制/output"
    pyfiles = []
    for root, _, files in os.walk(base):
        for f in files:
            if f.endswith(".py") and f not in ("_fix_paths.py",) and "common" not in root:
                pyfiles.append(os.path.join(root, f))
    total = 0
    for fp in pyfiles:
        with open(fp, "r", encoding="utf-8") as fh:
            orig = fh.read()
        if ROOT_LITERAL not in orig:
            continue
        new, changed = transform(orig)
        if not changed:
            continue
        # 确保 import os 存在
        if not re.search(r"^\s*import\s+os\b", new, re.M):
            new = "import os\n" + new
        # 注入 ROOT 定义（若尚无）
        if "ROOT = " not in new and "_project_root" not in new:
            # 放在首行之后（跳过 shebang/encoding）
            lines = new.split("\n")
            idx = 0
            if lines and lines[0].startswith("#!"):
                idx = 1
            if idx < len(lines) and lines[idx].startswith("#"):
                idx += 1
            lines.insert(idx, root_def().rstrip("\n"))
            new = "\n".join(lines)
        total += 1
        if not DRY_RUN:
            with open(fp, "w", encoding="utf-8") as fh:
                fh.write(new)
            print(f"[APPLIED] {fp}")
        else:
            print(f"[DRYRUN ] {fp}")
    print(f"总计待处理文件数: {total}  (DRY_RUN={DRY_RUN})")

if __name__ == "__main__":
    main()
