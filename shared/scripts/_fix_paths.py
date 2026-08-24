# -*- coding: utf-8 -*-
"""One-time tool: replace hardcoded project root paths in scripts with ROOT computed from __file__.
Performs string-level replacement only, without changing logic. Dry-run mode prints changes; after confirming, remove DRY_RUN=False to apply."""
import os, re, sys

ROOT_LITERAL = "C:/D/workbuddy/research/extracellular_GPX4_immunosuppression"
SUBDIRS = ["h2_bulk", "phase1_knowledge_gap_filling", "h1_pilot",
           "cgga_validation", "gse91061_validation", "immunotherapy_validation",
           "tcga_pancan", "comprehensive_lit_search", "pubmed_probe",
           "pubmed_probe2", "output"]

DRY_RUN = (len(sys.argv) > 1 and sys.argv[1] == "apply")

def root_def():
    # Search upward for the project root containing an 'output' subdirectory
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
    # 1) Full subdirectory path literals -> os.path.join(ROOT, "output", d)
    for d in SUBDIRS:
        full = ROOT_LITERAL + "/output/" + d
        for q in ('"', "'"):
            needle = q + full + q
            repl = 'os.path.join(ROOT, "output", "' + d + '")'
            if needle in text:
                text = text.replace(needle, repl)
                changed = True
    # 2) Remaining bare root literals -> ROOT
    for q in ('"', "'"):
        needle = q + ROOT_LITERAL + q
        if needle in text:
            text = text.replace(needle, 'ROOT')
            changed = True
    # 3) Forms like S + "/output/h2_bulk/..." already became ROOT + "..." in the previous step, which is valid.
    return text, changed

def main():
    base = "C:/D/workbuddy/research/extracellular_GPX4_immunosuppression/output"
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
        # Ensure import os exists
        if not re.search(r"^\s*import\s+os\b", new, re.M):
            new = "import os\n" + new
        # Inject ROOT definition (if not already present)
        if "ROOT = " not in new and "_project_root" not in new:
            # Place after the first line (skip shebang/encoding)
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
    print(f"Total files to process: {total}  (DRY_RUN={DRY_RUN})")

if __name__ == "__main__":
    main()
