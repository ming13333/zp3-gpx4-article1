# -*- coding: utf-8 -*-
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
"""Session-detached serial download daemon: melanoma resume -> GTEx -> Toil -> ALL_DONE
Launch: launch_detached.py spawns an independent process, unaffected by session interruption
Status check: download_status.txt  (real-time)
Log check: download_progress.log (appended)
All complete: ALL_DONE.txt
"""
import os, sys, time, requests, subprocess

BASE = ros.path.join(ROOT, "output", "phase1_knowledge_gap_filling")
SC_DIR = os.path.join(BASE, "sc_data")
LOG = os.path.join(BASE, "download_progress.log")
STATUS = os.path.join(BASE, "download_status.txt")

MEL_URL = "https://datasets.cellxgene.cziscience.com/1b76227b-c731-4807-9487-ad5e4d24e0d0.h5ad"
MEL_PATH = os.path.join(SC_DIR, "melanoma_myeloid.h5ad")
MEL_TARGET = 2606540453


def log(msg):
    try:
        with open(LOG, "a", encoding="utf-8") as f:
            f.write("[%s] %s\n" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg))
    except Exception:
        pass


def set_status(msg):
    try:
        with open(STATUS, "w", encoding="utf-8") as f:
            f.write(msg + "\n")
    except Exception:
        pass


def resume_download(url, path, target, label):
    existing = os.path.getsize(path) if os.path.exists(path) else 0
    if existing >= target:
        log("%s: already complete (%d)" % (label, existing))
        set_status(label + ": DONE")
        return True
    headers = {"Range": "bytes=%d-" % existing}
    log("%s: resuming from %d, target %d" % (label, existing, target))
    set_status("%s: resuming %d/%d (%.1f%%)" % (label, existing, target, existing * 100.0 / target))
    try:
        r = requests.get(url, headers=headers, stream=True, timeout=(30, 300))
        if r.status_code == 206:
            mode = "ab"
        else:
            mode = "wb"
            log("%s: server returned %d, downloading from scratch" % (label, r.status_code))
            existing = 0
        last_t = time.time()
        with open(path, mode) as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    f.write(chunk)
                    existing += len(chunk)
                    if time.time() - last_t >= 30:
                        last_t = time.time()
                        set_status("%s: resuming %d/%d (%.1f%%)" % (label, existing, target, existing * 100.0 / target))
        ok = os.path.getsize(path) >= target
        log("%s: %s, final size %d" % (label, "complete" if ok else "target not reached", os.path.getsize(path)))
        set_status(label + (": DONE" if ok else ": incomplete"))
        return ok
    except Exception as e:
        log("%s: exception %r" % (label, e))
        set_status(label + ": exception " + str(e))
        return False


def run_script(script_name, label):
    script = os.path.join(BASE, script_name)
    out_log = os.path.join(BASE, label + "_run.log")
    log("%s: starting %s" % (label, script_name))
    set_status(label + ": running")
    try:
        with open(out_log, "a", encoding="utf-8") as gl:
            subprocess.run([sys.executable, script], cwd=BASE, stdout=gl,
                           stderr=subprocess.STDOUT, timeout=4 * 3600)
        log("%s: script finished" % label)
        set_status(label + ": finished (see log)")
        return True
    except Exception as e:
        log("%s: exception %r" % (label, e))
        set_status(label + ": exception " + str(e))
        return False


def main():
    log("=== serial download daemon started (pid=%d) ===" % os.getpid())
    # 1. Melanoma resume download (CELLxGENE)
    mel_ok = resume_download(MEL_URL, MEL_PATH, MEL_TARGET, "Melanoma")
    # 2. GTEx extraction (toil.xenahubs.net)
    gtex_ok = run_script("gtex_brain_zp3_integrated.py", "GTEx")
    # 3. Toil isoforms (toil.xenahubs.net)
    toil_ok = run_script("toil_download_pure_python.py", "Toil")
    # All done
    with open(os.path.join(BASE, "ALL_DONE.txt"), "w", encoding="utf-8") as f:
        f.write("melanoma=%s\ngtex=%s\ntoil=%s\n" % (mel_ok, gtex_ok, toil_ok))
    log("=== All done: melanoma=%s gtex=%s toil=%s ===" % (mel_ok, gtex_ok, toil_ok))
    set_status("ALL_DONE")


if __name__ == "__main__":
    main()
