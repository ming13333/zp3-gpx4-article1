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
"""脱离会话的串行下载守护：黑色素瘤续传 → GTEx → Toil → ALL_DONE
启动方式: launch_detached.py 拉起到独立进程, 会话中断不影响
状态查询: download_status.txt  (实时)
日志查询: download_progress.log (追加)
全部完成: ALL_DONE.txt
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
        log("%s: 已完整 (%d)" % (label, existing))
        set_status(label + ": DONE")
        return True
    headers = {"Range": "bytes=%d-" % existing}
    log("%s: 从 %d 续传, 目标 %d" % (label, existing, target))
    set_status("%s: 续传中 %d/%d (%.1f%%)" % (label, existing, target, existing * 100.0 / target))
    try:
        r = requests.get(url, headers=headers, stream=True, timeout=(30, 300))
        if r.status_code == 206:
            mode = "ab"
        else:
            mode = "wb"
            log("%s: 服务器返回 %d, 从头下载" % (label, r.status_code))
            existing = 0
        last_t = time.time()
        with open(path, mode) as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    f.write(chunk)
                    existing += len(chunk)
                    if time.time() - last_t >= 30:
                        last_t = time.time()
                        set_status("%s: 续传中 %d/%d (%.1f%%)" % (label, existing, target, existing * 100.0 / target))
        ok = os.path.getsize(path) >= target
        log("%s: %s, 最终大小 %d" % (label, "完成" if ok else "未达目标", os.path.getsize(path)))
        set_status(label + (": DONE" if ok else ": 未完成"))
        return ok
    except Exception as e:
        log("%s: 异常 %r" % (label, e))
        set_status(label + ": 异常 " + str(e))
        return False


def run_script(script_name, label):
    script = os.path.join(BASE, script_name)
    out_log = os.path.join(BASE, label + "_run.log")
    log("%s: 启动 %s" % (label, script_name))
    set_status(label + ": 运行中")
    try:
        with open(out_log, "a", encoding="utf-8") as gl:
            subprocess.run([sys.executable, script], cwd=BASE, stdout=gl,
                           stderr=subprocess.STDOUT, timeout=4 * 3600)
        log("%s: 脚本结束" % label)
        set_status(label + ": 结束(见日志)")
        return True
    except Exception as e:
        log("%s: 异常 %r" % (label, e))
        set_status(label + ": 异常 " + str(e))
        return False


def main():
    log("=== 串行下载守护启动 (pid=%d) ===" % os.getpid())
    # 1. 黑色素瘤续传 (CELLxGENE)
    mel_ok = resume_download(MEL_URL, MEL_PATH, MEL_TARGET, "黑色素瘤")
    # 2. GTEx 提取 (toil.xenahubs.net)
    gtex_ok = run_script("gtex_brain_zp3_integrated.py", "GTEx")
    # 3. Toil 异构体 (toil.xenahubs.net)
    toil_ok = run_script("toil_download_pure_python.py", "Toil")
    # 全部完成
    with open(os.path.join(BASE, "ALL_DONE.txt"), "w", encoding="utf-8") as f:
        f.write("melanoma=%s\ngtex=%s\ntoil=%s\n" % (mel_ok, gtex_ok, toil_ok))
    log("=== 全部完成: melanoma=%s gtex=%s toil=%s ===" % (mel_ok, gtex_ok, toil_ok))
    set_status("ALL_DONE")


if __name__ == "__main__":
    main()
