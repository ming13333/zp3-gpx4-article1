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
import os
# -*- coding: utf-8 -*-
"""Launch resume_downloads.py detached from the session (DETACHED_PROCESS)"""
import subprocess, sys, os

BASE = ros.path.join(ROOT, "output", "phase1_knowledge_gap_filling")
DETACHED_PROCESS = 0x00000008
CREATE_NEW_PROCESS_GROUP = 0x00000200

p = subprocess.Popen(
    [sys.executable, os.path.join(BASE, "resume_downloads.py")],
    cwd=BASE,
    creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP,
    stdin=subprocess.DEVNULL,
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
    close_fds=True,
)
print("launched_pid=%d" % p.pid)
