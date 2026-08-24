# -*- coding: utf-8 -*-
"""Second verification: reverse-lookup candidate DOIs for suspicious entries to confirm title/year"""
import json, urllib.request, urllib.parse, time

UA = 'a2-ref-verify/1.0 (mailto:lm962272@gmail.com)'

def get_json(url, timeout=40):
    req = urllib.request.Request(url, headers={'User-Agent': UA})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode('utf-8'))

def check_doi(doi):
    try:
        d = get_json(f'https://api.crossref.org/works/{urllib.parse.quote(doi)}')
        m = d['message']
        title = ' '.join(m.get('title', []) or [])
        jr = '; '.join(m.get('container-title', []) or [])
        iss = m.get('issued', {}).get('date-parts', [[None]])[0][0]
        pp = m.get('published-print', {}).get('date-parts', [[None]])[0][0]
        po = m.get('published-online', {}).get('date-parts', [[None]])[0][0]
        vol = m.get('volume', ''); pg = m.get('page', '')
        return {'doi': doi, 'title': title, 'journal': jr, 'issued': iss,
                'print_year': pp, 'online_year': po, 'vol': vol, 'page': pg}
    except Exception as e:
        return {'doi': doi, 'error': str(e)}

candidates = {
    1:  ['10.1056/NEJMoa1003466', '10.1056/nejmx100063'],
    2:  ['10.1056/NEJMoa1200690', '10.1016/j.eururo.2014.12.052'],
    4:  ['10.1056/NEJMoa1302369', '10.1056/nejmc2501311'],
    10: ['10.1136/bmjebm-2023-112292'],
    23: ['10.1080/01616412.2025.2532039'],
    30: ['10.1038/nature25501', '10.3410/f.732660162.793559544'],
    31: ['10.1016/j.cell.2017.09.028', '10.3410/f.731996969.793570400'],
    32: ['10.1016/j.cell.2015.05.044', '10.1158/1538-7445.am2015-2972'],
}

expect = {
    1:  ('Improved Survival with Ipilimumab', 'Hodi'),
    2:  ('Safety, Activity, and Immune Correlates of Anti', 'Topalian'),
    4:  ('Nivolumab plus Ipilimumab in Advanced Melanoma', 'Wolchok'),
    10: ('Inverse publication reporting bias favouring null', 'Ioannidis'),
    23: ('Identification of potential biomarkers and therapeutic targets for cerebral venous', 'Song'),
    30: ('TGFbeta attenuates tumour response to PD-L1 blockade', 'Mariathasan'),
    31: ('Tumor and Microenvironment Evolution during Immunotherapy', 'Riaz'),
    32: ('Genomic Classification of Cutaneous Melanoma', 'Cancer Genome Atlas'),
}

for no, dois in candidates.items():
    print(f'--- Ref {no} ---')
    for doi in dois:
        r = check_doi(doi)
        t = r.get('title', '')
        exp = expect[no][0]
        ok = 'OK' if exp.lower() in t.lower() else 'MISMATCH'
        print(f"  [{ok}] {doi} | {t[:80]} | {r.get('journal','')[:40]} | issued={r.get('issued')} print={r.get('print_year')} online={r.get('online_year')} vol={r.get('vol')} pp={r.get('page')}")
        if r.get('error'):
            print(f"       ERROR: {r['error']}")
        time.sleep(0.2)
