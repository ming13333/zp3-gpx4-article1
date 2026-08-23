# -*- coding: utf-8 -*-
"""A2 参考文献 Crossref 批量核验脚本
用法: python verify_refs_a2.py <md_file> [--limit N] [--out report.csv]
"""
import re, sys, json, time, difflib, urllib.request, urllib.parse

UA = 'a2-ref-verify/1.0 (mailto:lm962272@gmail.com)'
BASE = 'https://api.crossref.org/works'

def get_json(url, params=None, timeout=40):
    if params:
        url = url + '?' + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={'User-Agent': UA})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode('utf-8'))

def parse_refs(md):
    refs = []
    lines = md.split('\n')
    in_refs = False
    for ln in lines:
        if re.match(r'^## References', ln):
            in_refs = True
            continue
        if in_refs:
            if re.match(r'^\d{2}\.\s', ln):
                refs.append(ln)
            elif ln.strip() == '':
                continue
            elif re.match(r'^#{1,3}\s', ln) and 'Reference' not in ln:
                break
    return refs

def extract(ref):
    ref = re.sub(r'^\d{2}\.\s*', '', ref).strip()
    star = ref.find('*')
    before = ref[:star] if star != -1 else ref
    # 标题: before 中最后一个 ".大写" 之后
    dots = list(re.finditer(r'\.\s[A-Z]', before))
    title = before[dots[-1].start()+2:].strip() if dots else before
    title = re.sub(r'\s+', ' ', title).strip()
    jm = re.search(r'\*([^*]+)\*', ref)
    journal = jm.group(1).strip() if jm else ''
    ym = re.search(r'\*\*(\d{4})\*\*', ref)
    year = ym.group(1) if ym else ''
    vm = re.search(r'\*\*(\d{4})\*\*,\s*\*(\d+)\*', ref)
    vol = vm.group(2) if vm else ''
    pm = re.search(r'\*\*(\d{4})\*\*,\s*\*(\d+)\*,\s*([\w\u2013\-]+(?:[.\w\u2013\-]+)*)', ref)
    pages = pm.group(3) if pm else ''
    return title, journal, year, vol, pages

def query(title):
    try:
        data = get_json(BASE, params={'query.bibliographic': title, 'rows': 20,
                        'select': 'DOI,title,author,issued,container-title,volume,page,type'})
        return data.get('message', {}).get('items', [])
    except Exception as e:
        return [{'__err__': str(e)}]

def year_of(it):
    try:
        p = it['issued']['date-parts'][0]
        return str(p[0]) if p and p[0] else ''
    except Exception:
        return ''

def best_match(title, year, items):
    best = None
    for it in items:
        if '__err__' in it:
            continue
        t = ' '.join(it.get('title', []) or []).strip()
        sim = difflib.SequenceMatcher(None, title.lower(), t.lower()).ratio()
        yr = year_of(it)
        score = sim + (0.15 if yr == year else 0.0)
        if best is None or score > best[0]:
            best = (score, sim, it, t, yr)
    return best

def main(md_path, limit=None, out='ref_verify_report.csv'):
    md = open(md_path, encoding='utf-8').read()
    refs = parse_refs(md)
    if limit:
        refs = refs[:limit]
    rows = []
    for i, ref in enumerate(refs, 1):
        title, journal, year, vol, pages = extract(ref)
        items = query(title)
        b = best_match(title, year, items)
        if b is None:
            rows.append({'no': i, 'status': 'NOT_FOUND', 'title': title, 'journal': journal,
                         'year': year, 'sim': '', 'doi': '', 'match_title': '', 'match_year': ''})
            print(f'[{i}] NOT_FOUND  {title[:70]}')
        else:
            score, sim, it, mt, myr = b
            doi = it.get('DOI', '')
            mj = '; '.join(it.get('container-title', []) or [])[:60]
            if sim >= 0.85 and myr == year:
                st = 'VERIFIED'
            elif sim >= 0.85 and myr != year and myr:
                st = 'YEAR_MISMATCH'
            elif sim >= 0.60:
                st = 'POSSIBLE'
            else:
                st = 'NOT_FOUND'
            rows.append({'no': i, 'status': st, 'title': title, 'journal': journal,
                         'year': year, 'sim': f'{sim:.2f}', 'doi': doi,
                         'match_title': mt, 'match_year': myr, 'match_journal': mj,
                         'match_vol': it.get('volume', ''), 'match_page': it.get('page', '')})
            print(f'[{i}] {st:<13} sim={sim:.2f} yr={year}/{myr} doi={doi}  {title[:50]}')
        time.sleep(0.35)
    import csv
    with open(out, 'w', newline='', encoding='utf-8-sig') as f:
        w = csv.DictWriter(f, fieldnames=['no','status','title','journal','year','sim','doi','match_title','match_year','match_journal','match_vol','match_page'])
        w.writeheader()
        w.writerows(rows)
    from collections import Counter
    print('\n=== SUMMARY ===')
    print(Counter(r['status'] for r in rows))
    print('report ->', out)

if __name__ == '__main__':
    md = sys.argv[1] if len(sys.argv) > 1 else 'Article2_A2_Cancers格式_v0.4.md'
    lim = int(sys.argv[sys.argv.index('--limit')+1]) if '--limit' in sys.argv else None
    out = sys.argv[sys.argv.index('--out')+1] if '--out' in sys.argv else 'ref_verify_report.csv'
    main(md, lim, out)
