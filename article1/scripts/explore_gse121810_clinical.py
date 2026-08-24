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
"""
Explore clinical information of GSE121810
"""
import os
import gzip
import re
import pandas as pd

OUT_DIR = os.path.join(ROOT, "article1", "results")
series_file = os.path.join(OUT_DIR, 'GSE121810_series_matrix.txt.gz')

print('Parsing series matrix...')
lines = []
with gzip.open(series_file, 'rt') as f:
    for line in f:
        if line.startswith('!series_matrix_table_begin'):
            break
        if line.startswith('!'):
            lines.append(line.rstrip())

print(f'Read {len(lines)} metadata lines')

# Extract all fields
fields = {}
for line in lines:
    m = re.match(r'!(\w+)\s+(.*)', line, re.S)
    if not m:
        continue
    key = m.group(1)
    cols = m.group(2).split('\t')
    cols = [c.strip('"') for c in cols]
    fields[key] = cols

print(f'Extracted {len(fields)} fields')

# View all fields
print('\nAll fields:')
for key, values in fields.items():
    print(f'{key}: {len(values)} values')
    if len(values) <= 5:
        print(f'  Values: {values}')
    else:
        print(f'  First 5 values: {values[:5]}')

# Focus on characteristic fields
print('\n=== Characteristic Field Details ===')
characteristic_fields = [k for k in fields.keys() if k.startswith('Sample_characteristics')]
for field in characteristic_fields:
    values = fields[field]
    print(f'\n{field} (length {len(values)}):')
    # Count unique values
    unique_vals = list(set(values))
    print(f'  Unique value count: {len(unique_vals)}')
    if len(unique_vals) <= 10:
        print(f'  Unique values: {unique_vals}')
    else:
        print(f'  Examples: {unique_vals[:10]}')

# Extract treatment group information
print('\n=== Treatment Group Information ===')
if 'Sample_characteristics_ch1' in fields:
    therapy_info = fields['Sample_characteristics_ch1']
    print(f'Treatment info: {therapy_info[:5]}...')
    
    # Count treatment types
    therapy_counts = pd.Series(therapy_info).value_counts()
    print(f'Treatment type distribution:')
    for therapy, count in therapy_counts.items():
        print(f'  {therapy}: {count} cases')

# Extract survival information
print('\n=== Survival Information ===')
survival_keywords = ['survival', 'os', 'overall', 'time', 'month', 'day', 'year', 'status', 'event', 'alive', 'dead', 'death', 'progression', 'response', 'outcome']
for key, values in fields.items():
    if any(keyword in key.lower() for keyword in survival_keywords):
        print(f'{key}: {values[:5]}...')

# Extract patient ID mapping
print('\n=== Patient ID Mapping ===')
if 'Sample_title' in fields:
    titles = fields['Sample_title']
    descriptions = fields.get('Sample_description', [])
    print(f'Sample titles (first 10):')
    for i, title in enumerate(titles[:10]):
        desc = descriptions[i] if i < len(descriptions) else 'N/A'
        print(f'  {title}: {desc}')

# Save all fields
all_fields_csv = os.path.join(OUT_DIR, 'gse121810_all_fields.csv')
all_fields_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in fields.items()]))
all_fields_df.to_csv(all_fields_csv, index=False)
print(f'\nAll fields saved: {all_fields_csv}')

# Try to extract treatment response information
print('\n=== Try to extract treatment response information ===')
# Extract from Sample_source_name_ch1
if 'Sample_source_name_ch1' in fields:
    source_names = fields['Sample_source_name_ch1']
    print(f'Sample sources (first 10):')
    for i, source in enumerate(source_names[:10]):
        print(f'  {i}: {source}')
    
    # Count source types
    source_counts = pd.Series(source_names).value_counts()
    print(f'Source type distribution:')
    for source, count in source_counts.items():
        print(f'  {source}: {count} cases')

# Extract from Sample_characteristics_ch2
if 'Sample_characteristics_ch2' in fields:
    ch2 = fields['Sample_characteristics_ch2']
    print(f'\nCharacteristic 2 (first 10):')
    for i, val in enumerate(ch2[:10]):
        print(f'  {i}: {val}')
    
    # Count unique values
    unique_ch2 = list(set(ch2))
    print(f'Number of unique values in characteristic 2: {len(unique_ch2)}')
    if len(unique_ch2) <= 10:
        print(f'Unique values of characteristic 2: {unique_ch2}')
