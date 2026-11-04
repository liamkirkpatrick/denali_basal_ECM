"""One-time core-name migration, preserving measurement CSV filenames."""
from pathlib import Path
import re, hashlib, json
ROOT = Path(__file__).resolve().parents[2]
raw_dirs = [ROOT/'data/clean', ROOT/'data/not_cleaned']
raw_hashes = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for d in raw_dirs for p in d.glob('*') if p.is_file()}
(ROOT/'agentic_workspace/core_rename/raw_sha256.json').write_text(json.dumps(raw_hashes, indent=2))
for p in ROOT.rglob('*'):
    if not p.is_file() or any(x in p.parts for x in ('.git','__pycache__','reference','core_rename')) or any(d in p.parents for d in raw_dirs): continue
    if p.suffix not in ('.py','.ipynb','.md','.csv','.json','.ndjson','.mjs','.txt','.yml'): continue
    s=p.read_text()
    # Keep literal source measurement filenames and glob expressions intact.
    protected=[]
    def protect(m):
        protected.append(m.group()); return f'__RAW_NAME_{len(protected)-1}__'
    t=re.sub(r'dic1-[^\s/\"\'\\]*?\.csv',protect,s,flags=re.I)
    t=re.sub('dic1','DEN13A',t,flags=re.I)
    for i,name in enumerate(protected): t=t.replace(f'__RAW_NAME_{i}__',name)
    if t!=s: p.write_text(t)
for p in sorted(ROOT.rglob('*'), key=lambda p:len(p.parts), reverse=True):
    if any(x in p.parts for x in ('.git','reference','core_rename')) or any(d in p.parents for d in raw_dirs): continue
    if re.search('dic1',p.name,re.I):
        target=p.with_name(re.sub('dic1','DEN13A',p.name,flags=re.I))
        assert not target.exists(),target
        p.rename(target)
