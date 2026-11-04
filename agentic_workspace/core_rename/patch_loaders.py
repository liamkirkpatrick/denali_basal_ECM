from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def edit(name, pairs):
    p=ROOT/name; s=p.read_text()
    for a,b in pairs:
        assert a in s,(name,a)
        s=s.replace(a,b)
    p.write_text(s)
edit('scripts/core_scripts/ECMclass.py',[
('from pathlib import Path','from pathlib import Path\nfrom core_names import canonical_core, source_core'),
("row = meta.loc[meta['core']==core]","core = canonical_core(core)\n        row = meta.loc[meta['core'].map(canonical_core)==core]"),
("self.core = row['core'].values[0]","self.core = canonical_core(row['core'].values[0])"),
("fname = self.core+'-'","fname = source_core(self.core)+'-'"),
("raw = pd.read_csv(Path(path_to_data) / fname)","self.source_filename = fname\n        raw = pd.read_csv(Path(path_to_data) / fname)")])
edit('scripts/core_scripts/button_edits_v3.py',[
('fname = f"{d.core}-{d.section}-{d.face}-{d.ACorDC}.csv"','fname = d.source_filename'),
("f'{Path(fname).stem}.png'","f'{d.core}-{d.section}-{d.face}-{d.ACorDC}.png'")])
edit('agentic_workspace/basic_plotting/build_notebook.py',[
('from ECMclass import ECM','from ECMclass import ECM\nfrom core_names import source_core'),
('f"{CORE}-{section}-{FACE}-{mode}.csv"','f"{source_core(CORE)}-{section}-{FACE}-{mode}.csv"')])
edit('agentic_workspace/dip_angles/feature_picker.py',[
('from ECMclass import ECM','from ECMclass import ECM\nfrom core_names import canonical_core, source_core'),
("item=ECM(self.core", "self.core = canonical_core(self.core)\n        item=ECM(self.core"),
("{args.core}-*-", "{source_core(args.core)}-*-")])
edit('agentic_workspace/dip_angles/dip_workflow.py',[
("FLAGS =", "import sys\nsys.path.insert(0, str(ROOT/'scripts/core_scripts'))\nfrom core_names import canonical_core, canonical_stem\n\nFLAGS ="),
('name = path.stem;', 'name = canonical_stem(path.stem);')])
edit('agentic_workspace/dip_angles/manual_workflow.py',[
('from feature_ranges import', 'from dip_workflow import canonical_core, canonical_stem\nfrom feature_ranges import'),
("selected=section_ranges(windows,core", "core=canonical_core(core)\n        selected=section_ranges(windows,core"),
("plot_section(path.stem+", "plot_section(canonical_stem(path.stem)+"),
("{path.stem}-manual-dip.png", "{canonical_stem(path.stem)}-manual-dip.png")])
