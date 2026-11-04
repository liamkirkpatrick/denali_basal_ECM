from pathlib import Path
import sys, os, json
os.environ['MPLBACKEND']='Agg'
os.environ['MPLCONFIGDIR']='/tmp/denali-rename-mpl'
ROOT=Path(__file__).resolve().parents[2]
os.chdir(ROOT)
sys.path.insert(0,str(ROOT/'agentic_workspace/dip_angles'))
import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager
for rel in ['scripts/plotting/basic_plots.ipynb','agentic_workspace/dip_angles/dip_angles.ipynb','agentic_workspace/dip_angles/manual_dip_angles.ipynb']:
    path=ROOT/rel
    nb=nbformat.read(path,as_version=4)
    km=KernelManager(kernel_name='python3')
    km.kernel_spec.argv[0]=sys.executable
    client=NotebookClient(nb,km=km,timeout=600,resources={'metadata':{'path':str(ROOT)}})
    client.execute()
    nbformat.write(nb,path)
    print('Refreshed',rel,flush=True)
import runpy
runpy.run_path(str(ROOT/'agentic_workspace/core_rename/refresh_saved_manual.py'))
