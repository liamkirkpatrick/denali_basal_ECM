# denali_basal_ECM

Summary
-------
This repository contains ECM data from the basal portion of the Denali Ice Core, processing scripts, and plotting utilities used to analyze and visualize the datasets produced by this analysis.

File structure
--------------
- **data/**  _Holds datasets, including pre and post-proccessing._
  - **data/not-cleaned/**  _Original unprocessed exports and raw text/CSV files. _
  - **data/clean/**  _Mirrors not-cleaned, but after running through button_edits (manual cleaning)._
  - **data/processed/**  _Contains processed datasets (none so far)._
    
- **scripts/**  _Python scripts and Jupyter notebooks used to process data and generate figures. Key subfolders:_
  - **scripts/core_scripts/**  _Core utilities and helper modules for data parsing, the ECM class implementation, and reusable functions (e.g., [ECMclass](/Users/Liam/Desktop/UW/ECM/COLDEX_ECM/COLDEX-ECM-data/scripts/core_scripts/ECMclass.py:1))._
  - **scripts/plotting/**  _Plot generation scripts and example notebooks. Includes notebooks like `basic_plots.ipynb` which creates a basic visual_
  - **scripts/processing/**  _Notebook-based helpers for data reformatting and one-line averaging._

- **figures/**  _Outputs from the plotting scripts organized by plot type and sample._
  - **figures/TBD** _to be filled in as figures created_

- **environment.yml**  _Conda environment specification (name: `COLDEX_ECM`). Use this to recreate the Python environment used for analysis._
- **.gitignore**  _Git ignore rules._
- **README.md**  _This file._
- **.vscode/**  _Editor settings and optional workspace configuration._

Notes
-----

If any of the above folders/files are missing or have different contents on your copy, please inspect the repository root for the exact layout.

Quick start / Instructions for use
---------------------------------
1. Clone the repository:

   git clone <repo-url>
   cd denali_basal_ECM

2. Create and activate the conda environment (recommended):

   conda env create -f environment.yml
   conda activate COLDEX_ECM


   Notes:
   - COLDEC_ECM isn't a typo - I'm copying an environment from some coldex work.
   - The environment file declares `name: COLDEX_ECM` and a prefix; adjust as needed if you prefer a different environment name or location.
   - If conda is not available, recreate the environment using your preferred Python package manager. The `environment.yml` lists the required packages and Python version (3.13.x).

3. Run processing scripts or notebooks:

   - To run the basic plotting script from the command line:
     python scripts/plotting/basic_plots.py

   - To open notebooks, launch JupyterLab or Jupyter Notebook from the activated environment so that kernels match:
     jupyter lab

   - Processing utilities are in `scripts/core_scripts/`

Core naming
-----------
The core is named **DEN13A** in metadata, analyses, figures, and processed output names.
Original measurement CSV filenames in `data/clean/` and `data/not_cleaned/` retain
the `dic1-` prefix. `scripts/core_scripts/core_names.py` maps between the public
core name and that source filename prefix; source paths in provenance tables
therefore still use the original filenames.
