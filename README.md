# evolution_of_hypermutation
This repo provides code and files to reproduce analyses for my paper on the evolution of hypermutated tumors.  This readme will explain the file structure, documentation and reproducibility instructions.  <br>
**Environment**: All jupyter notebooks and python scripts within this repository are run in python 2.7.  All R scripts are run using R version 3.5. Scripts not contained in this notebook and used to generate files are referred to in <code>generate_all_files_and_plots.ipynb</code> and are hosted on the juno computational cluster<br>
**Files** This repository contains all ipython, R and python scripts used to generate project plots except where otherwise specified in <code>generate_all_files_and_plots.ipynb</code>.  All files in the /files directory are hosted on google drive at __link___.  Download the directory at this __link__ and drag it to the files directory  <br>
<br>
**Directory Map** 
```
evolution_of_hypermutation
│   README.md
│   fileDirectory.txt (identifiers and paths to all project files)    
│   generate_all_files_and_plots.ipynb (ipython cells/written documentation of where files came from/how to create them)
│
└───scripts
│   └───figure1
|   |   make_figure_1.ipynb (all code to make figure 1 dataframes)
|   |   make_figure1_supplementary_figures.ipynb (all code to make figure S1 dataframes)
|   |   plotFigure1.R (plots all figure1 panels)
|   |   plotSupplementaryFiguresFig1.R (plots all Supplemental Figure 1 panels)
|   └───FIGURE1_PLOTTING_FILES
|   |   └───figurePdfs (PDFs of figures)
|.  
|   |   └───plotDataFiles (tsvs used to generate pdf plots

│   └───figure2 (same structure figure 1)
│   └───figure3 (same structure figure 1)
│   └───utilityScripts
│
└───files
    │   file021.txt
```
