### An easy-to-use wrapper function to run celltype identifiers in R

#### Supported celltype identifiers 
- __Marker-based identifiers (require HiCAT marker file)__
    1. HiCAT
    1. Garnett
    2. SCINA
    3. scSorter
    4. scType
    5. scCatch

- __Reference-based identifiers (require reference gene expression data with pre-annotation)__ 
    1. MarkerCount_Ref
    1. SingleR
    2. CaSTLe
    3. CHETAH
    4. scmap-cell
    5. scmap-cluster

#### Requirements & Installation
- __Python 3.8 or later__
- __R 4.0.2 or later__
- __Required python packages__: `numpy`, `pandas`, `scikit-learn`, `scipy`, `scikit-network`, __`MarkerCount`__ and __`HiCAT`__ (can be installed using `pip install <package name>`) (see git repo. combio-dku/MarkerCount and combio-dku/HiCAT)
- __Required R packages__: `igraph`, `scater`, `xgboost`, `SingleCellExperiment`, `dplyr`, `stringr`, `preprocessCore`, `Seurat`, `org.Hs.eg.db`, `scuttle`, `SingleR`, `CHETAH`, `scmap`, `SCINA`, `scSorter`, `garnett`, `scCATCH`, `reticulate`
- Once requirements are met, `CTIcollection` can be installed using the following command in R: 

   `devtools::install_github("combio-dku/CTIcollection")`

#### Using the package
See the jupyter notebook provided in this repo. (`CTIcollection_example.ipynb`)

#### Contact
Send email to syoon@dku.edu for any inquiry on the usages.

