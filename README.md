### An easy-to-use wrapper function to run celltype identifiers in R

- The package provides an easy-to-use wrapper function to run various celltype identifiers in R.
- It was used to obtain the results in the following papers.

    - "Hierarchical cell-type identifier accurately distinguishes immune-cell subtypes enabling precise profiling of tissue microenvironment with single-cell RNA-sequencing" Briefings in Bioinformatics. https://doi.org/10.1093/bib/bbad006, https://doi.org/10.1101/2022.07.27.501701.
    - "MarkerCount: A stable, count-based cell type identifier for single cell RNA-Seq experiments" Computational and Structural Biotechnology Journal, June 2022. https://doi.org/10.1016/j.csbj.2022.06.010.

#### Supported celltype identifiers 
- __Marker-based identifiers (require HiCAT marker file)__
    1. HiCAT        (https://github.com/combio-dku/HiCAT)
    1. Garnett      (https://cole-trapnell-lab.github.io/garnett/)
    2. SCINA        (https://github.com/jcao89757/SCINA)
    3. scSorter     (https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html)
    4. scType       (https://github.com/IanevskiAleksandr/sc-type)
    5. scCatch      (https://github.com/ZJUFanLab/scCATCH)

- __Reference-based identifiers (require reference gene expression data with pre-annotation)__ 
    1. MarkerCount_Ref  (https://github.com/combio-dku/MarkerCount)
    1. SingleR          (https://github.com/dviraran/SingleR)
    2. CaSTLe           (https://github.com/yuvallb/CaSTLe)
    3. CHETAH           (https://github.com/jdekanter/CHETAH)
    4. scmap-cell       (https://github.com/hemberg-lab/scmap)
    5. scmap-cluster    (https://github.com/hemberg-lab/scmap)

#### Requirements & Installation
- __Python 3.8 or later__
- __R 4.0.2 or later__
- __Required python packages__: `numpy`, `pandas`, `scikit-learn`, `scipy`, `scikit-network`, __`MarkerCount`__ and __`HiCAT`__ (can be installed using `pip install <package name>`) (see git repo. combio-dku/MarkerCount and combio-dku/HiCAT)
- __Required R packages__: `igraph`, `scater`, `xgboost`, `SingleCellExperiment`, `dplyr`, `stringr`, `preprocessCore`, `Seurat`, `org.Hs.eg.db`, `scuttle`, `SingleR`, `CHETAH`, `scmap`, `SCINA`, `scSorter`, `garnett`, `scCATCH`, `reticulate`
- Once requirements are met, `CTIcollection` can be installed using the following command in R: 

   1. run R
   2. run `devtools::install_github("combio-dku/RCTIcollection")` in R

#### Using the package
See the jupyter notebook provided in this repo. (`CTIcollection_example.ipynb`)

#### Contact
Send email to syoon@dku.edu for any inquiry on the usages.

