# FRAMEtags
This package automatically gates and quantifies flow cytometry or microscopy data captured from live cells barcoded with FRAME-tags (FTs) using openCyto pipeline.

***
# Installation
1. Download and install R v3.3.3 at https://cloud.r-project.org/
  
2. Open R
  
3. Paste and execute the following commands to install bioconductor (update packages if necessary)
```{r}  
source("https://bioconductor.org/biocLite.R")
biocLite()
```

4. Then the following commands to install FRAMEtags
```{r}  
library(devtools)
install_github("jmiguelj/FRAMEtags")
library("FRAMEtags")
```

***
# Required input files
(see "sample data" for example input folder/files)

You will need to prepare a folder on your computer containing the following files:

1. .fcs files to be analyzed
    + with at least 2 fluorescence, 1 forward scatter and 1 side scatter parameter
  
2. "Tag-Strain.csv" 
    + containing a list of the FTs present in the sample and their corresponding user-specified identities or strain names
  
3. "Tube-Treatment.csv" 
    + containing a list of the "Tube Name" or "Well IDs" and their corresponding user-specified sample treatment or experimental condition, each row corresponds to one input .fcs file and the "treatment" column is used mainly for labeling output statistics

***
# How to run
Into the R command line paste and execute the following command:
```{r} 
getFRAMEtagData()
```

then follow onscreen directions.

***
# Output files
In main input folder you should see:

1. "gating_template.csv""
    + can be modified and getFRAMEtagData() rerun to test modifications
  

2. "files-analyzed_log.csv"
    + includes list of .fcs files succesfully analyzed, getFRAMEtagData() will ignore these files
    + to reanalyze specific files, remove them from the log, or delete log entirely
  

3. "output" folder containing:

    i. "counts.csv"; final event counts for each FT and total parent population "nonMaxF"
  
    ii. "FT-gate.png"; visual overview of final FT gates for each .fcs
        + use this to determine any potential anomalous gating
  

    iii. "full-stats.csv"; event counts for every gate leading to and including the final FT gates
  
    iv. "percent.csv"; relative population fractions of each FT
  
  Leading number indicates processing round as listed in the log file.
  Output folder will also include detailed gating plots if user chooses to output them.
  
***  
  
  

