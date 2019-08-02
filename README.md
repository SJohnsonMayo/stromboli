# Stromboli: 

Statistically Robust Methods for Microbiome Analysis

# Installation:

Stromboli can be installed using the `devtools` package as follows:

```
devtools::install_github("SJohnsonMayo/stromboli")
```

# Motivation
The human microbiome, the collection of microbes associated with the human body, has received tremendous attention recently due to its important role in human health and disease.   Fueled by next-generation sequencing technologies, the human microbiome can now be studied by direct sequencing, generating a large amount of sequencing data.  However, due to its complex data characteristics, analysis of the microbiome sequencing data raises many statistical challenges. Here we introduce a new statistical pipeline MicrobiomeGPS for efficient and robust analysis of the microbiome data. MicrobiomeGPS is a web-based application, which takes the OTU count table, phylogenetic tree and sample mapping file as the input and generate s HTML reports of different levels of analysis including exploratory data analysis, alpha- and beta- diversity analysis, and differential abundance analysis at both the taxonomic and functional levels.  It also includes modules for predictive modeling, community subtype analysis and OTU network analysis.  Only the most robust statistical methods are included in our pipeline.  MicrobiomeGPS integrates covariate adjustment in each step and thus can address the confounding effects due to other technical, biological and clinical variables.   It could also analyze repeated measurement data by taking into account the within-subject correlation structure.  We believe that MicrobiomeGPS will be a useful tool for microbiome investigators to generate hypotheses, make discoveries and potentially translate the findings into clinical practice.


# Using stromboli
For a brief introduction to the functionality provided by stromboli, please see the `stromboli-overview.Rmd` vignette. 
