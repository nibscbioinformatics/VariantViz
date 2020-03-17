---
title: "README"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## SNP_inhouse_downstream script

This tool serves to describe a global overview of the variants found in a genomic sample.

The script R/SNP_inhouse_downstream.R has the following inputs:

  - VCF file path (annotated with the ClinVar database)
  
  - Sample name

And it has the following output:

  - A .pdf file (named "<sample_name>_SNP_inhouse_circos_plot.pdf") with a circular bar plot detailing the allele frequency of all mutations annotated by ClinVar. They are grouped into different clinical tiers
  
  - A .tsv file (named "<sample_name>_SNP_inhouse_table.tsv") with the previous described information, in a table

### Output image

By calling the tool this way:
Rscript /path/to/annotated_file_with_clinvar.vcf Example

We have the following output:

![](/home/AD/praposo/CSP/files/Example.circos.plot.png)
