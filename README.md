# medline-discoveries

## Overview

This repository contains the companion code for the paper "Mining impactful discoveries from the biomedical literature" **TODO link**.

This code implements the method described in the paper to extract "impactful discoveries" from [Medline](https://www.nlm.nih.gov/medline/index.html), using the MeSH descriptors to represent the concepts the articles. "Impactful discoveries" are the relations which have an exceptionally strong surge in the data at some point in time. The output is a list of relations together with their surge year (i.e. a list of discoveries). 

The main code is written in R and the documentation is provided as R Markdown documents for the sake of reproducibility. There are three parts in the documentation:

- The [data collection](https://erwanm.github.io/medline-discoveries/data-collection) process describes how the input data was created.
- The ["Detecting surges" main documentation](https://erwanm.github.io/medline-discoveries/1-detecting-surges.html) describes how to generate the output data (i.e. detecting surges) using the provided R implementation.
- The ["Analysis" part](https://erwanm.github.io/medline-discoveries/2-analysis.html) proposes some additional analysis of the results, including how to generate the graphs and tables found in the paper.  

The data can be downloaded from [https://zenodo.org/record/5888572](https://zenodo.org/record/5888572):

- The input data conists of the frequency by year for the Medline MeSH descriptors concepts and relations (individual and joint frequency)
- The output data is the extracted list of "impactful discoveries".

An [exploration tool](https://brainmend.adaptcentre.ie/) is also provided to visualize the output discoveries (code included in this repository).


## Acknowledgements

This work was conducted with the financial support of Irish Health Research Board (HRB) BRAINMend project (grant number HRB-JPI-JPND-2017-1) and with the support of the Science Foundation Ireland Research Centre ADAPT (grant number 13/RC/2106\_P2).

The interactive interfaces are implemented with [Shiny](https://shiny.rstudio.com) and the documentation uses [RMarkdown](https://rmarkdown.rstudio.com/). Graphs in the papers and documentation are produced with [ggplot2](https://ggplot2.tidyverse.org/), and the 'Medline discoveries' code heavily relies on the [data.table](https://rdatatable.gitlab.io/data.table/) library.
