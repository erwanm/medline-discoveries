# medline-discoveries

## Overview

This repository contains the companion code for the paper "Mining impactful discoveries from the biomedical literature" **TODO link**.

This code implements the method described in the paper to extract "impactful discoveries" from [Medline](https://www.nlm.nih.gov/medline/index.html), using the MeSH descriptors to represent the concepts the articles. "Impactful discoveries" are the relations which have an exceptionally strong surge in the data at some point in time. The output is a list of relations together with their surge year (i.e. a list of discoveries). 

The main code is written in R and the documentation is provided as R Markdown documents for the sake of reproducibility. There are three parts in the documentation:

- The [data collection](https://erwanm.github.io/medline-discoveries/data-collection) process describes how the input data was created.
- The ["Detecting surges" main documentation](https://erwanm.github.io/medline-discoveries/1-detecting-surges.html) describes how to generate the output data (i.e. detecting surges) using the provided R implementation.
- The ["Analysis" part](https://erwanm.github.io/medline-discoveries/2-analysis.html) proposes some additional analysis of the results, including how to generate the graphs and tables found in the paper.  

The data can be downloaded from  **TODO link**:

- The input data conists of the frequency by year for the Medline MeSH descriptors concepts and relations (individual and joint frequency)
- The output data is the extracted list of "impactful discoveries".

