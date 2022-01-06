# medline-discoveries: data collection

## Introduction

This directory contains the scripts and explanations (below) related to the data collection process, i.e. how the input data was generated.

This part is not required to reproduce the experiments in the paper. It is also possible to obtain the input data using some other sources and methods, as long as the input format is respected (see [details in the main documentation](https://erwanm.github.io/medline-discoveries/1-detecting-surges.html)).

## Extracting the MeSH descriptors by paper from Medline 

This part is a bit complex and relies on several other repositories:

- The [KD fork](https://github.com/erwanm/knowledgediscovery) can be used to extract the initial list
    - output `mesh-descriptors-by-pmid.tsv`
- The [KD Tools](https://github.com/erwanm/kd-data-tools) should be used to deduplicate the list
    - output `mesh-descriptors-by-pmid.deduplicated.tsv`
- Finally the [TDC Tools](https://github.com/erwanm/tdc-tools) can be used to obtain the DCM (Doc-Concept Matrix) format

Please follow the instructions provided in the [TDC Tools documentation](https://erwanm.github.io/tdc-tools/mesh-descriptors-by-pmid/). The final output should be the following two directories:

- `medline-mesh-descriptors/concept-freq`
- `medline-mesh-descriptors/joint-freq`

## 