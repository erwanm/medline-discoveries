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

## Preparing the input data

### Requirements

- This part requires some scripts from the [TDC Tools](https://github.com/erwanm/tdc-tools) repository. Replace `<my-path-to>` with the path in which `tdc-tools` is located in the command below:

```
export PATH=$PATH:<my-path-to>/tdc-tools/code/
```

- The UMLS metathesaurus data is required in order to map MeSH concepts ids to their corresponding term and group.
    - Access requires a free UTC account.
    - On the [UMLS download page](https://www.nlm.nih.gov/research/umls/licensedcontent/umlsknowledgesources.html) select "UMLS Metathesaurus Files" (no need for the "full release").
    - The [UMLS "semantic groups" file](https://lhncbc.nlm.nih.gov/ii/tools/MetaMap/documentation/SemanticTypesAndGroups.html) is also required. The 2021 version is provided for convenience in the `data/` directory.

- The generation of the ND subset requires a list of target concepts. The ND targets used for the paper experiments is provided in the `data/` directory.
    - The process to generate this list of targets is described in the [TDC Tools documentation](https://erwanm.github.io/tdc-tools/ND-use-case/#list-of-target-concepts).
- The process below requires 32G RAM.

### Generating the "full" variant of the input data

Replace `<my-path-to>` with the appropriate path in the command below:

```
data-collection/scripts/prepare-dataset.sh <my-path-to>/concept-freq/ <my-path-to>/joint-freq/ data/input
```

This process may take 15 to 30 mn. It should generate the files [described in the main documentation](https://erwanm.github.io/medline-discoveries/1-detecting-surges.html#input-data-format) in `data/input`: `indiv.full.min100`  `indiv.full.total`  `joint.full.min100`, as well as the `static` directory.


### Generating the ND subset of the input data

Replace `<path to umls>` with the path where the UMLS data is located (the one containing the `META` directory) in the command below:

```
data-collection/scripts/prepare-dataset-ND.sh -m data/input/ data/ND.mesh.targets <path to umls> data/SemGroups.txt
```

This process takes a few minutes. It creates the files `indiv.ND.min100` and `joint.ND.min100` in `data/input` as well as `indiv.min100.ND` `joint.min100.ND` in `data/input/static`.
