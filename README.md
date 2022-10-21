# CRC ATLAS

README Contents:

1. [Datasets Collected](https://github.com/saracardoso/CRC_ATLAS#1-datasets-collected)

2. [How this project is organized](https://github.com/saracardoso/CRC_ATLAS#2-how-this-project-is-organized)

3. [Where atlas files are available](https://github.com/saracardoso/CRC_ATLAS#3-where-atlas-files-are-available)


This project contains the scripts developed to construct a CRC atlas of tumour and normal matched tissue of patients. Also, samples from heatlhy donors were inserted. 

This atlas has a total of 163 810 cells, separated into 51 044 T-cells, 47 462 epithelial cells, 30 187 stromal
cells, 17 674 B-cells, and 17 443 myeloid cells. The number of cells for each cell-subtype identified is
summarised in the following table:

| A. Cell Groups          | B. Cell-types                        |                     | Number of cells | Total number of cells in A. |
|-------------------------|--------------------------------------|---------------------|-----------------|-----------------------------|
| Cancer cells            |                                      |                     |                 | 28 208                      |
| Normal Epithelial cells | Progenitor cells                     |                     | 2 818           | 19 254                      |
|                         | Secretory progenitors                |                     | 2 094           |                             |
|                         | Transit-amplifying (TA) cells        |                     | 3 862           |                             |
|                         | Tuft cells                           |                     | 191             |                             |
|                         | Enteroendocrine cells                |                     | 74              |                             |
|                         | Goblet cells                         | Not immature        | 391             |                             |
|                         |                                      | Immature            | 1 505           |                             |
|                         | Enterocyte/colonocyte cells          | Immature            | 3 926           |                             |
|                         |                                      | BEST4+              | 536             |                             |
|                         |                                      | Other               | 3 201           |                             |
|                         | Paneth-like cells                    |                     | 656             |                             |
| Stromal cells           | Fibroblasts                          |                     | 14 007          | 30 187                      |
|                         | CAFS                                 |                     | 4 914           |                             |
|                         | Myofibroblasts                       |                     | 693             |                             |
|                         | Pericytes                            |                     | 2 014           |                             |
|                         | Enteric glia cells                   |                     | 1 649           |                             |
|                         | VSMCs                                |                     | 791             |                             |
|                         | Endothelial cells                    | tip-like vascular   | 3 934           |                             |
|                         |                                      | stalk-like vascular | 1 792           |                             |
|                         |                                      | lymphatic           | 393             |                             |
| Myeloid cells           | Anti-inflammatory macro/mono lineage |                     | 4 963           | 17 443                      |
|                         |                                      | SPP1+               | 4 180           |                             |
|                         | Pro-inflammatory macro/mono lineage  |                     | 2 548           |                             |
|                         |                                      | FCN1+               | 1 639           |                             |
|                         | conventional DCs                     |                     | 2 072           |                             |
|                         | plasmacytoid DCs                     |                     | 335             |                             |
|                         | Mast cells                           |                     | 1 222           |                             |
|                         | Unkown                               |                     | 484             |                             |
| B-cells                 | Naive                                |                     | 2 872           | 17 674                      |
|                         | Memory                               |                     | 7 323           |                             |
|                         | Plasma cells                         | IgA+                | 5 624           |                             |
|                         |                                      | IgG+                | 685             |                             |
|                         | Proliferative                        |                     | 1 170           |                             |
| T-cells                 | CD4+                                 | Poliferative        | 377             | 51 044                      |
|                         |                                      | Naive               | 8 484           |                             |
|                         |                                      | Memory              | 13 484          |                             |
|                         |                                      | Regulatory          | 6 801           |                             |
|                         |                                      | Follicular          | 1 156           |                             |
|                         |                                      | IL22+               | 201             |                             |
|                         |                                      | IL17+               | 755             |                             |
|                         | CD8+                                 | Poliferative        | 439             |                             |
|                         |                                      | Naive               | 414             |                             |
|                         |                                      | Memory              | 6 534           |                             |
|                         |                                      | CXCL13+             | 1 977           |                             |
|                         |                                      | Cytotoxic           | 1 105           |                             |
|                         | Unconventional                       | $\gamma\delta$      | 2 317           |                             |
|                         |                                      | NKT                 | 1 759           |                             |
|                         |                                      | NK                  | 1 382           |                             |
|                         |                                      | CD8$\alpha\alpha$   | 2 488           |                             |
|                         |                                      | LTi cells           | 327             |                             |
|                         | Double-Negative                      |                     | 1 044           |                             |
| **TOTAL**               |                                      |                     |                 | **163 810**                 |



## 1. Datasets Collected

Raw counts from four publicly available datasets were used to construct the atlas. Three correspond to CRC studies, while one (*Colon_smillie*) has data from colon of patients with ulcerative colitis and healthy individuals. As such, the samples related to the healthy individuals were used. The next table shows an overview of the datasets collected.

|                          | *CRC_Qian*                        | *GSE132465*                       | *GSE144735*                       | *Colon_smillie*                    |
|--------------------------|-----------------------------------|-----------------------------------|-----------------------------------|------------------------------------|
| **Nº Patients**          | 7                                 | 23                                | 6                                 | 12                                 |
| **Nº Samples**           | 21                                | 33                                | 18                                | 48                                 |
| **Nº Cells (before QC)** | 44 684                            | 63 689                            | 27 414                            | 51 705                             |
| **Nº Cells (after QC)**  | 29 793                            | 57 804                            | 24 510                            | 50 811                             |
| **Technology**           | 3’ 10xGenomics                    | 3’ 10xGenomics                    | 3’ 10xGenomics                    | 3’ 10xGenomics                     |
| **Link to Article**      | [Link](https://doi.org/10.1038/s41422-020-0355-0) | [Link](https://doi.org/10.1038/s41588-020-0636-z) | [Link](https://doi.org/10.1038/s41588-020-0636-z) | [Link](https://doi.org/10.1016/j.cell.2019.06.029) |

Studies *CRC_Qian* and *GSE144735* have tumour (border and core) and normal matched mucosa samples for each patient. For study *GSE132465*, 10 of the 23 patients have tumour and normal matched mucosa samples, while the other 13 only have tumour samples. Regarding *Colon_smillie* dataset, each donor has a sample extracted from two different locations of the colon.



## 2. How this project is organized

The scripts developed to construct the CRC atlas are organized in this project as follows:

- __*0_initial_data_evaluation*__: quality control of each of the datasets used to construct the atlas;

- __*1_merge*__: merge and integration of the different datasets together into one;

- __*2_annotation*__: annotation of the cells into the various cell-types

- __*3_generate_atlas_files*__: script to generate the final atlas files


## 3. Where atlas files are available

The following files are available:

| File                                                                                                    | Format     | Size     | Observations                                                        |
|---------------------------------------------------------------------------------------------------------|------------|----------|---------------------------------------------------------------------|
| [CRCatlas.h5Seurat](https://drive.google.com/file/d/129dNBGHP9bo30P8Kj-m0jgMqe3e9-mgD/view?usp=sharing) | *h5Seurat* | 3.82GB   | **The full dataset**. It includes all annotations and metadata, with both *RNA* and *integrated* assays. |
| [CRCatlas.h5Seurat](https://drive.google.com/file/d/1C1kzMzEB2Txj-nuNhGMEcI_gnI6C3VxA/view?usp=sharing) | *h5Seurat* | 1.52GB   | It includes all annotations and metadata, but only the *RNA* assay. |
| [CRCatlas_PH.h5ad](https://drive.google.com/file/d/18aYKlknYvKTqX9WaCW1uQnMDUlXYAzk8/view?usp=sharing)  | *h5ad*     | 1.52GB   | It includes all annotations and metadata, but only the *RNA* assay. |
| [metadata.csv](https://drive.google.com/file/d/1ZxCG-DHYBJKzfMttvChYsKIzYqUm7NB4/view?usp=sharing)      | *CSV*      | 110.2 MB | Only the metadata                                                   |




<!-- ## 4. How to reference this atlas -->

