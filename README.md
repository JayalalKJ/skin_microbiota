# skin_microbiota

This repository contains the bioinformatics and microbial ecology analysis scripts used for 16S rRNA gene amplicon analysis of Atlantic salmon skin microbiota. The work was conducted as part of my PhD research in the Seafood Science Research Group at UiT – The Arctic University of Norway, under the supervision of Professor Edel O. Elvevoll, Professor Emeritus Bjarne Landfald, and Professor Karl-Erik Eilertsen.

**Diatom biomass as a functional feed ingredient: Effects on skin, gut, fillet quality and salmon lice resistance**

Authors: Hans Chr. Eilertsen, Jayalal K. Jayanthan, Anette Hustad, Dhivya Borra Thiyagarajan, Edel O. Elvevoll, Gunilla K. Eriksen, Jo H. Strømholt, John-Steinar Bergum, Espen Holst Hansen, Karl-Erik Eilertsen, Elisabeth Ytteborg, Gunhild Seljehaug Johansson, Stein Harris Olsen, Gerrit Timmerhaus, Mads Melingen, and Sten Siikavuopio.

## Project overview

The scripts in this repository were used to process and analyse 16S rRNA gene V3–V4 amplicon sequencing data generated from salmon skin microbiota samples. The workflow includes sequence quality control, contig assembly, filtering, chimera removal, taxonomic classification, OTU based community analysis, and downstream microbial ecology statistics and visualisation(phyloseq https://joey711.github.io/phyloseq/ AND microeco https://chiliubio.github.io/microeco_tutorial/).

The bioinformatics processing was based on the **mothur MiSeq SOP** workflow(https://mothur.org/wiki/miseq_sop/), adapted for the dataset used in this study. The mothur MiSeq SOP is described in:

Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD. 2013.  
**Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform.**  
*Applied and Environmental Microbiology* 79(17):5112–5120.

If using this repository or workflow, please cite the mothur MiSeq SOP and indicate the date the SOP was accessed.

## Open Science and Data Availability

This repository is made openly available in accordance with the open science principles of [UiT – The Arctic University of Norway](https://en.uit.no/) and the Data Management Plan established for this PhD project.

UiT supports the principle that research should be **“as open as possible, as closed as necessary”**, including open access to:

* Publications
* Research data
* Software and source code
* Research methods
* Algorithms and bioinformatics pipelines

The scripts and workflows provided in this repository are shared to improve the transparency, reproducibility, and re-use of the bioinformatics and statistical analyses performed in this study.

In accordance with the original PhD project agreement, the research data, code, software, algorithms, and bioinformatics pipelines generated during the project are stored and shared in a reproducible manner following publication.

### Sequencing Data

The 16S rRNA gene V3–V4 amplicon sequencing data generated in this study have been deposited in the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/home) under the following accession numbers:

* **BioProject accession:** [`PRJEB83870`](https://www.ebi.ac.uk/ena/browser/view/PRJEB83870)
* **Study accession:** [`ERP167439`](https://www.ebi.ac.uk/ena/browser/view/ERP167439)

# Atlantic Salmon Skin Mucus Microbiota

**README version:** 1.5

This repository contains the bioinformatics and microbial ecology analysis scripts used for 16S rRNA gene amplicon analysis of Atlantic salmon (*Salmo salar*) skin mucus microbiota.

The work was conducted as part of the PhD research of **Jayalal K. Jayanthan** in the Seafood Science Research Group at UiT – The Arctic University of Norway, under the supervision of Professor Edel O. Elvevoll, Professor Emeritus Bjarne Landfald, and Professor Karl-Erik Eilertsen.

---

## Related Manuscript

**Diatom biomass as a functional feed ingredient: Effects on skin, gut, fillet quality, and salmon lice resistance**

**Authors:**  
Hans Chr. Eilertsen, Jayalal K. Jayanthan, Anette Hustad, Dhivya Borra Thiyagarajan, Edel O. Elvevoll, Gunilla K. Eriksen, Jo H. Strømholt, John-Steinar Bergum, Espen Holst Hansen, Karl-Erik Eilertsen, Elisabeth Ytteborg, Gunhild Seljehaug Johansson, Stein Harris Olsen, Gerrit Timmerhaus, Mads Melingen, and Sten Siikavuopio.

---

## Project Overview

This repository provides the scripts and workflow used to process and analyse 16S rRNA gene V3–V4 amplicon sequencing data generated from Atlantic salmon skin mucus microbiota samples.

The analysis workflow includes:

- Sequence quality control
- Adapter and primer trimming
- Paired-end read merging
- Sequence filtering
- Chimera detection and removal
- Taxonomic classification
- Removal of non-bacterial lineages
- OTU-based microbial community analysis
- Contaminant assessment and removal
- Rarefaction analysis
- Alpha-diversity analysis
- Beta-diversity analysis
- Microbial community composition analysis
- Statistical analysis and visualisation

Downstream microbial ecology analyses were conducted in **R v4.5.2**, primarily using:

- [`phyloseq`](https://joey711.github.io/phyloseq/)
- [`microeco`](https://chiliubio.github.io/microeco_tutorial/)

---

## Experimental Design and Sampling

Atlantic salmon were allocated to nine tanks representing three dietary treatments, with three replicate tanks per diet:

| Diet code | Description |
|----------|-------------|
| `F0` | Control diet, 0% *Porosira glacialis* |
| `F1` | Diet containing 1% *Porosira glacialis* |
| `F3` | Diet containing 3% *Porosira glacialis*; approximately 2.2% inclusion in the later pellet formulation |

Skin mucus microbiota were sampled after the salmon lice challenge on **31 August 2023**, corresponding to sampling point **S6**.

Sampling was performed from one seawater grow-out tank per diet/tank group, giving:

- **3 sampled tanks**
- **12 randomly selected fish per tank**
- **36 skin mucus samples in total**

At the time of sampling:

- Water salinity: **34 ppt**
- Water temperature: **10 °C**

Sampling followed a cross-sectional strategy in accordance with Norwegian Food Safety Authority guidelines.

Because skin microbiota samples were collected from only one tank per diet/tank group, **diet and tank effects cannot be statistically separated**.

---

## Skin Mucus Sample Collection

Skin mucus was collected using sterile forensic swabs:

- **4N6 FLOQSwabs**, Copan, Italy

For each fish, mucus was collected from the lateral side of the fish, above or along the mid-lateral line.

The swab was gently rotated clockwise and counter-clockwise over the skin surface while applying light pressure to collect mucus and associated microbiota.

Three operators performed the sampling, while a fourth person recorded the sampling details.

After collection:

1. Swabs were immediately placed in **96% ethanol**.
2. Samples were kept on ice during transport.
3. Samples were stored at **−20 °C** until DNA extraction.

Six skin mucus samples did not yield successful sequencing data and were excluded from downstream microbiota analyses:

- `SM19`
- `SM61`
- `SM85`
- `SM97`
- `SM103`
- `SM106`

Thus, downstream skin microbiota analyses were performed on successfully sequenced samples only.

---

## Water Sample Collection

Water samples were also collected from each sampled tank.

For each tank:

- **60 mL** of water was collected approximately **5 cm below the surface**
- Sterile syringes were used
- Water was filtered through **0.22 µm Sterivex filter cartridges**
- Filter cartridges were used to collect microbial biomass

Sterivex filters were from:

- Merck Millipore, USA

---

## DNA Extraction and Sequencing

Genomic DNA was extracted using:

- **DNeasy PowerSoil Kit**, Qiagen

DNA extraction was performed according to the manufacturer’s instructions.

DNA quality was checked before samples were sent to BGI for 16S rRNA gene amplicon sequencing.

The V3–V4 region of the bacterial 16S rRNA gene was amplified using the following primers:

| Primer | Sequence |
|--------|----------|
| 338F | `5′-ACTCCTACGGGAGGCAGCAG-3′` |
| 806R | `5′-GGACTACHVGGGTWTCTAAT-3′` |

Adapter-linked fusion primers were used for library preparation.

PCR products were:

1. Checked for quality
2. Purified using **AMPure XP beads**
3. Measured and size-checked using an **Agilent 2100 Bioanalyzer**
4. Sequenced on the **DNBSEQ-G400 platform**, BGI

Sequencing generated:

- **2 × 300 bp paired-end reads**

---

## Bioinformatics Workflow

Sequence processing was based on the [mothur MiSeq SOP](https://mothur.org/wiki/miseq_sop/), adapted for this dataset.

The workflow included the following main steps:

1. Removal of sequencing adapters and primers using **Cutadapt**
2. Quality filtering of raw reads
3. Paired-end read merging using **FLASH**
4. Generation of clean paired-end reads
5. Processing using the mothur MiSeq workflow
6. Alignment against a tailored SILVA V3–V4 reference alignment
7. Chimera detection and removal
8. Taxonomic classification
9. Removal of non-bacterial lineages
10. OTU filtering based on prevalence and abundance
11. Downstream microbial ecology analysis in R

Paired-end reads were merged using FLASH with:

- Minimum overlap: **15 bp**
- Maximum mismatch ratio: **≤ 0.1**

A tailored 16S rRNA V3–V4 reference alignment was created using:

- **SILVA database release 138.2**
- Accessed: **11 July 2024**
- Coordinate reference: *Escherichia coli* 16S ribosomal RNA complete sequence, accession `J01859.1`

The tailored reference file was named:

```text
silva.v3.v4.fasta
## Funding and acknowledgements


This research was funded by Innovation Norway grant **2021/312146**, Norwegian Research Council grant **NFR 328654**, Skattefunn grant **15732**, and EU grant **IGNITION 101084651**.

Jayalal K. Jayanthan’s PhD research was funded by UiT – The Arctic University of Norway through the SECURE project, Cristin grant ID 2061344.

## Conflict of interest

The authors declare no conflict of interest.


##Contact

Jayalal K. Jayanthan
PhD Candidate, 2021–2026
Seafood Science Research Group
The Norwegian College of Fishery Science
Faculty of Biosciences, Fisheries and Economics
UiT – The Arctic University of Norway
Muninbakken 21, 9019 Tromsø, Norway
jayalal.p.kalathil@uit.no

## License

This repository is distributed under the MIT License.


