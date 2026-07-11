# Atlantic Salmon Skin Mucus Microbiota

> Reproducible bioinformatics and microbial ecology workflows for 16S rRNA gene amplicon analysis of Atlantic salmon (*Salmo salar*) skin mucus microbiota.

## Repository Scope

This repository contains the R scripts and supporting workflow used to process, quality-check, analyse, and visualise 16S rRNA gene amplicon data from Atlantic salmon skin mucus microbiota samples.

The repository supports the following analyses:

- Import of mothur-derived community data into R
- Assessment of extraction blanks and negative controls
- Identification and removal of potential contaminants
- Rarefaction and sequencing-depth assessment
- Alpha-diversity analysis
- Beta-diversity, ordination, PERMANOVA, and dispersion analysis
- Taxonomic composition analysis
- Comparison of skin mucus and tank-water microbial communities

The code was developed as part of the PhD research of **Jayalal K. Jayanthan** in the Seafood Science Research Group at UiT – The Arctic University of Norway.

## Important Interpretation of the Study Design

> **Study-design limitation**
>
> The feeding phase included three replicate tanks per dietary treatment. However, skin mucus microbiota samples collected at sampling point S6 were obtained from only one seawater tank per diet/tank group. Consequently, diet and tank effects are fully confounded at this sampling point.
>
> Individual fish represent biological subsamples within tanks and must not be interpreted as independent tank-level replicates of the dietary treatments. Differences among F0, F1, and F3 should therefore be described as **descriptive or exploratory diet/tank-associated patterns**, rather than as confirmatory causal effects of diet.

This limitation should be retained in all analyses, figures, captions, manuscripts, presentations, and secondary uses of the dataset.

## Associated Manuscript

**Diatom biomass as a functional feed ingredient: Effects on skin, gut, fillet quality, and salmon lice resistance**

**Authors:**  
Hans Chr. Eilertsen, Jayalal K. Jayanthan, Anette Hustad, Dhivya Borra Thiyagarajan, Edel O. Elvevoll, Gunilla K. Eriksen, Jo H. Strømholt, John-Steinar Bergum, Espen Holst Hansen, Karl-Erik Eilertsen, Elisabeth Ytteborg, Gunhild Seljehaug Johansson, Stein Harris Olsen, Gerrit Timmerhaus, Mads Melingen, and Sten Siikavuopio.

Add the final journal citation and DOI here after publication.

## Experimental Design and Sampling

Atlantic salmon were allocated to nine tanks during the replicated feeding phase, comprising three dietary treatments with three tanks per treatment.

| Diet code | Dietary treatment |
|---|---|
| `F0` | Control diet without *Porosira glacialis* |
| `F1` | Diet containing 1% *Porosira glacialis* |
| `F3` | Nominal 3% *Porosira glacialis* treatment; approximately 2.2% inclusion in the later pellet formulation |

Skin mucus microbiota samples were collected after the salmon lice challenge on **31 August 2023**, corresponding to sampling point **S6**.

At S6, sampling was conducted from one seawater tank per diet/tank group:

| Group | Samples collected | Samples retained after sequencing QC |
|---|---:|---:|
| `F0` | 12 | 11 |
| `F1` | 12 | 11 |
| `F3` | 12 | 8 |
| **Total** | **36** | **30** |

Environmental conditions at sampling were approximately:

- Salinity: **34 ppt**
- Water temperature: **10 °C**

The following samples did not yield usable sequencing data and were excluded from downstream microbiota analyses:

`SM19`, `SM61`, `SM85`, `SM97`, `SM103`, and `SM106`.

## Sample Collection

### Skin Mucus

Skin mucus was collected using sterile **4N6 FLOQSwabs** (Copan, Italy). For each fish, the lateral skin surface above or along the mid-lateral line was swabbed using gentle bidirectional rotation and light pressure.

After collection:

1. Swabs were placed immediately in 96% ethanol.
2. Samples were kept on ice during transport.
3. Samples were stored at −20 °C until DNA extraction.

### Tank Water

Tank-water samples were collected from each sampled tank. Approximately 60 mL of water was collected about 5 cm below the water surface using sterile syringes and filtered through 0.22 µm Sterivex filter cartridges (Merck Millipore, USA) to retain microbial biomass.

## DNA Extraction and Sequencing

Genomic DNA was extracted using the **DNeasy PowerSoil Kit** (Qiagen) according to the manufacturer’s instructions.

The bacterial 16S rRNA gene V3–V4 region was amplified using the following primers:

| Primer | Sequence |
|---|---|
| 338F | `5′-ACTCCTACGGGAGGCAGCAG-3′` |
| 806R | `5′-GGACTACHVGGGTWTCTAAT-3′` |

Adapter-linked fusion primers were used for library preparation. PCR products were purified using AMPure XP beads, quantified, and size-checked using an Agilent 2100 Bioanalyzer.

Sequencing was performed by BGI on the **DNBSEQ-G400** platform using **2 × 300 bp paired-end sequencing**.

> Before public release, verify the extraction-kit name, sequencing provider, sequencing platform, read length, primer sequences, and targeted variable region against the final laboratory and sequencing reports.

## Bioinformatics Workflow

Sequence processing was based on the [mothur MiSeq SOP](https://mothur.org/wiki/miseq_sop/), adapted to the characteristics of this dataset.

The main workflow comprised:

1. Adapter and primer removal using Cutadapt
2. Quality filtering of raw reads
3. Paired-end read merging using FLASH
4. Processing of cleaned reads using mothur
5. Alignment against a tailored SILVA V3–V4 reference alignment
6. Chimera detection and removal
7. Taxonomic classification
8. Removal of non-bacterial and other unwanted lineages
9. OTU construction and filtering
10. Assessment of extraction blanks and negative controls
11. Removal of potential contaminants
12. Rarefaction and alpha-diversity analysis
13. Beta-diversity, ordination, PERMANOVA, and PERMDISP analysis
14. Taxonomic composition analysis and visualisation in R

Paired-end reads were merged using FLASH with the following parameters:

- Minimum overlap: **15 bp**
- Maximum mismatch ratio: **0.1**

A tailored reference alignment was generated from:

- **SILVA release:** 138.2
- **Database access date:** 11 July 2024
- **Coordinate reference:** *Escherichia coli* 16S rRNA sequence, accession `J01859.1`
- **Reference filename:** `silva.v3.v4.fasta`

Downstream microbial ecology analyses were conducted primarily using:

- [`phyloseq`](https://joey711.github.io/phyloseq/)
- [`microeco`](https://chiliubio.github.io/microeco_tutorial/)
- [`vegan`](https://cran.r-project.org/package=vegan)

## Software and Reproducibility

The analysis was developed in **R 4.5.2**.

For a reproducible public release, the repository should include or report:

- mothur version
- R version
- R package versions
- Cutadapt version
- FLASH version
- SILVA release and access date
- Random seeds used in permutation-based analyses
- Filtering, prevalence, and abundance thresholds
- Rarefaction depth
- Operating system or container information

Recommended reproducibility files:

- `renv.lock`
- `sessionInfo.txt`
- `CITATION.cff`
- A tagged GitHub release corresponding to the submitted or published manuscript
- A Zenodo archive and DOI for the final release

## Repository Structure

| File | Purpose |
|---|---|
| `00_project_setup.R` | Defines project paths, loads packages, creates output directories, and sets shared analysis parameters |
| `01_import_mothur_to_microeco.R` | Imports mothur output and metadata into R and creates the analysis object |
| `02A_Blanks_visual_proofpack.R` | Produces diagnostic summaries and figures for blanks and negative controls |
| `02B_Negative_control_samples_Composition.R` | Examines taxonomic composition in negative-control samples |
| `02C_contaminant_removal_from_blanks.R` | Identifies and removes potential contaminants associated with blanks |
| `03_rarefaction_Skin_sample_for_alpha.R` | Evaluates sequencing depth and prepares skin samples for alpha-diversity analysis |
| `04_alpha_diversity_Skin_samples.R` | Calculates, tests, and visualises alpha-diversity metrics |
| `05_Beta_diversity_Skin_samples.R` | Performs beta-diversity, ordination, PERMANOVA, and dispersion analyses |
| `07_composition_analysis_Skin_water_samples.R` | Analyses and visualises skin mucus and water-community composition |
| `Run_all_pipeline.R` | Runs the principal scripts in the intended order |
| `README.md` | Describes the study, workflow, data availability, and repository use |
| `LICENSE` | Contains the software license |

The filename must be `README.md`, not `REAME.md`, so that GitHub renders it automatically on the repository landing page.

## Input Data

Raw sequencing reads are not stored directly in this GitHub repository. They are available from the European Nucleotide Archive under the accessions listed below.

Local input files required by the scripts may include:

- mothur OTU or shared files
- taxonomy files
- sample metadata
- blank and negative-control metadata
- environmental or treatment metadata
- reference taxonomy or alignment files

Expected paths, filenames, and metadata fields should be defined centrally in `00_project_setup.R`.

Sensitive, restricted, or personally identifiable information must not be committed to the repository.

## Running the Analysis

Clone the repository:

```bash
git clone <https://github.com/JayalalKJ/skin_microbiota.git>
cd skin_microbiota
```

Before running the workflow:

1. Review and update file paths in `00_project_setup.R`.
2. Confirm that all required input files are available.
3. Install the required R packages and command-line dependencies.
4. Confirm that sample identifiers and metadata columns match those expected by the scripts.
5. Set or record random seeds for analyses involving permutation or stochastic procedures.

Run the complete workflow from the command line:

```bash
Rscript Run_all_pipeline.R
```

Alternatively, run the numbered scripts sequentially from an R session.

## Expected Outputs

Depending on the configuration of the scripts, the workflow generates:

- Quality-control summaries
- Blank and negative-control diagnostic plots
- Contaminant tables
- Rarefaction curves
- Alpha-diversity tables and figures
- Beta-diversity distance matrices
- Ordination plots
- PERMANOVA and PERMDISP results
- Taxonomic abundance tables
- Community-composition figures
- Publication-ready summary tables and graphics

Generated outputs should be written to dedicated results directories and should not overwrite raw or intermediate input data.

## Data Availability

The 16S rRNA gene V3–V4 amplicon sequencing data generated in this study are deposited in the [European Nucleotide Archive](https://www.ebi.ac.uk/ena/browser/home):

- **BioProject accession:** [`PRJEB83870`](https://www.ebi.ac.uk/ena/browser/view/PRJEB83870)
- **Study accession:** [`ERP167439`](https://www.ebi.ac.uk/ena/browser/view/ERP167439)

## Open Science Statement

This repository is made openly available in accordance with the open science principles of [UiT – The Arctic University of Norway](https://en.uit.no/) and the Data Management Plan established for this PhD project.

The scripts are shared to support transparency, reproducibility, scrutiny, and responsible re-use of the bioinformatics and statistical workflow. The guiding principle is that research outputs should be **“as open as possible, as closed as necessary.”**

## Citation

When using this repository, please cite:

1. The associated research article once its final citation and DOI are available.
2. The archived release of this repository once a DOI has been assigned.
3. The principal software and workflow publications relevant to the reused analysis.

### Core Workflow References

Kozich, J. J., Westcott, S. L., Baxter, N. T., Highlander, S. K., and Schloss, P. D. (2013). Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform. *Applied and Environmental Microbiology*, 79(17), 5112–5120. <https://doi.org/10.1128/AEM.01043-13>

McMurdie, P. J., and Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. *PLOS ONE*, 8(4), e61217. <https://doi.org/10.1371/journal.pone.0061217>

Liu, C., Cui, Y., Li, X., and Yao, M. (2021). microeco: An R package for data mining in microbial community ecology. *FEMS Microbiology Ecology*, 97(2), fiaa255. <https://doi.org/10.1093/femsec/fiaa255>

Report the mothur version and the date on which the online MiSeq SOP was accessed.

## Funding and Acknowledgements

This research was supported by:

- Innovation Norway grant **2021/312146**
- Research Council of Norway grant **NFR 328654**
- SkatteFUNN grant **15732**
- European Union IGNITION grant **101084651**

Jayalal K. Jayanthan’s PhD research was funded by UiT – The Arctic University of Norway through the SECURE project, Cristin project ID **2061344**.

## Conflict of Interest

The authors declare no conflicts of interest.

## Contact

**Jayalal K. Jayanthan**  
PhD Candidate, 2021–2026  
Seafood Science Research Group  
The Norwegian College of Fishery Science  
Faculty of Biosciences, Fisheries and Economics  
UiT – The Arctic University of Norway  
Muninbakken 21  
9019 Tromsø, Norway  

Email: [jayalal.p.kalathil@uit.no](mailto:jayalal.p.kalathil@uit.no)

Questions, reproducibility problems, and bug reports should be submitted through the repository’s GitHub Issues page.

## License

The source code in this repository is distributed under the [MIT License](LICENSE).

The sequencing data are distributed through the European Nucleotide Archive and remain subject to the metadata, access conditions, and reuse requirements associated with the deposited study.
