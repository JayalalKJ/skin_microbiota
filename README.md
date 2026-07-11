# skin_microbiota

This repository contains the bioinformatics and microbial ecology analysis scripts used for 16S rRNA gene amplicon analysis of Atlantic salmon skin microbiota. The work was conducted as part of my PhD research in the Seafood Science Research Group at UiT – The Arctic University of Norway, under the supervision of Professor Edel O. Elvevoll, Professor Emeritus Bjarne Landfald, and Professor Karl-Erik Eilertsen.

**Diatom biomass as a functional feed ingredient: Effects on skin, gut, fillet quality and salmon lice resistance**

Authors: Hans Chr. Eilertsen, Jayalal K. Jayanthan, Anette Hustad, Dhivya Borra Thiyagarajan, Edel O. Elvevoll, Gunilla K. Eriksen, Jo H. Strømholt, John-Steinar Bergum, Espen Holst Hansen, Karl-Erik Eilertsen, Elisabeth Ytteborg, Gunhild Seljehaug Johansson, Stein Harris Olsen, Gerrit Timmerhaus, Mads Melingen, and Sten Siikavuopio.

## Project overview

The scripts in this repository were used to process and analyse 16S rRNA gene V3–V4 amplicon sequencing data generated from salmon skin microbiota samples. The workflow includes sequence quality control, contig assembly, filtering, chimera removal, taxonomic classification, OTU based community analysis, and downstream microbial ecology statistics and visualisation(phyloseq https://joey711.github.io/phyloseq/ AND microeco https://chiliubio.github.io/microeco_tutorial/).

The bioinformatics processing was based on the **mothur MiSeq SOP** workflow, adapted for the dataset used in this study. The mothur MiSeq SOP is described in:

Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD. 2013.  
**Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform.**  
*Applied and Environmental Microbiology* 79(17):5112–5120.

If using this repository or workflow, please cite the mothur MiSeq SOP and indicate the date the SOP was accessed.

## Open science statement

This repository is made openly available in line with UiT – The Arctic University of Norway’s open science principles and PhD research requirements. UiT supports the principle of making research **“as open as possible, as closed as necessary”**, including open access to publications, research data, software, source code, and methodology. By sharing these scripts, the aim is to improve transparency, reproducibility, and re-use of the bioinformatics and statistical workflow used in this study.

## Data availability

The 16S rRNA gene V3–V4 amplicon sequences generated in this study have been deposited in the European Nucleotide Archive, ENA, under BioProject accession:

**PRJEB83870**  
Study accession: **ERP167439**


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


