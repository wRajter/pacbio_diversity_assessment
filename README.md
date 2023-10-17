# PacBio Diversity Assessment of Eukaryotic Microbes (Protists)

This repository is dedicated to the exploration of eukaryotic microbial (protist) diversity leveraging the long-read PacBio sequencing method.

## Background

The raw sequences used in this project, representing full 18S rRNA, originate from the `vampyrella_2023` project. Detailed information about the samples and the pipeline from which they were derived can be found [here](https://github.com/wRajter/vampyrella_2023).


## Workflow Overview

1. **Protist Sequence Filtration**: Using the OTUs from the `vampyrella_2023` pipeline, we filter the dataset to retain only sequences belonging to protists.
   - Script: **<<add the name of the script here>>**

2. **Phylogenetic Placement**: Phylogenetic placement techniques are employed to get a panoramic view of the taxa present within the samples.
   - Scripts: **<<add the name of the scripts here>>**

3. **Targeted Taxa Exploration**: After a comprehensive evaluation of the taxonomic composition, we pinpoint specific taxa that are of particular interest for an in-depth exploration.

## Directory Structure:

- **notebooks:** This directory contains Jupyter notebooks that chronologically represent each major step in the analysis. They provide a step-by-step breakdown of the entire workflow, ensuring transparency, repeatability, and ease of understanding.

- **raw_data:** This is where all the input data files reside. Given the substantial size of sequencing data and other related files, they are maintained outside of the GitHub repository.

- **results:** Any significant output, be it in the form of tables, interim data files, or visualizations (figures), will be saved in this directory. This ensures a clear distinction between raw input and processed output.

- **scripts:** Supplementary code and utilities used within the Jupyter notebooks are housed here. These scripts may range from custom Python functions to Bash or R scripts required for specific tasks.


## Future Directions

This analysis serves as a foundational step in understanding the complexity of protist diversity in the sampled environments. The results and insights from this project will be instrumental in guiding further investigations and research into these eukaryotic microbes.
