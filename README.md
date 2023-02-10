# mycorrhizal_succession_iceCommunities

This repository contains the custom code and data to replicate the results published in Carteron et al. Mycorrhizal fungi dynamics and drivers after glacier retreat. Preprint: TBA. Published article: TBA.

For this study, we used environmental DNA metabarcoding and measurements of microhabitat characteristics to assess the evolution of mycorrhizal communities during soil development in 46 glacier forelands across the globe.

The scripts are ordered in this way:
- `import.R` (for accessing the data)
- `modeling.R` (for modeling)
- `figure.R` (for the figures)

Files available in the `data` folder:
- `data_filtered_Fung02.zip` Filtered sequence data from Fung02 marker in compressed format.
- `data_filtered_Sper01.zip` Filtered sequence data from Sper01 marker in compressed format.
- `env.myco.09.02.23.csv` Environmental data containing plot information, time since glacier retreat, soil chemistry (pH, C, N, P), temperature (meanT), productivity (NDVI), wetness (WTI) and plot coordinates. 
- `fung.phy.relax.rds` Phyloseq object of the fungal data.
- `sper.phy.relax.rds` Phyloseq object of the plant data.
- `full.table.09.02.23.csv` Table containing calculated diversity values with corresponding and environmental values. This file can be produce using the `import.R` script.  

Raw sequences from ITS and trnL amplification are deposited at https://doi.org/10.5281/zenodo.6620359  

Please cite original data if you used them as: Alessia Guerrieri, Aurélie Bonin, Ludovic Gielly, & Gentile Francesco Ficetola. (2022). Raw sequencing data for studying the colonization of soil communities after glacier retreat [Data set]. Zenodo. https://doi.org/10.5281/zenodo.6620359  

Related article: TBA

In this repository, there are in addition the `README.md` (this file) and `sessionInfo.md` for more information about R, package and OS versions.
