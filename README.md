<!---
This README uses Markdown syntax.
You can open it with R software or view it from Github.
--->

=====================

Provenance for this README
--------------------------

* File name: README.md
* Authors: Anne Chao
* Other contributors: Chun-Huo Chiu, Sebastian Vogel, Peter Kriegel, and Simon Thorn
* Date created: 2023-03-20

Dataset Attribution and Usage
--------------------------

* Dataset Title: Data from: Quantifying and estimating ecological network diversity based on incomplete sampling data.
* Persistent identifier: DOI: https://doi.org/10.5061/dryad.ksn02v78d Also available from Github https://github.com/AnneChao/PTB_iNEXT.link
* Dataset contributors: Anne Chao, Chun-Huo Chiu, Sebastian Vogel, Peter Kriegel, and Simon Thorn
* Suggested cittaions: 
  * Methodology citation:  
   > Chiu, C-H., Chao, A., Vogel, S., Kriegel, P. and Thorn, S. (2023). Network-diversity quantification and related statistical estimation: drawing on sampling models and methodologies from biodiversity research. Philosophical Transactions of the Royal Society B.
  * Dataset and code citation
   > Chao A, Chiu C-H, Vogel S, Kriegel P, and Thorn S. (2023). Data from: Quantifying and estimating ecological network diversity based on incomplete sampling data. Dryad Digital Repository. (doi: https://doi.org/10.5061/dryad.ksn02v78d)
   
   
Tree-beetle interaction data collection:
----------------------------------------
The experiment underlying the data was conducted by Vogel and colleagues in the Steigerwald forest located in northern Bavaria, Germany (49°32' N; 10°23' E). The Steigerwald forest covers around 16,500 ha of forests from colline up to sub-montane mountain range with a mean annual temperature around 8°C and an annual precipitation around 750 to 800 mm (forest climate station Ebrach). Temperatures higher than 8°C are reached on about 160 days per year. 
In 2015, six study sites (Plot A-F) were established. This dataset includes the data based on two subplots within each plot: (1) open habitat (sun-exposed on a forest meadow), and (2) closed habitat (canopy-shaded within a closed forest stand). Freshly cut logs and branch bundles of six tree species (Abies alba, Carpinus betulus, Fagus sylvatica, Pinus sylvestris, Populus tremula and Quercus petraea) were exposed on each subplot. 
Deadwood objects were sampled for saproxylic beetles using stem emergence traps. With this sampling technique it is possible to collect emerging individuals, developing in deadwood, while excluding other occasional visitors. Traps were filled with saturated saline solution as sampling fluid and emptied monthly between April and September in 2016, 2017, 2018, and 2019. Beetles were classified as saproxylic according to Schmidl and Bussler. 


Methodological Information
==========================

* Methods of data collection/generation: see the methodology paper by Chiu et al. (2023) for details


Data and File Overview
======================

Summary Metrics
---------------

* File count: 7
* Total file size: 7.72 MB
* Range of individual file sizes: 2 KB - 1430 KB

Naming Conventions
------------------

* File naming scheme: files with the format ".txt" or ".csv" denote data files. 


Table of Contents
-----------------

* Data_tree-beetle_interaction_frequency.csv
* Data_phylo_tree.txt
* Data_distance_matrix.txt
* R_code.R


Setup
-----

* Recommended software/tools: parentage- RStudio 2021.03.7; R version 4.0.0

* Required R packages: To run the R code, you must download the package "iNEXT.3D", "iNEXT.4steps", "iNEXT.beta3D", and "iNEXT.link" from Anne Chao's Github. 


File Details
===================

Details for: Data_tree-beetle_interaction_frequency.csv
---------------------------------------

* Description: a file containing tree-beetle interaction frequencies. The data used for examples are based on tree-beetle interaction frequency data collected from 2016 to 2019 in six study sites (Plot A to F) located in the Steigerwald forest, Germany. 

* Format(s): .csv

* Size(s): 305 KB

* Dimensions: 3891 rows x 9 columns

* Variables:
  * ID: unique identification code for each observation.
  * species: beetle species name (genius and species). Total number of species is 264.
  * plot: six study plots-- "A", "B", "C", "D", "E", and "F".
  * treatment: three treatments-- "G" (ground, canopy-shaded within closed habitat), "N" (net-shaded habitat), and "O" (sun-exposed open habitat).
  * tree_spec: "A", "E", "H", "K", "R", and "T", where T (Abies alba), H (Carpinus betulus), R (Fagus sylvatica), K (Pinus Sylvestris), A (Populus tremula) and E (Quercus petraea).
  * year: data collection year--2016, 2017, 2018 and 2019.
  * abundance: the tree-beetle interaction frequency. Range is from 1 to 1094.
  * tree_species: the same as the column "tree_spec".
  * inte: the tree-beetle interaction by the two columns "tree-species" and "species".
  
  
* Dataset citation:

  * Hagge J et al. 2021. What does a threatened saproxylic beetle look like? Modelling extinction risk using a new morphological trait database. Journal of Animal Ecology 90, 1934-1947.
  
 
Details for: Data_phylo_tree.txt
---------------------------------------

* Description: a file containing the phylogeny tree for 264 beetle species in newick format. 

* Format(s): .txt

* Size(s): 55.6 KB

* Dimensions: Phylogenetic tree with 1449 tips and 757 internal nodes. The tree height is 285 million years.

* Dataset citation:

  * Seibold S, Brandl R, Buse J, Hothorn T, Schmidl J, Thorn S, Muller J. 2015. Association of extinction risk of saproxylic beetles with ecological degradation of forests in Europe. Conservation Biology 29, 382-390.
  
  * Vogel S, Gossner MM, Mergner U, Muller J, Thorn S. 2020. Optimizing enrichment of deadwood for biodiversity by varying sun exposure and tree species: An experimental approach. Journal of Animal Ecology 57, 2075-2085.



Details for: Data_distance_matrix.txt
---------------------------------------

* Description: a file containing functional distance matrix between any two beetle species in the file "Data_tree-beetle_interaction_frequency.csv" based on the Gower distance. The distance ranges from zero and one.

* Format(s): .txt

* Size(s): 1.19 MB

* Dimensions: 265 rows x 265 columns

* Dataset citation:
  *	Chao A, Chiu, C-H, Villeger S, Sun I-F, Thorn S, Lin Y-C, Chiang JM, Sherwin WB. 2019. An attribute-diversity approach to functional diversity, functional beta diversity, and related (dis)similarity measures. Ecological Monographs 89, e01343. 


Details for: R_code.R
---------------------------------------

* Description: Main code for plotting all figures.

* Format(s): .R

* Size(s): 13.8 KB

* Dimensions: 356 lines

* Before using the R_code.R, you must download the package "iNEXT.3D", "iNEXT.4steps", "iNEXT.beta3D", and "iNEXT.link" from Anne Chao's Github.
* Because this code runs a random bootstrapping process in the background, with default = 100 replications/loops, to estimate standard error, the output (for se and 95% confidence intervals) will vary very slightly each time you enter the same data. Also, the bootstrap procedure is VERY time-consuming; please change the argument nboot = 100 (default) in any function to nboot = 0 (to skip the bootstrap process), or nboot = 20 (to roughly browse the resulting confidence intervals).   
   
  
- - -
END OF README
