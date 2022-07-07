General Information

This dataset includes all files necessary to perform all analyses and create all figures from Schumacher and Carlson 2022. 

This README is organized as followed:
	Brain Region Analyses -> all of the script files, data files, and any additional files necessary to perform all brain region analyses. Note that some files were downloaded from other sources and included for convenience of replication. I indicated in the file descriptions whether or not each file has been modified from the original source.
	Brain Size Analyses -> all of the script files, data files, and any additional files necessary to perform all brain size analyses. Note that some files were downloaded from other sources and included for convenience of replication. I indicated in the file descriptions whether or not each file has been modified from the original source.
	Synodontis Electrical Recordings -> the custom code necessary to record and plot the electrical recordings from Synodontis individuals
	The phylogenetic tree from Hughes et al 2018 used to make Figure 1
	Instructions to obtain raw micro-computed tomography scans

-----

Brain Region Analyses

The osetoglossiform and otophysan brain measurement data is available in SchumacherAndCarlson_SupplementaryFile1.csv, which is also a supplementary file for the paper. This contains all of the brain measurement data from this study obtained from micro-computed tomography scans.
Measured brain regions are: OB = olfactory bulbs, TEL = telencephalon, HB = hindbrain, OT = optic tectum, TS = torus semicircularis, CB = cerebellum, RoB = rest of brain.
Units are mm^3.
Other variables in this file include: 
	Species, 
	Family, 
	Specimen ID which indicates the individual, body mass in grams, total brain volume which is the sum of each of the region volumes, 
	Group which indicates the combined electrosensory phenotype and lineage for each individual (i.e. Pulse mormyroid, wave mormyroid, pulse gymnotiform, wave gymnotiform, Synodontis siluriform, non-electric siluriform, outgroup otophysan, and outgroup osteoglossiform), 
	Electroreceptor which indicates the type of electroreceptor that each individual has (i.e. tuberous and ampullary electroreceptors -> called tuberous, only ampullary electroreceptors, or none), 
	Electrogenic which indicates whether each individual is electrogenic or not, and 
	Lineage which indicates whether each individual is an osteoglossiform or otophysan.
Total rows = 113, total columns = 16

The maximum clade credibility tree estimated using beast is available in ConvMosEvol_tree.nex.
Note that to include data from unsequenced species, we used sequence data from the species in the same monophyletic genus with the shortest distance to the genus node.
Thus, these species names were changed in the .nex tree file.
	Mormyrus_tapirus = sequences from Mormyrus_rume
	Petrocephalus_tenuicauda = sequences from Petrocephalus_microphthalmus
	Campylomormyrus_sp = sequences from Campylomormyrus_numenius

Code for all brain region analyses is available in SchumacherAndCarlson_ConvMosEvol_Analysis.R.

phylo.fda.v0.2.R is the script from Motani and Schmitz 2011 for the pFDA, downloaded from github. This file has not been modified at all from the original source.

-----

Brain Mass Analyses

The compiled actinoptergyian brain and body mass data from Tsuboi et al 2018, Tsuboi 2021, and this study are available in SchumacherAndCarlson_SupplementaryFile2.csv, which is also a supplementary file for the paper. Units are g.
Variables from this file include: 
	Order as indicated in the original document from Tsuboi et al 2018, 
	Updated order which include classification changes since the original document, 
	Family, 
	Genus, 
	Species, 
	Subspecies, 
	Genus_species combined, 
	Sample size from the original sources, 
	Sex, 
	Method in which the brain mass was obtained (either measured as mass or converted from volume), 
	Body weight in grams of each individual, 
	Brain weight in grams of each individual, 
	Standard length in cm of each individual, 
	Total length in cm of each individual, 
	Age class of each individual, 
	Species name in the original article if it has since been reclassified, 
	Reference from which data from each individual originated, and 
	Any particular remarks about an individual.
Total rows = 4069, total columns = 18

The time calibrated tree from Rabosky et al 2018 is actinopt_12k_treePL.tre. This file has not been modified at all from the file included with Rabosky et al 2018.

Code for brain mass analyses are available in the following files:
	SchumacherAndCarlson_bayouGG.R: 1 global allometry across all taxa
	SchumacherAndCarlson_bayouSG.R: shifts allowed in intercept but not slope
	SchumacherAndCarlson_bayouSS.R: shifts allowed in both intercept and slope, figure, and PGLS comparisons
	SchumacherAndCarlson_bayouAF.R: fixed shifts at ampullary electroreceptor taxa
	SchumacherAndCarlson_bayouTF.R: fixed shifts at tuberous electroreceptor taxa
	SchumacherAndCarlson_bayouEF.R: fixed shifts at electrogenic taxa
	SchumacherAndCarlson_bayouOF.R: fixed shift at osteoglossiforms
-----	

Synodontis Electrical Recordings

Code for recording, extracting, and plotting the synodontid electrical recordings are available in the following files: SchumacherAndCarlson_SynodontisBehaviorScript.m, contrec1.m, contrec1.rcx, contrec1anal_ELS.m

-----

31Calibrations_Chronogram_newick.tre is the time calibrated phylogeny from Hughes et al 2018 used to make figure 1. This file has not been modified at all from the file included with Hughes et al 2018.

-----

The raw micro-computed tomography scans are too large to post (multiple TBs), but are available upon request. 
To request raw otophysan and/or osteoglossiform scans, contact the corresponding author. We ask that those who want access to the scan data send us an external hard drive, which we will upload all the data to and then return. 