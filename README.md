# SOMnmR
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7023316.svg)](https://doi.org/10.5281/zenodo.7023316)

A package for 13C NMR analysis (spinning sidebands, molecular mixing model, etc)

### What can be done?

The package mainly focuses on the integration of NMR regions. With this information it can either correct the spinning side bands,
give the integrals according to Bonanomi et al. (2011) or to the molecular mixing model (thereafter named MMM) according to Nelson et al. (2005). 

If the MMM option is chosen, then a fitting is made either constrained with the NC ration of the sample or not.

### What can't be done?

Background substraction, and phase correction is out of the scope of this package.

## Table of contents
- [How to cite](#how-to-cite)
- [How to install](#how-to-install)
  * [Step 1](#step-1)
  * [Step 2](#step-2)
  * [Step 3](#step-3)
- [Introduction](#introduction)
- [How to use](#how-to-use)
- [References](#references)
- [Contact](#contact)

## How to cite

This package is my effort of applying stuff that I learned with what I needed to make things faster during my Ph.D.
If you want me to be motivated to mantain this package and maybe add something or improve it. Show me the love/citations:

Just copy paste this:


Colocho Hurtarte, L. C. (2022). SOMnmR (Version 0.1.0) [Computer software]. https://doi.org/10.5281/zenodo.7023316


## How to install

### Step 1 (Install devtools)
To install a SOMnmR as a package, start by installing the devtools package. The best way to do this is from CRAN, by typing from your comand line:

```bash
# Install devtools
  install.packages("devtools")
```

### Step 2 (.....)
Install the package of interest from GitHub

Install SOMnmR from GitHub using the following code:

```bash
# load devtools
  library(devtools)

# load SOMnmR
  install_github("LuisCol8/SOMnmR")
```

 
### Step 3 (Profit)
Load the package
	
```bash
# load SOMnmR
  library(SOMnmR)

```

## Introduction

Hi My name is Luis Colocho.
I made my Ph.D. at the Chair of Soil Science of the Technical University of Munich, under Apl. Prof. Joerg Prietzel.
At our lab, a major equipment of use was the solid state NMR and for a manuscript (under review) I developed this package, which takes some functions
from the package developed by my college Florian Werner (check LCF if you do linear combination fitting of synchrotron data).

The package has some pre set regions like spinning side bands (which makes a spinning side band correction by subtracting the spinning sidebands from
carboxyl-C and O-N-alkyl-C region and adding the proportion to the aryl-C), Bonanomi,based on the reported regions of Bonanomi et al. (2011) and the MMM  according to Nelson et al. (2005).

A special thanks goes to Prof. Carol Aldair for encouraging me to finish.

## About the reference compounds of the Molecular mixing model 

The original paper from Nelson et al. (1999), used an amalgam of pure compounds as references

| Chemical shift region (ppm)   | Aminoacid | Hexose  | Lignin | Cutin  | Suberin | Char  | Chitin | Carbonyl  |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| 210 to 165 (Amide/Carboxyl) | 26.4 | 0.0  | 4.6  | 6.6  | 4.3  | 5.6  | 12.5  | 100.0  |
| 165 to 145 (Phenolic)  | 2.5  | 0.0  | 19.5  | 0.7 | 6.2  | 16.1  | 0.0  | 0.0  |
| 145 to 110 (Aromatic)  | 7.5  | 1.0  | 30.5  | 3.6  | 19.1  | 73.9  | 0.0  | 0.0  |
| 110 to 95 (Di-O-Alkyl)  | 0.0  | 15.7  | 8.6  | 0.0  | 0.0  | 4.3  | 12.5  | 0.0  |
| 95 to 60 (O-Alkyl)  | 2.1  | 79.0  | 12.5  | 9.0  | 6.8  | 0.0  | 50.0  | 0.0  |
| 60 to 45 (N-Alkyl/Methoxyl)  | 4.3  | 21.9  | 13.8 | 2.9  | 0.0  | 0.0  | 12.5  | 0.0  |
| 45 to -10 (Alkyl)  | 39.6 | 0.0  | 10.5  | 75.6  | 60.7  | 0.0  | 12.5  | 0.0  |
| Molar N:C   | 0.32  | 0.0  | 0.0  | 0.0  | 0.0  | 0.0  | 0.125  | 0.0  |

A subsequent paper from Baldock et al (2004), uses a modified version of this, adjusted for terrestrial ecosystems.

| Chemical shift region (ppm)   | Carbohydrate | Protein  | Lignin | Lipid  | Carbonyl | Char  |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| 210 to 165 (Amide/Carboxyl) | 0.0 | 28.3  | 4.6  | 5.6  | 100  | 3.9  |
| 165 to 145 (Phenolic)  | 0.0  | 1.3  | 19.5  | 0.0 | 0.0  | 15.2  |
| 145 to 110 (Aromatic)  | 0.0  | 8.9  | 30.5  | 0.0  | 0.0  | 72.0  |
| 110 to 95 (Di-O-Alkyl)  | 16.7  | 0.0  | 8.6  | 0.0  | 0.0  | 5.3  |
| 95 to 60 (O-Alkyl)  | 83.3  | 3.5  | 12.5  | 0.0  | 0.0  | 1.8  |
| 60 to 45 (N-Alkyl/Methoxyl)  | 0.0  | 22.6  | 13.8 | 0.5  | 0.0  | 1.7  |
| 45 to -10 (Alkyl)  | 0.0 | 35.4  | 10.5  | 94.4  | 0.0  | 0.0  |
| Molar N:C   | 0.0  | 0.275  | 0.0  | 0.0  | 0.0  | 0.0  |

with their molar elemental ratios (normalised against C)

| Element   | Carbohydrate | Protein  | Lignin | Lipid  | Carbonyl | Char  |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| C   | 1.0 | 1.0  | 1.0 | 1.0  | 1.0 | 1.0  |
| N   | 0.0 | 0.27  | 0.0 | 0.0  | 0.0 | 0.0  |
| H   | 1.67 | 1.1  | 1.24 | 1.94  | 1.0 | 0.45  |
| O   | 0.83 | 0.16  | 0.43 | 0.24  | 2.0 | 0.41  |

and acquatic ecosystems

| Chemical shift region (ppm)   | Carbohydrate | Protein  | Lignin | Lipid  | Carbonyl | Char  |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| 210 to 165 (Amide/Carboxyl) | 0.0 | 25.7  | 4.6  | 5.6  | 100  | 3.9  |
| 165 to 145 (Phenolic)  | 0.0  | 1.4  | 19.5  | 0.0 | 0.0  | 15.2  |
| 145 to 110 (Aromatic)  | 0.0  | 9.1  | 30.5  | 0.0  | 0.0  | 72.0  |
| 110 to 95 (Di-O-Alkyl)  | 16.7  | 0.0  | 8.6  | 0.0  | 0.0  | 5.3  |
| 95 to 60 (O-Alkyl)  | 83.3  | 2.3  | 12.5  | 0.0  | 0.0  | 1.8  |
| 60 to 45 (N-Alkyl/Methoxyl)  | 0.0  | 21.0  | 13.8 | 0.5  | 0.0  | 1.7  |
| 45 to -10 (Alkyl)  | 0.0 | 40.4  | 10.5  | 94.4  | 0.0  | 0.0  |
| Molar N:C   | 0.0  | 0.275  | 0.0  | 0.0  | 0.0  | 0.0  |

with their molar elemental ratios (normalised against C)

| Element   | Carbohydrate | Protein  | Lignin | Lipid  | Carbonyl | Char  |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| C   | 1.0 | 1.0  | 1.0 | 1.0  | 1.0 | 1.0  |
| N   | 0.0 | 0.27  | 0.0 | 0.0  | 0.0 | 0.0  |
| H   | 1.67 | 1.17  | 1.24 | 1.89  | 1.0 | 0.45  |
| O   | 0.83 | 0.12  | 0.43 | 0.11  | 2.0 | 0.41  |

suggests that two sets of compounds can be used, depending on the origin of the SOM of the sample, i.e. Terrestrial or Aquatic.

For this he considered the following composition for a terrestrial ecosystem

| Chemical shift region (ppm)   | Carbohydrate | Protein  | Lignin | Lipid  | Carbonyl | Char  |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| 210 to 165 (Amide/Carboxyl) | 0.0 | 26.4  | 4.6  | 6.6  | 100  | 5.7  |
| 165 to 145 (Phenolic)  | 0.0  | 2.5  | 19.5  | 0.7 | 0.0  | 16.1  |
| 145 to 110 (Aromatic)  | 1.0  | 7.5  | 30.5  | 3.6  | 0.0  | 73.9  |
| 110 to 95 (Di-O-Alkyl)  | 15.7  | 0.0  | 8.6  | 0.0  | 0.0  | 4.3  |
| 95 to 60 (O-Alkyl)  | 79.0  | 2.1  | 12.5  | 9.0  | 0.0  | 0.0  |
| 60 to 45 (N-Alkyl/Methoxyl)  | 4.3  | 21.9  | 13.8 | 4.5  | 0.0  | 0.0  |
| 45 to -10 (Alkyl)  | 0.0 | 39.6  | 10.5  | 75.6  | 0.0  | 0.0  |
| Molar N:C   | 0.0  | 0.32  | 0.0  | 0.0  | 0.0  | 0.0  |

While for an aquatic ecosystem, the composition of protein and lipid changes, as it is lead by phytoplancton

| Chemical shift region (ppm)   | Carbohydrate | Protein  | Lignin | Lipid  | Carbonyl | Char  |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| 210 to 165 (Amide/Carboxyl) | 0.0 | 30.4  | 4.6  | 5.6  | 100  | 5.7  |
| 165 to 145 (Phenolic)  | 0.0  | 1.0  | 19.5  | 0.0 | 0.0  | 16.1  |
| 145 to 110 (Aromatic)  | 1.0  | 4.5  | 30.5  | 0.0  | 0.0  | 73.9  |
| 110 to 95 (Di-O-Alkyl)  | 15.7  | 0.0  | 8.6  | 0.0  | 0.0  | 4.3  |
| 95 to 60 (O-Alkyl)  | 79.0  | 2.9  | 12.5  | 0.0  | 0.0  | 0.0  |
| 60 to 45 (N-Alkyl/Methoxyl)  | 4.3  | 24.7  | 13.8 | 0.5  | 0.0  | 0.0  |
| 45 to -10 (Alkyl)  | 0.0 | 36.6  | 10.5  | 94.4  | 0.0  | 0.0  |
| Molar N:C   | 0.0  | 0.266  | 0.0  | 0.0  | 0.0  | 0.0  |

However, a newer version of these standard compounds was calculated consiering, for cabohydrates a cellulose model, for protein based on extractable aminoacids from Friedel and Scheller (2002), the lignin is  assumed to be equivalent to  the average of those acquired for spruce and red alder lignin  presented by Wilson (1987), and the lipid is calculated from structures proposed by Kolattukudy (1980).

| Chemical shift region (ppm)   | Carbohydrate | Protein  | Lignin | Lipid  | Carbonyl | Char  |
| ------------- | ------------- | ------------- | ------------- | ------------- | ------------- | ------------- |
| 210 to 165 (Amide/Carboxyl) | 0.0 | 28.3  | 4.6  | 5.6  | 100  | 3.9  |
| 165 to 145 (Phenolic)  | 0.0  | 1.3  | 19.5  | 0.0 | 0.0  | 15.2  |
| 145 to 110 (Aromatic)  | 0.0  | 8.9  | 30.5  | 0.0  | 0.0  | 72.0  |
| 110 to 95 (Di-O-Alkyl)  | 16.7  | 0.0  | 8.6  | 0.0  | 0.0  | 5.3  |
| 95 to 60 (O-Alkyl)  | 83.3  | 3.5  | 12.5  | 0.0  | 0.0  | 1.8  |
| 60 to 45 (N-Alkyl/Methoxyl)  | 0.0  | 22.6  | 13.8 | 0.5  | 0.0  | 1.7  |
| 45 to -10 (Alkyl)  | 0.0 | 35.4  | 10.5  | 94.4  | 0.0  | 0.0  |
| Molar N:C   | 0.0  | 0.275  | 0.0  | 0.0  | 0.0  | 0.0  |

## How to use

As explained before the package is a simple wrap up of integration and fitting (in the case of the MMM). As an example I include a simple script that can be easily modified.

```bash
# Start by loading the necessary packages, which are cnma, minpack.lm, quadprog, pracma and ggplot2. 
# If you dont have them installedjust write install.packages("name_of_the_package_you_are_missing")

  library(cmna)
  library(minpack.lm)
  library(quadprog)
  library(pracma)
  (ggplot2)
  library(SOMnmR)

# set your working directory, I suggest you have a folder for where all the integrals are placed, 
# a folder for the C:N data and a folder for the plots

work.dir <- c("C:/Plant_NMR/")

# Load the folder containing the spectra, you can have them as a two column file,
# with the ppm and the intensity in each column

setwd(paste(work.dir,"csv",sep = "/"))
files <- list.files()
files <- files[grep(".csv", files)]

# Read all the spectra and make a variable containing them all. Types accepted: "Bruker", "coma", "tab"
spec <- read_raw_spec(files, filetype = "coma")

# If you are going to do the MMM fitting you need to create a file contaning the  Carbon and Nitrogen
# data of your spectra. With this function we will create an empty file which contains the columns 
# for C and N and the name of the files. Do not change the order as it is following the order of the loaded spectra

ncdata <- mk_nc_data(spec)

# save and fill with CN data
setwd(paste(work.dir,"CN",sep = "/"))
write.csv(ncdata,"./NC_plant.csv", row.names = FALSE)

# After you filled the file with the CN data load it into the script

ncdata <- nc_data("/NC/NC_plant.csv")

# Now you can calculate Molecular mixing model! Choose as NMR Meth either "MMMFixN" to fix the NC ratio of the fit
# or "MMM" to let it vary, you can also choose ecosys either "Terr" or "Aqua" which takes into account the different
# composition of terrestrial or aquatic ecosystems.

MMM <- region_calc(spec, NMRmeth = "MMMFixN", ecosys= "Terr", cndata = ncdata)

# Save your results using:

write.csv(MMM, "MMM_R.csv")

# If you want to check the integrals (for SSB calculated and corrected SSBcorr = TRUE,
# for SSB calculated SSBcorr = FALSE) of the MMM use:

MMMregions <- int_nmr(spec,NMRmeth = "MMM", SSBcorr = TRUE)

# You can save each result (change the number in x [[x]]) as follows: 

write.csv(unlist(MMMregions[[1]]$data$Integral), "MMMregions_name.csv")

# There is a simple function that makes a ggplot2 plot of your spectra. You can choose between different integration regions
# to be colored. Use "Bonanomi" or "MMM" or leave it blank to choose between the "Bonanomi", "MMM", or spinning side bands regions,
# respectively. 

plot_NMR(spec, NMRmeth = "MMM",file.output = TRUE,use.tiff = TRUE)



```

## References

Bonanomi, G., Incerti, G., Barile, E., Capodilupo, M., Antignani, V., Mingo, A., Lanzotti, V., Scala, F., & Mazzoleni, S. (2011). Phytotoxicity, not nitrogen immobilization, explains plant litter inhibitory effects: Evidence from solid-state 13C NMR spectroscopy. New Phytologist, 191(4), 1018–1030. https://doi.org/10.1111/j.1469-8137.2011.03765.x

Hockaday, W. C., Masiello, C. A., Randerson, J. T., Smernik, R. J., Baldock, J. A., Chadwick, O. A., & Harden, J. W. (2009). Measurement of soil carbon oxidation state and oxidative ratio by 13C nuclear magnetic resonance. Journal of Geophysical Research: Biogeosciences, 114(2), 1–14. https://doi.org/10.1029/2008JG000803

Kögel-Knabner, I. (2017). The macromolecular organic composition of plant and microbial residues as inputs to soil organic matter: Fourteen years on. Soil Biology and Biochemistry, 105, A3–A8. https://doi.org/10.1016/j.soilbio.2016.08.011

Nelson, P. N., & Baldock, J. A. (2005). Estimating the molecular composition of a diverse range of natural organic materials from solid-state 13 C NMR and elemental analyses. Biogeochemistry, 72(1), 1–34. https://doi.org/10.1007/s10533-004-0076-3

Wilson, M.A., 1987. NMR Techniques and Applications in Geochemistry and Soil Chemistry. Pergamon
Press, Oxford.

## Contact

If you have suggestions, complains and questions just send me an email to:

luis [dot] colocho [at] esrf [dot] fr
