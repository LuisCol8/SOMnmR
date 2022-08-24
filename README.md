# SOMnmR
A package for 13C NMR analysis (spinning sidebands, molecular mixing model, etc)

## Table of contents
* [How to install](##-How-to-install)
* [Introduction](##-Introduction)
* [How to use](##-How-to-use)

## How to install

### Step 1: 
To install a SOMnmR as a package, start by installing the devtools package. The best way to do this is from CRAN, by typing from your comand line:

```bash
# Install devtools
  install.packages("devtools")
```

### Step 2: 
Install the package of interest from GitHub

Install SOMnmR from GitHub using the following code:

```bash
# load devtools
  library(devtools)

# load SOMnmR
  install_github("LuisCol8/SOMnmR")
```

 
### Step 3: 
Load the package
	
```bash
# load SOMnmR
  library(SOMnmR)

```

## Introduction

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

# If you want to check the integrals of the MMM use:

MMMregions <- int_nmr(spec,NMRmeth = "MMM")

# You can save each result (change the number in x [[x]]) as follows: 

write.csv(unlist(MMMregions[[1]]$data$Integral), "MMMregions_name.csv")

# There is a simple function that makes a ggplot2 plot of your spectra. You can choose between different integration regions
# to be colored. Use "Bonanomi" or "MMM" or leave it blank to choose between the "Bonanomi", "MMM", or spinning side bands regions,
# respectively. 

plot_NMR(spec, NMRmeth = "MMMM",file.output = TRUE,use.tiff = TRUE)



```





