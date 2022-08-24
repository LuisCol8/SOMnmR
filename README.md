# SOMnmR
A package for 13C NMR analysis (spinning sidebands, molecular mixing model, etc)

## Table of contents
* [How to install](##How-to-install)
* [Introduction](##Introduction)
* [How to use](##How-to-use)

## How to install

### Step 1: 
To install a SOMnmR as a package, start by installing the devtools package. The best way to do this is from CRAN, by typing from your comand line:

```bash
# Install devtools
$ install.packages("devtools")
```

### Step 2: 
Install the package of interest from GitHub

Install SOMnmR from GitHub using the following code:

```bash
# load devtools
$ library(devtools)

# load SOMnmR
$ install_github("LuisCol8/SOMnmR")
```

 
### Step 3: 
Load the package
	
```bash
# load SOMnmR
$ library(SOMnmR)

```

## Introduction

## How to use

As explained before the package is a simple wrap up of integration and fitting (in the case of the MMM). As an example I include a simple script that can be easily modified.

```bash
# Start by loading the necessary packages, which are cnma, minpack.lm, quadprog, pracma and ggplot2. 
# If you dont have them installedjust write install.packages("name_of_the_package_you_are_missing")

$ library(cmna)
$ library(minpack.lm)
$ library(quadprog)
$ library(pracma)
$ (ggplot2)
$ library(SOMnmR)

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


```





