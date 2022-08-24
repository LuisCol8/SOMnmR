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
# Start by loading the necessary packages, which are cnma, minpack.lm, quadprog, pracma and ggplot2. If you dont have them installed
# just write install.packages("name_of_the_package_you_are_missing")

$ library(cmna)
$ library(minpack.lm)
$ library(quadprog)
$ library(pracma)
$ (ggplot2)
$ library(SOMnmR)


```





