# SOMnmR
A package for 13C NMR analysis (spinning sidebands, molecular mixing model, etc)

## Table of contents
* [How to install](##How-to-install)
* [Introduction](##Introduction)
* [How to use](##How-to-use)

## How to install

Step 1: Install the devtools package

To install a R package, start by installing the devtools package. The best way to do this is from CRAN, by typing:
1
	
install.packages("devtools")
Step 2: Install the package of interest from GitHub

Install the package of interest from GitHub using the following code, where you need to remember to list both the author and the name of the package (in GitHub jargon, the package is the repo, which is short for repository). In this example, we are installing the flipPlots package created by Displayr.
1
2
	
library(devtools)
install_github("Displayr/flipPlots")
Installing GitHub R packages

 
Step 3: Load the package
1
	
library(flipPlots)

## Introduction

## How to use
