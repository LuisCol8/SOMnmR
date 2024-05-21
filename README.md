# SOMnmR <a href='https://github.com/LuisCol8/SOMnmR'> <img src= https://github.com/LuisCol8/SOMnmR/assets/35764330/581949d2-386d-4f19-8c9c-5bac6f7b8bc9 align="right" height="300" /></a>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7023315.svg)](https://doi.org/10.5281/zenodo.10137768)

A package for 13C NMR analysis of soil and litter samples.

### What can be done?

The package mainly focuses on the integration of NMR regions. With this information it can either correct the spinning side bands,
give the integrals according to Bonanomi et al. (2011) or to the molecular mixing model (thereafter named MMM) according to Nelson et al. (2005). 

If the MMM option is chosen, then a fitting is made either constrained with the NC ratio of the sample or not.

### What can't be done?

Background substraction, and phase correction is out of the scope of this package.

## Table of contents
- [How to cite](#how-to-cite)
- [How to install](#how-to-install)
  * [Step 1](#step-1)
  * [Step 2](#step-2)
  * [Step 3](#step-3)
- [Introduction](#introduction)
- [Spinning Sidebands Correction](spinning-sidebands-correction.md)
- [Molecular Mixing Model](molecular-mixing-model.md)
- [Statistics](statistics.md)
- [How to use](how-to-use.md)
- [References](references.md)
- [Contact](#contact)

## How to cite

This package is my effort of applying stuff that I learned with what I needed to make things faster during my Ph.D.
If you want me to be motivated to mantain this package and maybe add something or improve it. Show me the love/citations:

Just copy paste this:


Colocho Hurtarte, L. C. (2022). SOMnmR (Version 0.2.0) [Computer software]. https://doi.org/10.5281/zenodo.7023315


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
At our lab, a major equipment of use was the solid state NMR and for a manuscript (under review) I developed this package for quick analysis of 13C NMR data.

The package can do the integration and spinning side bands correction for different sets of 13C NMR regions, namely "4region"(Alkyl, O-Alkyl, Aryl and Carboxyl), "Bonanomi", according to Bonanomi et al., and the Molecular mixing model regions and fitting, according to Nelson et al (2005).

A special thanks goes to Prof. Carol Aldair, and Gabriela Viyalba for encouraging me to finish.

## Contact

If you have suggestions, complains and questions just send me an email to:

luis [dot] colocho [at] diamond [dot] ac [dot] uk
