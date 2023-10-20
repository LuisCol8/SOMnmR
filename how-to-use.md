# How to use the package SOMnmR

[Back to the Main page](README.md)

## Basic example

Here we want to just load the data (NMR spectra) and do a simple integration, which will return the integrals and the SSB integrals.
This means that the output is not yet SSB corrected.

```bash
# load necessary packages
  library(minpack.lm)
  library(quadprog)
  library(pracma)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(IntervalSurgeon)
  library(SOMnmR)

# load data
## set working directory
work.dir <- c("C:/Documents/Data/Experiment_NMR/")

## go to directory containing all spectra
setwd(paste(work.dir,"NMR_integrals",sep = "/"))

## list all the files
files <- list.files()

## list all the files that end with .txt (or .csv, depending on how you saved the spectra)
files <- files[grep(".txt", files)]

# Load data, choose either  "Bruker" (standard Bruker file), "coma" (coma separated value), "tab" (tab separated value
spec <- read_raw_spec(files, filetype = "Bruker")

# Integrate the whole data set, NMRmeth can be "4region", "Bonanomi" and "MMM". I will chose 4 region.
# For this you will need to input the magnetic field of your NMR machine and the spinning frequency of the probe.
# In this example we used an NMR of 200 MHz and a spinning frequency of 6800 Hz.

Integralregions <- int_nmr(spec, NMRmeth = "4region", NMR_field = 200, NMR_rotation = 6800)

#To view the output of a single spectrum
View(Integralregions[[1]][["data"]][["Integral"]])

```
The output will be a list and in each list in the [["data"]][["Integral"]] path, you will find a table as follows

![image](https://github.com/LuisCol8/SOMnmR/assets/35764330/91030a35-df84-4fcb-9c94-f7bb886c8b08)

To understand this table check the [Spinning sidebands correction page](spinning-sidebands-correction.md)

## Correcting for the Spinning Side Bands (SSB)

We now know how to load a spectrum a perform a simple integration. There is a wraped function that performs the integration, and corrects for the SSBs. The example follows

```bash
# load necessary packages
  library(devtools)

# load data


# now use the region_calc funtion to perform both the integration and the SSB correction.
# For this you will need to input the magnetic field of your NMR machine and the spinning frequency of the probe.
# In this example we used an NMR of 200 MHz and a spinning frequency of 6800 Hz.



```

## Fitting the Molecular Mixing model (MMM)

In the case of The MMM, we can use the same region_calc funtion to perform the fitting. But for this, we also need to create a table of C and N data that matches the order of the files, and we also have to select which table of model compounds to use. The example follows:

```bash
# load necessary packages
  library(minpack.lm)
  library(quadprog)
  library(pracma)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(IntervalSurgeon)
  library(SOMnmR)

# load data
## set working directory
work.dir <- c("C:/Documents/Data/Experiment_NMR/")

## go to directory containing all spectra
setwd(paste(work.dir,"NMR_integrals",sep = "/"))

## list all the files
files <- list.files()

## list all the files that end with .txt (or .csv, depending on how you saved the spectra)
files <- files[grep(".txt", files)]

# Load data, choose either  "Bruker" (standard Bruker file), "coma" (coma separated value), "tab" (tab separated value
spec <- read_raw_spec(files, filetype = "Bruker")

# Integrate and correct for the SSB of the whole data set, NMRmeth can be "4region", "Bonanomi" and "MMM". I will chose 4 region.
# For this you will need to input the magnetic field of your NMR machine and the spinning frequency of the probe.
# In this example we used an NMR of 200 MHz and a spinning frequency of 6800 Hz.

IntegralSSBc <- region_calc(spec, NMRmeth = "4region", NMR_field = 200, NMR_rotation = 6800)

#To view the output of a single spectrum
View(Integralregions[[1]][["data"]][["Integral"]])

```

