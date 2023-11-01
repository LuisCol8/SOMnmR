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

# Integrate and correct for the SSB of the whole data set, NMRmeth can be "4region", "Bonanomi" and "MMM".
# I will chose 4 region.
# For this you will need to input the magnetic field of your NMR machine and the spinning frequency of the probe.
# In this example we used an NMR of 200 MHz and a spinning frequency of 6800 Hz.

IntegralSSBc <- region_calc(spec, NMRmeth = "4region", NMR_field = 200, NMR_rotation = 6800)

#To view the output of a single spectrum
View(IntegralSSBc[[1]])

#If you want to generate a table with all the results insead of having each in a list you can try this
Integrall <- NULL
for (i in 1:length(IntegralSSBc)) {
  Integral <- data.frame(IntegralSSBc[[i]])
  
  Integrall <- rbind(Integrall, Integral)
}
View(Integrall)

```

The output of region_calc  will be a list, you will find a table as follows

![image](https://github.com/LuisCol8/SOMnmR/assets/35764330/467b5fbd-ba28-4914-90b6-a5c6006d8ed0)


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

# We will now create an empty CN file, which we will save and open in excel and fill CN data.
#It is very important that you do not reorganize the Table, as it has to be in the same order as the input spectra.

ncdata <- mk_nc_data(spec)
write.csv(ncdata,"C:/Documents/Data/Experiment_NMR/NC/NC_table.csv", row.names = FALSE)

# Now we will load the NC table and use it for the MMM fitting
ncdata <- nc_data("C:/Documents/Data/Experiment_NMR/NC/NC_table.csv")

# We will now, integrate and correct for the SSB, and perform the fittin of the MMM for the whole data set
# For this the NMRMeth has to be "MMM", and I will choose to fix the NC data, using the FixNC = TRUE)
# You will need to input the magnetic field of your NMR machine and the spinning frequency of the probe.
# and the NCdata variable.
# In this example we used an NMR of 200 MHz and a spinning frequency of 6800 Hz.

MMM_TB_fix <- region_calc(spec, NMRmeth = "MMM", ecosys= "Terr_Baldock", cndata = ncdata, FixNC = TRUE)

#To view the output of a single spectrum
View(MMM_TB_fix)

```
The output of region_calc  will be a table, containing each sample and its succesive fittings, organized according to its lowest R-factor as follows

![image](https://github.com/LuisCol8/SOMnmR/assets/35764330/5e296275-a69b-47cc-9be8-bc4e5cdbbc31)

## Using a modified table for the MMM fitting

We can modifiy the input table as needed, ofcourse having in mind that the values should have physical meaning and be correctly reported in publications. The example follows:

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

# We will now create an empty CN file, which we will save and open in excel and fill CN data.
#It is very important that you do not reorganize the Table, as it has to be in the same order as the input spectra.

ncdata <- mk_nc_data(spec)
write.csv(ncdata,"C:/Documents/Data/Experiment_NMR/NC/NC_table.csv", row.names = FALSE)

# Now we will load the NC table and use it for the MMM fitting
ncdata <- nc_data("C:/Documents/Data/Experiment_NMR/NC/NC_table.csv")

# Now we will create a modified standard table. For this first we will load a table we want to modify.
# For instance the Terrestrial ecosystem table from Nelson et al.
modstd <-std_nmr(ecosys= "Terr_Nelson")

# We will now write the table as csv and open it in excel, and modify it.
# Currently, you can modify the values of the table, but not the names of the columns.
# If you add a new column or modify the name of the column it will not work.
write.csv(modstd[[1]],"C:/Documents/Data/Experiment_NMR/modstd.csv", row.names = FALSE)

# Now I will load this new modified table as a list.
modstd <- list(read.csv("C:/Users/zak69953/Downloads/Gabriela/modstd.csv"))

# I will now add the modified table as the table to be used for the Standards (mod_std = modstd)
# and will change the ecosys to "mod".
# We will now, integrate and correct for the SSB, and perform the fittin of the MMM for the whole data set
# For this the NMRMeth has to be "MMM", and I will choose to fix the NC data, using the FixNC = TRUE)
# You will need to input the magnetic field of your NMR machine and the spinning frequency of the probe.
# and the NCdata variable.
# In this example we used an NMR of 200 MHz and a spinning frequency of 6800 Hz.

MMM_TB_fix_mod <- region_calc(spec, NMRmeth = "MMM", ecosys= "mod", cndata = ncdata, FixNC = TRUE, mod_std = modstd,
 NMR_field = 200, NMR_rotation = 6800)

#To view the output of a single spectrum
View(MMM_TB_fix)

```

## Make a plot

There is a quick plot function that separates the regions of the spectrum.

```bash
# load necessary packages
  library(ggplot2)
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

# We will not change folder to one named "plot" and create the plots.
# ou can do the MMM regions by adding NMRmeth = "MMM", if it is empty it will do it as the "4region".

setwd(paste(work.dir,"plot",sep = "/"))

plot_NMR(spec, file.output = TRUE , use.tiff = FALSE)


```
