# How to use the package SOMnmR

[Back to the Main page](README.md)

## Basic example

Here we want to just load the data (NMR spectra) and do a simple integration, which will return the integrals and the SSB integrals.
This means that the output is not yet SSB corrected.

```bash
# load necessary packages
  library(devtools)

# load data

# integrate data
# Here we can choose between "4region", "Bonanomi" and "MMM". I will chose 4 region.



```

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
  library(devtools)

# load data

# integrate data
# Here we can choose between "4region", "Bonanomi" and "MMM". I will chose 4 region.



```

