library(cmna)
library(minpack.lm)
library(quadprog)
library(SOMnmR)
library(pracma)
library(ggplot2)
## set working directory

work.dir <- c("C:/Users/colochoh/Desktop/Plant_NMR/example")



## load standard files
setwd(paste(work.dir,"csv",sep = "/"))
files <- list.files()

files <- list.files()
files <- files[grep(".csv", files)]

# types Bruker, "coma", "tab"
spec <- read_raw_spec(files, filetype = "coma")

# create an empty CN file

ncdata <- mk_nc_data(spec)

# save and fill with CN data

write.csv(ncdata,"C:/Users/colochoh/Desktop/Plant_NMR/NC/NC_plant.csv", row.names = FALSE)

# load nc data
setwd(paste(work.dir,"nc",sep = "/"))
ncdata <- nc_data("C:/Users/colochoh/Desktop/Plant_NMR/NC/CNdata.csv")

# Calculate Molecular mixing model NMR Meth either MMMFixN or MMM, ecosys either Terr or Aqua
setwd(paste(work.dir,sep = "/"))

MMM <- region_calc2(spec, NMRmeth = "MMMFixN", ecosys= "Terr", cndata = ncdata)

write.csv(MMM, "MMM_R.csv")

###Obtain only integrals regions

MMMregions <- int_nmr(spec,NMRmeth = "MMM")

write.csv(unlist(MMMregions[[1]]$data$Integral), "MMMregions_test_trap.csv")

###plot
setwd(paste(work.dir,"plot",sep = "/"))

plot_NMR(spec, NMRmeth = "MMM",file.output = TRUE , use.tiff = FALSE)


work.dir <- c("N:/Documents/Data/carol_aldair/test_sample")


## load standard files
setwd(paste(work.dir,"csv",sep = "/"))
files <- list.files()

files <- list.files()
files <- files[grep(".csv", files)]

# types Bruker, "coma", "tab"
spec <- read_raw_spec(files, filetype = "coma")

MMMregions <- int_nmr(spec, NMRmeth = "MMM")

test <- data.frame(unlist(MMMregions[[1]]$data$Integral))

