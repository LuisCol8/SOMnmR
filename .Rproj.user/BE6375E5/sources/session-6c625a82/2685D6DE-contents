library(cmna)
library(minpack.lm)
library(quadprog)
library(SOMnmR)
library(pracma)
library(ggplot2)
## set working directory

# load devtools
library(devtools)
library(usethis)

# load SOMnmR
install_github("LuisCol8/SOMnmR")

work.dir <- c("N:/Documents/Data/Luis_Almeida/Antartica/Paper_Antarctic/NMR/NMR_batch/")


work.dir <- c("N:/Documents/Data/Jarosch_Klaus/NMR/")

work.dir <- c("N:/Documents/Data/Joerg/C-NEXAFS/NMR/")

## load standard files
setwd(paste(work.dir,"spectra",sep = "/"))
files <- list.files()

files <- list.files()
files <- files[grep(".txt", files)]

# types Bruker, "coma", "tab"
spec <- read_raw_spec(files, filetype = "Bruker")

# create an empty CN file

ncdata <- mk_nc_data(spec)

# save and fill with CN data

write.csv(ncdata,"N:/Documents/Data/Jarosch_Klaus/NMR/nc/NC_TPI.csv", row.names = FALSE)

# load nc data
setwd(paste(work.dir,"nc",sep = "/"))
ncdata <- nc_data("N:/Documents/Data/Jarosch_Klaus/NMR/nc/NC_TPI.csv")

# Calculate Molecular mixing model NMR Meth either MMMFixN or MMM, ecosys either Terr or Aqua
setwd(paste(work.dir,sep = "/"))

modstd <-std_nmr(ecosys= "Terr_Nelson")


MMM_TN_fix <- region_calc(spec, NMRmeth = "MMM", ecosys= "mod", cndata = ncdata,mod_std = modstd)

MMM_TN_fix <- region_calc(spec, NMRmeth = "MMMFixN", ecosys= "Terr_Nelson", cndata = ncdata, stats = TRUE)

MMM_TB_fix <- region_calc(spec, NMRmeth = "MMMFixN", ecosys= "Terr_Baldock", cndata = ncdata)

MMM_AN_fix <- region_calc(spec, NMRmeth = "MMMFixN", ecosys= "Aqua_Nelson", cndata = ncdata)

MMM_AB_fix <- region_calc(spec, NMRmeth = "MMMFixN", ecosys= "Aqua_Baldock", cndata = ncdata)

write.csv(MMM_TN_fix, "MMM_TN_fix.csv")
write.csv(MMM_TB_fix, "MMM_TB_fix.csv")
write.csv(MMM_AN_fix, "MMM_AN_fix.csv")
write.csv(MMM_AB_fix, "MMM_AB_fix.csv")


MMM_TN_free <- region_calc(spec, NMRmeth = "MMM", ecosys= "Terr_Nelson", cndata = ncdata)

MMM_TB_free <- region_calc(spec, NMRmeth = "MMM", ecosys= "Terr_Baldock", cndata = ncdata)

MMM_AN_free <- region_calc(spec, NMRmeth = "MMM", ecosys= "Aqua_Nelson", cndata = ncdata)

MMM_AB_free <- region_calc(spec, NMRmeth = "MMM", ecosys= "Aqua_Baldock", cndata = ncdata)

write.csv(MMM_TN_free, "MMM_TN_free.csv")
write.csv(MMM_TB_free, "MMM_TB_free.csv")
write.csv(MMM_AN_free, "MMM_AN_free.csv")
write.csv(MMM_AB_free, "MMM_AB_free.csv")

write.csv(MMM, "MMM_2.csv")

###Obtain only integrals regions

MMMregions_nc <- int_nmr(spec, NMRmeth = "MMM-SSB", SSBcorr = FALSE)

MMMregions <- int_nmr(spec, NMRmeth = "MMM-SSB", SSBcorr = TRUE)

Integrall <- NULL
for (i in 1:length(MMMregions)) {
  Integral <- c(unlist(MMMregions[[i]]$name),unlist(MMMregions[[i]]$data$Integral))

  Integrall <- rbind(Integrall, Integral)
}

write.csv(Integrall, "MMMregions_SSB_7regions.csv")

MMMregions <- int_nmr(spec)

integrall= do.call(rbind, MMMregions2)

write.csv(integrall, "MMMregions_SSB_4regions.csv")

Integrall <- NULL
for (i in 1:length(MMMregions)) {
  Integral <- c(MMMregions[[i]]$name,MMMregions[[i]]$data$Integral)

  Integrall <- rbind(Integrall, Integral)
  }

# NC extract
ncdataa <- NULL
for (i in 1:length(ncdata)) {
  ncdat <- c(ncdata[[i]]$name,ncdata[[i]]$NC)

  ncdataa <- rbind(ncdataa, ncdat)
}


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

setwd(paste(work.dir,"plot",sep = "/"))

plot_NMR(spec, NMRmeth = "Bonanomi",file.output = TRUE,use.tiff = TRUE)


Integrall <- NULL
for (i in 1:length(MMregions4)) {
  Integral <- c(MMregions4[[i]])
  
  Integrall <- rbind(Integrall, Integral)
}
write.csv(Integrall, "MMMregions_SSB_4regions.csv")
