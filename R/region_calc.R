#' Functional groups calculation
#'
#' This function loads, integates and calculates the functional group distribution from the raw spectra.
#' Output is a list with the raw data, integrals and corrected spectra.
#' @param raw.spec List of files already loaded with read_raw_spec()
#' @param batch_nmr Vector with file names, default
#' @param NMRmeth Regions to be integrated.
#' Default is spinning side bands, other methods available include: Bonanomi ("Bonanomi") and Molecular mixing model ("MMM").
#' @param ecosys Standards to be used for the MMM, can be Terrestrial("Terr") or Aquatic ("Aqua")
#' @param nc_constrain Boolean for the fitting of the MMM, TRUE for N/C constrained fit, FALS for non-constrained. Default is FALSE.
#' @param ncdata Initial correction and normalization parameters
#' @keywords
#' @export
#' @examples

region_calc <- function (batch_nmr = NULL, file = NULL, NMRmeth = NULL, ecosys=NULL,
                          cn_constrain=FALSE, cndata = NULL) {

  if (is.null(batch_nmr)) {

    if (is.null(file)) {

      stop("Please provide either a raw spectrum from the read_raw_spec function or a raw .xmu file")

    } else {

      ## read the raw spectra from a file list
      batch.nmr <- read_raw_spec(file = file)

    }

  } else {

    batch.nmr <- batch_nmr

  }
  if (is.null(NMRmeth)) {
    ## loop to process all samples
    NMR.end <- NULL

    batch.nmr <- int_nmr(batch.nmr)

    for (i in 1:length(batch.nmr)) {
      # file.name <- batch.nmr[[i]]$name
      Integral <- c(batch.nmr[[i]]$data$Integral)


      ##carboxyl C calculation
      carboxyl <-  setNames(data.frame(sum(2*sum(Integral$normalized.Int[30:33]), sum(Integral$normalized.Int[21:24]), -sum(Integral$normalized.Int[3:6]))), c("Carboxyl"))
      ##Aryl C calculation
      aryl <-  setNames(data.frame(sum(2*sum(Integral$normalized.Int[27:29]), sum(Integral$normalized.Int[18:20]), -sum(Integral$normalized.Int[1:2]))), c("Aryl"))
      ##O-Alkyl C calculation
      oalkyl <-  setNames(data.frame(sum(2*sum(Integral$normalized.Int[4:8]), sum(Integral$normalized.Int[13:17]), -sum(Integral$normalized.Int[31:33]))), c("O-Alkyl"))
      ##Alkyl C calculation
      alkyl <-  setNames(data.frame(sum(2*sum(Integral$normalized.Int[1:3]), sum(Integral$normalized.Int[10:12]), -sum(Integral$normalized.Int[28:30]))), c("Alkyl"))
      ##Put all together
      #NMR.end[[i]] <- list(file.name = file.name, data = data.frame(carboxyl, aryl, oalkyl, alkyl))
      NMR.end[[i]] <- data.frame(carboxyl, aryl, oalkyl, alkyl)
    }
  } else if (NMRmeth == "Bonanomi") {
    ## loop to process all samples
    NMR.end <- NULL

    for (i in 1:length(batch.nmr)) {
      file.name <- batch.nmr[[i]]$name
      sample <- list(batch.nmr[[i]])
      Integral <- int_nmr (raw.spec = sample, NMRmeth = NMRmeth)

      ##carboxyl C calculation
      carboxyl <-  setNames(data.frame(sum(2*sum(Integral[30:33,2]), sum(Integral[21:24,2]), -sum(Integral[3:6,2]))), c("Carboxyl"))
      ##Aryl C calculation
      aryl <-  setNames(data.frame(sum(2*sum(Integral[27:29,2]), sum(Integral[18:20,2]), -sum(Integral[1:2,2]))), c("Aryl"))
      ##O-Alkyl C calculation
      oalkyl <-  setNames(data.frame(sum(2*sum(Integral[4:8,2]), sum(Integral[13:17,2]), -sum(Integral[31:33,2]))), c("O-Alkyl"))
      ##Alkyl C calculation
      alkyl <-  setNames(data.frame(sum(2*sum(Integral[1:3,2]), sum(Integral[10:12,2]), -sum(Integral[28:30,2]))), c("Alkyl"))
      ##Put all together
      NMR.end[[i]] <- data.frame(file.name = file.name, data = data.frame(carboxyl, aryl, oalkyl, alkyl))
      write.table(NMR.end[[i]], "output.csv", append = TRUE, sep = ",", dec = ".",
                  col.names = TRUE, row.names = FALSE)
      sample <- NULL
    }
  } else if (NMRmeth == "MMM") {

    ## loop to process all samples
    NMR.end <- NULL
    raw.spec.end <- NULL
    batch.nmr <- int_nmr(batch.nmr, NMRmeth = "MMM-SSB")

    nmrmerge <- NULL

    for (i in 1:length(batch.nmr)) {

      raw.spec.end[[i]] <-  batch.nmr[[i]]
      NCval <- as.numeric(cndata[[i]]$NC)
      samplename <- batch.nmr[[i]]$name
      sampleraw.spec <- batch.nmr[[i]]$data$raw.spec
      sampleintegral <- as.data.frame(batch.nmr[[i]]$data$Integral)

      Alkyl <- ifelse(sampleintegral$normalized.Int[4] - sampleintegral$normalized.Int[15] < 0,0, sampleintegral$normalized.Int[4] - sampleintegral$normalized.Int[15])

      N_Alkyl_Methoxyl <- ifelse(sampleintegral$normalized.Int[5] < 0,0, sampleintegral$normalized.Int[5])

      O_Alkyl <- ifelse(sampleintegral$normalized.Int[6] < 0,0, sampleintegral$normalized.Int[6])

      Di_O_Alkyl <- ifelse(sampleintegral$normalized.Int[7] < 0,0, sampleintegral$normalized.Int[7] + sampleintegral$normalized.Int[1] + sampleintegral$normalized.Int[12])

      Aromatic <- ifelse(sampleintegral$normalized.Int[8] < 0,0, sampleintegral$normalized.Int[8] + sampleintegral$normalized.Int[2] + sampleintegral$normalized.Int[13])

      Phenolic <- ifelse(sampleintegral$normalized.Int[9] < 0,0, sampleintegral$normalized.Int[9] + sampleintegral$normalized.Int[3] + sampleintegral$normalized.Int[14])

      Amide_Carboxylic <- ifelse(sampleintegral$normalized.Int[10] < 0,0, sampleintegral$normalized.Int[10] + sampleintegral$normalized.Int[3] + 2*sampleintegral$normalized.Int[15])

      Ketone <- ifelse(sampleintegral$normalized.Int[11] < 0,0, sampleintegral$normalized.Int[11])

      Amide_to_Ketone <- c(Amide_Carboxylic + Ketone)

      sampleintegralend <- data.frame(Alkyl, N_Alkyl_Methoxyl, O_Alkyl, Di_O_Alkyl, Aromatic, Phenolic, Amide_to_Ketone)

      sampleintegralend <- t(sampleintegralend)

      sampleintegralend <- rbind(NCval,sampleintegralend)

      raw.spec.end[[i]] <- list("name" = samplename, "data" = list("raw.spec" = sampleraw.spec,"Integral" = sampleintegralend))

      if (ecosys == "Terr") {

        stdmat <- std_nmr(ecosys = "Terr")
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30, NMRmeth = "MMM")

        ## Elementar ratios (relative to C)

        Cmol <- as.numeric(NMR.end$Carbohydrates*1+NMR.end$Protein*1+NMR.end$Lignin*1+
                             NMR.end$Lipid*1+NMR.end$Carbonyl*1+NMR.end$Char*1)

        Nmol <- as.numeric(NMR.end$Carbohydrates*0+NMR.end$Protein*0.275+NMR.end$Lignin*0+
                             NMR.end$Lipid*0+NMR.end$Carbonyl*0+NMR.end$Char*0.004)

        Hmol <- as.numeric(NMR.end$Carbohydrates*1.667+NMR.end$Protein*1.101+NMR.end$Lignin*1.238+
                             NMR.end$Lipid*1.941+NMR.end$Carbonyl*1+NMR.end$Char*0.452)

        Omol <- as.numeric(NMR.end$Carbohydrates*0.833+NMR.end$Protein*0.155+NMR.end$Lignin*0.429+
                             NMR.end$Lipid*0.235+NMR.end$Carbonyl*2+NMR.end$Char*0.405)

        Cwgt <- as.numeric(NMR.end$Carbohydrates*1*12.0107+NMR.end$Protein*1*12.0107+NMR.end$Lignin*1*12.0107+
                             NMR.end$Lipid*1*12.0107+NMR.end$Carbonyl*1*12.0107+NMR.end$Char*1*12.0107)

        Nwgt <- as.numeric(NMR.end$Carbohydrates*0*14.0067+NMR.end$Protein*0.275*14.0067+NMR.end$Lignin*0*14.0067+
                             NMR.end$Lipid*0*14.0067+NMR.end$Carbonyl*0*14.0067+NMR.end$Char*0.004*14.0067)

        Hwgt <- as.numeric(NMR.end$Carbohydrates*1.667*1.00794+NMR.end$Protein*1.101*1.00794+NMR.end$Lignin*1.238*1.00794+
                             NMR.end$Lipid*1.941*1.00794+NMR.end$Carbonyl*1*1.00794+NMR.end$Char*0.452*1.00794)

        Owgt <- as.numeric(NMR.end$Carbohydrates*0.833*15.994+NMR.end$Protein*0.155*15.994+NMR.end$Lignin*0.429*15.994+
                             NMR.end$Lipid*0.235*15.994+NMR.end$Carbonyl*2*15.994+NMR.end$Char*0.405*15.994)

        swgt <- c(Cwgt + Nwgt + Hwgt +Owgt)

        NOSC <- as.numeric(4+((2*Omol+3*Nmol-1*Hmol-4*Cmol)/Cmol))

        ## Final result

        NMR.end <- cbind(NMR.end, Cmol, Nmol, Hmol, Omol, Cwgt/swgt, Nwgt/swgt, Hwgt/swgt, Owgt/swgt, NOSC)

      } else if (ecosys == "Aqua") {

        stdmat <- std_nmr(ecosys = "Aqua")
        NMR.end <- fit_athena(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30,  NMRmeth = "MMM")

        ## Elementar ratios (relative to C)

        Cmol <- as.numeric(NMR.end$Carbohydrates*1+NMR.end$Protein*1+NMR.end$Lignin*1+
                             NMR.end$Lipid*1+NMR.end$Carbonyl*1+NMR.end$Char*1)

        Nmol <- as.numeric(NMR.end$Carbohydrates*0+NMR.end$Protein*0.275+NMR.end$Lignin*0+
                             NMR.end$Lipid*0+NMR.end$Carbonyl*0+NMR.end$Char*0.004)

        Hmol <- as.numeric(NMR.end$Carbohydrates*1.667+NMR.end$Protein*1.101+NMR.end$Lignin*1.238+
                             NMR.end$Lipid*1.941+NMR.end$Carbonyl*1+NMR.end$Char*0.452)

        Omol <- as.numeric(NMR.end$Carbohydrates*0.833+NMR.end$Protein*0.155+NMR.end$Lignin*0.429+
                             NMR.end$Lipid*0.235+NMR.end$Carbonyl*2+NMR.end$Char*0.405)

        Cwgt <- as.numeric(NMR.end$Carbohydrates*1*12.0107+NMR.end$Protein*1*12.0107+NMR.end$Lignin*1*12.0107+
                             NMR.end$Lipid*1*12.0107+NMR.end$Carbonyl*1*12.0107+NMR.end$Char*1*12.0107)

        Nwgt <- as.numeric(NMR.end$Carbohydrates*0*14.0067+NMR.end$Protein*0.275*14.0067+NMR.end$Lignin*0*14.0067+
                             NMR.end$Lipid*0*14.0067+NMR.end$Carbonyl*0*14.0067+NMR.end$Char*0.004*14.0067)

        Hwgt <- as.numeric(NMR.end$Carbohydrates*1.667*1.00794+NMR.end$Protein*1.101*1.00794+NMR.end$Lignin*1.238*1.00794+
                             NMR.end$Lipid*1.941*1.00794+NMR.end$Carbonyl*1*1.00794+NMR.end$Char*0.452*1.00794)

        Owgt <- as.numeric(NMR.end$Carbohydrates*0.833*15.994+NMR.end$Protein*0.155*15.994+NMR.end$Lignin*0.429*15.994+
                             NMR.end$Lipid*0.235*15.994+NMR.end$Carbonyl*2*15.994+NMR.end$Char*0.405*15.994)

        swgt <- c(Cwgt + Nwgt + Hwgt +Owgt)

        NOSC <- as.numeric(4+((2*Omol+3*Nmol-1*Hmol-4*Cmol)/Cmol))

        ## Final result

        NMR.end <- cbind(NMR.end, Cmol, Nmol, Hmol, Omol, Cwgt/swgt, Nwgt/swgt, Hwgt/swgt, Owgt/swgt, NOSC)
      }
    }
    ## return the corrected spectra list
  } else if (NMRmeth == "MMMFixN") {

    ## loop to process all samples
    NMR.end <- NULL
    raw.spec.end <- NULL
    batch.nmr <- int_nmr(batch.nmr, NMRmeth = "MMM-SSB")

    nmrmerge <- NULL

    for (i in 1:length(batch.nmr)) {

      raw.spec.end[[i]] <-  batch.nmr[[i]]
      NCval <- as.numeric(cndata[[i]]$NC)
      samplename <- batch.nmr[[i]]$name
      sampleraw.spec <- batch.nmr[[i]]$data$raw.spec
      sampleintegral <- as.data.frame(batch.nmr[[i]]$data$Integral)

      Alkyl <- ifelse(sampleintegral$normalized.Int[4] - sampleintegral$normalized.Int[15] < 0,0, sampleintegral$normalized.Int[4] - sampleintegral$normalized.Int[15])

      N_Alkyl_Methoxyl <- ifelse(sampleintegral$normalized.Int[5] < 0,0, sampleintegral$normalized.Int[5])

      O_Alkyl <- ifelse(sampleintegral$normalized.Int[6] < 0,0, sampleintegral$normalized.Int[6])

      Di_O_Alkyl <- ifelse(sampleintegral$normalized.Int[7] < 0,0, sampleintegral$normalized.Int[7] + sampleintegral$normalized.Int[1] + sampleintegral$normalized.Int[12])

      Aromatic <- ifelse(sampleintegral$normalized.Int[8] < 0,0, sampleintegral$normalized.Int[8] + sampleintegral$normalized.Int[2] + sampleintegral$normalized.Int[13])

      Phenolic <- ifelse(sampleintegral$normalized.Int[9] < 0,0, sampleintegral$normalized.Int[9] + sampleintegral$normalized.Int[3] + sampleintegral$normalized.Int[14])

      Amide_Carboxylic <- ifelse(sampleintegral$normalized.Int[10] < 0,0, sampleintegral$normalized.Int[10] + sampleintegral$normalized.Int[3] + 2*sampleintegral$normalized.Int[15])

      Ketone <- ifelse(sampleintegral$normalized.Int[11] < 0,0, sampleintegral$normalized.Int[11])

      Amide_to_Ketone <- c(Amide_Carboxylic + Ketone)

      sampleintegralend <- data.frame(Alkyl, N_Alkyl_Methoxyl, O_Alkyl, Di_O_Alkyl, Aromatic, Phenolic, Amide_to_Ketone)

      sampleintegralend <- t(sampleintegralend)

      sampleintegralend <- rbind(NCval,sampleintegralend)

      raw.spec.end[[i]] <- list("name" = samplename, "data" = list("raw.spec" = sampleraw.spec,"Integral" = sampleintegralend))

      if (ecosys == "Terr") {

        stdmat <- std_nmr(ecosys = "Terr")
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30, NMRmeth = "MMMFixN")

        ## Elementar ratios (relative to C)

        Cmol <- as.numeric(NMR.end$Carbohydrates*1+NMR.end$Protein*1+NMR.end$Lignin*1+
                             NMR.end$Lipid*1+NMR.end$Carbonyl*1+NMR.end$Char*1)

        Nmol <- as.numeric(NMR.end$Carbohydrates*0+NMR.end$Protein*0.275+NMR.end$Lignin*0+
                             NMR.end$Lipid*0+NMR.end$Carbonyl*0+NMR.end$Char*0.004)

        Hmol <- as.numeric(NMR.end$Carbohydrates*1.667+NMR.end$Protein*1.101+NMR.end$Lignin*1.238+
                             NMR.end$Lipid*1.941+NMR.end$Carbonyl*1+NMR.end$Char*0.452)

        Omol <- as.numeric(NMR.end$Carbohydrates*0.833+NMR.end$Protein*0.155+NMR.end$Lignin*0.429+
                             NMR.end$Lipid*0.235+NMR.end$Carbonyl*2+NMR.end$Char*0.405)

        Cwgt <- as.numeric(NMR.end$Carbohydrates*1*12.0107+NMR.end$Protein*1*12.0107+NMR.end$Lignin*1*12.0107+
                             NMR.end$Lipid*1*12.0107+NMR.end$Carbonyl*1*12.0107+NMR.end$Char*1*12.0107)

        Nwgt <- as.numeric(NMR.end$Carbohydrates*0*14.0067+NMR.end$Protein*0.275*14.0067+NMR.end$Lignin*0*14.0067+
                             NMR.end$Lipid*0*14.0067+NMR.end$Carbonyl*0*14.0067+NMR.end$Char*0.004*14.0067)

        Hwgt <- as.numeric(NMR.end$Carbohydrates*1.667*1.00794+NMR.end$Protein*1.101*1.00794+NMR.end$Lignin*1.238*1.00794+
                             NMR.end$Lipid*1.941*1.00794+NMR.end$Carbonyl*1*1.00794+NMR.end$Char*0.452*1.00794)

        Owgt <- as.numeric(NMR.end$Carbohydrates*0.833*15.994+NMR.end$Protein*0.155*15.994+NMR.end$Lignin*0.429*15.994+
                             NMR.end$Lipid*0.235*15.994+NMR.end$Carbonyl*2*15.994+NMR.end$Char*0.405*15.994)

        swgt <- c(Cwgt + Nwgt + Hwgt +Owgt)

        NOSC <- as.numeric(4+((2*Omol+3*Nmol-1*Hmol-4*Cmol)/Cmol))

        ## Final result

        NMR.end <- cbind(NMR.end, Cmol, Nmol, Hmol, Omol, Cwgt/swgt, Nwgt/swgt, Hwgt/swgt, Owgt/swgt, NOSC)

      } else if (ecosys == "Aqua") {

        stdmat <- std_nmr(ecosys = "Aqua")
        NMR.end <- fit_athena(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30,  NMRmeth = "MMMFixN")

        ## Elementar ratios (relative to C)

        Cmol <- as.numeric(NMR.end$Carbohydrates*1+NMR.end$Protein*1+NMR.end$Lignin*1+
                             NMR.end$Lipid*1+NMR.end$Carbonyl*1+NMR.end$Char*1)

        Nmol <- as.numeric(NMR.end$Carbohydrates*0+NMR.end$Protein*0.275+NMR.end$Lignin*0+
                             NMR.end$Lipid*0+NMR.end$Carbonyl*0+NMR.end$Char*0.004)

        Hmol <- as.numeric(NMR.end$Carbohydrates*1.667+NMR.end$Protein*1.101+NMR.end$Lignin*1.238+
                             NMR.end$Lipid*1.941+NMR.end$Carbonyl*1+NMR.end$Char*0.452)

        Omol <- as.numeric(NMR.end$Carbohydrates*0.833+NMR.end$Protein*0.155+NMR.end$Lignin*0.429+
                             NMR.end$Lipid*0.235+NMR.end$Carbonyl*2+NMR.end$Char*0.405)

        Cwgt <- as.numeric(NMR.end$Carbohydrates*1*12.0107+NMR.end$Protein*1*12.0107+NMR.end$Lignin*1*12.0107+
                             NMR.end$Lipid*1*12.0107+NMR.end$Carbonyl*1*12.0107+NMR.end$Char*1*12.0107)

        Nwgt <- as.numeric(NMR.end$Carbohydrates*0*14.0067+NMR.end$Protein*0.275*14.0067+NMR.end$Lignin*0*14.0067+
                             NMR.end$Lipid*0*14.0067+NMR.end$Carbonyl*0*14.0067+NMR.end$Char*0.004*14.0067)

        Hwgt <- as.numeric(NMR.end$Carbohydrates*1.667*1.00794+NMR.end$Protein*1.101*1.00794+NMR.end$Lignin*1.238*1.00794+
                             NMR.end$Lipid*1.941*1.00794+NMR.end$Carbonyl*1*1.00794+NMR.end$Char*0.452*1.00794)

        Owgt <- as.numeric(NMR.end$Carbohydrates*0.833*15.994+NMR.end$Protein*0.155*15.994+NMR.end$Lignin*0.429*15.994+
                             NMR.end$Lipid*0.235*15.994+NMR.end$Carbonyl*2*15.994+NMR.end$Char*0.405*15.994)

        swgt <- c(Cwgt + Nwgt + Hwgt +Owgt)

        NOSC <- as.numeric(4+((2*Omol+3*Nmol-1*Hmol-4*Cmol)/Cmol))

        ## Final result

        NMR.end <- cbind(NMR.end, Cmol, Nmol, Hmol, Omol, Cwgt/swgt, Nwgt/swgt, Hwgt/swgt, Owgt/swgt, NOSC)
      }
    }
  }
  citation  <- setNames(data.frame(matrix(ncol = 1, nrow= nrow(NMR.end))), c("Plz cite this work as an incentive to its curation"))
  NMR.end <- cbind(NMR.end,citation)
  return(NMR.end)
}
