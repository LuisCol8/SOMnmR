#' Functional groups calculation
#'
#' This function loads, integates and calculates the functional group distribution from the raw spectra.
#' Output is a list with the raw data, integrals and corrected spectra.
#' @param raw.spec List of files already loaded with read_raw_spec()
#' @param batch_nmr Vector with file names, default
#' @param NMRmeth Regions to be integrated.
#' Default is spinning side bands, other methods available include: Bonanomi ("Bonanomi") and Molecular mixing model ("MMM").
#' @param ecosys Standards to be used for the MMM, can be Terrestrial("Terr_Nelson" or "Terr_Baldock") or Aquatic ("Aqua_Nelson" or "Aqua_Baldock")
#' @param ncdata Initial correction and normalization parameters
#' @keywords
#' @export
#' @examples

region_calc <- function (batch_nmr = NULL, file = NULL, NMRmeth = NULL, ecosys=NULL,
                          cndata = NULL, mod_std = NULL, stats = FALSE) {

  if (is.null(batch_nmr)) {

    if (is.null(file)) {

      stop("Please provide either a raw spectrum from the read_raw_spec function or a raw NMR file")

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
      file.name <- batch.nmr[[i]]$name
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
      integral.end <- data.frame(carboxyl, aryl, oalkyl, alkyl)
      
      norm <- sum(integral.end)
      normalized.Int <- (integral.end/norm)*100
      integral.end <- data.frame(normalized.Int)
      #integral.end <- data.frame(file.name,integral.end)
      NMR.end[[i]] <- data.frame(file.name, integral.end)
      
    }
  } else if (NMRmeth == "Bonanomi") {
    ## loop to process all samples
    NMR.end <- NULL

    for (i in 1:length(batch.nmr)) {
      file.name <- batch.nmr[[i]]$name
      sample <- list(batch.nmr[[i]])
      Integral <- int_nmr (raw.spec = sample, NMRmeth = NMRmeth)

      ##carboxyl C calculation
      carboxyl <-  setNames(data.frame(sum(sum(Integral[30:33,2]), sum(Integral[21:24,2]), -2*sum(Integral[3:6,2]))), c("Carboxyl"))
      ##Aryl C calculation
      aryl <-  setNames(data.frame(sum(sum(Integral[27:29,2]), sum(Integral[18:20,2]), -2*sum(Integral[1:2,2]))), c("Aryl"))
      ##O-Alkyl C calculation
      oalkyl <-  setNames(data.frame(sum(sum(Integral[4:8,2]), sum(Integral[13:17,2]), -2*sum(Integral[31:33,2]))), c("O-Alkyl"))
      ##Alkyl C calculation
      alkyl <-  setNames(data.frame(sum(sum(Integral[1:3,2]), sum(Integral[10:12,2]), -2*sum(Integral[28:30,2]))), c("Alkyl"))
      ##Put all together
      NMR.end[[i]] <- data.frame(file.name = file.name, data = data.frame(carboxyl, aryl, oalkyl, alkyl))

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
      
      ##Alklyl C calculation
      Alkyl <-  setNames(data.frame(sum(2*sum(sampleintegral[1:3,1]), sum(sampleintegral[10:12,1]), -sum(sampleintegral[28:30,1]))), c("Alkyl"))
      
      ##N_Alkyl_Methoxyl C calculation
      N_Alkyl_Methoxyl <-  setNames(data.frame(sum(2*sum(sampleintegral[4:5,1]), sum(sampleintegral[13:14,1]), -sum(sampleintegral[31:32,1]))), c("N_Alkyl_Methoxyl"))
      
      ##O-Alkyl C calculation
      O_Alkyl <-  setNames(data.frame(sum(2*sum(sampleintegral[6:7,1]), sum(sampleintegral[15:16,1]), -sum(sampleintegral[33:33,1]))), c("O-Alkyl"))
      
      ##Di_O_Alkyl C calculation
      Di_O_Alkyl <-  setNames(data.frame(sum(2*sum(sampleintegral[8:8,1]), sum(sampleintegral[17:17,1]))), c("Di_O_Alkyl"))
      
      ##Aromatic C calculation
      Aromatic <-  setNames(data.frame(sum(2*sum(sampleintegral[27:28,1]), sum(sampleintegral[18:19,1]), -sum(sampleintegral[1:1,1]))), c("Aromatic"))
      
      ##Phenolic C calculation
      Phenolic <-  setNames(data.frame(sum(2*sum(sampleintegral[29:29,1]), sum(sampleintegral[20:20,1]), -sum(sampleintegral[2:2,1]))), c("Phenolic"))
      
      ##Amide_Carboxylic C calculation
      Amide_Carboxylic <-  setNames(data.frame(sum(2*sum(sampleintegral[30:32,1]), sum(sampleintegral[21:23,1]), -sum(sampleintegral[3:5,1]))), c("Amide_Carboxylic"))
      
      ##Ketone C calculation
      Ketone <-  setNames(data.frame(sum(2*sum(sampleintegral[33:33,1]), sum(sampleintegral[24:24,1]), -sum(sampleintegral[6:6,1]))), c("Ketone"))
      
      ##Put all together
      Amide_to_Ketone <- c(Amide_Carboxylic + Ketone)

      sampleintegraljoin <- data.frame(Alkyl, N_Alkyl_Methoxyl, O_Alkyl, Di_O_Alkyl, Aromatic, Phenolic, Amide_to_Ketone)

      sampleintegraljoin <- t(sampleintegraljoin)

      sampleintegralend <- rbind(NCval,sampleintegraljoin)

      raw.spec.end[[i]] <- list("name" = samplename, "data" = list("raw.spec" = sampleraw.spec,"Integral" = sampleintegralend))

      if (ecosys == "Terr_Nelson") {

        stdmat <- std_nmr(ecosys = "Terr_Nelson")
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30, NMRmeth = "MMM")

        ## Elementar ratios (relative to C)

        Cmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]+NMR.end$Protein*stdmat[[1]][9,1]+NMR.end$Lignin*stdmat[[1]][9,3]+
                             NMR.end$Lipid*stdmat[[1]][9,4]+NMR.end$Carbonyl*stdmat[[1]][9,5]+NMR.end$Char*stdmat[[1]][9,6])

        Nmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]+NMR.end$Protein*stdmat[[1]][10,1]+NMR.end$Lignin*stdmat[[1]][10,3]+
                             NMR.end$Lipid*stdmat[[1]][10,4]+NMR.end$Carbonyl*stdmat[[1]][10,5]+NMR.end$Char*stdmat[[1]][10,6])

        Hmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]+NMR.end$Protein*stdmat[[1]][11,1]+NMR.end$Lignin*stdmat[[1]][11,3]+
                             NMR.end$Lipid*stdmat[[1]][11,4]+NMR.end$Carbonyl*stdmat[[1]][11,5]+NMR.end$Char*stdmat[[1]][11,6])

        Omol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]+NMR.end$Protein*stdmat[[1]][12,1]+NMR.end$Lignin*stdmat[[1]][12,3]+
                             NMR.end$Lipid*stdmat[[1]][12,4]+NMR.end$Carbonyl*stdmat[[1]][12,5]+NMR.end$Char*stdmat[[1]][12,6])

        Cwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]*12.0107+NMR.end$Protein*stdmat[[1]][9,1]*12.0107+NMR.end$Lignin*stdmat[[1]][9,3]*12.0107+
                             NMR.end$Lipid*stdmat[[1]][9,4]*12.0107+NMR.end$Carbonyl*stdmat[[1]][9,5]*12.0107+NMR.end$Char*stdmat[[1]][9,6]*12.0107)

        Nwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]*14.0067+NMR.end$Protein*stdmat[[1]][10,1]*14.0067+NMR.end$Lignin*stdmat[[1]][10,3]*14.0067+
                             NMR.end$Lipid*stdmat[[1]][10,4]*14.0067+NMR.end$Carbonyl*stdmat[[1]][10,5]*14.0067+NMR.end$Char*stdmat[[1]][10,6]*14.0067)

        Hwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]*1.00794+NMR.end$Protein*stdmat[[1]][11,1]*1.00794+NMR.end$Lignin*stdmat[[1]][11,3]*1.00794+
                             NMR.end$Lipid*stdmat[[1]][11,4]*1.00794+NMR.end$Carbonyl*stdmat[[1]][11,5]*1.00794+NMR.end$Char*stdmat[[1]][11,6]*1.00794)

        Owgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]*15.994+NMR.end$Protein*stdmat[[1]][12,1]*15.994+NMR.end$Lignin*stdmat[[1]][12,3]*15.994+
                             NMR.end$Lipid*stdmat[[1]][12,4]*15.994+NMR.end$Carbonyl*stdmat[[1]][12,5]*15.994+NMR.end$Char*stdmat[[1]][12,6]*15.994)

        swgt <- c(Cwgt + Nwgt + Hwgt +Owgt)

        NOSC <- as.numeric(4+((2*Omol+3*Nmol-1*Hmol-4*Cmol)/Cmol))

        ## Final result

        NMR.end <- cbind(NMR.end, Cmol, Nmol, Hmol, Omol, Cwgt/swgt, Nwgt/swgt, Hwgt/swgt, Owgt/swgt, NOSC)
        
        ## Back calculated NMR results
        
        ##Alklyl C calculation
        Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][1,2]+NMR.end$Protein*stdmat[[1]][1,1]+NMR.end$Lignin*stdmat[[1]][1,3]+NMR.end$Lipid*stdmat[[1]][1,4]+
                               NMR.end$Carbonyl*stdmat[[1]][1,5]+NMR.end$Char*stdmat[[1]][1,6], c("Alkyl"))
        
        ##N_Alkyl_Methoxyl C calculation
        N_Alkyl_Methoxyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][2,2]+NMR.end$Protein*stdmat[[1]][2,1]+NMR.end$Lignin*stdmat[[1]][2,3]+NMR.end$Lipid*stdmat[[1]][2,4]+
                                          NMR.end$Carbonyl*stdmat[[1]][2,5]+NMR.end$Char*stdmat[[1]][2,6], c("N_Alkyl_Methoxyl"))
        
        ##O-Alkyl C calculation
        O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][3,2]+NMR.end$Protein*stdmat[[1]][3,1]+NMR.end$Lignin*stdmat[[1]][3,3]+NMR.end$Lipid*stdmat[[1]][3,4]+
                                 NMR.end$Carbonyl*stdmat[[1]][3,5]+NMR.end$Char*stdmat[[1]][3,6], c("O-Alkyl"))
        
        ##Di_O_Alkyl C calculation
        Di_O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][4,2]+NMR.end$Protein*stdmat[[1]][4,1]+NMR.end$Lignin*stdmat[[1]][4,3]+NMR.end$Lipid*stdmat[[1]][4,4]+
                                    NMR.end$Carbonyl*stdmat[[1]][4,5]+NMR.end$Char*stdmat[[1]][4,6], c("Di_O_Alkyl"))
        
        ##Aromatic C calculation
        Aromatic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][5,2]+NMR.end$Protein*stdmat[[1]][5,1]+NMR.end$Lignin*stdmat[[1]][5,3]+NMR.end$Lipid*stdmat[[1]][5,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][5,5]+NMR.end$Char*stdmat[[1]][5,6], c("Aromatic"))
        
        ##Phenolic C calculation
        Phenolic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][6,2]+NMR.end$Protein*stdmat[[1]][6,1]+NMR.end$Lignin*stdmat[[1]][6,3]+NMR.end$Lipid*stdmat[[1]][6,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][6,5]+NMR.end$Char*stdmat[[1]][6,6], c("Phenolic"))
        
        ##Amide_Carboxylic C calculation
        Amide_to_Ketone_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][7,2]+NMR.end$Protein*stdmat[[1]][7,1]+NMR.end$Lignin*stdmat[[1]][7,3]+NMR.end$Lipid*stdmat[[1]][7,4]+
                                         NMR.end$Carbonyl*stdmat[[1]][7,5]+NMR.end$Char*stdmat[[1]][7,6], c("Amide_to_Ketone"))
        
        sum_m <- setNames(Alkyl_m + N_Alkyl_Methoxyl_m + O_Alkyl_m + Di_O_Alkyl_m + Aromatic_m + Phenolic_m + Amide_to_Ketone_m, c("Sum"))
        
        sum_c <- sum(sampleintegraljoin)
        
        sampleintegraljoin <-rbind(sampleintegraljoin, sum_c)
        
        sample_stats <- data.frame(Alkyl_m, N_Alkyl_Methoxyl_m, O_Alkyl_m, Di_O_Alkyl_m, Aromatic_m, Phenolic_m, Amide_to_Ketone_m, sum_m)
        
        nmrrest <- NULL
        for (i in 1:nrow(sample_stats)) {
          nmrrestt <- c(sampleintegraljoin)
          nmrrest <- rbind(nmrrest, nmrrestt)
        }
        colnames(nmrrest) <-  c("Alkyl", "N_Alkyl_Methoxyl", "O-Alkyl", "Di_O_Alkyl", "Aromatic", "Phenolic", "Amide_to_Ketone", "Sum")
        
        sample_stats <- cbind(sample_stats,nmrrest)
        ssq_sample <- data.frame((sample_stats$Alkyl-sample_stats$Alkyl_m)^2, (sample_stats$N_Alkyl_Methoxyl-sample_stats$N_Alkyl_Methoxyl_m)^2,
                                 (sample_stats$`O-Alkyl` -sample_stats$O_Alkyl_m)^2, (sample_stats$Di_O_Alkyl-sample_stats$Di_O_Alkyl_m)^2, 
                                 (sample_stats$Aromatic -sample_stats$Aromatic_m)^2, (sample_stats$Phenolic -sample_stats$Phenolic_m)^2,
                                 (sample_stats$Amide_to_Ketone -sample_stats$Amide_to_Ketone_m)^2, (sample_stats$Sum -sample_stats$sum_m)^2)
        
        colnames(ssq_sample) <-  c("Alkyl_ssq", "N_Alkyl_Methoxyl_ssq", "O-Alkyl_ssq", "Di_O_Alkyl_ssq", "Aromatic_ssq", "Phenolic_ssq", "Amide_to_Ketone_ssq", "Sum_ssq")
        
        sample_stats <- cbind(sample_stats,ssq_sample)

      } else if (ecosys == "Terr_Baldock") {

        stdmat <- std_nmr(ecosys = "Terr_Baldock")
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30, NMRmeth = "MMM")

        ## Elementar ratios (relative to C)

        Cmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]+NMR.end$Protein*stdmat[[1]][9,1]+NMR.end$Lignin*stdmat[[1]][9,3]+
                             NMR.end$Lipid*stdmat[[1]][9,4]+NMR.end$Carbonyl*stdmat[[1]][9,5]+NMR.end$Char*stdmat[[1]][9,6])

        Nmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]+NMR.end$Protein*stdmat[[1]][10,1]+NMR.end$Lignin*stdmat[[1]][10,3]+
                             NMR.end$Lipid*stdmat[[1]][10,4]+NMR.end$Carbonyl*stdmat[[1]][10,5]+NMR.end$Char*stdmat[[1]][10,6])

        Hmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]+NMR.end$Protein*stdmat[[1]][11,1]+NMR.end$Lignin*stdmat[[1]][11,3]+
                             NMR.end$Lipid*stdmat[[1]][11,4]+NMR.end$Carbonyl*stdmat[[1]][11,5]+NMR.end$Char*stdmat[[1]][11,6])

        Omol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]+NMR.end$Protein*stdmat[[1]][12,1]+NMR.end$Lignin*stdmat[[1]][12,3]+
                             NMR.end$Lipid*stdmat[[1]][12,4]+NMR.end$Carbonyl*stdmat[[1]][12,5]+NMR.end$Char*stdmat[[1]][12,6])

        Cwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]*12.0107+NMR.end$Protein*stdmat[[1]][9,1]*12.0107+NMR.end$Lignin*stdmat[[1]][9,3]*12.0107+
                             NMR.end$Lipid*stdmat[[1]][9,4]*12.0107+NMR.end$Carbonyl*stdmat[[1]][9,5]*12.0107+NMR.end$Char*stdmat[[1]][9,6]*12.0107)

        Nwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]*14.0067+NMR.end$Protein*stdmat[[1]][10,1]*14.0067+NMR.end$Lignin*stdmat[[1]][10,3]*14.0067+
                             NMR.end$Lipid*stdmat[[1]][10,4]*14.0067+NMR.end$Carbonyl*stdmat[[1]][10,5]*14.0067+NMR.end$Char*stdmat[[1]][10,6]*14.0067)

        Hwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]*1.00794+NMR.end$Protein*stdmat[[1]][11,1]*1.00794+NMR.end$Lignin*stdmat[[1]][11,3]*1.00794+
                             NMR.end$Lipid*stdmat[[1]][11,4]*1.00794+NMR.end$Carbonyl*stdmat[[1]][11,5]*1.00794+NMR.end$Char*stdmat[[1]][11,6]*1.00794)

        Owgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]*15.994+NMR.end$Protein*stdmat[[1]][12,1]*15.994+NMR.end$Lignin*stdmat[[1]][12,3]*15.994+
                             NMR.end$Lipid*stdmat[[1]][12,4]*15.994+NMR.end$Carbonyl*stdmat[[1]][12,5]*15.994+NMR.end$Char*stdmat[[1]][12,6]*15.994)

        swgt <- c(Cwgt + Nwgt + Hwgt +Owgt)

        NOSC <- as.numeric(4+((2*Omol+3*Nmol-1*Hmol-4*Cmol)/Cmol))

        ## Final result

        NMR.end <- cbind(NMR.end, Cmol, Nmol, Hmol, Omol, Cwgt/swgt, Nwgt/swgt, Hwgt/swgt, Owgt/swgt, NOSC)

        ## Back calculated NMR results
        
        ##Alklyl C calculation
        Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][1,2]+NMR.end$Protein*stdmat[[1]][1,1]+NMR.end$Lignin*stdmat[[1]][1,3]+NMR.end$Lipid*stdmat[[1]][1,4]+
                               NMR.end$Carbonyl*stdmat[[1]][1,5]+NMR.end$Char*stdmat[[1]][1,6], c("Alkyl"))
        
        ##N_Alkyl_Methoxyl C calculation
        N_Alkyl_Methoxyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][2,2]+NMR.end$Protein*stdmat[[1]][2,1]+NMR.end$Lignin*stdmat[[1]][2,3]+NMR.end$Lipid*stdmat[[1]][2,4]+
                                          NMR.end$Carbonyl*stdmat[[1]][2,5]+NMR.end$Char*stdmat[[1]][2,6], c("N_Alkyl_Methoxyl"))
        
        ##O-Alkyl C calculation
        O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][3,2]+NMR.end$Protein*stdmat[[1]][3,1]+NMR.end$Lignin*stdmat[[1]][3,3]+NMR.end$Lipid*stdmat[[1]][3,4]+
                                 NMR.end$Carbonyl*stdmat[[1]][3,5]+NMR.end$Char*stdmat[[1]][3,6], c("O-Alkyl"))
        
        ##Di_O_Alkyl C calculation
        Di_O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][4,2]+NMR.end$Protein*stdmat[[1]][4,1]+NMR.end$Lignin*stdmat[[1]][4,3]+NMR.end$Lipid*stdmat[[1]][4,4]+
                                    NMR.end$Carbonyl*stdmat[[1]][4,5]+NMR.end$Char*stdmat[[1]][4,6], c("Di_O_Alkyl"))
        
        ##Aromatic C calculation
        Aromatic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][5,2]+NMR.end$Protein*stdmat[[1]][5,1]+NMR.end$Lignin*stdmat[[1]][5,3]+NMR.end$Lipid*stdmat[[1]][5,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][5,5]+NMR.end$Char*stdmat[[1]][5,6], c("Aromatic"))
        
        ##Phenolic C calculation
        Phenolic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][6,2]+NMR.end$Protein*stdmat[[1]][6,1]+NMR.end$Lignin*stdmat[[1]][6,3]+NMR.end$Lipid*stdmat[[1]][6,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][6,5]+NMR.end$Char*stdmat[[1]][6,6], c("Phenolic"))
        
        ##Amide_Carboxylic C calculation
        Amide_to_Ketone_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][7,2]+NMR.end$Protein*stdmat[[1]][7,1]+NMR.end$Lignin*stdmat[[1]][7,3]+NMR.end$Lipid*stdmat[[1]][7,4]+
                                         NMR.end$Carbonyl*stdmat[[1]][7,5]+NMR.end$Char*stdmat[[1]][7,6], c("Amide_to_Ketone"))
        
        sum_m <- setNames(Alkyl_m + N_Alkyl_Methoxyl_m + O_Alkyl_m + Di_O_Alkyl_m + Aromatic_m + Phenolic_m + Amide_to_Ketone_m, c("Sum"))
        
        sum_c <- sum(sampleintegraljoin)
        
        sampleintegraljoin <-rbind(sampleintegraljoin, sum_c)
        
        sample_stats <- data.frame(Alkyl_m, N_Alkyl_Methoxyl_m, O_Alkyl_m, Di_O_Alkyl_m, Aromatic_m, Phenolic_m, Amide_to_Ketone_m, sum_m)
        
        nmrrest <- NULL
        for (i in 1:nrow(sample_stats)) {
          nmrrestt <- c(sampleintegraljoin)
          nmrrest <- rbind(nmrrest, nmrrestt)
        }
        colnames(nmrrest) <-  c("Alkyl", "N_Alkyl_Methoxyl", "O-Alkyl", "Di_O_Alkyl", "Aromatic", "Phenolic", "Amide_to_Ketone", "Sum")
        
        sample_stats <- cbind(sample_stats,nmrrest)
        ssq_sample <- data.frame((sample_stats$Alkyl-sample_stats$Alkyl_m)^2, (sample_stats$N_Alkyl_Methoxyl-sample_stats$N_Alkyl_Methoxyl_m)^2,
                                 (sample_stats$`O-Alkyl` -sample_stats$O_Alkyl_m)^2, (sample_stats$Di_O_Alkyl-sample_stats$Di_O_Alkyl_m)^2, 
                                 (sample_stats$Aromatic -sample_stats$Aromatic_m)^2, (sample_stats$Phenolic -sample_stats$Phenolic_m)^2,
                                 (sample_stats$Amide_to_Ketone -sample_stats$Amide_to_Ketone_m)^2, (sample_stats$Sum -sample_stats$sum_m)^2)
        
        colnames(ssq_sample) <-  c("Alkyl_ssq", "N_Alkyl_Methoxyl_ssq", "O-Alkyl_ssq", "Di_O_Alkyl_ssq", "Aromatic_ssq", "Phenolic_ssq", "Amide_to_Ketone_ssq", "Sum_ssq")
        
        sample_stats <- cbind(sample_stats,ssq_sample)        
      } else if (ecosys == "Aqua_Nelson") {

        stdmat <- std_nmr(ecosys = "Aqua_Nelson")
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30,  NMRmeth = "MMM")

        ## Elementar ratios (relative to C)

        Cmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]+NMR.end$Protein*stdmat[[1]][9,1]+NMR.end$Lignin*stdmat[[1]][9,3]+
                             NMR.end$Lipid*stdmat[[1]][9,4]+NMR.end$Carbonyl*stdmat[[1]][9,5]+NMR.end$Char*stdmat[[1]][9,6])

        Nmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]+NMR.end$Protein*stdmat[[1]][10,1]+NMR.end$Lignin*stdmat[[1]][10,3]+
                             NMR.end$Lipid*stdmat[[1]][10,4]+NMR.end$Carbonyl*stdmat[[1]][10,5]+NMR.end$Char*stdmat[[1]][10,6])

        Hmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]+NMR.end$Protein*stdmat[[1]][11,1]+NMR.end$Lignin*stdmat[[1]][11,3]+
                             NMR.end$Lipid*stdmat[[1]][11,4]+NMR.end$Carbonyl*stdmat[[1]][11,5]+NMR.end$Char*stdmat[[1]][11,6])

        Omol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]+NMR.end$Protein*stdmat[[1]][12,1]+NMR.end$Lignin*stdmat[[1]][12,3]+
                             NMR.end$Lipid*stdmat[[1]][12,4]+NMR.end$Carbonyl*stdmat[[1]][12,5]+NMR.end$Char*stdmat[[1]][12,6])

        Cwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]*12.0107+NMR.end$Protein*stdmat[[1]][9,1]*12.0107+NMR.end$Lignin*stdmat[[1]][9,3]*12.0107+
                             NMR.end$Lipid*stdmat[[1]][9,4]*12.0107+NMR.end$Carbonyl*stdmat[[1]][9,5]*12.0107+NMR.end$Char*stdmat[[1]][9,6]*12.0107)

        Nwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]*14.0067+NMR.end$Protein*stdmat[[1]][10,1]*14.0067+NMR.end$Lignin*stdmat[[1]][10,3]*14.0067+
                             NMR.end$Lipid*stdmat[[1]][10,4]*14.0067+NMR.end$Carbonyl*stdmat[[1]][10,5]*14.0067+NMR.end$Char*stdmat[[1]][10,6]*14.0067)

        Hwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]*1.00794+NMR.end$Protein*stdmat[[1]][11,1]*1.00794+NMR.end$Lignin*stdmat[[1]][11,3]*1.00794+
                             NMR.end$Lipid*stdmat[[1]][11,4]*1.00794+NMR.end$Carbonyl*stdmat[[1]][11,5]*1.00794+NMR.end$Char*stdmat[[1]][11,6]*1.00794)

        Owgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]*15.994+NMR.end$Protein*stdmat[[1]][12,1]*15.994+NMR.end$Lignin*stdmat[[1]][12,3]*15.994+
                             NMR.end$Lipid*stdmat[[1]][12,4]*15.994+NMR.end$Carbonyl*stdmat[[1]][12,5]*15.994+NMR.end$Char*stdmat[[1]][12,6]*15.994)

        swgt <- c(Cwgt + Nwgt + Hwgt +Owgt)

        NOSC <- as.numeric(4+((2*Omol+3*Nmol-1*Hmol-4*Cmol)/Cmol))

        ## Final result

        NMR.end <- cbind(NMR.end, Cmol, Nmol, Hmol, Omol, Cwgt/swgt, Nwgt/swgt, Hwgt/swgt, Owgt/swgt, NOSC)
        
        ## Back calculated NMR results
        
        ##Alklyl C calculation
        Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][1,2]+NMR.end$Protein*stdmat[[1]][1,1]+NMR.end$Lignin*stdmat[[1]][1,3]+NMR.end$Lipid*stdmat[[1]][1,4]+
                               NMR.end$Carbonyl*stdmat[[1]][1,5]+NMR.end$Char*stdmat[[1]][1,6], c("Alkyl"))
        
        ##N_Alkyl_Methoxyl C calculation
        N_Alkyl_Methoxyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][2,2]+NMR.end$Protein*stdmat[[1]][2,1]+NMR.end$Lignin*stdmat[[1]][2,3]+NMR.end$Lipid*stdmat[[1]][2,4]+
                                          NMR.end$Carbonyl*stdmat[[1]][2,5]+NMR.end$Char*stdmat[[1]][2,6], c("N_Alkyl_Methoxyl"))
        
        ##O-Alkyl C calculation
        O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][3,2]+NMR.end$Protein*stdmat[[1]][3,1]+NMR.end$Lignin*stdmat[[1]][3,3]+NMR.end$Lipid*stdmat[[1]][3,4]+
                                 NMR.end$Carbonyl*stdmat[[1]][3,5]+NMR.end$Char*stdmat[[1]][3,6], c("O-Alkyl"))
        
        ##Di_O_Alkyl C calculation
        Di_O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][4,2]+NMR.end$Protein*stdmat[[1]][4,1]+NMR.end$Lignin*stdmat[[1]][4,3]+NMR.end$Lipid*stdmat[[1]][4,4]+
                                    NMR.end$Carbonyl*stdmat[[1]][4,5]+NMR.end$Char*stdmat[[1]][4,6], c("Di_O_Alkyl"))
        
        ##Aromatic C calculation
        Aromatic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][5,2]+NMR.end$Protein*stdmat[[1]][5,1]+NMR.end$Lignin*stdmat[[1]][5,3]+NMR.end$Lipid*stdmat[[1]][5,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][5,5]+NMR.end$Char*stdmat[[1]][5,6], c("Aromatic"))
        
        ##Phenolic C calculation
        Phenolic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][6,2]+NMR.end$Protein*stdmat[[1]][6,1]+NMR.end$Lignin*stdmat[[1]][6,3]+NMR.end$Lipid*stdmat[[1]][6,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][6,5]+NMR.end$Char*stdmat[[1]][6,6], c("Phenolic"))
        
        ##Amide_Carboxylic C calculation
        Amide_to_Ketone_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][7,2]+NMR.end$Protein*stdmat[[1]][7,1]+NMR.end$Lignin*stdmat[[1]][7,3]+NMR.end$Lipid*stdmat[[1]][7,4]+
                                         NMR.end$Carbonyl*stdmat[[1]][7,5]+NMR.end$Char*stdmat[[1]][7,6], c("Amide_to_Ketone"))
        
        sum_m <- setNames(Alkyl_m + N_Alkyl_Methoxyl_m + O_Alkyl_m + Di_O_Alkyl_m + Aromatic_m + Phenolic_m + Amide_to_Ketone_m, c("Sum"))
        
        sum_c <- sum(sampleintegraljoin)
        
        sampleintegraljoin <-rbind(sampleintegraljoin, sum_c)
        
        sample_stats <- data.frame(Alkyl_m, N_Alkyl_Methoxyl_m, O_Alkyl_m, Di_O_Alkyl_m, Aromatic_m, Phenolic_m, Amide_to_Ketone_m, sum_m)
        
        nmrrest <- NULL
        for (i in 1:nrow(sample_stats)) {
          nmrrestt <- c(sampleintegraljoin)
          nmrrest <- rbind(nmrrest, nmrrestt)
        }
        colnames(nmrrest) <-  c("Alkyl", "N_Alkyl_Methoxyl", "O-Alkyl", "Di_O_Alkyl", "Aromatic", "Phenolic", "Amide_to_Ketone", "Sum")
        
        sample_stats <- cbind(sample_stats,nmrrest)
        ssq_sample <- data.frame((sample_stats$Alkyl-sample_stats$Alkyl_m)^2, (sample_stats$N_Alkyl_Methoxyl-sample_stats$N_Alkyl_Methoxyl_m)^2,
                                 (sample_stats$`O-Alkyl` -sample_stats$O_Alkyl_m)^2, (sample_stats$Di_O_Alkyl-sample_stats$Di_O_Alkyl_m)^2, 
                                 (sample_stats$Aromatic -sample_stats$Aromatic_m)^2, (sample_stats$Phenolic -sample_stats$Phenolic_m)^2,
                                 (sample_stats$Amide_to_Ketone -sample_stats$Amide_to_Ketone_m)^2, (sample_stats$Sum -sample_stats$sum_m)^2)
        
        colnames(ssq_sample) <-  c("Alkyl_ssq", "N_Alkyl_Methoxyl_ssq", "O-Alkyl_ssq", "Di_O_Alkyl_ssq", "Aromatic_ssq", "Phenolic_ssq", "Amide_to_Ketone_ssq", "Sum_ssq")
        
        sample_stats <- cbind(sample_stats,ssq_sample)
        
      } else if (ecosys == "Aqua_Baldock") {

        stdmat <- std_nmr(ecosys = "Aqua_Baldock")
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30,  NMRmeth = "MMM")

        ## Elementar ratios (relative to C)

        Cmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]+NMR.end$Protein*stdmat[[1]][9,1]+NMR.end$Lignin*stdmat[[1]][9,3]+
                             NMR.end$Lipid*stdmat[[1]][9,4]+NMR.end$Carbonyl*stdmat[[1]][9,5]+NMR.end$Char*stdmat[[1]][9,6])

        Nmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]+NMR.end$Protein*stdmat[[1]][10,1]+NMR.end$Lignin*stdmat[[1]][10,3]+
                             NMR.end$Lipid*stdmat[[1]][10,4]+NMR.end$Carbonyl*stdmat[[1]][10,5]+NMR.end$Char*stdmat[[1]][10,6])

        Hmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]+NMR.end$Protein*stdmat[[1]][11,1]+NMR.end$Lignin*stdmat[[1]][11,3]+
                             NMR.end$Lipid*stdmat[[1]][11,4]+NMR.end$Carbonyl*stdmat[[1]][11,5]+NMR.end$Char*stdmat[[1]][11,6])

        Omol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]+NMR.end$Protein*stdmat[[1]][12,1]+NMR.end$Lignin*stdmat[[1]][12,3]+
                             NMR.end$Lipid*stdmat[[1]][12,4]+NMR.end$Carbonyl*stdmat[[1]][12,5]+NMR.end$Char*stdmat[[1]][12,6])

        Cwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]*12.0107+NMR.end$Protein*stdmat[[1]][9,1]*12.0107+NMR.end$Lignin*stdmat[[1]][9,3]*12.0107+
                             NMR.end$Lipid*stdmat[[1]][9,4]*12.0107+NMR.end$Carbonyl*stdmat[[1]][9,5]*12.0107+NMR.end$Char*stdmat[[1]][9,6]*12.0107)

        Nwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]*14.0067+NMR.end$Protein*stdmat[[1]][10,1]*14.0067+NMR.end$Lignin*stdmat[[1]][10,3]*14.0067+
                             NMR.end$Lipid*stdmat[[1]][10,4]*14.0067+NMR.end$Carbonyl*stdmat[[1]][10,5]*14.0067+NMR.end$Char*stdmat[[1]][10,6]*14.0067)

        Hwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]*1.00794+NMR.end$Protein*stdmat[[1]][11,1]*1.00794+NMR.end$Lignin*stdmat[[1]][11,3]*1.00794+
                             NMR.end$Lipid*stdmat[[1]][11,4]*1.00794+NMR.end$Carbonyl*stdmat[[1]][11,5]*1.00794+NMR.end$Char*stdmat[[1]][11,6]*1.00794)

        Owgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]*15.994+NMR.end$Protein*stdmat[[1]][12,1]*15.994+NMR.end$Lignin*stdmat[[1]][12,3]*15.994+
                             NMR.end$Lipid*stdmat[[1]][12,4]*15.994+NMR.end$Carbonyl*stdmat[[1]][12,5]*15.994+NMR.end$Char*stdmat[[1]][12,6]*15.994)

        swgt <- c(Cwgt + Nwgt + Hwgt +Owgt)

        NOSC <- as.numeric(4+((2*Omol+3*Nmol-1*Hmol-4*Cmol)/Cmol))
        
        ## Back calculated NMR results
        
        ##Alklyl C calculation
        Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][1,2]+NMR.end$Protein*stdmat[[1]][1,1]+NMR.end$Lignin*stdmat[[1]][1,3]+NMR.end$Lipid*stdmat[[1]][1,4]+
                               NMR.end$Carbonyl*stdmat[[1]][1,5]+NMR.end$Char*stdmat[[1]][1,6], c("Alkyl"))
        
        ##N_Alkyl_Methoxyl C calculation
        N_Alkyl_Methoxyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][2,2]+NMR.end$Protein*stdmat[[1]][2,1]+NMR.end$Lignin*stdmat[[1]][2,3]+NMR.end$Lipid*stdmat[[1]][2,4]+
                                          NMR.end$Carbonyl*stdmat[[1]][2,5]+NMR.end$Char*stdmat[[1]][2,6], c("N_Alkyl_Methoxyl"))
        
        ##O-Alkyl C calculation
        O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][3,2]+NMR.end$Protein*stdmat[[1]][3,1]+NMR.end$Lignin*stdmat[[1]][3,3]+NMR.end$Lipid*stdmat[[1]][3,4]+
                                 NMR.end$Carbonyl*stdmat[[1]][3,5]+NMR.end$Char*stdmat[[1]][3,6], c("O-Alkyl"))
        
        ##Di_O_Alkyl C calculation
        Di_O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][4,2]+NMR.end$Protein*stdmat[[1]][4,1]+NMR.end$Lignin*stdmat[[1]][4,3]+NMR.end$Lipid*stdmat[[1]][4,4]+
                                    NMR.end$Carbonyl*stdmat[[1]][4,5]+NMR.end$Char*stdmat[[1]][4,6], c("Di_O_Alkyl"))
        
        ##Aromatic C calculation
        Aromatic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][5,2]+NMR.end$Protein*stdmat[[1]][5,1]+NMR.end$Lignin*stdmat[[1]][5,3]+NMR.end$Lipid*stdmat[[1]][5,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][5,5]+NMR.end$Char*stdmat[[1]][5,6], c("Aromatic"))
        
        ##Phenolic C calculation
        Phenolic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][6,2]+NMR.end$Protein*stdmat[[1]][6,1]+NMR.end$Lignin*stdmat[[1]][6,3]+NMR.end$Lipid*stdmat[[1]][6,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][6,5]+NMR.end$Char*stdmat[[1]][6,6], c("Phenolic"))
        
        ##Amide_Carboxylic C calculation
        Amide_to_Ketone_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][7,2]+NMR.end$Protein*stdmat[[1]][7,1]+NMR.end$Lignin*stdmat[[1]][7,3]+NMR.end$Lipid*stdmat[[1]][7,4]+
                                         NMR.end$Carbonyl*stdmat[[1]][7,5]+NMR.end$Char*stdmat[[1]][7,6], c("Amide_to_Ketone"))
        
        sum_m <- setNames(Alkyl_m + N_Alkyl_Methoxyl_m + O_Alkyl_m + Di_O_Alkyl_m + Aromatic_m + Phenolic_m + Amide_to_Ketone_m, c("Sum"))
        
        sum_c <- sum(sampleintegraljoin)
        
        sampleintegraljoin <-rbind(sampleintegraljoin, sum_c)
        
        sample_stats <- data.frame(Alkyl_m, N_Alkyl_Methoxyl_m, O_Alkyl_m, Di_O_Alkyl_m, Aromatic_m, Phenolic_m, Amide_to_Ketone_m, sum_m)
        
        nmrrest <- NULL
        for (i in 1:nrow(sample_stats)) {
          nmrrestt <- c(sampleintegraljoin)
          nmrrest <- rbind(nmrrest, nmrrestt)
        }
        colnames(nmrrest) <-  c("Alkyl", "N_Alkyl_Methoxyl", "O-Alkyl", "Di_O_Alkyl", "Aromatic", "Phenolic", "Amide_to_Ketone", "Sum")
        
        sample_stats <- cbind(sample_stats,nmrrest)
        ssq_sample <- data.frame((sample_stats$Alkyl-sample_stats$Alkyl_m)^2, (sample_stats$N_Alkyl_Methoxyl-sample_stats$N_Alkyl_Methoxyl_m)^2,
                                 (sample_stats$`O-Alkyl` -sample_stats$O_Alkyl_m)^2, (sample_stats$Di_O_Alkyl-sample_stats$Di_O_Alkyl_m)^2, 
                                 (sample_stats$Aromatic -sample_stats$Aromatic_m)^2, (sample_stats$Phenolic -sample_stats$Phenolic_m)^2,
                                 (sample_stats$Amide_to_Ketone -sample_stats$Amide_to_Ketone_m)^2, (sample_stats$Sum -sample_stats$sum_m)^2)
        
        colnames(ssq_sample) <-  c("Alkyl_ssq", "N_Alkyl_Methoxyl_ssq", "O-Alkyl_ssq", "Di_O_Alkyl_ssq", "Aromatic_ssq", "Phenolic_ssq", "Amide_to_Ketone_ssq", "Sum_ssq")
        
        sample_stats <- cbind(sample_stats,ssq_sample)
        
        ## Final result

        NMR.end <- cbind(NMR.end, Cmol, Nmol, Hmol, Omol, Cwgt/swgt, Nwgt/swgt, Hwgt/swgt, Owgt/swgt, NOSC)

      } else if (ecosys == "mod") {

        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = mod_std, amoSTD = 6, best.fits = 30,  NMRmeth = "MMM")

        ## Elementar ratios (relative to C)

        Cmol <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][9,2]+NMR.end$Protein*mod_std[[1]][9,1]+NMR.end$Lignin*mod_std[[1]][9,3]+
                             NMR.end$Lipid*mod_std[[1]][9,4]+NMR.end$Carbonyl*mod_std[[1]][9,5]+NMR.end$Char*mod_std[[1]][9,6])

        Nmol <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][10,2]+NMR.end$Protein*mod_std[[1]][10,1]+NMR.end$Lignin*mod_std[[1]][10,3]+
                             NMR.end$Lipid*mod_std[[1]][10,4]+NMR.end$Carbonyl*mod_std[[1]][10,5]+NMR.end$Char*mod_std[[1]][10,6])

        Hmol <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][11,2]+NMR.end$Protein*mod_std[[1]][11,1]+NMR.end$Lignin*mod_std[[1]][11,3]+
                             NMR.end$Lipid*mod_std[[1]][11,4]+NMR.end$Carbonyl*mod_std[[1]][11,5]+NMR.end$Char*mod_std[[1]][11,6])

        Omol <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][12,2]+NMR.end$Protein*mod_std[[1]][12,1]+NMR.end$Lignin*mod_std[[1]][12,3]+
                             NMR.end$Lipid*mod_std[[1]][12,4]+NMR.end$Carbonyl*mod_std[[1]][12,5]+NMR.end$Char*mod_std[[1]][12,6])

        Cwgt <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][9,2]*12.0107+NMR.end$Protein*mod_std[[1]][9,1]*12.0107+NMR.end$Lignin*mod_std[[1]][9,3]*12.0107+
                             NMR.end$Lipid*mod_std[[1]][9,4]*12.0107+NMR.end$Carbonyl*mod_std[[1]][9,5]*12.0107+NMR.end$Char*mod_std[[1]][9,6]*12.0107)

        Nwgt <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][10,2]*14.0067+NMR.end$Protein*mod_std[[1]][10,1]*14.0067+NMR.end$Lignin*mod_std[[1]][10,3]*14.0067+
                             NMR.end$Lipid*mod_std[[1]][10,4]*14.0067+NMR.end$Carbonyl*mod_std[[1]][10,5]*14.0067+NMR.end$Char*mod_std[[1]][10,6]*14.0067)

        Hwgt <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][11,2]*1.00794+NMR.end$Protein*mod_std[[1]][11,1]*1.00794+NMR.end$Lignin*mod_std[[1]][11,3]*1.00794+
                             NMR.end$Lipid*mod_std[[1]][11,4]*1.00794+NMR.end$Carbonyl*mod_std[[1]][11,5]*1.00794+NMR.end$Char*mod_std[[1]][11,6]*1.00794)

        Owgt <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][12,2]*15.994+NMR.end$Protein*mod_std[[1]][12,1]*15.994+NMR.end$Lignin*mod_std[[1]][12,3]*15.994+
                             NMR.end$Lipid*mod_std[[1]][12,4]*15.994+NMR.end$Carbonyl*mod_std[[1]][12,5]*15.994+NMR.end$Char*mod_std[[1]][12,6]*15.994)

        swgt <- c(Cwgt + Nwgt + Hwgt +Owgt)

        NOSC <- as.numeric(4+((2*Omol+3*Nmol-1*Hmol-4*Cmol)/Cmol))

        ## Final result

        NMR.end <- cbind(NMR.end, Cmol, Nmol, Hmol, Omol, Cwgt/swgt, Nwgt/swgt, Hwgt/swgt, Owgt/swgt, NOSC)
        
        ## Back calculated NMR results
        
        ##Alklyl C calculation
        Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][1,2]+NMR.end$Protein*stdmat[[1]][1,1]+NMR.end$Lignin*stdmat[[1]][1,3]+NMR.end$Lipid*stdmat[[1]][1,4]+
                               NMR.end$Carbonyl*stdmat[[1]][1,5]+NMR.end$Char*stdmat[[1]][1,6], c("Alkyl"))
        
        ##N_Alkyl_Methoxyl C calculation
        N_Alkyl_Methoxyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][2,2]+NMR.end$Protein*stdmat[[1]][2,1]+NMR.end$Lignin*stdmat[[1]][2,3]+NMR.end$Lipid*stdmat[[1]][2,4]+
                                          NMR.end$Carbonyl*stdmat[[1]][2,5]+NMR.end$Char*stdmat[[1]][2,6], c("N_Alkyl_Methoxyl"))
        
        ##O-Alkyl C calculation
        O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][3,2]+NMR.end$Protein*stdmat[[1]][3,1]+NMR.end$Lignin*stdmat[[1]][3,3]+NMR.end$Lipid*stdmat[[1]][3,4]+
                                 NMR.end$Carbonyl*stdmat[[1]][3,5]+NMR.end$Char*stdmat[[1]][3,6], c("O-Alkyl"))
        
        ##Di_O_Alkyl C calculation
        Di_O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][4,2]+NMR.end$Protein*stdmat[[1]][4,1]+NMR.end$Lignin*stdmat[[1]][4,3]+NMR.end$Lipid*stdmat[[1]][4,4]+
                                    NMR.end$Carbonyl*stdmat[[1]][4,5]+NMR.end$Char*stdmat[[1]][4,6], c("Di_O_Alkyl"))
        
        ##Aromatic C calculation
        Aromatic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][5,2]+NMR.end$Protein*stdmat[[1]][5,1]+NMR.end$Lignin*stdmat[[1]][5,3]+NMR.end$Lipid*stdmat[[1]][5,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][5,5]+NMR.end$Char*stdmat[[1]][5,6], c("Aromatic"))
        
        ##Phenolic C calculation
        Phenolic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][6,2]+NMR.end$Protein*stdmat[[1]][6,1]+NMR.end$Lignin*stdmat[[1]][6,3]+NMR.end$Lipid*stdmat[[1]][6,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][6,5]+NMR.end$Char*stdmat[[1]][6,6], c("Phenolic"))
        
        ##Amide_Carboxylic C calculation
        Amide_to_Ketone_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][7,2]+NMR.end$Protein*stdmat[[1]][7,1]+NMR.end$Lignin*stdmat[[1]][7,3]+NMR.end$Lipid*stdmat[[1]][7,4]+
                                         NMR.end$Carbonyl*stdmat[[1]][7,5]+NMR.end$Char*stdmat[[1]][7,6], c("Amide_to_Ketone"))
        
        sum_m <- setNames(Alkyl_m + N_Alkyl_Methoxyl_m + O_Alkyl_m + Di_O_Alkyl_m + Aromatic_m + Phenolic_m + Amide_to_Ketone_m, c("Sum"))
        
        sum_c <- sum(sampleintegraljoin)
        
        sampleintegraljoin <-rbind(sampleintegraljoin, sum_c)
        
        sample_stats <- data.frame(Alkyl_m, N_Alkyl_Methoxyl_m, O_Alkyl_m, Di_O_Alkyl_m, Aromatic_m, Phenolic_m, Amide_to_Ketone_m, sum_m)
        
        nmrrest <- NULL
        for (i in 1:nrow(sample_stats)) {
          nmrrestt <- c(sampleintegraljoin)
          nmrrest <- rbind(nmrrest, nmrrestt)
        }
        colnames(nmrrest) <-  c("Alkyl", "N_Alkyl_Methoxyl", "O-Alkyl", "Di_O_Alkyl", "Aromatic", "Phenolic", "Amide_to_Ketone", "Sum")
        
        sample_stats <- cbind(sample_stats,nmrrest)
        ssq_sample <- data.frame((sample_stats$Alkyl-sample_stats$Alkyl_m)^2, (sample_stats$N_Alkyl_Methoxyl-sample_stats$N_Alkyl_Methoxyl_m)^2,
                                 (sample_stats$`O-Alkyl` -sample_stats$O_Alkyl_m)^2, (sample_stats$Di_O_Alkyl-sample_stats$Di_O_Alkyl_m)^2, 
                                 (sample_stats$Aromatic -sample_stats$Aromatic_m)^2, (sample_stats$Phenolic -sample_stats$Phenolic_m)^2,
                                 (sample_stats$Amide_to_Ketone -sample_stats$Amide_to_Ketone_m)^2, (sample_stats$Sum -sample_stats$sum_m)^2)
        
        colnames(ssq_sample) <-  c("Alkyl_ssq", "N_Alkyl_Methoxyl_ssq", "O-Alkyl_ssq", "Di_O_Alkyl_ssq", "Aromatic_ssq", "Phenolic_ssq", "Amide_to_Ketone_ssq", "Sum_ssq")
        
        sample_stats <- cbind(sample_stats,ssq_sample)
        
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

      ##Alklyl C calculation
      Alkyl <-  setNames(data.frame(sum(2*sum(sampleintegral[1:3,1]), sum(sampleintegral[10:12,1]), -sum(sampleintegral[28:30,1]))), c("Alkyl"))
      
      ##N_Alkyl_Methoxyl C calculation
      N_Alkyl_Methoxyl <-  setNames(data.frame(sum(2*sum(sampleintegral[4:5,1]), sum(sampleintegral[13:14,1]), -sum(sampleintegral[31:32,1]))), c("N_Alkyl_Methoxyl"))
      
      ##O-Alkyl C calculation
      O_Alkyl <-  setNames(data.frame(sum(2*sum(sampleintegral[6:7,1]), sum(sampleintegral[15:16,1]), -sum(sampleintegral[33:33,1]))), c("O-Alkyl"))
      
      ##Di_O_Alkyl C calculation
      Di_O_Alkyl <-  setNames(data.frame(sum(2*sum(sampleintegral[8:8,1]), sum(sampleintegral[17:17,1]))), c("Di_O_Alkyl"))
      
      ##Aromatic C calculation
      Aromatic <-  setNames(data.frame(sum(2*sum(sampleintegral[27:28,1]), sum(sampleintegral[18:19,1]), -sum(sampleintegral[1:1,1]))), c("Aromatic"))
      
      ##Phenolic C calculation
      Phenolic <-  setNames(data.frame(sum(2*sum(sampleintegral[29:29,1]), sum(sampleintegral[20:20,1]), -sum(sampleintegral[2:2,1]))), c("Phenolic"))
      
      ##Amide_Carboxylic C calculation
      Amide_Carboxylic <-  setNames(data.frame(sum(2*sum(sampleintegral[30:32,1]), sum(sampleintegral[21:23,1]), -sum(sampleintegral[3:5,1]))), c("Amide_Carboxylic"))
      
      ##Ketone C calculation
      Ketone <-  setNames(data.frame(sum(2*sum(sampleintegral[33:33,1]), sum(sampleintegral[24:24,1]), -sum(sampleintegral[6:6,1]))), c("Ketone"))
      
      ##Put all together
      Amide_to_Ketone <- c(Amide_Carboxylic + Ketone)
      
      sampleintegraljoin <- data.frame(Alkyl, N_Alkyl_Methoxyl, O_Alkyl, Di_O_Alkyl, Aromatic, Phenolic, Amide_to_Ketone)
      
      sampleintegraljoin <- t(sampleintegraljoin)
      
      sampleintegralend <- rbind(NCval,sampleintegraljoin)
      
      raw.spec.end[[i]] <- list("name" = samplename, "data" = list("raw.spec" = sampleraw.spec,"Integral" = sampleintegralend))
      
      if (ecosys == "Terr_Nelson") {

        stdmat <- std_nmr(ecosys = "Terr_Nelson")
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30, NMRmeth = "MMMFixN")

        ## Elementar ratios (relative to C)

        Cmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]+NMR.end$Protein*stdmat[[1]][9,1]+NMR.end$Lignin*stdmat[[1]][9,3]+
                             NMR.end$Lipid*stdmat[[1]][9,4]+NMR.end$Carbonyl*stdmat[[1]][9,5]+NMR.end$Char*stdmat[[1]][9,6])

        Nmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]+NMR.end$Protein*stdmat[[1]][10,1]+NMR.end$Lignin*stdmat[[1]][10,3]+
                             NMR.end$Lipid*stdmat[[1]][10,4]+NMR.end$Carbonyl*stdmat[[1]][10,5]+NMR.end$Char*stdmat[[1]][10,6])

        Hmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]+NMR.end$Protein*stdmat[[1]][11,1]+NMR.end$Lignin*stdmat[[1]][11,3]+
                             NMR.end$Lipid*stdmat[[1]][11,4]+NMR.end$Carbonyl*stdmat[[1]][11,5]+NMR.end$Char*stdmat[[1]][11,6])

        Omol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]+NMR.end$Protein*stdmat[[1]][12,1]+NMR.end$Lignin*stdmat[[1]][12,3]+
                             NMR.end$Lipid*stdmat[[1]][12,4]+NMR.end$Carbonyl*stdmat[[1]][12,5]+NMR.end$Char*stdmat[[1]][12,6])

        Cwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]*12.0107+NMR.end$Protein*stdmat[[1]][9,1]*12.0107+NMR.end$Lignin*stdmat[[1]][9,3]*12.0107+
                             NMR.end$Lipid*stdmat[[1]][9,4]*12.0107+NMR.end$Carbonyl*stdmat[[1]][9,5]*12.0107+NMR.end$Char*stdmat[[1]][9,6]*12.0107)

        Nwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]*14.0067+NMR.end$Protein*stdmat[[1]][10,1]*14.0067+NMR.end$Lignin*stdmat[[1]][10,3]*14.0067+
                             NMR.end$Lipid*stdmat[[1]][10,4]*14.0067+NMR.end$Carbonyl*stdmat[[1]][10,5]*14.0067+NMR.end$Char*stdmat[[1]][10,6]*14.0067)

        Hwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]*1.00794+NMR.end$Protein*stdmat[[1]][11,1]*1.00794+NMR.end$Lignin*stdmat[[1]][11,3]*1.00794+
                             NMR.end$Lipid*stdmat[[1]][11,4]*1.00794+NMR.end$Carbonyl*stdmat[[1]][11,5]*1.00794+NMR.end$Char*stdmat[[1]][11,6]*1.00794)

        Owgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]*15.994+NMR.end$Protein*stdmat[[1]][12,1]*15.994+NMR.end$Lignin*stdmat[[1]][12,3]*15.994+
                             NMR.end$Lipid*stdmat[[1]][12,4]*15.994+NMR.end$Carbonyl*stdmat[[1]][12,5]*15.994+NMR.end$Char*stdmat[[1]][12,6]*15.994)

        swgt <- c(Cwgt + Nwgt + Hwgt +Owgt)

        NOSC <- as.numeric(4+((2*Omol+3*Nmol-1*Hmol-4*Cmol)/Cmol))

        ## Final result

        NMR.end <- cbind(NMR.end, Cmol, Nmol, Hmol, Omol, Cwgt/swgt, Nwgt/swgt, Hwgt/swgt, Owgt/swgt, NOSC)
        
        ## Back calculated NMR results
        
        ##Alklyl C calculation
        Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][1,2]+NMR.end$Protein*stdmat[[1]][1,1]+NMR.end$Lignin*stdmat[[1]][1,3]+NMR.end$Lipid*stdmat[[1]][1,4]+
                               NMR.end$Carbonyl*stdmat[[1]][1,5]+NMR.end$Char*stdmat[[1]][1,6], c("Alkyl"))
        
        ##N_Alkyl_Methoxyl C calculation
        N_Alkyl_Methoxyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][2,2]+NMR.end$Protein*stdmat[[1]][2,1]+NMR.end$Lignin*stdmat[[1]][2,3]+NMR.end$Lipid*stdmat[[1]][2,4]+
                                          NMR.end$Carbonyl*stdmat[[1]][2,5]+NMR.end$Char*stdmat[[1]][2,6], c("N_Alkyl_Methoxyl"))
        
        ##O-Alkyl C calculation
        O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][3,2]+NMR.end$Protein*stdmat[[1]][3,1]+NMR.end$Lignin*stdmat[[1]][3,3]+NMR.end$Lipid*stdmat[[1]][3,4]+
                                 NMR.end$Carbonyl*stdmat[[1]][3,5]+NMR.end$Char*stdmat[[1]][3,6], c("O-Alkyl"))
        
        ##Di_O_Alkyl C calculation
        Di_O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][4,2]+NMR.end$Protein*stdmat[[1]][4,1]+NMR.end$Lignin*stdmat[[1]][4,3]+NMR.end$Lipid*stdmat[[1]][4,4]+
                                    NMR.end$Carbonyl*stdmat[[1]][4,5]+NMR.end$Char*stdmat[[1]][4,6], c("Di_O_Alkyl"))
        
        ##Aromatic C calculation
        Aromatic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][5,2]+NMR.end$Protein*stdmat[[1]][5,1]+NMR.end$Lignin*stdmat[[1]][5,3]+NMR.end$Lipid*stdmat[[1]][5,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][5,5]+NMR.end$Char*stdmat[[1]][5,6], c("Aromatic"))
        
        ##Phenolic C calculation
        Phenolic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][6,2]+NMR.end$Protein*stdmat[[1]][6,1]+NMR.end$Lignin*stdmat[[1]][6,3]+NMR.end$Lipid*stdmat[[1]][6,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][6,5]+NMR.end$Char*stdmat[[1]][6,6], c("Phenolic"))
        
        ##Amide_Carboxylic C calculation
        Amide_to_Ketone_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][7,2]+NMR.end$Protein*stdmat[[1]][7,1]+NMR.end$Lignin*stdmat[[1]][7,3]+NMR.end$Lipid*stdmat[[1]][7,4]+
                                         NMR.end$Carbonyl*stdmat[[1]][7,5]+NMR.end$Char*stdmat[[1]][7,6], c("Amide_to_Ketone"))
        
        sum_m <- setNames(Alkyl_m + N_Alkyl_Methoxyl_m + O_Alkyl_m + Di_O_Alkyl_m + Aromatic_m + Phenolic_m + Amide_to_Ketone_m, c("Sum"))
        
        sum_c <- sum(sampleintegraljoin)
        
        sampleintegraljoin <-rbind(sampleintegraljoin, sum_c)
        
        sample_stats <- data.frame(Alkyl_m, N_Alkyl_Methoxyl_m, O_Alkyl_m, Di_O_Alkyl_m, Aromatic_m, Phenolic_m, Amide_to_Ketone_m, sum_m)
        
        nmrrest <- NULL
        for (i in 1:nrow(sample_stats)) {
          nmrrestt <- c(sampleintegraljoin)
          nmrrest <- rbind(nmrrest, nmrrestt)
        }
        colnames(nmrrest) <-  c("Alkyl", "N_Alkyl_Methoxyl", "O-Alkyl", "Di_O_Alkyl", "Aromatic", "Phenolic", "Amide_to_Ketone", "Sum")
        
        sample_stats <- cbind(sample_stats,nmrrest)
        ssq_sample <- data.frame((sample_stats$Alkyl-sample_stats$Alkyl_m)^2, (sample_stats$N_Alkyl_Methoxyl-sample_stats$N_Alkyl_Methoxyl_m)^2,
                                 (sample_stats$`O-Alkyl` -sample_stats$O_Alkyl_m)^2, (sample_stats$Di_O_Alkyl-sample_stats$Di_O_Alkyl_m)^2, 
                                 (sample_stats$Aromatic -sample_stats$Aromatic_m)^2, (sample_stats$Phenolic -sample_stats$Phenolic_m)^2,
                                 (sample_stats$Amide_to_Ketone -sample_stats$Amide_to_Ketone_m)^2, (sample_stats$Sum -sample_stats$sum_m)^2)
        
        colnames(ssq_sample) <-  c("Alkyl_ssq", "N_Alkyl_Methoxyl_ssq", "O-Alkyl_ssq", "Di_O_Alkyl_ssq", "Aromatic_ssq", "Phenolic_ssq", "Amide_to_Ketone_ssq", "Sum_ssq")
        
        sample_stats <- cbind(sample_stats,ssq_sample)
        
      } else if (ecosys == "Terr_Baldock") {

        stdmat <- std_nmr(ecosys = "Terr_Baldock")
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30,  NMRmeth = "MMMFixN")

        ## Elementar ratios (relative to C)

        Cmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]+NMR.end$Protein*stdmat[[1]][9,1]+NMR.end$Lignin*stdmat[[1]][9,3]+
                             NMR.end$Lipid*stdmat[[1]][9,4]+NMR.end$Carbonyl*stdmat[[1]][9,5]+NMR.end$Char*stdmat[[1]][9,6])

        Nmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]+NMR.end$Protein*stdmat[[1]][10,1]+NMR.end$Lignin*stdmat[[1]][10,3]+
                             NMR.end$Lipid*stdmat[[1]][10,4]+NMR.end$Carbonyl*stdmat[[1]][10,5]+NMR.end$Char*stdmat[[1]][10,6])

        Hmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]+NMR.end$Protein*stdmat[[1]][11,1]+NMR.end$Lignin*stdmat[[1]][11,3]+
                             NMR.end$Lipid*stdmat[[1]][11,4]+NMR.end$Carbonyl*stdmat[[1]][11,5]+NMR.end$Char*stdmat[[1]][11,6])

        Omol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]+NMR.end$Protein*stdmat[[1]][12,1]+NMR.end$Lignin*stdmat[[1]][12,3]+
                             NMR.end$Lipid*stdmat[[1]][12,4]+NMR.end$Carbonyl*stdmat[[1]][12,5]+NMR.end$Char*stdmat[[1]][12,6])

        Cwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]*12.0107+NMR.end$Protein*stdmat[[1]][9,1]*12.0107+NMR.end$Lignin*stdmat[[1]][9,3]*12.0107+
                             NMR.end$Lipid*stdmat[[1]][9,4]*12.0107+NMR.end$Carbonyl*stdmat[[1]][9,5]*12.0107+NMR.end$Char*stdmat[[1]][9,6]*12.0107)

        Nwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]*14.0067+NMR.end$Protein*stdmat[[1]][10,1]*14.0067+NMR.end$Lignin*stdmat[[1]][10,3]*14.0067+
                             NMR.end$Lipid*stdmat[[1]][10,4]*14.0067+NMR.end$Carbonyl*stdmat[[1]][10,5]*14.0067+NMR.end$Char*stdmat[[1]][10,6]*14.0067)

        Hwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]*1.00794+NMR.end$Protein*stdmat[[1]][11,1]*1.00794+NMR.end$Lignin*stdmat[[1]][11,3]*1.00794+
                             NMR.end$Lipid*stdmat[[1]][11,4]*1.00794+NMR.end$Carbonyl*stdmat[[1]][11,5]*1.00794+NMR.end$Char*stdmat[[1]][11,6]*1.00794)

        Owgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]*15.994+NMR.end$Protein*stdmat[[1]][12,1]*15.994+NMR.end$Lignin*stdmat[[1]][12,3]*15.994+
                             NMR.end$Lipid*stdmat[[1]][12,4]*15.994+NMR.end$Carbonyl*stdmat[[1]][12,5]*15.994+NMR.end$Char*stdmat[[1]][12,6]*15.994)

        swgt <- c(Cwgt + Nwgt + Hwgt +Owgt)

        NOSC <- as.numeric(4+((2*Omol+3*Nmol-1*Hmol-4*Cmol)/Cmol))

        ## Final result

        NMR.end <- cbind(NMR.end, Cmol, Nmol, Hmol, Omol, Cwgt/swgt, Nwgt/swgt, Hwgt/swgt, Owgt/swgt, NOSC)
        
        ## Back calculated NMR results
        
        ##Alklyl C calculation
        Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][1,2]+NMR.end$Protein*stdmat[[1]][1,1]+NMR.end$Lignin*stdmat[[1]][1,3]+NMR.end$Lipid*stdmat[[1]][1,4]+
                               NMR.end$Carbonyl*stdmat[[1]][1,5]+NMR.end$Char*stdmat[[1]][1,6], c("Alkyl"))
        
        ##N_Alkyl_Methoxyl C calculation
        N_Alkyl_Methoxyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][2,2]+NMR.end$Protein*stdmat[[1]][2,1]+NMR.end$Lignin*stdmat[[1]][2,3]+NMR.end$Lipid*stdmat[[1]][2,4]+
                                          NMR.end$Carbonyl*stdmat[[1]][2,5]+NMR.end$Char*stdmat[[1]][2,6], c("N_Alkyl_Methoxyl"))
        
        ##O-Alkyl C calculation
        O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][3,2]+NMR.end$Protein*stdmat[[1]][3,1]+NMR.end$Lignin*stdmat[[1]][3,3]+NMR.end$Lipid*stdmat[[1]][3,4]+
                                 NMR.end$Carbonyl*stdmat[[1]][3,5]+NMR.end$Char*stdmat[[1]][3,6], c("O-Alkyl"))
        
        ##Di_O_Alkyl C calculation
        Di_O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][4,2]+NMR.end$Protein*stdmat[[1]][4,1]+NMR.end$Lignin*stdmat[[1]][4,3]+NMR.end$Lipid*stdmat[[1]][4,4]+
                                    NMR.end$Carbonyl*stdmat[[1]][4,5]+NMR.end$Char*stdmat[[1]][4,6], c("Di_O_Alkyl"))
        
        ##Aromatic C calculation
        Aromatic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][5,2]+NMR.end$Protein*stdmat[[1]][5,1]+NMR.end$Lignin*stdmat[[1]][5,3]+NMR.end$Lipid*stdmat[[1]][5,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][5,5]+NMR.end$Char*stdmat[[1]][5,6], c("Aromatic"))
        
        ##Phenolic C calculation
        Phenolic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][6,2]+NMR.end$Protein*stdmat[[1]][6,1]+NMR.end$Lignin*stdmat[[1]][6,3]+NMR.end$Lipid*stdmat[[1]][6,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][6,5]+NMR.end$Char*stdmat[[1]][6,6], c("Phenolic"))
        
        ##Amide_Carboxylic C calculation
        Amide_to_Ketone_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][7,2]+NMR.end$Protein*stdmat[[1]][7,1]+NMR.end$Lignin*stdmat[[1]][7,3]+NMR.end$Lipid*stdmat[[1]][7,4]+
                                         NMR.end$Carbonyl*stdmat[[1]][7,5]+NMR.end$Char*stdmat[[1]][7,6], c("Amide_to_Ketone"))
        
        sum_m <- setNames(Alkyl_m + N_Alkyl_Methoxyl_m + O_Alkyl_m + Di_O_Alkyl_m + Aromatic_m + Phenolic_m + Amide_to_Ketone_m, c("Sum"))
        
        sum_c <- sum(sampleintegraljoin)
        
        sampleintegraljoin <-rbind(sampleintegraljoin, sum_c)
        
        sample_stats <- data.frame(Alkyl_m, N_Alkyl_Methoxyl_m, O_Alkyl_m, Di_O_Alkyl_m, Aromatic_m, Phenolic_m, Amide_to_Ketone_m, sum_m)
        
        nmrrest <- NULL
        for (i in 1:nrow(sample_stats)) {
          nmrrestt <- c(sampleintegraljoin)
          nmrrest <- rbind(nmrrest, nmrrestt)
        }
        colnames(nmrrest) <-  c("Alkyl", "N_Alkyl_Methoxyl", "O-Alkyl", "Di_O_Alkyl", "Aromatic", "Phenolic", "Amide_to_Ketone", "Sum")
        
        sample_stats <- cbind(sample_stats,nmrrest)
        ssq_sample <- data.frame((sample_stats$Alkyl-sample_stats$Alkyl_m)^2, (sample_stats$N_Alkyl_Methoxyl-sample_stats$N_Alkyl_Methoxyl_m)^2,
                                 (sample_stats$`O-Alkyl` -sample_stats$O_Alkyl_m)^2, (sample_stats$Di_O_Alkyl-sample_stats$Di_O_Alkyl_m)^2, 
                                 (sample_stats$Aromatic -sample_stats$Aromatic_m)^2, (sample_stats$Phenolic -sample_stats$Phenolic_m)^2,
                                 (sample_stats$Amide_to_Ketone -sample_stats$Amide_to_Ketone_m)^2, (sample_stats$Sum -sample_stats$sum_m)^2)
        
        colnames(ssq_sample) <-  c("Alkyl_ssq", "N_Alkyl_Methoxyl_ssq", "O-Alkyl_ssq", "Di_O_Alkyl_ssq", "Aromatic_ssq", "Phenolic_ssq", "Amide_to_Ketone_ssq", "Sum_ssq")
        
        sample_stats <- cbind(sample_stats,ssq_sample)
        
      } else if (ecosys == "Aqua_Nelson") {

        stdmat <- std_nmr(ecosys = "Aqua_Nelson")
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30,  NMRmeth = "MMMFixN")

        ## Elementar ratios (relative to C)

        Cmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]+NMR.end$Protein*stdmat[[1]][9,1]+NMR.end$Lignin*stdmat[[1]][9,3]+
                             NMR.end$Lipid*stdmat[[1]][9,4]+NMR.end$Carbonyl*stdmat[[1]][9,5]+NMR.end$Char*stdmat[[1]][9,6])

        Nmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]+NMR.end$Protein*stdmat[[1]][10,1]+NMR.end$Lignin*stdmat[[1]][10,3]+
                             NMR.end$Lipid*stdmat[[1]][10,4]+NMR.end$Carbonyl*stdmat[[1]][10,5]+NMR.end$Char*stdmat[[1]][10,6])

        Hmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]+NMR.end$Protein*stdmat[[1]][11,1]+NMR.end$Lignin*stdmat[[1]][11,3]+
                             NMR.end$Lipid*stdmat[[1]][11,4]+NMR.end$Carbonyl*stdmat[[1]][11,5]+NMR.end$Char*stdmat[[1]][11,6])

        Omol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]+NMR.end$Protein*stdmat[[1]][12,1]+NMR.end$Lignin*stdmat[[1]][12,3]+
                             NMR.end$Lipid*stdmat[[1]][12,4]+NMR.end$Carbonyl*stdmat[[1]][12,5]+NMR.end$Char*stdmat[[1]][12,6])

        Cwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]*12.0107+NMR.end$Protein*stdmat[[1]][9,1]*12.0107+NMR.end$Lignin*stdmat[[1]][9,3]*12.0107+
                             NMR.end$Lipid*stdmat[[1]][9,4]*12.0107+NMR.end$Carbonyl*stdmat[[1]][9,5]*12.0107+NMR.end$Char*stdmat[[1]][9,6]*12.0107)

        Nwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]*14.0067+NMR.end$Protein*stdmat[[1]][10,1]*14.0067+NMR.end$Lignin*stdmat[[1]][10,3]*14.0067+
                             NMR.end$Lipid*stdmat[[1]][10,4]*14.0067+NMR.end$Carbonyl*stdmat[[1]][10,5]*14.0067+NMR.end$Char*stdmat[[1]][10,6]*14.0067)

        Hwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]*1.00794+NMR.end$Protein*stdmat[[1]][11,1]*1.00794+NMR.end$Lignin*stdmat[[1]][11,3]*1.00794+
                             NMR.end$Lipid*stdmat[[1]][11,4]*1.00794+NMR.end$Carbonyl*stdmat[[1]][11,5]*1.00794+NMR.end$Char*stdmat[[1]][11,6]*1.00794)

        Owgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]*15.994+NMR.end$Protein*stdmat[[1]][12,1]*15.994+NMR.end$Lignin*stdmat[[1]][12,3]*15.994+
                             NMR.end$Lipid*stdmat[[1]][12,4]*15.994+NMR.end$Carbonyl*stdmat[[1]][12,5]*15.994+NMR.end$Char*stdmat[[1]][12,6]*15.994)

        swgt <- c(Cwgt + Nwgt + Hwgt +Owgt)

        NOSC <- as.numeric(4+((2*Omol+3*Nmol-1*Hmol-4*Cmol)/Cmol))

        ## Final result

        NMR.end <- cbind(NMR.end, Cmol, Nmol, Hmol, Omol, Cwgt/swgt, Nwgt/swgt, Hwgt/swgt, Owgt/swgt, NOSC)
        
        ## Back calculated NMR results
        
        ##Alklyl C calculation
        Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][1,2]+NMR.end$Protein*stdmat[[1]][1,1]+NMR.end$Lignin*stdmat[[1]][1,3]+NMR.end$Lipid*stdmat[[1]][1,4]+
                               NMR.end$Carbonyl*stdmat[[1]][1,5]+NMR.end$Char*stdmat[[1]][1,6], c("Alkyl"))
        
        ##N_Alkyl_Methoxyl C calculation
        N_Alkyl_Methoxyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][2,2]+NMR.end$Protein*stdmat[[1]][2,1]+NMR.end$Lignin*stdmat[[1]][2,3]+NMR.end$Lipid*stdmat[[1]][2,4]+
                                          NMR.end$Carbonyl*stdmat[[1]][2,5]+NMR.end$Char*stdmat[[1]][2,6], c("N_Alkyl_Methoxyl"))
        
        ##O-Alkyl C calculation
        O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][3,2]+NMR.end$Protein*stdmat[[1]][3,1]+NMR.end$Lignin*stdmat[[1]][3,3]+NMR.end$Lipid*stdmat[[1]][3,4]+
                                 NMR.end$Carbonyl*stdmat[[1]][3,5]+NMR.end$Char*stdmat[[1]][3,6], c("O-Alkyl"))
        
        ##Di_O_Alkyl C calculation
        Di_O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][4,2]+NMR.end$Protein*stdmat[[1]][4,1]+NMR.end$Lignin*stdmat[[1]][4,3]+NMR.end$Lipid*stdmat[[1]][4,4]+
                                    NMR.end$Carbonyl*stdmat[[1]][4,5]+NMR.end$Char*stdmat[[1]][4,6], c("Di_O_Alkyl"))
        
        ##Aromatic C calculation
        Aromatic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][5,2]+NMR.end$Protein*stdmat[[1]][5,1]+NMR.end$Lignin*stdmat[[1]][5,3]+NMR.end$Lipid*stdmat[[1]][5,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][5,5]+NMR.end$Char*stdmat[[1]][5,6], c("Aromatic"))
        
        ##Phenolic C calculation
        Phenolic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][6,2]+NMR.end$Protein*stdmat[[1]][6,1]+NMR.end$Lignin*stdmat[[1]][6,3]+NMR.end$Lipid*stdmat[[1]][6,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][6,5]+NMR.end$Char*stdmat[[1]][6,6], c("Phenolic"))
        
        ##Amide_Carboxylic C calculation
        Amide_to_Ketone_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][7,2]+NMR.end$Protein*stdmat[[1]][7,1]+NMR.end$Lignin*stdmat[[1]][7,3]+NMR.end$Lipid*stdmat[[1]][7,4]+
                                         NMR.end$Carbonyl*stdmat[[1]][7,5]+NMR.end$Char*stdmat[[1]][7,6], c("Amide_to_Ketone"))
        
        sum_m <- setNames(Alkyl_m + N_Alkyl_Methoxyl_m + O_Alkyl_m + Di_O_Alkyl_m + Aromatic_m + Phenolic_m + Amide_to_Ketone_m, c("Sum"))
        
        sum_c <- sum(sampleintegraljoin)
        
        sampleintegraljoin <-rbind(sampleintegraljoin, sum_c)
        
        sample_stats <- data.frame(Alkyl_m, N_Alkyl_Methoxyl_m, O_Alkyl_m, Di_O_Alkyl_m, Aromatic_m, Phenolic_m, Amide_to_Ketone_m, sum_m)
        
        nmrrest <- NULL
        for (i in 1:nrow(sample_stats)) {
          nmrrestt <- c(sampleintegraljoin)
          nmrrest <- rbind(nmrrest, nmrrestt)
        }
        colnames(nmrrest) <-  c("Alkyl", "N_Alkyl_Methoxyl", "O-Alkyl", "Di_O_Alkyl", "Aromatic", "Phenolic", "Amide_to_Ketone", "Sum")
        
        sample_stats <- cbind(sample_stats,nmrrest)
        ssq_sample <- data.frame((sample_stats$Alkyl-sample_stats$Alkyl_m)^2, (sample_stats$N_Alkyl_Methoxyl-sample_stats$N_Alkyl_Methoxyl_m)^2,
                                 (sample_stats$`O-Alkyl` -sample_stats$O_Alkyl_m)^2, (sample_stats$Di_O_Alkyl-sample_stats$Di_O_Alkyl_m)^2, 
                                 (sample_stats$Aromatic -sample_stats$Aromatic_m)^2, (sample_stats$Phenolic -sample_stats$Phenolic_m)^2,
                                 (sample_stats$Amide_to_Ketone -sample_stats$Amide_to_Ketone_m)^2, (sample_stats$Sum -sample_stats$sum_m)^2)
        
        colnames(ssq_sample) <-  c("Alkyl_ssq", "N_Alkyl_Methoxyl_ssq", "O-Alkyl_ssq", "Di_O_Alkyl_ssq", "Aromatic_ssq", "Phenolic_ssq", "Amide_to_Ketone_ssq", "Sum_ssq")
        
        sample_stats <- cbind(sample_stats,ssq_sample)
        
      } else if (ecosys == "Aqua_Baldock") {

        stdmat <- std_nmr(ecosys = "Aqua_Baldock")
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30,  NMRmeth = "MMMFixN")

        ## Elementar ratios (relative to C)

        Cmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]+NMR.end$Protein*stdmat[[1]][9,1]+NMR.end$Lignin*stdmat[[1]][9,3]+
                             NMR.end$Lipid*stdmat[[1]][9,4]+NMR.end$Carbonyl*stdmat[[1]][9,5]+NMR.end$Char*stdmat[[1]][9,6])

        Nmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]+NMR.end$Protein*stdmat[[1]][10,1]+NMR.end$Lignin*stdmat[[1]][10,3]+
                             NMR.end$Lipid*stdmat[[1]][10,4]+NMR.end$Carbonyl*stdmat[[1]][10,5]+NMR.end$Char*stdmat[[1]][10,6])

        Hmol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]+NMR.end$Protein*stdmat[[1]][11,1]+NMR.end$Lignin*stdmat[[1]][11,3]+
                             NMR.end$Lipid*stdmat[[1]][11,4]+NMR.end$Carbonyl*stdmat[[1]][11,5]+NMR.end$Char*stdmat[[1]][11,6])

        Omol <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]+NMR.end$Protein*stdmat[[1]][12,1]+NMR.end$Lignin*stdmat[[1]][12,3]+
                             NMR.end$Lipid*stdmat[[1]][12,4]+NMR.end$Carbonyl*stdmat[[1]][12,5]+NMR.end$Char*stdmat[[1]][12,6])

        Cwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][9,2]*12.0107+NMR.end$Protein*stdmat[[1]][9,1]*12.0107+NMR.end$Lignin*stdmat[[1]][9,3]*12.0107+
                             NMR.end$Lipid*stdmat[[1]][9,4]*12.0107+NMR.end$Carbonyl*stdmat[[1]][9,5]*12.0107+NMR.end$Char*stdmat[[1]][9,6]*12.0107)

        Nwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][10,2]*14.0067+NMR.end$Protein*stdmat[[1]][10,1]*14.0067+NMR.end$Lignin*stdmat[[1]][10,3]*14.0067+
                             NMR.end$Lipid*stdmat[[1]][10,4]*14.0067+NMR.end$Carbonyl*stdmat[[1]][10,5]*14.0067+NMR.end$Char*stdmat[[1]][10,6]*14.0067)

        Hwgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][11,2]*1.00794+NMR.end$Protein*stdmat[[1]][11,1]*1.00794+NMR.end$Lignin*stdmat[[1]][11,3]*1.00794+
                             NMR.end$Lipid*stdmat[[1]][11,4]*1.00794+NMR.end$Carbonyl*stdmat[[1]][11,5]*1.00794+NMR.end$Char*stdmat[[1]][11,6]*1.00794)

        Owgt <- as.numeric(NMR.end$Carbohydrates*stdmat[[1]][12,2]*15.994+NMR.end$Protein*stdmat[[1]][12,1]*15.994+NMR.end$Lignin*stdmat[[1]][12,3]*15.994+
                             NMR.end$Lipid*stdmat[[1]][12,4]*15.994+NMR.end$Carbonyl*stdmat[[1]][12,5]*15.994+NMR.end$Char*stdmat[[1]][12,6]*15.994)

        swgt <- c(Cwgt + Nwgt + Hwgt +Owgt)

        NOSC <- as.numeric(4+((2*Omol+3*Nmol-1*Hmol-4*Cmol)/Cmol))

        ## Final result

        NMR.end <- cbind(NMR.end, Cmol, Nmol, Hmol, Omol, Cwgt/swgt, Nwgt/swgt, Hwgt/swgt, Owgt/swgt, NOSC)
        
        ## Back calculated NMR results
        
        ##Alklyl C calculation
        Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][1,2]+NMR.end$Protein*stdmat[[1]][1,1]+NMR.end$Lignin*stdmat[[1]][1,3]+NMR.end$Lipid*stdmat[[1]][1,4]+
                               NMR.end$Carbonyl*stdmat[[1]][1,5]+NMR.end$Char*stdmat[[1]][1,6], c("Alkyl"))
        
        ##N_Alkyl_Methoxyl C calculation
        N_Alkyl_Methoxyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][2,2]+NMR.end$Protein*stdmat[[1]][2,1]+NMR.end$Lignin*stdmat[[1]][2,3]+NMR.end$Lipid*stdmat[[1]][2,4]+
                                          NMR.end$Carbonyl*stdmat[[1]][2,5]+NMR.end$Char*stdmat[[1]][2,6], c("N_Alkyl_Methoxyl"))
        
        ##O-Alkyl C calculation
        O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][3,2]+NMR.end$Protein*stdmat[[1]][3,1]+NMR.end$Lignin*stdmat[[1]][3,3]+NMR.end$Lipid*stdmat[[1]][3,4]+
                                 NMR.end$Carbonyl*stdmat[[1]][3,5]+NMR.end$Char*stdmat[[1]][3,6], c("O-Alkyl"))
        
        ##Di_O_Alkyl C calculation
        Di_O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][4,2]+NMR.end$Protein*stdmat[[1]][4,1]+NMR.end$Lignin*stdmat[[1]][4,3]+NMR.end$Lipid*stdmat[[1]][4,4]+
                                    NMR.end$Carbonyl*stdmat[[1]][4,5]+NMR.end$Char*stdmat[[1]][4,6], c("Di_O_Alkyl"))
        
        ##Aromatic C calculation
        Aromatic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][5,2]+NMR.end$Protein*stdmat[[1]][5,1]+NMR.end$Lignin*stdmat[[1]][5,3]+NMR.end$Lipid*stdmat[[1]][5,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][5,5]+NMR.end$Char*stdmat[[1]][5,6], c("Aromatic"))
        
        ##Phenolic C calculation
        Phenolic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][6,2]+NMR.end$Protein*stdmat[[1]][6,1]+NMR.end$Lignin*stdmat[[1]][6,3]+NMR.end$Lipid*stdmat[[1]][6,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][6,5]+NMR.end$Char*stdmat[[1]][6,6], c("Phenolic"))
        
        ##Amide_Carboxylic C calculation
        Amide_to_Ketone_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][7,2]+NMR.end$Protein*stdmat[[1]][7,1]+NMR.end$Lignin*stdmat[[1]][7,3]+NMR.end$Lipid*stdmat[[1]][7,4]+
                                         NMR.end$Carbonyl*stdmat[[1]][7,5]+NMR.end$Char*stdmat[[1]][7,6], c("Amide_to_Ketone"))
        
        sum_m <- setNames(Alkyl_m + N_Alkyl_Methoxyl_m + O_Alkyl_m + Di_O_Alkyl_m + Aromatic_m + Phenolic_m + Amide_to_Ketone_m, c("Sum"))
        
        sum_c <- sum(sampleintegraljoin)
        
        sampleintegraljoin <-rbind(sampleintegraljoin, sum_c)
        
        sample_stats <- data.frame(Alkyl_m, N_Alkyl_Methoxyl_m, O_Alkyl_m, Di_O_Alkyl_m, Aromatic_m, Phenolic_m, Amide_to_Ketone_m, sum_m)
        
        nmrrest <- NULL
        for (i in 1:nrow(sample_stats)) {
          nmrrestt <- c(sampleintegraljoin)
          nmrrest <- rbind(nmrrest, nmrrestt)
        }
        colnames(nmrrest) <-  c("Alkyl", "N_Alkyl_Methoxyl", "O-Alkyl", "Di_O_Alkyl", "Aromatic", "Phenolic", "Amide_to_Ketone", "Sum")
        
        sample_stats <- cbind(sample_stats,nmrrest)
        ssq_sample <- data.frame((sample_stats$Alkyl-sample_stats$Alkyl_m)^2, (sample_stats$N_Alkyl_Methoxyl-sample_stats$N_Alkyl_Methoxyl_m)^2,
                                 (sample_stats$`O-Alkyl` -sample_stats$O_Alkyl_m)^2, (sample_stats$Di_O_Alkyl-sample_stats$Di_O_Alkyl_m)^2, 
                                 (sample_stats$Aromatic -sample_stats$Aromatic_m)^2, (sample_stats$Phenolic -sample_stats$Phenolic_m)^2,
                                 (sample_stats$Amide_to_Ketone -sample_stats$Amide_to_Ketone_m)^2, (sample_stats$Sum -sample_stats$sum_m)^2)
        
        colnames(ssq_sample) <-  c("Alkyl_ssq", "N_Alkyl_Methoxyl_ssq", "O-Alkyl_ssq", "Di_O_Alkyl_ssq", "Aromatic_ssq", "Phenolic_ssq", "Amide_to_Ketone_ssq", "Sum_ssq")
        
        sample_stats <- cbind(sample_stats,ssq_sample)
        
      } else if (ecosys == "mod") {

        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = mod_std, amoSTD = 6, best.fits = 30,  NMRmeth = "MMMFixN")

        ## Elementar ratios (relative to C)

        Cmol <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][9,2]+NMR.end$Protein*mod_std[[1]][9,1]+NMR.end$Lignin*mod_std[[1]][9,3]+
                             NMR.end$Lipid*mod_std[[1]][9,4]+NMR.end$Carbonyl*mod_std[[1]][9,5]+NMR.end$Char*mod_std[[1]][9,6])

        Nmol <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][10,2]+NMR.end$Protein*mod_std[[1]][10,1]+NMR.end$Lignin*mod_std[[1]][10,3]+
                             NMR.end$Lipid*mod_std[[1]][10,4]+NMR.end$Carbonyl*mod_std[[1]][10,5]+NMR.end$Char*mod_std[[1]][10,6])

        Hmol <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][11,2]+NMR.end$Protein*mod_std[[1]][11,1]+NMR.end$Lignin*mod_std[[1]][11,3]+
                             NMR.end$Lipid*mod_std[[1]][11,4]+NMR.end$Carbonyl*mod_std[[1]][11,5]+NMR.end$Char*mod_std[[1]][11,6])

        Omol <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][12,2]+NMR.end$Protein*mod_std[[1]][12,1]+NMR.end$Lignin*mod_std[[1]][12,3]+
                             NMR.end$Lipid*mod_std[[1]][12,4]+NMR.end$Carbonyl*mod_std[[1]][12,5]+NMR.end$Char*mod_std[[1]][12,6])

        Cwgt <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][9,2]*12.0107+NMR.end$Protein*mod_std[[1]][9,1]*12.0107+NMR.end$Lignin*mod_std[[1]][9,3]*12.0107+
                             NMR.end$Lipid*mod_std[[1]][9,4]*12.0107+NMR.end$Carbonyl*mod_std[[1]][9,5]*12.0107+NMR.end$Char*mod_std[[1]][9,6]*12.0107)

        Nwgt <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][10,2]*14.0067+NMR.end$Protein*mod_std[[1]][10,1]*14.0067+NMR.end$Lignin*mod_std[[1]][10,3]*14.0067+
                             NMR.end$Lipid*mod_std[[1]][10,4]*14.0067+NMR.end$Carbonyl*mod_std[[1]][10,5]*14.0067+NMR.end$Char*mod_std[[1]][10,6]*14.0067)

        Hwgt <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][11,2]*1.00794+NMR.end$Protein*mod_std[[1]][11,1]*1.00794+NMR.end$Lignin*mod_std[[1]][11,3]*1.00794+
                             NMR.end$Lipid*mod_std[[1]][11,4]*1.00794+NMR.end$Carbonyl*mod_std[[1]][11,5]*1.00794+NMR.end$Char*mod_std[[1]][11,6]*1.00794)

        Owgt <- as.numeric(NMR.end$Carbohydrates*mod_std[[1]][12,2]*15.994+NMR.end$Protein*mod_std[[1]][12,1]*15.994+NMR.end$Lignin*mod_std[[1]][12,3]*15.994+
                             NMR.end$Lipid*mod_std[[1]][12,4]*15.994+NMR.end$Carbonyl*mod_std[[1]][12,5]*15.994+NMR.end$Char*mod_std[[1]][12,6]*15.994)

        swgt <- c(Cwgt + Nwgt + Hwgt +Owgt)

        NOSC <- as.numeric(4+((2*Omol+3*Nmol-1*Hmol-4*Cmol)/Cmol))

        ## Final result

        NMR.end <- cbind(NMR.end, Cmol, Nmol, Hmol, Omol, Cwgt/swgt, Nwgt/swgt, Hwgt/swgt, Owgt/swgt, NOSC)
        
        ## Back calculated NMR results
        
        ##Alklyl C calculation
        Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][1,2]+NMR.end$Protein*stdmat[[1]][1,1]+NMR.end$Lignin*stdmat[[1]][1,3]+NMR.end$Lipid*stdmat[[1]][1,4]+
                               NMR.end$Carbonyl*stdmat[[1]][1,5]+NMR.end$Char*stdmat[[1]][1,6], c("Alkyl"))
        
        ##N_Alkyl_Methoxyl C calculation
        N_Alkyl_Methoxyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][2,2]+NMR.end$Protein*stdmat[[1]][2,1]+NMR.end$Lignin*stdmat[[1]][2,3]+NMR.end$Lipid*stdmat[[1]][2,4]+
                                          NMR.end$Carbonyl*stdmat[[1]][2,5]+NMR.end$Char*stdmat[[1]][2,6], c("N_Alkyl_Methoxyl"))
        
        ##O-Alkyl C calculation
        O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][3,2]+NMR.end$Protein*stdmat[[1]][3,1]+NMR.end$Lignin*stdmat[[1]][3,3]+NMR.end$Lipid*stdmat[[1]][3,4]+
                                 NMR.end$Carbonyl*stdmat[[1]][3,5]+NMR.end$Char*stdmat[[1]][3,6], c("O-Alkyl"))
        
        ##Di_O_Alkyl C calculation
        Di_O_Alkyl_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][4,2]+NMR.end$Protein*stdmat[[1]][4,1]+NMR.end$Lignin*stdmat[[1]][4,3]+NMR.end$Lipid*stdmat[[1]][4,4]+
                                    NMR.end$Carbonyl*stdmat[[1]][4,5]+NMR.end$Char*stdmat[[1]][4,6], c("Di_O_Alkyl"))
        
        ##Aromatic C calculation
        Aromatic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][5,2]+NMR.end$Protein*stdmat[[1]][5,1]+NMR.end$Lignin*stdmat[[1]][5,3]+NMR.end$Lipid*stdmat[[1]][5,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][5,5]+NMR.end$Char*stdmat[[1]][5,6], c("Aromatic"))
        
        ##Phenolic C calculation
        Phenolic_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][6,2]+NMR.end$Protein*stdmat[[1]][6,1]+NMR.end$Lignin*stdmat[[1]][6,3]+NMR.end$Lipid*stdmat[[1]][6,4]+
                                  NMR.end$Carbonyl*stdmat[[1]][6,5]+NMR.end$Char*stdmat[[1]][6,6], c("Phenolic"))
        
        ##Amide_Carboxylic C calculation
        Amide_to_Ketone_m <-  setNames(NMR.end$Carbohydrates*stdmat[[1]][7,2]+NMR.end$Protein*stdmat[[1]][7,1]+NMR.end$Lignin*stdmat[[1]][7,3]+NMR.end$Lipid*stdmat[[1]][7,4]+
                                         NMR.end$Carbonyl*stdmat[[1]][7,5]+NMR.end$Char*stdmat[[1]][7,6], c("Amide_to_Ketone"))
        
        sum_m <- setNames(Alkyl_m + N_Alkyl_Methoxyl_m + O_Alkyl_m + Di_O_Alkyl_m + Aromatic_m + Phenolic_m + Amide_to_Ketone_m, c("Sum"))
        
        sum_c <- sum(sampleintegraljoin)
        
        sampleintegraljoin <-rbind(sampleintegraljoin, sum_c)
        
        sample_stats <- data.frame(Alkyl_m, N_Alkyl_Methoxyl_m, O_Alkyl_m, Di_O_Alkyl_m, Aromatic_m, Phenolic_m, Amide_to_Ketone_m, sum_m)
        
        nmrrest <- NULL
        for (i in 1:nrow(sample_stats)) {
          nmrrestt <- c(sampleintegraljoin)
          nmrrest <- rbind(nmrrest, nmrrestt)
        }
        colnames(nmrrest) <-  c("Alkyl", "N_Alkyl_Methoxyl", "O-Alkyl", "Di_O_Alkyl", "Aromatic", "Phenolic", "Amide_to_Ketone", "Sum")
        
        sample_stats <- cbind(sample_stats,nmrrest)
        ssq_sample <- data.frame((sample_stats$Alkyl-sample_stats$Alkyl_m)^2, (sample_stats$N_Alkyl_Methoxyl-sample_stats$N_Alkyl_Methoxyl_m)^2,
                                 (sample_stats$`O-Alkyl` -sample_stats$O_Alkyl_m)^2, (sample_stats$Di_O_Alkyl-sample_stats$Di_O_Alkyl_m)^2, 
                                 (sample_stats$Aromatic -sample_stats$Aromatic_m)^2, (sample_stats$Phenolic -sample_stats$Phenolic_m)^2,
                                 (sample_stats$Amide_to_Ketone -sample_stats$Amide_to_Ketone_m)^2, (sample_stats$Sum -sample_stats$sum_m)^2)
        
        colnames(ssq_sample) <-  c("Alkyl_ssq", "N_Alkyl_Methoxyl_ssq", "O-Alkyl_ssq", "Di_O_Alkyl_ssq", "Aromatic_ssq", "Phenolic_ssq", "Amide_to_Ketone_ssq", "Sum_ssq")
        
        sample_stats <- cbind(sample_stats,ssq_sample)        
      }
    }
  }
  #citation  <- setNames(data.frame(matrix(ncol = 1, nrow= nrow(NMR.end))), c("Plz cite this work as an incentive to its curation"))
  #NMR.end <- cbind(NMR.end,citation)
  
  if (stats == FALSE) {
    
    return(NMR.end)
    
  } else if (stats == TRUE) {
    
    return(sample_stats)
    }
}
