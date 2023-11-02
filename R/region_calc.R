#' Functional groups calculation
#'
#' This function loads, integrates and calculates the functional group distribution from the raw spectra.
#' Produces also the molecular mixing model fitting if NC data is provided.
#' Output is a list with the raw data, integrals and corrected spectra.
#' @param raw.spec List of files already loaded with read_raw_spec()
#' @param batch_nmr Vector with file names, default
#' @param NMRmeth Regions to be integrated.
#' @param FixNC TRUE or FALSE, for fixing or not the NC ratio on the sample fitting.
#' Default is spinning side bands, other methods available include: Bonanomi ("Bonanomi") and Molecular mixing model ("MMM").
#' @param ecosys Standards to be used for the MMM, can be Terrestrial("Terr_Nelson" or "Terr_Baldock") or Aquatic ("Aqua_Nelson" or "Aqua_Baldock")
#' @param ncdata Initial correction and normalization parameters
#' @param NMR_field Magnetic field of the NMR
#' @param NMR_rotation Rotation frequency of the sample probe in the NMR
#' @param mod_std File containing a modified NMR table
#' @keywords fitting, Molecular Mixing model, Bonanomi
#' @export
#' @examples


region_calc <- function (batch_nmr = NULL, file = NULL, NMRmeth = NULL, FixNC,
                         NMR_field=NULL, NMR_rotation=NULL, ecosys=NULL,
                         cndata = NULL, mod_std = NULL) {

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

    stop("Please choose an preset region model composition by typing 'MMM' for Molecular mixing model, 'Bonanomi' or '4region'")

    ## Start processing for an Integration of 4 regions
  } else if (NMRmeth == "4region") {

    ## Integral function
      Integral <- int_nmr (raw.spec = batch.nmr, NMRmeth = NMRmeth, NMR_field = NMR_field, NMR_rotation = NMR_rotation)

      ## Variable of end result
      NMR.end <- NULL

      ## Start loop to process all integrals
      for (i in 1:length(Integral)) {

        ## Name of sample
        name <- Integral[[i]]$name

        ## paste which sample is used now
        print(paste("Integrating Sample: ", name, ", start date: ", Sys.time(), sep = ""))

        ## create progress bar for the standards combinations
        pb <- txtProgressBar(min = 1, max = length(batch.nmr), style = 3)

        setTxtProgressBar(pb, i)

        ## Removed messages from tidy
        ## Sum of the Central integrals, without Spinning Sidebands
        suppressMessages(Integral[[i]]$data$Integral %>%
          filter(is.na(Component) == FALSE) %>%
          group_by(Component) %>%
          summarise(Integral = sum(Integral), n = n())) -> sums_Integral

        ## Removed messages from tidy
        ## Sum of the Spinning Sidebands, discriminating each side
        suppressMessages(Integral[[i]]$data$Integral %>%
          filter(is.na(Component) == TRUE) %>%
          group_by(Component, Component_ssb, ssb_ofset) %>%
          summarise(Integral = sum(Integral), n = n())) -> double_ssb

        ## Removed messages from tidy
        ## Chooses the max value between the sum of the SSB on each side, and multiplies by 2
        suppressMessages(double_ssb %>%
          group_by(Component_ssb, n) %>%
          summarize(Integral = 2*max(Integral))) ->double_ssb

        ## Removed messages from tidy
        ## Renames the SSB from Component_ssb to Component
        suppressMessages(double_ssb <- double_ssb %>%
          group_by(Component_ssb) %>%
          filter(n == max(n)) %>%
          rename(Component = Component_ssb) %>%
          ungroup())

        ## Removed messages from tidy
        ## takes the Integral table and selects the regions in which the main integral is overlaping with the spinning sidebands
        suppressMessages(Integral[[i]]$data$Integral %>%
          filter(is.na(Component) == FALSE, is.na(sbb_index) == FALSE) %>%
          group_by(Component, sbb_index) %>%
          summarise(n = n())) -> ossb

        ## Removed messages from tidy
        ## takes the Integral table and selects the regions in the spinning sidebands is not overlaping and multiplies it by -1
        suppressMessages(Integral[[i]]$data$Integral %>%
          filter(is.na(Component) == TRUE) %>%
          group_by(sbb_index) %>%
          summarise(Integral = sum(Integral)*-1,n = n())) -> ossb2

        ## takes the Integral table and selects the regions in the spinning sidebands is not overlaping and multiplies it by -1
        matching_rows <- NULL
        for (j in 1:nrow(ossb)) {
          index <- ossb$sbb_index[j]
          result <- ossb2 %>%
            filter(sbb_index == index)
            matching_rows <- bind_rows(matching_rows, result)
        }

        joined_df <- merge(ossb, matching_rows, by = "sbb_index")

        ## Creates a table with sums_Integral: the sum of the main integrals, double_ssb: the ssb area multiplied by 2, and joined_df: the equivalent to the overalpping ssb to be substracted
        table_merged <- rbind(sums_Integral[1:2], double_ssb[, c('Component', 'Integral')], joined_df[, c('Component', 'Integral')])

        ## Creates the final sum of integrals
        table_merged %>%
          group_by(Component) %>%
          summarise(Integral = sum(Integral)) -> final_integral

        ## Normalizing the sum of integrals to 100
          ## Total value of the sum
        norm <- sum(data.frame(final_integral$Integral))

        ## Normalized values
        normalized.Int <- t(final_integral$Integral/norm)*100

        ## Giving back the names to the columns
        colnames(normalized.Int) <- as.character(final_integral$Component)

        ##Final result
        NMR.end[[i]] <- data.frame(cbind(name, normalized.Int))
        close(pb)
      }
    ## Start of the processing for regions according to Bonanomi et al.
    } else if (NMRmeth == "Bonanomi") {

    ## Integral function
    Integral <- int_nmr (raw.spec = batch.nmr, NMRmeth = NMRmeth, NMR_field = NMR_field, NMR_rotation = NMR_rotation)

    ## Variable of end result
    NMR.end <- NULL

    ## Start loop to process all integrals
    for (i in 1:length(Integral)) {

      ## Name of sample
      name <- Integral[[i]]$name

      ## paste which sample is used now
      print(paste("Integrating Sample: ", name, ", start date: ", Sys.time(), sep = ""))

      ## create progress bar for the standards combinations
      pb <- txtProgressBar(min = 1, max = length(batch.nmr), style = 3)

      setTxtProgressBar(pb, i)

      ## Removed messages from tidy
      ## Sum of the Central integrals, without Spinning Sidebands
      suppressMessages(Integral[[i]]$data$Integral %>%
                         filter(is.na(Component) == FALSE) %>%
                         group_by(Component) %>%
                         summarise(Integral = sum(Integral), n = n())) -> sums_Integral

      ## Removed messages from tidy
      ## Sum of the Spinning Sidebands, discriminating each side
      suppressMessages(Integral[[i]]$data$Integral %>%
                         filter(is.na(Component) == TRUE) %>%
                         group_by(Component, Component_ssb, ssb_ofset) %>%
                         summarise(Integral = sum(Integral), n = n())) -> double_ssb

      ## Removed messages from tidy
      ## Chooses the max value between the sum of the SSB on each side, and multiplies by 2
      suppressMessages(double_ssb %>%
                         group_by(Component_ssb, n) %>%
                         summarize(Integral = 2*max(Integral))) ->double_ssb

      ## Removed messages from tidy
      ## Renames the SSB from Component_ssb to Component
      suppressMessages(double_ssb <- double_ssb %>%
                         group_by(Component_ssb) %>%
                         filter(n == max(n)) %>%
                         rename(Component = Component_ssb) %>%
                         ungroup())

      ## Removed messages from tidy
      ## takes the Integral table and selects the regions in which the main integral is overlaping with the spinning sidebands
      suppressMessages(Integral[[i]]$data$Integral %>%
                         filter(is.na(Component) == FALSE, is.na(sbb_index) == FALSE) %>%
                         group_by(Component, sbb_index) %>%
                         summarise(n = n())) -> ossb

      ## Removed messages from tidy
      ## takes the Integral table and selects the regions in the spinning sidebands is not overlaping and multiplies it by -1
      suppressMessages(Integral[[i]]$data$Integral %>%
                         filter(is.na(Component) == TRUE) %>%
                         group_by(sbb_index) %>%
                         summarise(Integral = sum(Integral)*-1,n = n())) -> ossb2

      ## takes the Integral table and selects the regions in the spinning sidebands is not overlaping and multiplies it by -1
      matching_rows <- NULL
      for (j in 1:nrow(ossb)) {
        index <- ossb$sbb_index[j]
        result <- ossb2 %>%
          filter(sbb_index == index)
        matching_rows <- bind_rows(matching_rows, result)
      }

      joined_df <- merge(ossb, matching_rows, by = "sbb_index")

      ## Creates a table with sums_Integral: the sum of the main integrals, double_ssb: the ssb area multiplied by 2, and joined_df: the equivalent to the overalpping ssb to be substracted
      table_merged <- rbind(sums_Integral[1:2], double_ssb[, c('Component', 'Integral')], joined_df[, c('Component', 'Integral')])

      ## Creates the final sum of integrals
      table_merged %>%
        group_by(Component) %>%
        summarise(Integral = sum(Integral)) -> final_integral

      ## Normalizing the sum of integrals to 100
      ## Total value of the sum
      norm <- sum(data.frame(final_integral$Integral))

      ## Normalized values
      normalized.Int <- t(final_integral$Integral/norm)*100

      ## Giving back the names to the columns
      colnames(normalized.Int) <- as.character(final_integral$Component)

      ## Final result
      NMR.end[[i]] <- data.frame(cbind(name, normalized.Int))
      close(pb)

    }
  } else if (NMRmeth == "MMM") {

    ## loop to process all samples

    Integral <- int_nmr (raw.spec = batch.nmr, NMRmeth = NMRmeth, NMR_field = NMR_field, NMR_rotation = NMR_rotation)

    ###sum of raw integrals
    NMR.end <- NULL
    raw.spec.end <- NULL

    for (i in 1:length(Integral)) {

      NCval <- c(as.numeric(cndata[[i]]$NC))
      name <- Integral[[i]]$name

      ## paste which sample is used now
      print(paste("Integrating Sample: ", name, ", start date: ", Sys.time(), sep = ""))

      ## create progress bar for the standards combinations
      pb <- txtProgressBar(min = 1, max = length(batch.nmr), style = 3)

      setTxtProgressBar(pb, i)

      suppressMessages(Integral[[i]]$data$Integral %>%
                         filter(is.na(Component) == FALSE) %>%
                         group_by(Component) %>%
                         summarise(Integral = sum(Integral), n = n())) -> sums_Integral

      ###sum of double spinning sidebands

      suppressMessages(Integral[[i]]$data$Integral %>%
                         filter(is.na(Component) == TRUE) %>%
                         group_by(Component, Component_ssb, ssb_ofset) %>%
                         summarise(Integral = sum(Integral), n = n())) -> double_ssb

      suppressMessages(double_ssb %>%
                         group_by(Component_ssb, n) %>%
                         summarize(Integral = 2*max(Integral))) ->double_ssb

      suppressMessages(double_ssb <- double_ssb %>%
                         group_by(Component_ssb) %>%
                         filter(n == max(n)) %>%
                         rename(Component = Component_ssb) %>%
                         ungroup())

      ###sum of overlaping spinning sidebands

      suppressMessages(Integral[[i]]$data$Integral %>%
                         filter(is.na(Component) == FALSE, is.na(sbb_index) == FALSE) %>%
                         group_by(Component, sbb_index) %>%
                         summarise(n = n())) -> ossb

      suppressMessages(Integral[[i]]$data$Integral %>%
                         filter(is.na(Component) == TRUE) %>%
                         group_by(sbb_index) %>%
                         summarise(Integral = sum(Integral)*-1,n = n())) -> ossb2

      matching_rows <- NULL
      for (j in 1:nrow(ossb)) {
        index <- ossb$sbb_index[j]
        #print(index)
        result <- ossb2 %>%
          filter(sbb_index == index)
        matching_rows <- bind_rows(matching_rows, result)
      }

      joined_df <- merge(ossb, matching_rows, by = "sbb_index")

      table_merged <- rbind(sums_Integral[1:2], double_ssb[, c('Component', 'Integral')], joined_df[, c('Component', 'Integral')])

      table_merged %>%
        group_by(Component) %>%
        summarise(Integral = sum(Integral)) -> final_integral

      ###finished calculation of integrals

      if (!is.null(ecosys)){

        norm <- sum(data.frame(final_integral$Integral))
        normalized.Int <- (final_integral$Integral/norm)*100
        normalized.Int <- data.frame(final_integral$Component, normalized.Int)
        #print(normalized.Int)
        #normalized.Int <- c(normalized.Int[1:7], sum(normalized.Int[8:9]))
        #normalized.Int <- data.frame(final_integral$Component,normalized.Int)

        normalized.Int <-setNames(normalized.Int,c("Component", "Integral"))

        #rownames(normalized.Int) <- as.character(final_integral$Component)

        int_NMR <- NMR_table(NMRmeth = NMRmeth)

        normalized.Int <- merge(int_NMR, normalized.Int, by =  'Component')

        normalized.Int %>%
          group_by(ID) %>%
          summarise(From = From, To = To, Integral = sum(Integral))  %>%
          arrange(From) -> normalized.Int

        #normalized.Int <- data.frame(normalized.Int$Component, normalized.Int$Integral)

        normalized.Int <- rbind(NCval, normalized.Int[,4])

        raw.spec.end[[i]] <- list("name" = name, "data" = list("Integral" = data.frame(normalized.Int[1:8,])))
        close(pb)

      } else if (is.null(ecosys))  {

        norm <- sum(data.frame(final_integral$Integral))
        normalized.Int <- t(final_integral$Integral/norm)*100
        colnames(normalized.Int) <- as.character(final_integral$Component)
        #final_integral <- data.frame(final_integral$Component, normalized.Int)
        #integral.end <- data.frame(file.name,integral.end)
        NMR.end[[i]] <- data.frame(cbind(name, normalized.Int))
        close(pb)

      }
    }
    if (is.null(ecosys)) {

      return(NMR.end)
    }

    else if (ecosys == "Terr_Nelson") {

      stdmat <- std_nmr(ecosys = "Terr_Nelson")
      if (FixNC == TRUE) {
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30, NMRmeth = NMRmeth, FixNC = TRUE)

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


      } else if (FixNC == FALSE) {
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30, NMRmeth = NMRmeth, FixNC = FALSE)

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

      }
    }
    else if (ecosys == "Terr_Baldock") {

      stdmat <- std_nmr(ecosys = "Terr_Baldock")
      if (FixNC == TRUE) {
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30, NMRmeth = NMRmeth, FixNC = TRUE)

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


      } else if (FixNC == FALSE) {
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30, NMRmeth = NMRmeth, FixNC = FALSE)

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

      }
    }

    else if (ecosys == "Aqua_Nelson") {

      stdmat <- std_nmr(ecosys = "Aqua_Nelson")
      if (FixNC == TRUE) {
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30, NMRmeth = NMRmeth, FixNC = TRUE)

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


      } else if (FixNC == FALSE) {
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30, NMRmeth = NMRmeth, FixNC = FALSE)

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

      }
    }
    else if (ecosys == "Aqua_Baldock") {

      stdmat <- std_nmr(ecosys = "Aqua_Baldock")
      if (FixNC == TRUE) {

        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30, NMRmeth = NMRmeth, FixNC = TRUE)

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


      } else if (FixNC == FALSE) {


        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30, NMRmeth = NMRmeth, FixNC = FALSE)

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

      }
    }

    else if (ecosys == "mod") {

      if (FixNC == TRUE) {
        print("here")

        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = mod_std, amoSTD = 6, best.fits = 30, NMRmeth = NMRmeth, FixNC = TRUE)

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


      } else if (FixNC == FALSE) {
        NMR.end <- fit_LCF(all.samples = raw.spec.end, all.standards = stdmat, amoSTD = 6, best.fits = 30, NMRmeth = NMRmeth, FixNC = FALSE)

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
      }
    }
  }

  return(NMR.end)
}
