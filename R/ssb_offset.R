#' Spinnning side bands ofset calculation function
#' This function calculates the spinning side band ofset for a given 13C NMR table.
#' The function returns the 13C NMR integration table to be used in the int_nmr function.
#' @param NMRmeth Regions to be integrated.
#' @param NMR_field Magnetic field of the NMR
#' @param NMR_rotation Rotation frequency of the sample probe in the NMR
#' @keywords normalization integration
#' @export
#' @import data.table
#' @examples
#' see_ofset <- ssb_offset (NMRmeth='4region', NMR_field = 200, NMR_rotation = 6800)

ssb_offset <- function(NMRmeth=NULL, NMR_field=NULL, NMR_rotation=NULL) {

  NMR_field = NMR_field

  NMR_rotation = NMR_rotation

  C_resonance =  c(NMR_field/4)

  min_MAS_rate = 1/2 *(0.0002)* C_resonance

  sb_ofset = (NMR_rotation/C_resonance)

  if (is.null(NMRmeth)) {

    stop("Please choose an preset region model composition by typing 'MMM' for Molecular mixing model, 'Bonanomi' or '4region'")

  } else if (!is.null(NMRmeth)) {

  int_NMR <- NMR_table(NMRmeth = NMRmeth)

  # comvert to data.table format
  setDT(int_NMR)
  # create high and low intervals
  From <- NULL
  To <- NULL
  Component <- NULL
  Component_ssb <- NULL
  i.Component <- NULL
  i.Component_ssb <- NULL
  ssb_ofset <- NULL
  i.ssb_ofset <- NULL
  `.` <- list

  int_NMR_high <- copy(int_NMR)[, `:=`(From = From + sb_ofset, To = To + sb_ofset, Component_ssb = paste0(Component), ssb_ofset =  "High")]
  int_NMR_low <- copy(int_NMR)[, `:=`(From = From - sb_ofset, To = To - sb_ofset, Component_ssb = paste0(Component), ssb_ofset =  "Low")]
  # create one data.table of all intervals
  all_int <- rbindlist(list(int_NMR_low, int_NMR, int_NMR_high), fill = TRUE)
  # create a new data.table with all the non-overlapping intervals
  final <- data.table(sections(breaks(as.matrix(all_int[, 2:3]))))
  setnames(final, c("From", "To"))
  # perform overlap joins
  final[all_int[is.na(Component_ssb), ], Component := i.Component, on = .(From < To, To > From )]
  final[all_int[!is.na(Component_ssb), ], Component_ssb := i.Component_ssb, on = .(From < To, To > From)]
  final[all_int[!is.na(Component_ssb), ], ssb_ofset := i.ssb_ofset, on = .(From < To, To > From )]

  int_NMR2 <- final[complete.cases(final$Component),]

  int_NMR2 <- int_NMR2[,-4:-5]

  setDT(int_NMR2)

  int_NMR_high2 <- copy(int_NMR2)[, `:=`(From = From + sb_ofset, To = To + sb_ofset, Component_ssb = paste0(Component), ssb_ofset =  "High")]
  int_NMR_low2 <- copy(int_NMR2)[, `:=`(From = From - sb_ofset, To = To - sb_ofset, Component_ssb = paste0(Component), ssb_ofset =  "Low")]

  all_int2 <- rbindlist(list(int_NMR_low2,int_NMR2, int_NMR_high2), fill = TRUE)

  # create a new data.table with all the non-overlapping intervals
  table_all <- data.table(sections(breaks(as.matrix(all_int2[, 1:2]))))
  setnames(table_all, c("From", "To"))
  # perform overlap joins
  table_all[all_int2[is.na(Component_ssb), ], Component := i.Component, on = .(From < To, To > From )]
  table_all[all_int2[!is.na(Component_ssb), ], Component_ssb := i.Component_ssb, on = .(From < To, To > From )]
  table_all[all_int2[!is.na(Component_ssb), ], ssb_ofset := i.ssb_ofset, on = .(From < To, To > From )]
  }

  table_all <- table_all %>%
    mutate(Component_index = ifelse(!is.na(Component), cumsum(!is.na(Component)), NA))

  table_all <- table_all %>%
    group_by(ssb_ofset) %>%
    mutate(
    sbb_index = ifelse(!is.na(Component_ssb), cumsum(!is.na(Component_ssb)), NA)) %>%
    ungroup()

  table_all <- table_all[, c("From", "To", "Component", "Component_index", "Component_ssb", "sbb_index", "ssb_ofset")]

  return(table_all)
}
