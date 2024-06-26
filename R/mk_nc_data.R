#' Create .csv file for CN data
#'
#' This function copies the spectra read using the read_spec function and creates
#' a .csv file with a colum with their names and two empty columns where the user must add the C and N values.
#' Thereafter the file is read with the function nc_data
#'
#' @param raw.spec The uploaded spectra read using the read_spec function
#' @returns A data frame with three columns, one containing the names extracted from the raw.spec, and two columns to be filled manually with the carbon and nitrogen values.
#' @keywords CN file
#' @export
#' @importFrom utils read.table
#' @examples
#' ## any .txt file as output from BRUKER

mk_nc_data <- function(raw.spec) {

    name.end  <- NULL

    nc.end <- setNames(data.frame(matrix(ncol = 2, nrow = length(raw.spec))), c("C", "N"))

      for (i in 1:length(raw.spec)) {
        name.end[[i]] <- raw.spec[[i]]$name

      }

    name <- data.frame(unlist(name.end))

    name <- setNames(name, c("name"))

    nc.end <- cbind(name,nc.end)

    return(nc.end)
    ## close function
}

