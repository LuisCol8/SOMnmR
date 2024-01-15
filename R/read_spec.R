#' Read spectra
#'
#' This function reads the raw file, Bruker, tab separated or coma separated
#' extracts the spectra and returns a list with name, the raw spectrum
#' @param file The raw  file
#' @param filetype The raw  file type "Bruker", .csv ("tab"), csv ("coma")
#' @keywords integration
#' @export
#' @importFrom utils read.table
#' @examples
#' ## any .txt file as output from BRUKER

read_raw_spec <- function (file = NULL, filetype = NULL) {

  if (is.null(filetype)) {

    stop("Please provide a file type (Bruker, .csv (tab), csv (coma)")

  } else if (filetype == "Bruker") {
    ## create a dummy vector to be filled
    raw.spec.end <- NULL

    ## loop over all files
    for (i in 1:length(file)) {

      ## extract file name
      file.name <- strsplit(file[i],"\\.")[[1]][1]

      ## extract file header
      file.head <- readLines(file[i], n = length(grep("#", readLines(file[i]))))

      ## read the table data
      raw.spec <- read.table(file[i], header = FALSE, skip = 1, sep = ",")

      ## set the column names accordingly
      colnames(raw.spec) <- c("Index", "raw.intensity","raw.cycle","ppm")

      ## create a list with name, element, edge and the data of the spectrum
      raw.spec.end[[i]] <- list("name" = file.name, "data" = list("raw.spec" = raw.spec))

      ## close file loop
    }

    ## .csv with only two columns

    } else if (filetype == "coma") {
      ## create a dummy vector to be filled
      raw.spec.end <- NULL

      ## loop over all files
      for (i in 1:length(file)) {

        ## extract file name
        file.name <- strsplit(file[i],"\\.")[[1]][1]

        ## extract file header
        file.head <- readLines(file[i], n = length(grep("#", readLines(file[i]))))

        ## read the table data
        raw.spec <- read.table(file[i], header = FALSE, skip = 0, sep = ",")

        ## set the column names accordingly
        colnames(raw.spec) <- c("ppm", "raw.intensity")

        ## create a list with name, element, edge and the data of the spectrum
        raw.spec.end[[i]] <- list("name" = file.name, "data" = list("raw.spec" = raw.spec))

        ## close file loop
        }
    } else if (filetype == "tab") {
      ## create a dummy vector to be filled
      raw.spec.end <- NULL

      ## loop over all files
      for (i in 1:length(file)) {

        ## extract file name
        file.name <- strsplit(file[i],"\\.")[[1]][1]

        ## extract file header
        file.head <- readLines(file[i], n = length(grep("#", readLines(file[i]))))

        ## read the table data
        raw.spec <- read.table(file[i], header = FALSE, skip = 0, sep = "\t")

        ## set the column names accordingly
        colnames(raw.spec) <- c("ppm", "raw.intensity")

        ## create a list with name, element, edge and the data of the spectrum
        raw.spec.end[[i]] <- list("name" = file.name, "data" = list("raw.spec" = raw.spec))

        ## close file loop
      }
    }

  ## create and return a list with E zero and the raw spectrum
  return(raw.spec.end)

  ## close function
}
