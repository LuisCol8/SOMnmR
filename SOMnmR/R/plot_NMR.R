#' NMR Plotting Function
#'
#' This function allows you to plot the 13C-NMR spectra using marking different integration regions.
#' The loaded Raw spectra are intensity normalized and plotted with the chosen integration regions, either spinning side bands (default),
#' the Bonanomi("Bonanomi") regions or the Molecular Mixing Model regions("MMM").
#' The function returns the plots as images either tiff or png, normalized and flattened spectrum
#' @param raw.spec loaded NMR spectra
#' @param NMRmeth Regions to be integrated.
#' Default is spinning side bands, other methods available include: Bonanomi ("Bonanomi") and Molecular mixing model ("MMM").
#' @param use.tiff Logical, default to FALSE (use png)
#' @param set.plot.ymax Set maximum of plot y axis, defaults to NULL
#' @param file.output Logical, default to FALSE
#' @keywords normalization, integration
#' @export
#' @importFrom ggplot2
#' @examples

plot_NMR <- function (raw.spec, NMRmeth = NULL,  use.tiff = NULL,
                      set.plot.ymax = NULL, file.output = NULL) {

  ## check and assign FALSE if optional parameters are not set
  if(is.null(use.tiff)) {use.tiff <- FALSE}
  if(is.null(file.output)) {file.output <- FALSE}

  for (i in 1:length(raw.spec)) {

    sample.name <- raw.spec[[i]]$name
    corr.spec  <- raw.spec[[i]]$data$raw.spec
    corr.spec$raw.intensity <- (corr.spec$raw.intensity/max(corr.spec$raw.intensity))*100


    ## create list of new sample for fitting function
    new.spec <- list("name" = raw.spec[[i]]$name, "data" = list("corr.spec" = corr.spec))

    if (file.output == TRUE) {

      output.spec <- data.frame(cbind(ppm = corr.spec$ppm, norm.intensity = corr.spec$raw.intensity))
      write.csv2(output.spec, paste("output.", sample.name, "csv", sep = "."), row.names = FALSE)

    } else if (file.output == FALSE) {

      output.spec <- data.frame(cbind(ppm = corr.spec$ppm, norm.intensity = corr.spec$raw.intensity))
    }

    ## create image file name with high resolution (png or tiff)
    if (file.output == TRUE) {
      if (use.tiff == TRUE) {
        tiff(filename = paste(sample.name, "normalized.tiff", sep = "."),
             width = 3200, height = 3200, units = "px", res = 800,
             compression = c("none"))
      } else {
        png(filename = paste(sample.name, "normalized.png", sep = "."),
            width = 3200, height = 3200, units = "px", res = 800)
      }

      ## create header name
      head.name <- paste(sample.name, "normalized", sep = " ")

    } else {

      ## simply use sample name as header name
      head.name <- paste(sample.name)
    }

    ## set squared plots (pty), margins (mar), and size of plot (cex)
    if (file.output == TRUE) {
      par(pty="s",
          mar=c(5,3,2,0)+0.1, # c(bottom, left, top, right)
          cex = 0.5
      )
    }

    if (is.null(set.plot.ymax)) {
      plot.ymax <- max(output.spec$norm.intensity)+0.5
    } else {
      plot.ymax <- set.plot.ymax
    }

    if (is.null(NMRmeth)) {

      raincloud_theme = theme(
        text = element_text(size = 16),
        axis.title.x=element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.position = "bottom",
        plot.title = element_text(lineheight=.8, face="bold", size = 16),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))


      ## plot the coordinate system
      plt  <- ggplot(data = output.spec, mapping = aes(x = ppm, y = norm.intensity)) +

        ## area of carboxyl-C
        geom_area(mapping = aes(x = ifelse(ppm >= 295 & ppm <= 355, ppm,0)), fill = "#1F78B4", alpha = 0.8, size = 0.6) +
        geom_area(mapping = aes(x = ifelse(ppm>160 & ppm< 220 , ppm, 0)), fill = "#1F78B4" , alpha = 0.8, size = 0.6) +
        geom_area(mapping = aes(x = ifelse(ppm>-110 & ppm< -50 , ppm, 0)), fill = "#1F78B4", alpha = 0.8, size = 0.6) +

        ## area of Aryl-C
        geom_area(mapping = aes(x = ifelse(ppm>245 & ppm< 295 , ppm, 0)), fill = "#33A02C", alpha = 0.8, size = 0.6) +
        geom_area(mapping = aes(x = ifelse(ppm>110 & ppm< 160 , ppm, 0)), fill = "#33A02C", alpha = 0.8, size = 0.6) +
        geom_area(mapping = aes(x = ifelse(ppm>-145 & ppm< -110 , ppm, 0)), fill = "#33A02C", alpha = 0.8, size = 0.6) +

        ## area of O Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm>-90 & ppm< -25 , ppm, 0)), fill = "#E31A1C", alpha = 0.6, size = 0.6) +
        geom_area(mapping = aes(x = ifelse(ppm>45 & ppm< 110 , ppm, 0)), fill = "#E31A1C", alpha = 0.6, size = 0.6) +
        geom_area(mapping = aes(x = ifelse(ppm>315 & ppm< 355 , ppm, 0)), fill = "#E31A1C", alpha = 0.6, size = 0.6) +

        ## area of Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm>-145 & ppm< -90 , ppm, 0)), fill = "#FF7F00", alpha = 0.6, size = 0.6) +
        geom_area(mapping = aes(x = ifelse(ppm>-10 & ppm< 45 , ppm, 0)), fill = "#FF7F00", alpha = 0.6, size = 0.6) +
        geom_area(mapping = aes(x = ifelse(ppm>260 & ppm< 315 , ppm, 0)), fill = "#FF7F00", alpha = 0.6, size = 0.6) +

        ## plot the NMR spectrum
        geom_line(size = 1)+

        ## Axis title
        xlab("Chemical shift (ppm)") +
        ylab("Intensity")+

        ## create line y = 0
        ylim(-5, plot.ymax)+
        xlim(-200, 400)+
        raincloud_theme

    } else if (NMRmeth == "Bonanomi") {

      raincloud_theme = theme(
        text = element_text(size = 16),
        axis.title.x=element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.position = "bottom",
        plot.title = element_text(lineheight=.8, face="bold", size = 16),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

      ## plot the coordinate system
      plt  <- ggplot(data = output.spec, mapping = aes(x = ppm, y = norm.intensity)) +


        ## area of Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm >= 0 & ppm <= 45, ppm,0)), fill = "#A6CEE3", alpha = 0.6, size = 0.5) +

        ## area of Methoxyl-C
        geom_area(mapping = aes(x = ifelse(ppm>45 & ppm< 60 , ppm, 0)), fill = "#1F78B4", alpha = 0.8, size = 0.5) +

        ## area of O-Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm>-60 & ppm< 90 , ppm, 0)), fill = "#33A02C", alpha = 0.8, size = 0.5) +

        ## area of Di-O-Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm>90 & ppm< 110 , ppm, 0)), fill = "#E31A1C", alpha = 0.8, size = 0.5) +

        ## area of H and C substituted aromatic -C
        geom_area(mapping = aes(x = ifelse(ppm>110 & ppm< 140 , ppm, 0)), fill = "#B2DF8A", alpha = 0.8, size = 0.5) +

        ## area of O substituted aromatic -C
        geom_area(mapping = aes(x = ifelse(ppm>-140 & ppm< 160 , ppm, 0)), fill = "#6A3D9A", alpha = 0.8, size = 0.5) +

        ## area of carboxyl-C
        geom_area(mapping = aes(x = ifelse(ppm>-160 & ppm< 190 , ppm, 0)), fill = "#33A02C" , alpha = 0.8, size = 0.5) +

        ## plot the NMR spectrum
        geom_line(size = 1)+

        ## Axis title
        xlab("Chemical shift (ppm)") +
        ylab("Intensity")+

        ## create line y = 0
        ylim(-5, plot.ymax)+
        xlim(-200, 400)+
        raincloud_theme


    } else if (NMRmeth == "MMM") {

      raincloud_theme = theme(
        text = element_text(size = 16),
        axis.title.x=element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.position = "bottom",
        plot.title = element_text(lineheight=.8, face="bold", size = 16),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

      ## plot the coordinate system
      plt  <- ggplot(data = output.spec, mapping = aes(x = ppm, y = norm.intensity)) +

        ## area of Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm >= 0 & ppm <= 45, ppm,0)), fill = "#A6CEE3", alpha = 0.6, size = 0.5) +

        ## area of Methoxyl-C
        geom_area(mapping = aes(x = ifelse(ppm>45 & ppm< 60 , ppm, 0)), fill = "#1F78B4", alpha = 0.6, size = 0.5) +

        ## area of O-Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm>-60 & ppm< 95 , ppm, 0)), fill = "#B2DF8A", alpha = 0.6, size = 0.5) +

        ## area of Di-O-Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm>95 & ppm< 110 , ppm, 0)), fill = "#33A02C", alpha = 0.6, size = 0.5) +

        ## area of H and C substituted aromatic -C
        geom_area(mapping = aes(x = ifelse(ppm>110 & ppm< 145 , ppm, 0)), fill = "#FB9A99", alpha = 0.6, size = 0.5) +

        ## area of O substituted aromatic -C
        geom_area(mapping = aes(x = ifelse(ppm>-145 & ppm< 165 , ppm, 0)), fill = "#E31A1C", alpha = 0.6, size = 0.5) +

        ## area of carboxyl-C
        geom_area(mapping = aes(x = ifelse(ppm>-165 & ppm< 215 , ppm, 0)), fill = "#FDBF6F", alpha = 0.6, size = 0.5) +

        ## plot the NMR spectrum
        geom_line(size = 1)+

        ## Axis title
        xlab("Chemical shift (ppm)") +
        ylab("Intensity")+

        ## create line y = 0
        ylim(-5, plot.ymax)+
        xlim(-200, 400)+
        raincloud_theme


    }

    ## close image file
    if (file.output == TRUE) {
      dev.off()
    } else {
      print(plt)
    }
  }
  #return(plt)
  ## close function
}
