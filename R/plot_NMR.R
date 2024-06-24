#' NMR Plotting Function
#'
#' This function allows you to plot the 13C-NMR spectra using marking different integration regions.
#' The loaded Raw spectra are intensity normalized and plotted with the chosen integration regions, either spinning side bands (default),
#' the Bonanomi("Bonanomi") regions or the Molecular Mixing Model regions("MMM").
#' The function returns the plots as images either tiff or png, normalized and flattened spectrum
#' @param raw.spec loaded NMR spectra
#' @param NMRmeth  Regions to be integrated, methods available include: "4region", "Bonanomi", "Smernik" and Molecular mixing model ("MMM").
#' @param use.tiff Logical, default to FALSE (use png)
#' @param set.plot.ymax Set maximum of plot y axis, defaults to NULL
#' @param file.output Logical, default to FALSE
#' @returns A plot of the NMR spectrum and a  csv file of the data plotted.
#' @keywords normalization integration
#' @export
#' @import ggplot2
#' @examples
#' \donttest{
#' library(ggplot2)
#' data("GarciaF200")
#' plot_NMR(GarciaF200, NMRmeth = "MMM", file.output = FALSE, use.tiff = FALSE)
#' }

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

    ppm <- NULL
    norm.intensity <- NULL

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
             width = 3200, height = 3200, units = "px", res = 400,
             compression = c("none"))
      } else {
        png(filename = paste(sample.name, "normalized.png", sep = "."),
            width = 3200, height = 3200, units = "px", res = 400)
      }

      ## create header name
      head.name <- paste(sample.name, "normalized", sep = " ")

    } else {

      ## simply use sample name as header name
      head.name <- paste(sample.name)
    }

    ## preserving old parameters
    oldpar <- par(no.readonly = TRUE)

    # Ensure the original graphical parameters are restored on exit
    on.exit(par(oldpar))


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
        legend.position = c(0.8,0.85),
        plot.title = element_text(lineheight=.8, face="bold", size = 16),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

      colors  <- c("Alkyl-C" = "#CC79A7", "O/N-Alkyl-C" = "#0072B2", "Aryl-C" = "#E69F00", "Carboxyl-C" = "#009E73")


      ## plot the coordinate system
      plt <- ggplot(data = output.spec, mapping = aes(x = ppm, y = norm.intensity)) +

        scale_fill_manual(breaks = c("Alkyl-C", "O/N-Alkyl-C", "Aryl-C", "Carboxyl-C"),values=colors, labels = c("Alkyl-C", "O/N-Alkyl-C", "Aryl-C", "Carboxyl-C"))+

        ## area of Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm>0 & ppm< 45 , ppm, NA), fill = "Alkyl-C"), alpha = 0.8, size = 0.6) +

        ## area of O Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm>45 & ppm< 110 , ppm, NA), fill = "O/N-Alkyl-C"), alpha = 0.8, size = 0.6) +

        ## area of Aryl-C
        geom_area(mapping = aes(x = ifelse(ppm>110 & ppm< 160 , ppm, NA), fill = "Aryl-C"), alpha = 0.8, size = 0.6) +

        ## area of carboxyl-C
        geom_area(mapping = aes(x = ifelse(ppm>160 & ppm< 220 , ppm, NA), fill = "Carboxyl-C"), alpha = 0.8, size = 0.6) +

        ## plot the NMR spectrum
        geom_line(size = 0.25)+

        ## Axis title
        xlab("Chemical shift (ppm)") +
        ylab("Intensity")+

        ## create line y = 0
        ylim(-5, plot.ymax)+
        scale_x_reverse(limits=c(400, -200))+
        raincloud_theme

      print(plt)

    }  else if (NMRmeth == "Smernik") {

      raincloud_theme = theme(
        text = element_text(size = 16),
        axis.title.x=element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.position = c(0.8,0.85),
        plot.title = element_text(lineheight=.8, face="bold", size = 16),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

      colors  <- c("Alkyl-C" = "#CC79A7", "O/N-Alkyl-C" = "#0072B2", "Aryl-C" = "#E69F00", "Carboxyl-C" = "#009E73")


      ## plot the coordinate system
      plt <- ggplot(data = output.spec, mapping = aes(x = ppm, y = norm.intensity)) +

        scale_fill_manual(breaks = c("Alkyl-C", "O/N-Alkyl-C", "Aryl-C", "Carboxyl-C"),values=colors, labels = c("Alkyl-C", "O/N-Alkyl-C", "Aryl-C", "Carboxyl-C"))+

        ## area of Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm>0 & ppm< 45 , ppm, NA), fill = "Alkyl-C"), alpha = 0.8, size = 0.6) +

        ## area of O Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm>45 & ppm< 110 , ppm, NA), fill = "O/N-Alkyl-C"), alpha = 0.8, size = 0.6) +

        ## area of Aryl-C
        geom_area(mapping = aes(x = ifelse(ppm>110 & ppm< 165 , ppm, NA), fill = "Aryl-C"), alpha = 0.8, size = 0.6) +

        ## area of carboxyl-C
        geom_area(mapping = aes(x = ifelse(ppm>165 & ppm< 185 , ppm, NA), fill = "Carboxyl-C"), alpha = 0.8, size = 0.6) +

        ## plot the NMR spectrum
        geom_line(size = 0.25)+

        ## Axis title
        xlab("Chemical shift (ppm)") +
        ylab("Intensity")+

        ## create line y = 0
        ylim(-5, plot.ymax)+
        scale_x_reverse(limits=c(400, -200))+
        raincloud_theme

      print(plt)

    } else if (NMRmeth == "Bonanomi") {

      raincloud_theme = theme(
        text = element_text(size = 16),
        axis.title.x=element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.position = c(0.8,0.85),
        plot.title = element_text(lineheight=.8, face="bold", size = 16),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

      colors  <- c("Alkyl-C" = "#CC79A7", "Methoxyl/N-Alkyl-C"= "#D55E00", "O-Alkyl-C" = "#0072B2",  "Di-O-Alkyl-C" = "#F0E442",
                   "Aromatic-C" = "#56B4E9", "Phenolic-C" = "#E69F00", "Carbonyl and amide-C" = "#999999", "Ketone-C" = "#009E73")

      ## plot the coordinate system
      plt  <- ggplot(data = output.spec, mapping = aes(x = ppm, y = norm.intensity)) +

        scale_fill_manual(breaks = c("Alkyl-C", "Methoxyl/N-Alkyl-C", "O-Alkyl-C", "Di-O-Alkyl-C", "Aromatic-C",
                                     "Phenolic-C", "Carbonyl and amide-C", "Ketone-C"),
                          values=colors,
                          labels = c("Alkyl-C", "Methoxyl/N-Alkyl-C", "O-Alkyl-C", "Di-O-Alkyl-C", "Aromatic-C",
                                     "Phenolic-C", "Carbonyl and amide-C","Ketone-C"))+

        ## area of Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm >= 0 & ppm <= 45, ppm,NA), y = norm.intensity, fill = "Alkyl-C"), alpha = 0.8, size = 0.5) +

        ## area of Methoxyl-C
        geom_area(mapping = aes(x = ifelse(ppm> 46 & ppm< 60 , ppm, NA), y = norm.intensity, fill = "Methoxyl/N-Alkyl-C"), alpha = 0.8, size = 0.5) +

        ## area of O-Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm> 61 & ppm< 90 , ppm, NA), y = norm.intensity, fill = "O-Alkyl-C"), alpha = 0.8, size = 0.5) +

        ## area of Di-O-Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm> 91 & ppm< 110 , ppm, NA), y = norm.intensity, fill = "Di-O-Alkyl-C"), alpha = 0.8, size = 0.5) +

        ## area of H and C substituted aromatic -C
        geom_area(mapping = aes(x = ifelse(ppm> 111 & ppm< 140 , ppm, NA), y = norm.intensity, fill = "Aromatic-C"), alpha = 0.8, size = 0.5) +

        ## area of O substituted aromatic -C
        geom_area(mapping = aes(x = ifelse(ppm> 141 & ppm< 160 , ppm, NA), y = norm.intensity, fill = "Phenolic-C"), alpha = 0.8, size = 0.5) +

        ## area of carboxyl-C
        geom_area(mapping = aes(x = ifelse(ppm> 161 & ppm< 190 , ppm, NA), y = norm.intensity, fill = "Carbonyl and amide-C"), alpha = 0.8, size = 0.5) +

        ## plot the NMR spectrum
        geom_line(size = 0.25)+

        ## Axis title
        xlab("Chemical shift (ppm)") +
        ylab("Intensity")+

        ## create line y = 0
        ylim(-5, plot.ymax)+
        scale_x_reverse(limits=c(400, -200))+
        raincloud_theme

      print(plt)


    } else if (NMRmeth == "MMM") {

      raincloud_theme = theme(
        text = element_text(size = 16),
        axis.title.x=element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 0, vjust = 0.5),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        legend.position = c(0.8,0.85),
        plot.title = element_text(lineheight=.8, face="bold", size = 16),
        panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
        axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))

      colors  <- c("Alkyl-C" = "#CC79A7", "Methoxyl/N-Alkyl-C"= "#D55E00", "O-Alkyl-C" = "#0072B2",  "Di-O-Alkyl-C" = "#F0E442",
                   "Aromatic-C" = "#56B4E9", "Phenolic-C" = "#E69F00", "Carbonyl and amide-C" = "#999999", "Ketone-C" = "#009E73")

      ## plot the coordinate system
      plt  <- ggplot(data = output.spec, mapping = aes(x = ppm, y = norm.intensity)) +

        scale_fill_manual(breaks = c("Alkyl-C", "Methoxyl/N-Alkyl-C", "O-Alkyl-C", "Di-O-Alkyl-C", "Aromatic-C",
                                     "Phenolic-C", "Carbonyl and amide-C", "Ketone-C"),
                          values=colors,
                          labels = c("Alkyl-C", "Methoxyl/N-Alkyl-C", "O-Alkyl-C", "Di-O-Alkyl-C", "Aromatic-C",
                                     "Phenolic-C", "Carbonyl and amide-C","Ketone-C"))+
        ## area of Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm >= 0 & ppm <= 45, ppm,NA), y = norm.intensity, fill = "Alkyl-C"), alpha = 0.8, size = 0.5) +

        ## area of Methoxyl-C
        geom_area(mapping = aes(x = ifelse(ppm> 45 & ppm< 60 , ppm, NA), y = norm.intensity, fill = "Methoxyl/N-Alkyl-C"), alpha = 0.8, size = 0.5) +

        ## area of O-Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm> 60 & ppm< 95 , ppm, NA), y = norm.intensity, fill = "O-Alkyl-C"), alpha = 0.8, size = 0.5) +

        ## area of Di-O-Alkyl-C
        geom_area(mapping = aes(x = ifelse(ppm> 95 & ppm< 110 , ppm, NA), y = norm.intensity, fill = "Di-O-Alkyl-C"), alpha = 0.8, size = 0.5) +

        ## area of H and C substituted aromatic -C
        geom_area(mapping = aes(x = ifelse(ppm> 110 & ppm< 145 , ppm, NA), y = norm.intensity, fill = "Aromatic-C"), alpha = 0.8, size = 0.5) +

        ## area of O substituted aromatic -C
        geom_area(mapping = aes(x = ifelse(ppm> 145 & ppm< 165 , ppm, NA), y = norm.intensity, fill = "Phenolic-C"), alpha = 0.8, size = 0.5) +

        ## area of carboxyl-C
        geom_area(mapping = aes(x = ifelse(ppm> 165 & ppm< 190 , ppm, NA), y = norm.intensity, fill = "Carbonyl and amide-C"), alpha = 0.8, size = 0.5) +

        ## area of carboxyl-C
        geom_area(mapping = aes(x = ifelse(ppm> 190 & ppm< 210 , ppm, NA), y = norm.intensity, fill = "Ketone-C"), alpha = 0.8, size = 0.5) +

        ## plot the NMR spectrum
        geom_line(size = 0.25)+

        ## Axis title
        xlab("Chemical shift (ppm)") +
        ylab("Intensity")+

        ## create line y = 0
        ylim(-5, plot.ymax)+
        scale_x_reverse(limits=c(400, -200))+
        raincloud_theme

      print(plt)


    }

    ## close image file
    if (file.output == TRUE) {
      dev.off()
    } else {
      print(plt)
    }
  }
  return(plt)
  ## close function
}
