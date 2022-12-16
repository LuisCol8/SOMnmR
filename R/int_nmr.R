#' Integration function
#'
#' This function allows you to integrate the 13C-NMR spectra using diferent integration regions.
#' The loaded Raw spectra can be integrated using the spinning side bands regions(default),
#' the Bonanomi("Bonanomi") regions or the Molecular Mixing Model regions("MMM").
#' The function returns the corrected, normalized and flattened spectrum
#' @param raw.spec Raw spectrum
#' @param NMRmeth Regions to be integrated.
#' Default is spinning side bands, other methods available include: Bonanomi ("Bonanomi") and Molecular mixing model ("MMM").
#' @keywords normalization, integration
#' @export
#' @importFrom cmna simp
#' @examples

int_nmr <- function(raw.spec, NMRmeth=NULL, SSBcorr=FALSE) {

  raw.spec.end <- NULL

  if (is.null(NMRmeth)) {

    raw.spec.end <- NULL

    for (i in 1:length(raw.spec)) {
      Integral <- NULL
      name <- raw.spec[[i]]$name
      raw.spec.end[[i]] <- raw.spec[[i]]
      spectrum <- raw.spec[[i]]$data$raw.spec

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.1 <- which(abs(spectrum[[c("ppm")]]-(-145)) == min(abs(spectrum[[c("ppm")]]-(-145))))
      Int.max.1 <- which(abs(spectrum[[c("ppm")]]-(-130)) == min(abs(spectrum[[c("ppm")]]-(-130))))
      Int.x.1 <- c(spectrum[[c("ppm")]][(Int.min.1:Int.max.1)])
      Int.y.1 <- c(spectrum[[c("raw.intensity")]][(Int.min.1:Int.max.1)])
      Int1 <- data.frame(Int.x.1,Int.y.1)
      Int_model <- function (Int.x.1) {(Int.y.1)}
      Integral <- append(Integral,trap(Int_model, -145,-130, m = 240))
      #print(length((Int.x.1)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.2 <- which(abs(spectrum[[c("ppm")]]-(-130)) == min(abs(spectrum[[c("ppm")]]-(-130))))
      Int.max.2 <- which(abs(spectrum[[c("ppm")]]-(-110)) == min(abs(spectrum[[c("ppm")]]-(-110))))
      Int.x.2 <- c(spectrum[[c("ppm")]][(Int.min.2:Int.max.2)])
      Int.y.2 <- c(spectrum[[c("raw.intensity")]][(Int.min.2:Int.max.2)])
      Int2 <- data.frame(Int.x.2,Int.y.2)
      Int_model <- function (Int.x.2) {(Int.y.2)}
      Integral <- append(Integral,trap(f = Int_model, -130,-110, m = 300))
      #print(length((Int.x.2)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.3 <- which(abs(spectrum[[c("ppm")]]-(-110)) == min(abs(spectrum[[c("ppm")]]-(-110))))
      Int.max.3 <- which(abs(spectrum[[c("ppm")]]-(-90)) == min(abs(spectrum[[c("ppm")]]-(-90))))
      Int.x.3 <- c(spectrum[[c("ppm")]][(Int.min.3:Int.max.3)])
      Int.y.3 <- c(spectrum[[c("raw.intensity")]][(Int.min.3:Int.max.3)])
      Int3 <- data.frame(Int.x.3,Int.y.3)
      Int_model <- function (Int.x.3) {(Int.y.3)}
      Integral <- append(Integral,trap(f = Int_model, -110,-90, m = 300))
      #print(length((Int.x.3)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.4 <- which(abs(spectrum[[c("ppm")]]-(-90)) == min(abs(spectrum[[c("ppm")]]-(-90))))
      Int.max.4 <- which(abs(spectrum[[c("ppm")]]-(-86)) == min(abs(spectrum[[c("ppm")]]-(-86))))
      Int.x.4 <- c(spectrum[[c("ppm")]][(Int.min.4:Int.max.4)])
      Int.y.4 <- c(spectrum[[c("raw.intensity")]][(Int.min.4:Int.max.4)])
      Int4 <- data.frame(Int.x.4,Int.y.4)
      Int_model <- function (Int.x.4) {(Int.y.4)}
      Integral <- append(Integral,trap(f = Int_model, -90,-86, m = 66))
      #print(length((Int.x.4)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.5 <- which(abs(spectrum[[c("ppm")]]-(-86)) == min(abs(spectrum[[c("ppm")]]-(-86))))
      Int.max.5 <- which(abs(spectrum[[c("ppm")]]-(-75)) == min(abs(spectrum[[c("ppm")]]-(-75))))
      Int.x.5 <- c(spectrum[[c("ppm")]][(Int.min.5:Int.max.5)])
      Int.y.5 <- c(spectrum[[c("raw.intensity")]][(Int.min.5:Int.max.5)])
      Int5 <- data.frame(Int.x.5,Int.y.5)
      Int_model <- function (Int.x.5) {(Int.y.5)}
      Integral <- append(Integral,trap(f = Int_model, -86,-75, m = 180))
      #print(length((Int.x.5)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.6 <- which(abs(spectrum[[c("ppm")]]-(-75)) == min(abs(spectrum[[c("ppm")]]-(-75))))
      Int.max.6 <- which(abs(spectrum[[c("ppm")]]-(-50)) == min(abs(spectrum[[c("ppm")]]-(-50))))
      Int.x.6 <- c(spectrum[[c("ppm")]][(Int.min.6:Int.max.6)])
      Int.y.6 <- c(spectrum[[c("raw.intensity")]][(Int.min.6:Int.max.6)])
      Int6 <- data.frame(Int.x.6,Int.y.6)
      Int_model <- function (Int.x.6) {(Int.y.6)}
      Integral <- append(Integral,trap(f = Int_model, -75,-50, m = 400))
      #print(length((Int.x.6)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.7 <- which(abs(spectrum[[c("ppm")]]-(-50)) == min(abs(spectrum[[c("ppm")]]-(-50))))
      Int.max.7 <- which(abs(spectrum[[c("ppm")]]-(-45)) == min(abs(spectrum[[c("ppm")]]-(-45))))
      Int.x.7 <- c(spectrum[[c("ppm")]][(Int.min.7:Int.max.7)])
      Int.y.7 <- c(spectrum[[c("raw.intensity")]][(Int.min.7:Int.max.7)])
      Int7 <- data.frame(Int.x.7,Int.y.7)
      Int_model <- function (Int.x.7) {(Int.y.7)}
      Integral <- append(Integral,trap(f = Int_model, -50,-45, m = 80))
      #print(length((Int.x.7)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.8 <- which(abs(spectrum[[c("ppm")]]-(-45)) == min(abs(spectrum[[c("ppm")]]-(-45))))
      Int.max.8 <- which(abs(spectrum[[c("ppm")]]-(-25)) == min(abs(spectrum[[c("ppm")]]-(-25))))
      Int.x.8 <- c(spectrum[[c("ppm")]][(Int.min.8:Int.max.8)])
      Int.y.8 <- c(spectrum[[c("raw.intensity")]][(Int.min.8:Int.max.8)])
      Int8 <- data.frame(Int.x.8,Int.y.8)
      Int_model <- function (Int.x.8) {(Int.y.8)}
      Integral <- append(Integral,trap(f = Int_model, -45,-25, m = 300))
      #print(length((Int.x.8)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.9 <- which(abs(spectrum[[c("ppm")]]-(-25)) == min(abs(spectrum[[c("ppm")]]-(-25))))
      Int.max.9 <- which(abs(spectrum[[c("ppm")]]-(-10)) == min(abs(spectrum[[c("ppm")]]-(-10))))
      Int.x.9 <- c(spectrum[[c("ppm")]][(Int.min.9:Int.max.9)])
      Int.y.9 <- c(spectrum[[c("raw.intensity")]][(Int.min.9:Int.max.9)])
      Int9 <- data.frame(Int.x.9,Int.y.9)
      Int_model <- function (Int.x.9) {(Int.y.9)}
      Integral <- append(Integral,trap(f = Int_model, -25,-10, m = 240))
      #print(length((Int.x.9)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.10 <- which(abs(spectrum[[c("ppm")]]-(-10)) == min(abs(spectrum[[c("ppm")]]-(-10))))
      Int.max.10 <- which(abs(spectrum[[c("ppm")]]-(5)) == min(abs(spectrum[[c("ppm")]]-(5))))
      Int.x.10 <- c(spectrum[[c("ppm")]][(Int.min.10:Int.max.10)])
      Int.y.10 <- c(spectrum[[c("raw.intensity")]][(Int.min.10:Int.max.10)])
      Int10 <- data.frame(Int.x.10,Int.y.10)
      Int_model <- function (Int.x.10) {(Int.y.10)}
      Integral <- append(Integral,trap(f = Int_model, -10,5, m = 240))
      #print(length((Int.x.10)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.11 <- which(abs(spectrum[[c("ppm")]]-(5)) == min(abs(spectrum[[c("ppm")]]-(5))))
      Int.max.11 <- which(abs(spectrum[[c("ppm")]]-(25)) == min(abs(spectrum[[c("ppm")]]-(25))))
      Int.x.11 <- c(spectrum[[c("ppm")]][(Int.min.11:Int.max.11)])
      Int.y.11 <- c(spectrum[[c("raw.intensity")]][(Int.min.11:Int.max.11)])
      Int11 <- data.frame(Int.x.11,Int.y.11)
      Int_model <- function (Int.x.11) {(Int.y.11)}
      Integral <- append(Integral,trap(f = Int_model, 5,25, m = 300))
      #print(length((Int.x.11)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.12 <- which(abs(spectrum[[c("ppm")]]-(25)) == min(abs(spectrum[[c("ppm")]]-(25))))
      Int.max.12 <- which(abs(spectrum[[c("ppm")]]-(45)) == min(abs(spectrum[[c("ppm")]]-(45))))
      Int.x.12 <- c(spectrum[[c("ppm")]][(Int.min.12:Int.max.12)])
      Int.y.12 <- c(spectrum[[c("raw.intensity")]][(Int.min.12:Int.max.12)])
      Int12 <- data.frame(Int.x.12,Int.y.12)
      Int_model <- function (Int.x.12) {(Int.y.12)}
      Integral <- append(Integral,trap(f = Int_model, 25,45, m = 300))
      #print(length((Int.x.12)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.13 <- which(abs(spectrum[[c("ppm")]]-(45)) == min(abs(spectrum[[c("ppm")]]-(45))))
      Int.max.13 <- which(abs(spectrum[[c("ppm")]]-(49)) == min(abs(spectrum[[c("ppm")]]-(49))))
      Int.x.13 <- c(spectrum[[c("ppm")]][(Int.min.13:Int.max.13)])
      Int.y.13 <- c(spectrum[[c("raw.intensity")]][(Int.min.13:Int.max.13)])
      Int13 <- data.frame(Int.x.13,Int.y.13)
      Int_model <- function (Int.x.13) {(Int.y.13)}
      Integral <- append(Integral,trap(f = Int_model, 45,49, m = 66))
      #print(length((Int.x.13)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.14 <- which(abs(spectrum[[c("ppm")]]-(49)) == min(abs(spectrum[[c("ppm")]]-(49))))
      Int.max.14 <- which(abs(spectrum[[c("ppm")]]-(60)) == min(abs(spectrum[[c("ppm")]]-(60))))
      Int.x.14 <- c(spectrum[[c("ppm")]][(Int.min.14:Int.max.14)])
      Int.y.14 <- c(spectrum[[c("raw.intensity")]][(Int.min.14:Int.max.14)])
      Int14 <- data.frame(Int.x.14,Int.y.14)
      Int_model <- function (Int.x.14) {(Int.y.14)}
      Integral <- append(Integral,trap(f = Int_model, 49,60, m = 180))
      #print(length((Int.x.14)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.15 <- which(abs(spectrum[[c("ppm")]]-(60)) == min(abs(spectrum[[c("ppm")]]-(60))))
      Int.max.15 <- which(abs(spectrum[[c("ppm")]]-(85)) == min(abs(spectrum[[c("ppm")]]-(85))))
      Int.x.15 <- c(spectrum[[c("ppm")]][(Int.min.15:Int.max.15)])
      Int.y.15 <- c(spectrum[[c("raw.intensity")]][(Int.min.15:Int.max.15)])
      Int15 <- data.frame(Int.x.15,Int.y.15)
      Int_model <- function (Int.x.15) {(Int.y.15)}
      Integral <- append(Integral,trap(f = Int_model, 60,85, m = 400))
      #print(length((Int.x.15)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.16 <- which(abs(spectrum[[c("ppm")]]-(85)) == min(abs(spectrum[[c("ppm")]]-(85))))
      Int.max.16 <- which(abs(spectrum[[c("ppm")]]-(90)) == min(abs(spectrum[[c("ppm")]]-(90))))
      Int.x.16 <- c(spectrum[[c("ppm")]][(Int.min.16:Int.max.16)])
      Int.y.16 <- c(spectrum[[c("raw.intensity")]][(Int.min.16:Int.max.16)])
      Int16 <- data.frame(Int.x.16,Int.y.16)
      Int_model <- function (Int.x.16) {(Int.y.16)}
      Integral <- append(Integral,trap(f = Int_model, 85,90, m = 80))
      #print(length((Int.x.16)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.17 <- which(abs(spectrum[[c("ppm")]]-(90)) == min(abs(spectrum[[c("ppm")]]-(90))))
      Int.max.17 <- which(abs(spectrum[[c("ppm")]]-(110)) == min(abs(spectrum[[c("ppm")]]-(110))))
      Int.x.17 <- c(spectrum[[c("ppm")]][(Int.min.17:Int.max.17)])
      Int.y.17 <- c(spectrum[[c("raw.intensity")]][(Int.min.17:Int.max.17)])
      Int17 <- data.frame(Int.x.17,Int.y.17)
      Int_model <- function (Int.x.17) {(Int.y.17)}
      Integral <- append(Integral,trap(f = Int_model, 90,110, m = 300))
      #print(length((Int.x.17)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.18 <- which(abs(spectrum[[c("ppm")]]-(110)) == min(abs(spectrum[[c("ppm")]]-(110))))
      Int.max.18 <- which(abs(spectrum[[c("ppm")]]-(125)) == min(abs(spectrum[[c("ppm")]]-(125))))
      Int.x.18 <- c(spectrum[[c("ppm")]][(Int.min.18:Int.max.18)])
      Int.y.18 <- c(spectrum[[c("raw.intensity")]][(Int.min.18:Int.max.18)])
      Int18 <- data.frame(Int.x.18,Int.y.18)
      Int_model <- function (Int.x.18) {(Int.y.18)}
      Integral <- append(Integral,trap(f = Int_model, 110,125, m = 240))
      #print(length((Int.x.18)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.19 <- which(abs(spectrum[[c("ppm")]]-(125)) == min(abs(spectrum[[c("ppm")]]-(125))))
      Int.max.19 <- which(abs(spectrum[[c("ppm")]]-(140)) == min(abs(spectrum[[c("ppm")]]-(140))))
      Int.x.19 <- c(spectrum[[c("ppm")]][(Int.min.19:Int.max.19)])
      Int.y.19 <- c(spectrum[[c("raw.intensity")]][(Int.min.19:Int.max.19)])
      Int19 <- data.frame(Int.x.19,Int.y.19)
      Int_model <- function (Int.x.19) {(Int.y.19)}
      Integral <- append(Integral,trap(f = Int_model, 125,140, m = 240))
      #print(length((Int.x.19)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.20 <- which(abs(spectrum[[c("ppm")]]-(140)) == min(abs(spectrum[[c("ppm")]]-(140))))
      Int.max.20 <- which(abs(spectrum[[c("ppm")]]-(160)) == min(abs(spectrum[[c("ppm")]]-(160))))
      Int.x.20 <- c(spectrum[[c("ppm")]][(Int.min.20:Int.max.20)])
      Int.y.20 <- c(spectrum[[c("raw.intensity")]][(Int.min.20:Int.max.20)])
      Int20 <- data.frame(Int.x.20,Int.y.20)
      Int_model <- function (Int.x.20) {(Int.y.20)}
      Integral <- append(Integral,trap(f = Int_model, 140,160, m = 300))
      #print(length((Int.x.20)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.21 <- which(abs(spectrum[[c("ppm")]]-(160)) == min(abs(spectrum[[c("ppm")]]-(160))))
      Int.max.21 <- which(abs(spectrum[[c("ppm")]]-(180)) == min(abs(spectrum[[c("ppm")]]-(180))))
      Int.x.21 <- c(spectrum[[c("ppm")]][(Int.min.21:Int.max.21)])
      Int.y.21 <- c(spectrum[[c("raw.intensity")]][(Int.min.21:Int.max.21)])
      Int21 <- data.frame(Int.x.21,Int.y.21)
      Int_model <- function (Int.x.21) {(Int.y.21)}
      Integral <- append(Integral,trap(f = Int_model, 160, 180, m = 300))
      #print(length((Int.x.21)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.22 <- which(abs(spectrum[[c("ppm")]]-(180)) == min(abs(spectrum[[c("ppm")]]-(180))))
      Int.max.22 <- which(abs(spectrum[[c("ppm")]]-(185)) == min(abs(spectrum[[c("ppm")]]-(185))))
      Int.x.22 <- c(spectrum[[c("ppm")]][(Int.min.22:Int.max.22)])
      Int.y.22<- c(spectrum[[c("raw.intensity")]][(Int.min.22:Int.max.22)])
      Int22 <- data.frame(Int.x.22,Int.y.22)
      Int_model <- function (Int.x.22) {(Int.y.22)}
      Integral <- append(Integral,trap(f = Int_model, 180,185, m = 80))
      #print(length((Int.x.22)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.23 <- which(abs(spectrum[[c("ppm")]]-(185)) == min(abs(spectrum[[c("ppm")]]-(185))))
      Int.max.23 <- which(abs(spectrum[[c("ppm")]]-(195)) == min(abs(spectrum[[c("ppm")]]-(195))))
      Int.x.23 <- c(spectrum[[c("ppm")]][(Int.min.23:Int.max.23)])
      Int.y.23 <- c(spectrum[[c("raw.intensity")]][(Int.min.23:Int.max.23)])
      Int23 <- data.frame(Int.x.23,Int.y.23)
      Int_model <- function (Int.x.23) {(Int.y.23)}
      Integral <- append(Integral,trap(f = Int_model, 185,195, m = 160))
      #print(length((Int.x.23)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.24 <- which(abs(spectrum[[c("ppm")]]-(195)) == min(abs(spectrum[[c("ppm")]]-(195))))
      Int.max.24 <- which(abs(spectrum[[c("ppm")]]-(220)) == min(abs(spectrum[[c("ppm")]]-(220))))
      Int.x.24 <- c(spectrum[[c("ppm")]][(Int.min.24:Int.max.24)])
      Int.y.24 <- c(spectrum[[c("raw.intensity")]][(Int.min.24:Int.max.24)])
      Int24 <- data.frame(Int.x.24,Int.y.24)
      Int_model <- function (Int.x.24) {(Int.y.24)}
      Integral <- append(Integral,trap(f = Int_model, 195,220, m = 400))
      #print(length((Int.x.24)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.25 <- which(abs(spectrum[[c("ppm")]]-(220)) == min(abs(spectrum[[c("ppm")]]-(220))))
      Int.max.25 <- which(abs(spectrum[[c("ppm")]]-(225)) == min(abs(spectrum[[c("ppm")]]-(225))))
      Int.x.25 <- c(spectrum[[c("ppm")]][(Int.min.25:Int.max.25)])
      Int.y.25 <- c(spectrum[[c("raw.intensity")]][(Int.min.25:Int.max.25)])
      Int25 <- data.frame(Int.x.25,Int.y.25)
      Int_model <- function (Int.x.25) {(Int.y.25)}
      Integral <- append(Integral,trap(f = Int_model, 220,225, m = 80))
      #print(length((Int.x.25)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.26 <- which(abs(spectrum[[c("ppm")]]-(225)) == min(abs(spectrum[[c("ppm")]]-(225))))
      Int.max.26 <- which(abs(spectrum[[c("ppm")]]-(245)) == min(abs(spectrum[[c("ppm")]]-(245))))
      Int.x.26 <- c(spectrum[[c("ppm")]][(Int.min.26:Int.max.26)])
      Int.y.26 <- c(spectrum[[c("raw.intensity")]][(Int.min.26:Int.max.26)])
      Int26 <- data.frame(Int.x.26,Int.y.26)
      Int_model <- function (Int.x.26) {(Int.y.26)}
      Integral <- append(Integral,trap(f = Int_model, 225,245, m = 300))
      #print(length((Int.x.26)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.27 <- which(abs(spectrum[[c("ppm")]]-(245)) == min(abs(spectrum[[c("ppm")]]-(245))))
      Int.max.27 <- which(abs(spectrum[[c("ppm")]]-(260)) == min(abs(spectrum[[c("ppm")]]-(260))))
      Int.x.27 <- c(spectrum[[c("ppm")]][(Int.min.27:Int.max.27)])
      Int.y.27 <- c(spectrum[[c("raw.intensity")]][(Int.min.27:Int.max.27)])
      Int27 <- data.frame(Int.x.27,Int.y.27)
      Int_model <- function (Int.x.27) {(Int.y.27)}
      Integral <- append(Integral,trap(f = Int_model, 245,260, m = 240))
      #print(length((Int.x.27)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.28 <- which(abs(spectrum[[c("ppm")]]-(260)) == min(abs(spectrum[[c("ppm")]]-(260))))
      Int.max.28 <- which(abs(spectrum[[c("ppm")]]-(275)) == min(abs(spectrum[[c("ppm")]]-(275))))
      Int.x.28 <- c(spectrum[[c("ppm")]][(Int.min.28:Int.max.28)])
      Int.y.28 <- c(spectrum[[c("raw.intensity")]][(Int.min.28:Int.max.28)])
      Int28 <- data.frame(Int.x.28,Int.y.28)
      Int_model <- function (Int.x.28) {(Int.y.28)}
      Integral <- append(Integral,trap(f = Int_model, 260,275, m = 240))
      #print(length((Int.x.28)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.29 <- which(abs(spectrum[[c("ppm")]]-(275)) == min(abs(spectrum[[c("ppm")]]-(275))))
      Int.max.29 <- which(abs(spectrum[[c("ppm")]]-(290)) == min(abs(spectrum[[c("ppm")]]-(290))))
      Int.x.29 <- c(spectrum[[c("ppm")]][(Int.min.29:Int.max.29)])
      Int.y.29 <- c(spectrum[[c("raw.intensity")]][(Int.min.29:Int.max.29)])
      Int29 <- data.frame(Int.x.29,Int.y.29)
      Int_model <- function (Int.x.29) {(Int.y.29)}
      Integral <- append(Integral,trap(f = Int_model, 275,295, m = 240))
      #print(length((Int.x.29)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.30 <- which(abs(spectrum[[c("ppm")]]-(295)) == min(abs(spectrum[[c("ppm")]]-(295))))
      Int.max.30 <- which(abs(spectrum[[c("ppm")]]-(315)) == min(abs(spectrum[[c("ppm")]]-(315))))
      Int.x.30 <- c(spectrum[[c("ppm")]][(Int.min.30:Int.max.30)])
      Int.y.30 <- c(spectrum[[c("raw.intensity")]][(Int.min.30:Int.max.30)])
      Int30 <- data.frame(Int.x.30,Int.y.30)
      Int_model <- function (Int.x.30) {(Int.y.30)}
      Integral <- append(Integral,trap(f = Int_model, 295,315, m = 300))
      #print(length((Int.x.30)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.31 <- which(abs(spectrum[[c("ppm")]]-(315)) == min(abs(spectrum[[c("ppm")]]-(315))))
      Int.max.31 <- which(abs(spectrum[[c("ppm")]]-(320)) == min(abs(spectrum[[c("ppm")]]-(320))))
      Int.x.31 <- c(spectrum[[c("ppm")]][(Int.min.31:Int.max.31)])
      Int.y.31 <- c(spectrum[[c("raw.intensity")]][(Int.min.31:Int.max.31)])
      Int31 <- data.frame(Int.x.31,Int.y.31)
      Int_model <- function (Int.x.31) {(Int.y.31)}
      Integral <- append(Integral,trap(f = Int_model, 315,320, m = 80))
      #print(length((Int.x.31)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.32 <- which(abs(spectrum[[c("ppm")]]-(320)) == min(abs(spectrum[[c("ppm")]]-(320))))
      Int.max.32 <- which(abs(spectrum[[c("ppm")]]-(330)) == min(abs(spectrum[[c("ppm")]]-(330))))
      Int.x.32 <- c(spectrum[[c("ppm")]][(Int.min.32:Int.max.32)])
      Int.y.32 <- c(spectrum[[c("raw.intensity")]][(Int.min.32:Int.max.32)])
      Int32 <- data.frame(Int.x.32,Int.y.32)
      Int_model <- function (Int.x.32) {(Int.y.32)}
      Integral <- append(Integral,trap(f = Int_model, 320,330, m = 160))
      #print(length((Int.x.32)))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.33 <- which(abs(spectrum[[c("ppm")]]-(330)) == min(abs(spectrum[[c("ppm")]]-(330))))
      Int.max.33 <- which(abs(spectrum[[c("ppm")]]-(355)) == min(abs(spectrum[[c("ppm")]]-(355))))
      Int.x.33 <- c(spectrum[[c("ppm")]][(Int.min.33:Int.max.33)])
      Int.y.33 <- c(spectrum[[c("raw.intensity")]][(Int.min.33:Int.max.33)])
      Int33 <- data.frame(Int.x.33,Int.y.33)
      Int_model <- function (Int.x.33) {(Int.y.33)}
      Integral <- append(Integral,trap(f = Int_model, 330,355, m = 400))
      #print(length((Int.x.33)))

      norm <- sum(Integral[1:33])
      normalized.Int <- (Integral/norm)*100
      Integral <- data.frame(normalized.Int)
      raw.spec.end[[i]] <- list("name" = name, "data" = list("raw.spec" = spectrum,"Integral" = Integral))
    }
  }  else if (NMRmeth == "MMM-SSB") {

    raw.spec.end <- NULL

    for (i in 1:length(raw.spec)) {
      Integral <- NULL
      name <- raw.spec[[i]]$name
      raw.spec.end[[i]] <- raw.spec[[i]]
      spectrum <- raw.spec[[i]]$data$raw.spec

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.1 <- which(abs(spectrum[[c("ppm")]]-(-70)) == min(abs(spectrum[[c("ppm")]]-(-70))))
      Int.max.1 <- which(abs(spectrum[[c("ppm")]]-(-50)) == min(abs(spectrum[[c("ppm")]]-(-50))))
      Int.x.1 <- c(spectrum[[c("ppm")]][(Int.min.1:Int.max.1)])
      Int.y.1 <- c(spectrum[[c("raw.intensity")]][(Int.min.1:Int.max.1)])
      Integral <- append(Integral,trapz(Int.x.1,Int.y.1))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.2 <- which(abs(spectrum[[c("ppm")]]-(-50)) == min(abs(spectrum[[c("ppm")]]-(-50))))
      Int.max.2 <- which(abs(spectrum[[c("ppm")]]-(-15)) == min(abs(spectrum[[c("ppm")]]-(-15))))
      Int.x.2 <- c(spectrum[[c("ppm")]][(Int.min.2:Int.max.2)])
      Int.y.2 <- c(spectrum[[c("raw.intensity")]][(Int.min.2:Int.max.2)])
      Integral <- append(Integral,trapz(Int.x.2,Int.y.2))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.3 <- which(abs(spectrum[[c("ppm")]]-(-15)) == min(abs(spectrum[[c("ppm")]]-(-15))))
      Int.max.3 <- which(abs(spectrum[[c("ppm")]]-(0)) == min(abs(spectrum[[c("ppm")]]-(0))))
      Int.x.3 <- c(spectrum[[c("ppm")]][(Int.min.3:Int.max.3)])
      Int.y.3 <- c(spectrum[[c("raw.intensity")]][(Int.min.3:Int.max.3)])
      Integral <- append(Integral,trapz(Int.x.3,Int.y.3))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.4 <- which(abs(spectrum[[c("ppm")]]-(0)) == min(abs(spectrum[[c("ppm")]]-(0))))
      Int.max.4 <- which(abs(spectrum[[c("ppm")]]-(45)) == min(abs(spectrum[[c("ppm")]]-(45))))
      Int.x.4 <- c(spectrum[[c("ppm")]][(Int.min.4:Int.max.4)])
      Int.y.4 <- c(spectrum[[c("raw.intensity")]][(Int.min.4:Int.max.4)])
      Integral <- append(Integral,trapz(Int.x.4,Int.y.4))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.5 <- which(abs(spectrum[[c("ppm")]]-(45)) == min(abs(spectrum[[c("ppm")]]-(45))))
      Int.max.5 <- which(abs(spectrum[[c("ppm")]]-(60)) == min(abs(spectrum[[c("ppm")]]-(60))))
      Int.x.5 <- c(spectrum[[c("ppm")]][(Int.min.5:Int.max.5)])
      Int.y.5 <- c(spectrum[[c("raw.intensity")]][(Int.min.5:Int.max.5)])
      Integral <- append(Integral,trapz(Int.x.5,Int.y.5))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.6 <- which(abs(spectrum[[c("ppm")]]-(60)) == min(abs(spectrum[[c("ppm")]]-(60))))
      Int.max.6 <- which(abs(spectrum[[c("ppm")]]-(95)) == min(abs(spectrum[[c("ppm")]]-(95))))
      Int.x.6 <- c(spectrum[[c("ppm")]][(Int.min.6:Int.max.6)])
      Int.y.6 <- c(spectrum[[c("raw.intensity")]][(Int.min.6:Int.max.6)])
      Integral <- append(Integral,trapz(Int.x.6,Int.y.6))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.7 <- which(abs(spectrum[[c("ppm")]]-(95)) == min(abs(spectrum[[c("ppm")]]-(95))))
      Int.max.7 <- which(abs(spectrum[[c("ppm")]]-(110)) == min(abs(spectrum[[c("ppm")]]-(110))))
      Int.x.7 <- c(spectrum[[c("ppm")]][(Int.min.7:Int.max.7)])
      Int.y.7 <- c(spectrum[[c("raw.intensity")]][(Int.min.7:Int.max.7)])
      Integral <- append(Integral,trapz(Int.x.7,Int.y.7))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.8 <- which(abs(spectrum[[c("ppm")]]-(110)) == min(abs(spectrum[[c("ppm")]]-(110))))
      Int.max.8 <- which(abs(spectrum[[c("ppm")]]-(145)) == min(abs(spectrum[[c("ppm")]]-(145))))
      Int.x.8 <- c(spectrum[[c("ppm")]][(Int.min.8:Int.max.8)])
      Int.y.8 <- c(spectrum[[c("raw.intensity")]][(Int.min.8:Int.max.8)])
      Integral <- append(Integral,trapz(Int.x.8,Int.y.8))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.9 <- which(abs(spectrum[[c("ppm")]]-(145)) == min(abs(spectrum[[c("ppm")]]-(145))))
      Int.max.9 <- which(abs(spectrum[[c("ppm")]]-(165)) == min(abs(spectrum[[c("ppm")]]-(165))))
      Int.x.9 <- c(spectrum[[c("ppm")]][(Int.min.9:Int.max.9)])
      Int.y.9 <- c(spectrum[[c("raw.intensity")]][(Int.min.9:Int.max.9)])
      Integral <- append(Integral,trapz(Int.x.9,Int.y.9))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.10 <- which(abs(spectrum[[c("ppm")]]-(165)) == min(abs(spectrum[[c("ppm")]]-(165))))
      Int.max.10 <- which(abs(spectrum[[c("ppm")]]-(190)) == min(abs(spectrum[[c("ppm")]]-(190))))
      Int.x.10 <- c(spectrum[[c("ppm")]][(Int.min.10:Int.max.10)])
      Int.y.10 <- c(spectrum[[c("raw.intensity")]][(Int.min.10:Int.max.10)])
      Integral <- append(Integral,trapz(Int.x.10,Int.y.10))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.11 <- which(abs(spectrum[[c("ppm")]]-(190)) == min(abs(spectrum[[c("ppm")]]-(190))))
      Int.max.11 <- which(abs(spectrum[[c("ppm")]]-(215)) == min(abs(spectrum[[c("ppm")]]-(215))))
      Int.x.11 <- c(spectrum[[c("ppm")]][(Int.min.11:Int.max.11)])
      Int.y.11 <- c(spectrum[[c("raw.intensity")]][(Int.min.11:Int.max.11)])
      Integral <- append(Integral,trapz(Int.x.11,Int.y.11))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.12 <- which(abs(spectrum[[c("ppm")]]-(250)) == min(abs(spectrum[[c("ppm")]]-(250))))
      Int.max.12 <- which(abs(spectrum[[c("ppm")]]-(270)) == min(abs(spectrum[[c("ppm")]]-(270))))
      Int.x.12 <- c(spectrum[[c("ppm")]][(Int.min.12:Int.max.12)])
      Int.y.12 <- c(spectrum[[c("raw.intensity")]][(Int.min.12:Int.max.12)])
      Integral <- append(Integral,trapz(Int.x.12,Int.y.12))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.13 <- which(abs(spectrum[[c("ppm")]]-(270)) == min(abs(spectrum[[c("ppm")]]-(270))))
      Int.max.13 <- which(abs(spectrum[[c("ppm")]]-(305)) == min(abs(spectrum[[c("ppm")]]-(305))))
      Int.x.13 <- c(spectrum[[c("ppm")]][(Int.min.13:Int.max.13)])
      Int.y.13 <- c(spectrum[[c("raw.intensity")]][(Int.min.13:Int.max.13)])
      Integral <- append(Integral,trapz(Int.x.13,Int.y.13))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.14 <- which(abs(spectrum[[c("ppm")]]-(305)) == min(abs(spectrum[[c("ppm")]]-(305))))
      Int.max.14 <- which(abs(spectrum[[c("ppm")]]-(325)) == min(abs(spectrum[[c("ppm")]]-(325))))
      Int.x.14 <- c(spectrum[[c("ppm")]][(Int.min.14:Int.max.14)])
      Int.y.14 <- c(spectrum[[c("raw.intensity")]][(Int.min.14:Int.max.14)])
      Integral <- append(Integral,trapz(Int.x.14,Int.y.14))

      ## Extract ppm (x) and the intesity (y) for predifined intervals
      Int.min.15 <- which(abs(spectrum[[c("ppm")]]-(325)) == min(abs(spectrum[[c("ppm")]]-(325))))
      Int.max.15 <- which(abs(spectrum[[c("ppm")]]-(350)) == min(abs(spectrum[[c("ppm")]]-(350))))
      Int.x.15 <- c(spectrum[[c("ppm")]][(Int.min.15:Int.max.15)])
      Int.y.15 <- c(spectrum[[c("raw.intensity")]][(Int.min.15:Int.max.15)])
      Integral <- append(Integral,trapz(Int.x.15,Int.y.15))

      norm <- sum(Integral[1:14])
      normalized.Int <- (Integral/norm)*100
      Integral <- data.frame(normalized.Int)

      if (SSBcorr == TRUE) {

        Alkyl <- ifelse(Integral[4,1] - Integral[15,1] < 0,0, Integral[4,1] - Integral[15,1])

        N_Alkyl_Methoxyl <- ifelse(Integral[5,1] < 0,0, Integral[5,1])

        O_Alkyl <- ifelse(Integral[6,1] < 0,0, Integral[6,1])

        Di_O_Alkyl <- ifelse(Integral[7,1] < 0,0, Integral[7,1] + Integral[1,1] + Integral[12,1])

        Aromatic <- ifelse(Integral[8,1] < 0,0, Integral[8,1] + Integral[2,1] + Integral[13,1])

        Phenolic <- ifelse(Integral[9,1] < 0,0, Integral[9,1] + Integral[3,1] + Integral[14,1])

        Amide_Carboxylic <- ifelse(Integral[10,1] < 0,0, Integral[10,1] + Integral[3,1] + 2*Integral[15,1])

        Ketone <- ifelse(Integral[11,1] < 0,0, Integral[11,1])

        Amide_to_Ketone <- c(Amide_Carboxylic + Ketone)

        Integral <- data.frame(Alkyl, N_Alkyl_Methoxyl, O_Alkyl, Di_O_Alkyl, Aromatic, Phenolic, Amide_to_Ketone)

        Integral <- Integral/sum(Integral)

        raw.spec.end[[i]] <- list("name" = name, "data" = list("raw.spec" = spectrum, "Integral" = Integral))

      } else if (SSBcorr == FALSE) {
        raw.spec.end[[i]] <- list("name" = name, "data" = list("raw.spec" = spectrum, "Integral" = Integral))
      }
    }

  }
  return(raw.spec.end)
}
}
