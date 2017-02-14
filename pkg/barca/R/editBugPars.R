#' convenience function for setting hyperparameters and point density
#' @param blines character vector of bugs text
#' @param beta_a numeric beta hyperparameter value
#' @param beta_b numeric beta hyperparameter value
#' @param ini_x initial value of sequence of concentrations for pointwise posterior generation
#' @param del_x gap value for sequence of concentrations for pointwise posterior generation
#' @param npts number of concentration values to survey
#' @export
editBugPars = function(blines, beta_a = 25, beta_b = 2500,
   ini_x = .1, del_x = .1, npts=20 ) {
   blines = gsub("%%BETA_A%%", beta_a, blines)
   blines = gsub("%%BETA_B%%", beta_b, blines)
   blines = gsub("%%INI_X%%", ini_x, blines)
   blines = gsub("%%DEL_X%%", del_x, blines)
   blines = gsub("%%N_PTS%%", npts, blines)
   blines
}

