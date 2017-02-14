#' convert CSV in a specific format to a jagsdata file
#' @import rjags
#' @importFrom utils read.csv
#' @param type character string in \code{c("TypeIII", "TypeIa", "TypeV")}
#' @param package name of package in which the file will be found under \code{csv/} in installed folder
#' @param dodump logical indicating that the base::dump function will be used to create a .jagsdata file, with prefix given by the value of \code{type}
#' @export
csv2jags = function(type, package="barca",
   dodump = TRUE) {
 if (!(type %in% c("TypeIII", "TypeIa", "TypeV"))) stop("need Type*")
 targ = system.file(paste0("/csv/", type, ".csv"), package=package)
 sspec = read.csv(targ, stringsAsFactors=FALSE)
#         Subject Serotype Case.Maternal.Ab Control.Maternal.Ab
 sspec = sspec[-nrow(sspec),] # drop median record
 N = nrow(sspec)
 A = sspec$Control.Maternal.Ab # general Ab record
 A[seq(1,N,4)] = sspec[seq(1,N,4), 3]
 str = rep(1:(N/4), each=4)
 Nstrat = N/4
 Nsubg = 2
 d = rep(c(1,0,0,0), Nstrat)
 if (dodump) dump(c("N", "Nstrat", "Nsubg", "str", "d", "A"), 
   paste0(type, ".jagsdata"))
 else list(N=N, Nstrat=Nstrat, Nsubg=Nsubg, str=str, d=d, A=A)
}
