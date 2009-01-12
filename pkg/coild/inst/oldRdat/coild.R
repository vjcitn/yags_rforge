demod <- matrix(c(1,2,3,NA,NA,
                  2,2,NA,NA,NA,
                  1,2,3,4,6,
                  2,NA,NA,NA,NA,
                  3,3,3,NA,NA,
                  2,2,NA,NA,NA,
                  1,2,3,4,3,
                  2,NA,NA,NA,NA,
                  3,3,3,5,NA,
                  2,2,NA,NA,NA,
                  1,2,3,4,2,
                  2,NA,NA,NA,NA,
                  3,3,3,6,NA,
                  2,2,NA,NA,NA,
                  1,2,3,4,7,
                  2,NA,NA,NA,NA,
                  3,3,3,5,NA), nc=5, byrow=T)

demod.id <- 1:nrow(demod)
demod.times <- t(matrix(rep(1:ncol(demod),nrow(demod)),nc=nrow(demod)))

