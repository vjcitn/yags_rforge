
library(geeni)
data(bigd)



source("engine.R")

library(foreach)
library(doMC)
registerDoMC(cores=8)
library(ff)

choppeddd = function( fmla, src )
   parols( fmla, src, family=gaussian(), binit=c(0,0) )

ww = .dfbychunk$new("dfbychunk", data=bigd[1:1e6,], chunksize=125000L)
dd8c = list()
dd8c$M1M = list()
dd8c$M1M[[1]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M1M[[2]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M1M[[3]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M1M[[4]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M1M[[5]] = unix.time( choppeddd( isT~GCcon, ww ) )

ww = .dfbychunk$new("dfbychunk", data=bigd[1:2e6,], chunksize=250000L)
dd8c$M2M = list()
dd8c$M2M[[1]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M2M[[2]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M2M[[3]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M2M[[4]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M2M[[5]] = unix.time( choppeddd( isT~GCcon, ww ) )

ww = .dfbychunk$new("dfbychunk", data=bigd[1:3e6,], chunksize=375000L)
dd8c$M3M = list()
dd8c$M3M[[1]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M3M[[2]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M3M[[3]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M3M[[4]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M3M[[5]] = unix.time( choppeddd( isT~GCcon, ww ) )

ww = .dfbychunk$new("dfbychunk", data=bigd[1:4e6,], chunksize=500000L)
dd8c$M4M = list()
dd8c$M4M[[1]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M4M[[2]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M4M[[3]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M4M[[4]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M4M[[5]] = unix.time( choppeddd( isT~GCcon, ww ) )

ww = .dfbychunk$new("dfbychunk", data=bigd, chunksize=625000L)
dd8c$M5M = list()
dd8c$M5M[[1]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M5M[[2]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M5M[[3]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M5M[[4]] = unix.time( choppeddd( isT~GCcon, ww ) )
dd8c$M5M[[5]] = unix.time( choppeddd( isT~GCcon, ww ) )

save(dd8c, file="dd8c.rda")
