

library(geeni)
data(bigd)



source("engine.R")


library(biglm)


choppedlm = function(fmla, data, chunks) {
  ai = function(x) as.integer(as.hi(x))
  m = biglm(fmla, data[ai(chunks[[1]]),])
  if (length(chunks)==1) return(m)
  for (j in 2:length(chunks)) m = update(m, bigd[ai(chunks[[j]]),])
  m
}

library(ff)

c1m = chunk(from=1, to=5e6, by=1e6, max=5e6)
c100k = chunk(from=1, to=5e6, by=1e5, max=5e6)

qr1m = list()
qr1m$M1M = list()
qr1m$M1M[[1]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1] ) )
qr1m$M1M[[2]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1] ) )
qr1m$M1M[[3]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1] ) )
qr1m$M1M[[4]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1] ) )
qr1m$M1M[[5]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1] ) )

qr1m$M2M = list()
qr1m$M2M[[1]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1:2] ) )
qr1m$M2M[[2]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1:2] ) )
qr1m$M2M[[3]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1:2] ) )
qr1m$M2M[[4]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1:2] ) )
qr1m$M2M[[5]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1:2] ) )

qr1m$M3M = list()
qr1m$M3M[[1]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1:3] ) )
qr1m$M3M[[2]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1:3] ) )
qr1m$M3M[[3]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1:3] ) )
qr1m$M3M[[4]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1:3] ) )
qr1m$M3M[[5]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1:3] ) )

qr1m$M4M = list()
qr1m$M4M[[1]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1:4] ) )
qr1m$M4M[[2]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1:4] ) )
qr1m$M4M[[3]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1:4] ) )
qr1m$M4M[[4]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1:4] ) )
qr1m$M4M[[5]] = unix.time( choppedlm( isT~GCcon, bigd, c1m[1:4] ) )

qr1m$M5M = list()
qr1m$M5M[[1]] = unix.time( choppedlm( isT~GCcon, bigd, c1m ) )
qr1m$M5M[[2]] = unix.time( choppedlm( isT~GCcon, bigd, c1m ) )
qr1m$M5M[[3]] = unix.time( choppedlm( isT~GCcon, bigd, c1m ) )
qr1m$M5M[[4]] = unix.time( choppedlm( isT~GCcon, bigd, c1m ) )
qr1m$M5M[[5]] = unix.time( choppedlm( isT~GCcon, bigd, c1m ) )
