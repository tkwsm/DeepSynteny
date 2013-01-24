#!/usr/bin/ruby

inprefix = ARGV.shift
xsize = ( ARGV.shift ).to_i
ysize = ( ARGV.shift ).to_i
outprefix = File.new( ARGV.shift, "w+" )

outprefix.print "t2 <- read.table(\"#{inprefix}.plot\")\n"
outprefix.print "t2x <- read.table(\"#{inprefix}.x\")\n"
outprefix.print "t2y <- read.table(\"#{inprefix}.y\")\n"
#outprefix.print "plot( t2[,3], t2[,4])\n"
outprefix.print "plot(t2[,3], t2[,4], xlim=c(0,100000000), ylim=c(0,100000000), pch=20)\n"
outprefix.print "abline(v=0)\n"
outprefix.print "abline(h=0)\n"
outprefix.print "for(i in 1:#{xsize}){abline(v=t2x[i,4],lwd=0.3,col='gray')}\n"
outprefix.print "for(i in 1:#{ysize}){abline(h=t2y[i,4],lwd=0.3,col='gray')}\n"

