#!/usr/bin/ruby

inprefix = ARGV.shift
xsize = ( ARGV.shift ).to_i
ysize = ( ARGV.shift ).to_i
outprefix = File.new( ARGV.shift, "w+" )

outprefix.print "t1 <- read.table(\"#{inprefix}.plot\")\n"
outprefix.print "tx <- read.table(\"#{inprefix}.x\")\n"
outprefix.print "ty <- read.table(\"#{inprefix}.y\")\n"
outprefix.print "data <- read.table(\"#{inprefix}.clst\", header=TRUE, row.names=1)\n"
outprefix.print "data.r <- read.table(\"#{inprefix}.r.clst\", header=TRUE, row.names=1)\n"
#outprefix.print "plot( t1[,3], t1[,4])\n"
outprefix.print "plot(t1[,3], t1[,4], xlim=c(0,30000000), ylim=c(0,30000000))\n"
outprefix.print "abline(v=0)\n"
outprefix.print "abline(h=0)\n"
outprefix.print "for(i in 1:#{xsize}){ abline(v=tx[i, 4]) }\n"
outprefix.print "for(i in 1:#{ysize}){ abline(h=ty[i, 4]) }\n"

outprefix.print "d <- dist( data, method=\"euclidean\") \n"
#outprefix.print "d <- dist( data, method=\"maximum\") \n"
#outprefix.print "d <- dist( data, method=\"manhattan\") \n"

#outprefix.print "d.hc <- hclust( d, method=\"ward\" ) \n"
outprefix.print "d.hc <- hclust( d, method=\"complete\" ) \n"

outprefix.print "plot( d.hc, main=\"Clustering1\", hang=-1 ) \n"
outprefix.print "cutcol <- cutree( d.hc, k=20 ) \n"
outprefix.print "write.table( cutcol, file=\"#{inprefix}.cutcol\", col.names=TRUE) \n"

outprefix.print "dr <- dist( data.r, method=\"euclidean\") \n"
#outprefix.print "dr <- dist( data.r, method=\"maximum\") \n"
#outprefix.print "dr <- dist( data.r, method=\"manhattan\") \n"

#outprefix.print "dr.hc <- hclust( dr, method=\"ward\" ) \n"
outprefix.print "dr.hc <- hclust( dr, method=\"complete\" ) \n"

outprefix.print "plot( dr.hc, main=\"Clustering1\", hang=-1 ) \n"
outprefix.print "cutcol.r <- cutree( dr.hc, k=20 ) \n"
outprefix.print "write.table( cutcol.r, file=\"#{inprefix}.r.cutcol\", col.names=TRUE) \n"

