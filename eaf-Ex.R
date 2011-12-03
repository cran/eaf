pkgname <- "eaf"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('eaf')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("HybridGA")
### * HybridGA

flush(stderr()); flush(stdout())

### Name: HybridGA
### Title: Results of Hybrid GA on vanzyl and Richmond water networks
### Aliases: HybridGA
### Keywords: datasets

### ** Examples

data(HybridGA)
print(HybridGA$vanzyl)
print(HybridGA$richmond)



cleanEx()
nameEx("SPEA2minstoptimeRichmond")
### * SPEA2minstoptimeRichmond

flush(stderr()); flush(stdout())

### Name: SPEA2minstoptimeRichmond
### Title: Results of SPEA2 when minimising electrical cost and maximising
###   the minimum idle time of pumps on Richmond water network.
### Aliases: SPEA2minstoptimeRichmond
### Keywords: datasets

### ** Examples

data(HybridGA)
data(SPEA2minstoptimeRichmond)
SPEA2minstoptimeRichmond[,2] <- SPEA2minstoptimeRichmond[,2] / 60
eafplot (SPEA2minstoptimeRichmond, xlab = expression(C[E]),
         ylab = "Minimum idle time (minutes)",
         las = 1, log = "y", ymaximise = TRUE)



cleanEx()
nameEx("SPEA2relativeRichmond")
### * SPEA2relativeRichmond

flush(stderr()); flush(stdout())

### Name: SPEA2relativeRichmond
### Title: Results of SPEA2 with relative time-controlled triggers on
###   Richmond water network.
### Aliases: SPEA2relativeRichmond
### Keywords: datasets

### ** Examples

data(HybridGA)
data(SPEA2relativeRichmond)
eafplot (SPEA2relativeRichmond, percentiles = c(25, 50, 75),
         xlab = expression(C[E]), ylab = "Total switches",
         xlim = c(90, 140), ylim = c(0, 25),
         extra.points = HybridGA$richmond, extra.lty = "dashed",
         extra.legend = "Hybrid GA")



cleanEx()
nameEx("SPEA2relativeVanzyl")
### * SPEA2relativeVanzyl

flush(stderr()); flush(stdout())

### Name: SPEA2relativeVanzyl
### Title: Results of SPEA2 with relative time-controlled triggers on
###   Vanzyl's water network.
### Aliases: SPEA2relativeVanzyl
### Keywords: datasets

### ** Examples

data(HybridGA)
data(SPEA2relativeVanzyl)
eafplot(SPEA2relativeVanzyl, percentiles = c(25, 50, 75), 
        xlab = expression(C[E]), ylab = "Total switches", xlim = c(320, 400),
        extra.points = HybridGA$vanzyl, extra.legend = "Hybrid GA")



cleanEx()
nameEx("eaf-package")
### * eaf-package

flush(stderr()); flush(stdout())

### Name: eaf-package
### Title: Plots of the Empirical Attainment Function
### Aliases: eaf-package eaf
### Keywords: package graphs multivariate optimize Time-quality algorithm
###   profile Empirical attainment function

### ** Examples

data(gcp2x2)
tabucol<-subset(gcp2x2, alg!="TSinN1")
tabucol$alg<-tabucol$alg[drop=TRUE]
eafplot(time+best~run,data=tabucol,subset=tabucol$inst=="DSJC500.5")

eafplot(time+best~run|inst,groups=alg,data=gcp2x2)
eafplot(time+best~run|inst,groups=alg,data=gcp2x2,
	percentiles=c(0,50,100),include.extremes=TRUE,
	cex=1.4, lty=c(2,1,2),lwd=c(2,2,2),
        col=c("black","blue","grey50"))

A1<-read.data.sets(file.path(system.file(package="eaf"),"extdata","ALG_1_dat"))
A2<-read.data.sets(file.path(system.file(package="eaf"),"extdata","ALG_2_dat"))
eafplot(A1,A2, percentiles=c(50))
eafplot(list(A1=A1, A2=A2), percentiles=c(50))
eafdiffplot(A1, A2)
## Not run: dev.copy2pdf(file="eaf.pdf", onefile=TRUE, width=5, height=4)



cleanEx()
nameEx("eafdiffplot")
### * eafdiffplot

flush(stderr()); flush(stdout())

### Name: eafdiffplot
### Title: Empirical attainment function differences
### Aliases: eafdiffplot
### Keywords: graphs

### ** Examples

A1<-read.data.sets(file.path(system.file(package="eaf"),"extdata","ALG_1_dat"))
A2<-read.data.sets(file.path(system.file(package="eaf"),"extdata","ALG_2_dat"))
eafdiffplot(A1,A2)
eafdiffplot(A1,A2, full.eaf=TRUE)



cleanEx()
nameEx("eafplot")
### * eafplot

flush(stderr()); flush(stdout())

### Name: eafplot
### Title: Plot the Empirical Attainment Function for two objectives
### Aliases: eafplot eafplot.default eafplot.data.frame eafplot.formula
###   eafplot.list
### Keywords: graphs

### ** Examples

data(gcp2x2)
tabucol<-subset(gcp2x2, alg!="TSinN1")
tabucol$alg<-tabucol$alg[drop=TRUE]
eafplot(time+best~run,data=tabucol,subset=tabucol$inst=="DSJC500.5")

eafplot(time+best~run|inst,groups=alg,data=gcp2x2)
eafplot(time+best~run|inst,groups=alg,data=gcp2x2,
	percentiles=c(0,50,100),include.extremes=TRUE,
	cex=1.4, lty=c(2,1,2),lwd=c(2,2,2),
        col=c("black","blue","grey50"))

A1 <- read.data.sets(file.path(system.file(package = "eaf"), "extdata", "ALG_1_dat"))
A2 <- read.data.sets(file.path(system.file(package = "eaf"), "extdata", "ALG_2_dat"))
eafplot(A1, A2, percentiles = c(50))
eafplot(list(A1 = A1, A2 = A2), percentiles = c(50))
## Not run: dev.copy2pdf(file = "eaf.pdf", onefile = TRUE, width = 5, height = 4)

## Using extra.points
data(HybridGA)
data(SPEA2relativeVanzyl)
eafplot(SPEA2relativeVanzyl, percentiles = c(25, 50, 75), 
        xlab = expression(C[E]), ylab = "Total switches", xlim = c(320, 400),
        extra.points = HybridGA$vanzyl, extra.legend = "Hybrid GA")

data(SPEA2relativeRichmond)
eafplot (SPEA2relativeRichmond, percentiles = c(25, 50, 75),
         xlab = expression(C[E]), ylab = "Total switches",
         xlim = c(90, 140), ylim = c(0, 25),
         extra.points = HybridGA$richmond, extra.lty = "dashed",
         extra.legend = "Hybrid GA")

data(SPEA2minstoptimeRichmond)
SPEA2minstoptimeRichmond[,2] <- SPEA2minstoptimeRichmond[,2] / 60
eafplot (SPEA2minstoptimeRichmond, xlab = expression(C[E]),
         ylab = "Minimum idle time (minutes)",
         las = 1, log = "y", ymaximise = TRUE, main = "SPEA2 (Richmond)")



cleanEx()
nameEx("gcp2x2")
### * gcp2x2

flush(stderr()); flush(stdout())

### Name: gcp2x2
### Title: Metaheuristics for solving the Graph Vertex Coloring Problem
### Aliases: gcp2x2
### Keywords: datasets

### ** Examples

data(gcp2x2)
## maybe str(gcp2x2)



cleanEx()
nameEx("read.data.sets")
### * read.data.sets

flush(stderr()); flush(stdout())

### Name: read.data.sets
### Title: Read several data.frame sets
### Aliases: read.data.sets
### Keywords: file

### ** Examples

A1<-read.data.sets(file.path(system.file(package="eaf"),"extdata","ALG_1_dat"))
str(A1)
A2<-read.data.sets(file.path(system.file(package="eaf"),"extdata","ALG_2_dat"))
str(A2)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
