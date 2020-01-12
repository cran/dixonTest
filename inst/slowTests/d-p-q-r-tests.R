####	d|ensity
####	p|robability (cumulative)
####	q|uantile
####	r|andom number generation
####
####	Functions for  ``d/p/q/r''
require(dixonTest)

.ptime <- proc.time()
F <- FALSE
T <- TRUE

options(warn = 2)
##      ======== No warnings, unless explicitly asserted via
assertWarning <- tools::assertWarning

as.nan <- function(x) { x[is.na(x) & !is.nan(x)] <- NaN ; x }
###-- these are identical in ./arith-true.R ["fixme": use source(..)]
opt.conformance <- 0
Meps <- .Machine $ double.eps
xMax <- .Machine $ double.xmax
options(rErr.eps = 1e-30)
rErr <- function(approx, true, eps = getOption("rErr.eps", 1e-30))
{
    ifelse(Mod(true) >= eps,
	   1 - approx / true, # relative error
	   true - approx)     # absolute error (e.g. when true=0)
}
## Numerical equality: Here want "rel.error" almost always:
All.eq <- function(x,y) {
    all.equal.numeric(x,y, tolerance = 64*.Machine$double.eps,
                      scale = max(0, mean(abs(x), na.rm=TRUE)))
}
if(!interactive())
    set.seed(123)

## The prefixes of ALL the PDQ & R functions
PDQRinteg <- c("dixon")
PDQR <- c(PDQRinteg)
##PQonly <- c("tukey")

###--- Discrete Distributions --- Consistency Checks  pZZ = cumsum(dZZ)

##for(pre in PDQRinteg) { n <- paste("d",pre,sep=""); cat(n,": "); str(get(n))}

### must safely stop
i <- c( 2, 0, 2, 1, 4)
j <- c( 1, 1, 3, 0, 3)
n <- c(31, 7, 7, 7, 7)

options(try.outFile = stdout())
for (k in 1:5) {
    try(qdixon(
            p = 0.1,
            n = n[k],
            i = i[k],
            j = j[k]
        ),
        silent = FALSE)
}


###-------- Continuous Distributions ----------


###===== Random numbers -- first, just output:

set.seed(123)
# .Random.seed <- c(0L, 17292L, 29447L, 24113L)
n <- 20
## for(pre in PDQR) { n <- paste("r",pre,sep=""); cat(n,": "); str(get(n))}
(Rdixon	  <- rdixon   (n, i = 3, j = 2) )

(Pdixon	  <- pdixon    (Rdixon, n = n, i = 3, j = 2) )

ddixon	 (Rdixon, n = n, i = 3, j = 2)

## Check q*(p*(.)) = identity
All.eq(Rdixon,	  qdixon	   (Pdixon, n = n, i = 3, j = 2))

## Same with "upper tail":
All.eq(Rdixon,	  qdixon	  (1- Pdixon, n = n, i = 3, j = 2, lower=F))

## Check q*(p* ( log ), log) = identity
All.eq(Rdixon,	  qdixon	   (log(Pdixon), n = n, i = 3, j = 2,  log=T))

## same q*(p* (log) log) with upper tail:
All.eq(Rdixon,	  qdixon	   (log1p(-Pdixon), n = n, i = 3, j = 2, lower=F, log=T))

## Check log( upper.tail ):
All.eq(log1p(-Pdixon),	  pdixon	   (Rdixon, n = n, i = 3, j = 2, lower=F, log=T))

## Check data with test.data from McBane (2006)
xd <- read.table(file = "test1.ref.output.txt", header = F, skip = 2)
r <- xd$V1
All.eq(xd$V3, (p <- pdixon(r, n = 7, i = 2, j = 1, lower.tail = F)))
All.eq(xd$V2, (x <- ddixon(r, n = 7, i = 2, j = 1)))

xd <- read.table(file = "test2.ref.output.txt", header = T, skip = 1)
r <- sapply(xd$n, function(n) qdixon(0.1, n, i = 2, j = 1, lower.tail=T))
All.eq(xd$Rcrit , r)

xd <- read.table(file = "test4.ref.output.txt", header = F, skip = 152)
alpha <- c(0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40,
           0.50, 0.60, 0.70, 0.80, 0.90, 0.95)
for (n in 6:30) {
   r <- qdixon(p = alpha, n = n, i = 3, j = 2)
   o <- unlist(xd[n-5,-1])
   cat( All.eq(o, r), "\n" )
}

cat("Time elapsed: ", proc.time() - .ptime,"\n")
