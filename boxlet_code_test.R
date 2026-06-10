# if (!requireNamespace("mgcv", quietly = TRUE)) install.packages("mgcv")
# if (!requireNamespace("akima", quietly = TRUE)) install.packages("akima")
# 
library(mgcv)
library(akima)

mgcv.noclus.names <- paste("mgcv.noclus.sim", 1:16, sep="")
mgcv.clus.names <- paste("mgcv.clus.sim", 1:16, sep="")

rectangle.area2d.func <- function(xlow, xhi, ylow, yhi){
  xdiffs <- xhi - xlow
  ydiffs <- yhi - ylow
  griddiffs <- expand.grid(x=xdiffs, y=ydiffs)
  griddiffs$x * griddiffs$y
}

jeff.covariates.plot <- function(covar = jeff.covariates){
  image(unique(covar$x), unique(covar$y), matrix(covar$habitat, nrow=100, byrow=T))
}

draw.boxes.func <- function(uvals.x, uvals.y, w, h, col=1)
{
  k <- length(uvals.x)
  m <- length(uvals.y)
  for(xline in 1:k){
    for(yrow in 1:m){
      xpoly <- c(rep(uvals.x[xline]-w, 2), rep(uvals.x[xline]+w, 2))
      ypoly.tmp <- c(uvals.y[yrow]-h, uvals.y[yrow]+h)
      ypoly <- c(ypoly.tmp, rev(ypoly.tmp))
      polygon(xpoly, ypoly, border=col)
    }
  }
}

striplet.bounds.func <- function(w, k, s, min.width=1e-8){
  line.int <- (1 - 2*w) / k
  bvec <- seq(from = w, to = w + line.int, by = w/s)
  B <- length(bvec)
  striplet.pos <- numeric(0)
  for(b in 1:B) striplet.pos <- c(striplet.pos, bvec[b] + (0:(k-1)) * line.int)
  for(b in 1:B) striplet.pos <- c(striplet.pos, bvec[b] + (0:(k-1)) * line.int - w)
  for(b in 1:B) striplet.pos <- c(striplet.pos, bvec[b] + (0:(k-1)) * line.int + w)
  striplet.pos <- sort(unique(striplet.pos))
  striplet.diff <- diff(striplet.pos)
  ndiff <- length(striplet.diff)
  remove.elts <- (1:ndiff)[striplet.diff < min.width]
  if(length(remove.elts)>0) striplet.pos <- striplet.pos[-remove.elts]
  Nstriplets <- length(striplet.pos) - 1
  cat("Number of striplets: ", Nstriplets, "\n")
  cat("Summary of striplet widths: \n", summary(diff(striplet.pos)), "\n\n")
  striplet.pos
}

image.grid.func <- function(grid, zname){
  x.pos <- unique(grid$x)
  y.pos <- unique(grid$y)
  if(any(grid[,c("x", "y")] !=expand.grid(x=x.pos, y=y.pos))) stop("Grid supplied in wrong alignment")
  zmat <- matrix(grid[,zname], nrow=length(x.pos))
  image(x.pos, y.pos, zmat, col=heat.colors(20), xlab="x", ylab="y")
  title(main=zname)
}

persp.grid.func <- function(grid, zname, theta=330, phi=40, expand=0.3, ltheta=-120, shade=0.75){
  x.pos <- unique(grid$x)
  y.pos <- unique(grid$y)
  if(any(grid[,c("x", "y")] !=expand.grid(x=x.pos, y=y.pos))) stop("Grid supplied in wrong alignment")
  zmat <- matrix(grid[,zname], nrow=length(x.pos))
  persp(x.pos, y.pos, zmat, theta=theta, phi=phi, expand=expand, ltheta=ltheta, shade=shade, xlab="x", ylab="y")
  title(main=zname)
}

mgcv.boxlet.func <- function(whkmNpars, cluster, Nsurv=1000, v=1, sx=max(1, round(w*100)),
                             sy=max(1, round(h*100)), min.width=1e-8, est.df=T, dfval=missing(), use.musum=F,
                             plotit=F, covariates=jeff.covariates, Nmax=NULL)
{
  if(v!=1) stop("I don't think it's properly encoded for v != 1 : better go read the code and check.")
  
  whkmNpars <- unlist(whkmNpars)
  w <- whkmNpars[1]
  h <- whkmNpars[2]
  k <- whkmNpars[3]
  m <- whkmNpars[4]
  N <- whkmNpars[5]
  
  print(c(w, h, k, m, N))
  
  line.int <- (1 - 2*w) / k
  row.int <- (v - 2*h) / m
  
  print(c(sx=sx, sy=sy))
  
  bvec.x <- seq(from = w, to = w + line.int, by = w/sx)
  bvec.y <- seq(from = h, to = h + row.int, by = h/sy)
  
  Bx <- length(bvec.x)
  By <- length(bvec.y)
  B <- Bx * By
  
  striplet.bounds.x <- striplet.bounds.func(w=w, k=k, s=sx, min.width=min.width)
  n.striplet.x <- length(striplet.bounds.x) - 1
  
  striplet.bounds.y <- striplet.bounds.func(w=h, k=m, s=sy, min.width=min.width)
  n.striplet.y <- length(striplet.bounds.y) - 1
  
  J <- n.striplet.x * n.striplet.y
  
  cat("Number of boxlets tesselating region, J = ", J, "\n")
  cat("Number of grid start-points, B = ", B, "\n")
  
  midvec.xpoints <- (striplet.bounds.x[-1] + striplet.bounds.x[-(n.striplet.x+1)])/2
  midvec.ypoints <- (striplet.bounds.y[-1] + striplet.bounds.y[-(n.striplet.y+1)])/2
  midvec <- expand.grid(x=midvec.xpoints, y=midvec.ypoints)
  
  boxlet.area.vec <- rectangle.area2d.func(xlow=striplet.bounds.x[-(n.striplet.x+1)], xhi=striplet.bounds.x[-1],
                                           ylow=striplet.bounds.y[-(n.striplet.y+1)], yhi=striplet.bounds.y[-1])
  
  gbmat <- matrix(0, nrow=J, ncol=B)
  
  b.ind <- 0
  for(b.indy in (1:By)){
    b.y <- bvec.y[b.indy]
    rows.grid <- seq(b.y, length=m, by=row.int)
    
    for(b.indx in (1:Bx)){
      b.x <- bvec.x[b.indx]
      lines.grid <- seq(b.x, length=k, by = line.int)
      b.ind <- b.ind + 1
      
      min.dist.x <- sapply(midvec$x, function(xval) min(abs(xval - lines.grid)))
      min.dist.y <- sapply(midvec$y, function(yval) min(abs(yval - rows.grid)))
      
      boxlets.included.b <- (1:J)[min.dist.x < w & min.dist.y < h]
      gbmat[boxlets.included.b, b.ind] <- 1
    }
  }
  
  if(is.null(Nmax)) Nmax <- round(qpois(1-1e-10, N) + 0.25*N)
  survest.func <- function(survey_idx){
    print(survey_idx)
    
    jeffpops <- one.2d.ipp.func(Npoints=N, cluster=cluster, covariates=covariates, Nmax=Nmax)
    xvals <- jeffpops$x
    yvals <- jeffpops$y
    if(any(is.na(xvals))) stop("Wrong number of points: quitting this survey")
    
    surv.bval.x <- runif(1, w, w + line.int)
    surv.bval.y <- runif(1, h, h + row.int)
    
    uvals.x <- surv.bval.x + seq(from = 0, by = line.int, length = k)
    uvals.y <- surv.bval.y + seq(from = 0, by =row.int, length = m)
    
    area.surv.vec <- rectangle.area2d.func(xlow=uvals.x-w, xhi=uvals.x+w, ylow=uvals.y-h, yhi=uvals.y+h)
    area.surv <- sum(area.surv.vec)
    
    boxcount.mat <- matrix(-1, ncol=k, nrow=m)
    for(xline in 1:k){
      for(yrow in 1:m){
        boxcount.mat[yrow, xline] <- length((1:N)[(abs(xvals - uvals.x[xline])
                                                   <= w) & (abs(yvals - uvals.y[yrow]) <= h)])
      }
    }
    
    ntot <- sum(boxcount.mat)
    Dhat <- ntot/area.surv
    
    a.offset <- area.surv.vec
    nspotted <- as.vector(t(boxcount.mat))
    
    if(Nsurv==1 | plotit==T){
      par(mfrow=c(3,2))
      jeff.covariates.plot(covar=covariates)
      plot(xvals, yvals, pch=16, xlim=c(0,1), ylim=c(0,1))
      draw.boxes.func(uvals.x=uvals.x, uvals.y=uvals.y, w=w, h=h, col=2)
    }
    
    box.dat <- expand.grid(x=uvals.x, y=uvals.y)
    box.dat$n <-  nspotted
    box.dat$loga.offset <- log(a.offset)
    
    if(est.df==T) box.gam <- gam(n~offset(loga.offset) + s(x, y), family=poisson(link=log), data=box.dat)
    else box.gam <- gam(n~offset(loga.offset) + s(x, y, fx=TRUE, k=dfval+1),
                        family=poisson(link=log), data=box.dat)
    
    if(Nsurv==1 | plotit==T) print(box.gam)
    
    box.terms <- predict.gam(box.gam, type="terms")
    box.lambda <- exp(box.terms + attributes(box.terms)$constant)
    
    sur.x <- midvec.xpoints
    sur.y <- midvec.ypoints
    
    sur.x[sur.x < min(uvals.x)] <- min(uvals.x)
    sur.x[sur.x > max(uvals.x)] <- max(uvals.x)
    
    sur.y[sur.y < min(uvals.y)] <- min(uvals.y)
    sur.y[sur.y > max(uvals.y)] <- max(uvals.y)
    
    boxlet.lambda.mat <- interp(box.dat$x, box.dat$y, box.lambda, xo=sur.x, yo=sur.y)$z
    boxlet.lambda <- as.vector(boxlet.lambda.mat)
    
    muvec <- boxlet.lambda * boxlet.area.vec
    musum <- sum(muvec)
    
    mu.times.g.mat <- muvec * gbmat
    Abvec <- apply(mu.times.g.mat, 2, sum)
    
    if(use.musum) var.n.striplet <- 1/B * sum(Abvec + Abvec^2*(1-1/musum))- (1/B * sum(Abvec))^2
    else{
      Nhat <- v * Dhat
      var.n.striplet <- 1/B * sum(Abvec + Abvec^2*(1-1/Nhat))- (1/B * sum(Abvec))^2
    }
    varD.striplet <- var.n.striplet/area.surv^2
    
    box.dat$lambda <- box.lambda
    boxlet.dat <- midvec
    boxlet.dat$lambda <- boxlet.lambda
    boxlet.dat$muvec <- muvec
    
    if(Nsurv==1 | plotit==T){
      cat("muvec:\n")
      print(summary(muvec))
      image.grid.func(box.dat, zname="lambda")
      persp.grid.func(box.dat, zname="lambda")
      image.grid.func(boxlet.dat, zname="lambda")
      persp.grid.func(boxlet.dat, zname="lambda")
    }
    
    nbar <- mean(boxcount.mat)
    nvec <- as.vector(boxcount.mat)
    varD.P1 <- (1/area.surv^2) * k*m/(k*m-1) * sum((nvec-nbar)^2)
    
    strat.sum <- 0
    for(xline in 1:(k-1)){
      for(yrow in 1:(m-1)){
        xs <- c(xline, xline+1)
        ys <- c(yrow, yrow+1)
        strat.mat <- boxcount.mat[ys, xs]
        strat.mean <- mean(strat.mat)
        strat.vec <- as.vector(strat.mat)
        strat.var <- 1/3 * sum((strat.vec - strat.mean)^2)
        strat.sum <- strat.sum + strat.var
      }
    }
    varD.O1 <- (1/area.surv^2) * k * m * strat.sum / ((k-1)*(m-1))
    
    boxcount.extend <- cbind(NA, boxcount.mat, NA)
    boxcount.extend <- rbind(NA, boxcount.extend, NA)
    
    delta.sum.hvd <- 0
    delta.nsq.sum.hvd <- 0
    delta.sum.hv <- 0
    delta.nsq.sum.hv <- 0
    
    for(xline in 1:k){
      for(yrow in 1:m){
        n.centre <- boxcount.mat[yrow, xline]
        xs <- seq(xline-1, xline+1)+1
        ys <- seq(yrow-1, yrow+1)+1
        mini.mat <- boxcount.extend[ys, xs]
        mini.vec.hvd <- as.vector(mini.mat)
        mini.vec.hvd <- mini.vec.hvd[!is.na(mini.vec.hvd)]
        delta.sum.hvd <- delta.sum.hvd + length(mini.vec.hvd) - 1
        delta.nsq.sum.hvd <- delta.nsq.sum.hvd + sum((mini.vec.hvd - n.centre)^2)
        mini.mat.hv <- mini.mat
        mini.mat.hv[1, 1] <- NA
        mini.mat.hv[3, 1] <- NA
        mini.mat.hv[1, 3] <- NA
        mini.mat.hv[3, 3] <- NA
        mini.vec.hv <- as.vector(mini.mat.hv)
        mini.vec.hv <- mini.vec.hv[!is.na(mini.vec.hv)]
        delta.sum.hv <- delta.sum.hv + length(mini.vec.hv) - 1
        delta.nsq.sum.hv <- delta.nsq.sum.hv + sum((mini.vec.hv - n.centre)^2)
      }
    }
    varD.P4.hvd <- (1/area.surv^2) * k * m * delta.nsq.sum.hvd / (2 * delta.sum.hvd)
    varD.P4.hv <- (1/area.surv^2) * k * m * delta.nsq.sum.hv / (2 * delta.sum.hv)
    
    surv.result <- c(n = ntot, musum=musum, Dhat = Dhat,
                     varD.striplet=varD.striplet, varD.P1 = varD.P1, varD.O1=varD.O1,
                     varD.P4.hvd=varD.P4.hvd, varD.P4.hv=varD.P4.hv,
                     b.x=surv.bval.x, b.y=surv.bval.y)
    if(plotit==T && interactive()) readline()
    surv.result
  }
  all.results <- lapply(1:Nsurv, function(survey) try(survest.func(survey)))
  all.succ <- all.results[sapply(all.results, function(res) !inherits(res, "try-error"))]
  all.succ.mat <- matrix(unlist(all.succ), ncol=length(all.succ[[1]]), byrow=T)
  all.succ.df <- as.data.frame(all.succ.mat)
  names(all.succ.df) <- c("n", "musum", "Dhat", "striplet", "P1", "O1", "P4.hvd", "P4.hv", "b.x", "b.y")
  attributes(all.succ.df)$inputs <- whkmNpars
  all.succ.df
}

oneres.boxplot <- function(res, sd.plot=T)
{
  trueval <- var(res$Dhat)
  comps <- list(res$striplet, res$P1, res$O1, res$P4.hvd, res$P4.hv)
  y.hi <- max(c(unlist(comps), trueval), na.rm=TRUE)
  y.lo <- min(c(unlist(comps), trueval), na.rm=TRUE) - 1
  
  if(sd.plot == T) {
    comps <- lapply(comps, sqrt)
    y.hi <- sqrt(y.hi)
    y.lo <- sqrt(pmax(0, y.lo + 1)) - 0.4
    trueval <- sqrt(trueval)
  }
  
  bp <- boxplot(comps, names = c("striplet", "P1", "O1", "P4.hvd", "P4.hv"), ylim = c(y.lo, y.hi),
                pars=list(medlty="blank", boxfill="light blue", outwex = 0.3, boxwex = 0.6,
                          outpch=NA, outlty="solid", staplewex=0.6, staplelwd=2), cex.axis=1.2)
  abline(h = trueval, col = 2, lwd = 4)
  for(i in 1:length(comps))
    lines(c(i - 0.4, i + 0.4), rep(mean(comps[[i]], na.rm=TRUE), 2), lwd = 2)
  
  print(list(true=trueval, means=lapply(comps, mean, na.rm=TRUE)))
}

boxlet.simgrid.func <- function(plotit=F){
  jeff.grid <- expand.grid(En = c(75, 400), P=c(0.04, 0.25), K=c(10, 20), M=c(5, 16))
  ngrid <- nrow(jeff.grid)
  rachel.grid <- NULL
  for(i in 1:ngrid){
    k <- jeff.grid$K[i]
    m <- jeff.grid$M[i]
    P <- jeff.grid$P[i]
    En <- jeff.grid$En[i]
    if(P==0.04){
      w <- 0.1/k
      h <- 0.1/m
    }
    if(P==0.25){
      w <- 0.5*1.25/(2*k)
      h <- 0.5/(1.25*2*m)
    }
    N <- round(En / P)
    rachel.grid <- rbind(rachel.grid, c(w, h, k, m, N))
  }
  rachel.grid <- as.data.frame(rachel.grid)
  names(rachel.grid) <- c("w", "h", "k", "m", "N")
  rachel.grid <- rachel.grid[order(rachel.grid$N),]
  rownames(rachel.grid) <- 1:ngrid
  rachel.grid
}

jeff.covariates <- expand.grid(x = 1:100, y = 1:100)
jeff.covariates$habitat <- as.numeric(jeff.covariates$x > 50) * 2 + as.numeric(jeff.covariates$y > 50) + 1

one.2d.ipp.func <- function(Npoints, cluster, covariates=jeff.covariates, Nmax=NULL) {
  if (cluster == FALSE) {
    pop.x <- runif(Npoints, 0, 1)
    pop.y <- runif(Npoints, 0, 1)
  } else {
    n_centers <- 4
    cx <- c(0.25, 0.75, 0.25, 0.75)
    cy <- c(0.25, 0.25, 0.75, 0.75)
    alloc <- sample(1:n_centers, Npoints, replace = TRUE)
    pop.x <- cx[alloc] + rnorm(Npoints, 0, 0.1)
    pop.y <- cy[alloc] + rnorm(Npoints, 0, 0.1)
    pop.x <- pmax(0, pmin(1, pop.x))
    pop.y <- pmax(0, pmin(1, pop.y))
  }
  cat("Points generated = ", length(pop.x), "\n")
  list(x=pop.x, y=pop.y)
}

set.seed(123)
boxlet.simgrid <- boxlet.simgrid.func(plotit = FALSE)

single_survey_demo <- mgcv.boxlet.func(
  whkmNpars = boxlet.simgrid[1, ], 
  cluster = TRUE, 
  Nsurv = 1, 
  plotit = TRUE
)

simulation_batch <- mgcv.boxlet.func(
  whkmNpars = boxlet.simgrid[1, ], 
  cluster = TRUE, 
  Nsurv = 20, 
  plotit = FALSE
)

dev.new()
oneres.boxplot(simulation_batch, sd.plot = TRUE)