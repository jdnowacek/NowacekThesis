# For boxlets / quadrat samples, we don't estimate detection.  These are straightforward counts
# per quadrat.


# NOTE: I don't think that the boxlet method is properly coded for the case where the vertical height
# of the region is not v=1.  It seems to be partially encoded, but (e.g.) the boxlet boundaries don't seem
# to be done.  Have to go through the code carefully to fix this scenario.

# library(DSpat)
library(mgcv)
library(akima)  # for interp

mgcv.noclus.names <- paste("mgcv.noclus.sim", 1:16, sep="")
mgcv.clus.names <- paste("mgcv.clus.sim", 1:16, sep="")


######################################################################################


rectangle.area2d.func <- function(xlow, xhi, ylow, yhi){
	# rectangle.area2d.func 5/5/09
	# For rectangular regions, calculates the area of a set of boxes or boxlets where the lower x-boundaries
	# are in xlow, the upper x-boundaries in xhi, and similarly for ylow and yhi.
	# The x- and y- values in the arguments are not cycled: they are only given once.
	# For example, if xlow and xhi had length 6, and ylow and yhi had length 3, the total number of
	# boxlets would be 6*3 = 18, with arrangement:
	# 13 14 15  16 17 18
	#   7   8   9  10 11 12
	#   1   2   3   4   5   6
	#
	# If calling the function to find areas of full boxes, rather than boxlets, use:
	# 	rectangle.area2d.func(xlow=uvals.x-w, xhi=uvals.x+w, ylow=uvals.y-h, yhi=uvals.y+h)
	#
	#----------------------------------------------------------------------------------

	xdiffs <- xhi - xlow
	ydiffs <- yhi - ylow
	griddiffs <- expand.grid(x=xdiffs, y=ydiffs)
	griddiffs$x * griddiffs$y
}



###################################################################


jeff.covariates.plot <- function(covar = jeff.covariates){
	# jeff.covariates.plot 21/4/09
	image(unique(covar$x), unique(covar$y), matrix(covar$habitat, nrow=100, byrow=T))
}


##############################################################


draw.boxes.func <- function(uvals.x, uvals.y, w, h, col=1)
{
        ## draw.boxes.func 1/5/09
        ## Draws the grid of selected quadrats on a 1x1 square.

        k <- length(uvals.x)
        m <- length(uvals.y)
        for(xline in 1:k){
                for(yrow in 1:m){
                        ## lines(rep(uvals.x[xline]-w, 2), c(uvals.y[yrow]-h, uvals.y[yrow]+h), col=col)
                        ## lines(rep(uvals.x[xline]+w, 2), c(uvals.y[yrow]-h, uvals.y[yrow]+h), col=col)
                        ## lines(c(uvals.x[xline]-w, uvals.x[xline]+w), rep(uvals.y[yrow]-h, 2), col=col)
                        ## lines(c(uvals.x[xline]-w, uvals.x[xline]+w), rep(uvals.y[yrow]+h, 2), col=col)
                        xpoly <- c(rep(uvals.x[xline]-w, 2), rep(uvals.x[xline]+w, 2))
                        ypoly.tmp <- c(uvals.y[yrow]-h, uvals.y[yrow]+h)
                        ypoly <- c(ypoly.tmp, rev(ypoly.tmp))
                        polygon(xpoly, ypoly, border=col)
                }
        }


}



###################################################################

striplet.bounds.func <- function(w, k, s, min.width=1e-8){
	# striplet.bounds.func 3/4/09
	# Creates striplet boundaries so that every full strip in any systematic sample from the
	# grid of all possible systematic samples is covered by an exact number of striplets.
	# Striplets are not of the same sizes, unless (I think) line.int / (w/s) is an integer.
	#
	# However, b-values denoting the start of the B systematic sampling grids are all
	# equally spaced.
	#
	# w = full strip half-width
	# k = number of lines in one systematic sample
	# s = ideal number of striplets per half-strip (an integer).  The actual number of striplets
	# per half-strip will not be less than s, but might be more.  Note that w/s stripulates the
	# maximum allowed striplet width.
	#
	# min.width is the minimum striplet width accepted.  Note that the survey area is always
	# scaled to be on a scale from 0 to 1.  Sometimes this will be metres or km, sometimes it
	# might be thousands of km.  So min.width = 1e-8 has to be interpreted in the light of the
	# scale of the survey area.  A smaller min.width might be appropriate at times.
	#
	# EXAMPLE:
	# To match the 2006 settings:
	# striplet.bounds.func(w=0.02, k=20, s=40)
	# Gives a total of 2000 striplets, each of width 0.0005.
	#----------------------------------------------------------------------------------

	line.int <- (1 - 2*w) / k

	bvec <- seq(from = w, to = w + line.int, by = w/s)
	# bvec is the vector of of starting points for the grid of systematic samples.
	# Each systematic sample is bvec[i] + (0:(k-1)) * line.int  for some i.

	# Note that bvec might not finish exactly on w + line.int : it will stop at the nearest place below it.

	B <- length(bvec)
	# B is the number of starting values in the grid of systematic samples.

	striplet.pos <- numeric(0)

	# striplet.pos is the vector of striplet locations between 0 and 1.
	# First ensure that every transect line in the grid of systematic samples lands on a striplet boundary:
	for(b in 1:B) striplet.pos <- c(striplet.pos, bvec[b] + (0:(k-1)) * line.int)

	# Next ensure the LH border of every full strip lands on a striplet boundary:
	for(b in 1:B) striplet.pos <- c(striplet.pos, bvec[b] + (0:(k-1)) * line.int - w)

	# Finally ensure the RH border of every full strip lands on a striplet boundary:
	for(b in 1:B) striplet.pos <- c(striplet.pos, bvec[b] + (0:(k-1)) * line.int + w)

	# The vector striplet.pos now contains everything wanted.  Get into order and
	# weed out the duplicates:
	striplet.pos <- sort(unique(striplet.pos))

	# There will still be some "duplicates" less than min.width apart from each other.
	# Element i of striplet.diff corresponds to striplet.pos[i+1] - striplet.pos[i]:
	striplet.diff <- diff(striplet.pos)

	ndiff <- length(striplet.diff)
	remove.elts <- (1:ndiff)[striplet.diff < min.width]
	# e.g. if the i'th elt of striplet.diff is too small, can remove the i'th elt of striplet.pos.
	if(length(remove.elts)>0) striplet.pos <- striplet.pos[-remove.elts]

	# Plotting commands for checking: comment out when not required.
	# plot(-1, -1, xlim=c(0, 1), ylim=c(0, 1))
	# abline(v=striplet.pos, col=4)
	# for(b in 1:B) points(bvec[b] + (0:(k-1)) * line.int - w, rep(0.8, k))
	# for(b in 1:B) points(bvec[b] + (0:(k-1)) * line.int, rep(0.6, k))
	# for(b in 1:B) points(bvec[b] + (0:(k-1)) * line.int + w, rep(0.4, k))

	Nstriplets <- length(striplet.pos) - 1
	cat("Number of striplets: ", Nstriplets, "\n")
	cat("Summary of striplet widths: \n", summary(diff(striplet.pos)), "\n\n")
	striplet.pos
}


#################################################################################



image.grid.func <- function(grid, zname){
	# image.grid.func 6/5/09
	# Takes a grid which MUST be in the following format (although the numbers can be different):
	# 	x	y	zname
	#	1	1	:
	#	2	1	:
	#	3	1	:
	#	1	2	:
	#	2	2	:
	#	3	2	:
	#	etc
	# Plots an image with heat.colors so that the orientation is correct and plots as:
	# 	|
	#          y	|
	#	|___________
	#		x
	#
	# EXAMPLE:
	# delgrid <- expand.grid(x=1:3, y=10:14)
	# delgrid$n <- 1:15
	# image.grid.func(delgrid, "n")

	x.pos <- unique(grid$x)
	y.pos <- unique(grid$y)
	if(any(grid[,c("x", "y")] !=expand.grid(x=x.pos, y=y.pos))) stop("Grid supplied in wrong alignment")
	zmat <- matrix(grid[,zname], nrow=length(x.pos))
	image(x.pos, y.pos, zmat, col=heat.colors(20), xlab="x", ylab="y")
	title(main=zname)
}

#####################################################################################

persp.grid.func <- function(grid, zname, theta=330, phi=40, expand=0.3, ltheta=-120, shade=0.75){
	# persp.grid.func 6/5/09
	# Takes a grid which MUST be in the following format (although the numbers can be different):
	# 	x	y	zname
	#	1	1	:
	#	2	1	:
	#	3	1	:
	#	1	2	:
	#	2	2	:
	#	3	2	:
	#	etc
	# Plots a perspective plot with shading in the following alignment:
	# 	|
	#          y	|
	#	|___________
	#		x
	#
	# EXAMPLE:
	# delgrid <- expand.grid(x=1:3, y=10:14)
	# delgrid$n <- 1:15
	# persp.grid.func(delgrid, "n")

	x.pos <- unique(grid$x)
	y.pos <- unique(grid$y)
	if(any(grid[,c("x", "y")] !=expand.grid(x=x.pos, y=y.pos))) stop("Grid supplied in wrong alignment")
	zmat <- matrix(grid[,zname], nrow=length(x.pos))
	persp(x.pos, y.pos, zmat, theta=theta, phi=phi, expand=expand, ltheta=ltheta, shade=shade, xlab="x", ylab="y")
	title(main=zname)
}

####################################################################################

one.2d.ipp.func <- function(Npoints, cluster, covariates=jeff.covariates,
                            Nmax=round(qpois(1-1e-10, Npoints) + 0.25*Npoints))
{
	# one.2d.ipp.func 6/5/09
	# Like generate.2d.ipp.func in other directories, but generates only one population.
	#
	# ipp = inhomogeneous Poisson process - i.e. the intensity depends upon covariates and river.
	# This function assumes dependence on BOTH covariates AND river: need another function if
	# not interested in river.
	#
	# cluster = T or F.
	# If cluster = T then an overdispersed IPP will be generated, with clustering of objects,
	# using the same settings as Jeff used in his paper.  The cluster settings are hard-coded
	# into this function: model="gauss", cor.par=c(0.25, -25/log(0.05)).
	#
	# If cluster = F then no clustering takes place: the objects are independently distributed
	# according to the specified IPP.
	# In Jeff's paper, clustering and overdispersion seem to relate to the same thing.
	#
	# covariates is a pre-defined set of covariates:
	# jeff.covariates <- simCovariates(hab.range=30, probs=c(1/3, 2/3), river.loc=50)
	#
	# USAGE:
	#...
	#----------------------


	# Simulate Npoints realisations of points: deliberately overshoot on EN to pick
	# off the exact number of points desired.
	if(cluster==F){
		pts <- simPts(covariates=covariates, int.formula = ~factor(habitat)+river,
			int.par=c(0, 1, 2, -1), EN=Nmax)
	}
	else if(cluster==T){
		pts <- simPts(covariates=covariates, int.formula = ~factor(habitat)+river,
			int.par=c(0, 1, 2, -1), EN=Nmax, model="gauss",
			cor.par=c(0.25, -25/log(0.05)))
	}
	# Now get it into an x-vector of the correct length, and rescale so that the
	# coordinates are between 0 and 1:
	# Return both xvals and yvals
	pop.x <- pts$x[1:Npoints]/100
	pop.y <- pts$y[1:Npoints]/100
                 cat("Points generated = ", length(pts$x), "\n")

                 # Enough to check for NA's in pop.x only: pop.y will be the same.
	if(any(is.na(pop.x))) print("NA")

	list(x=pop.x, y=pop.y)
}



#################################################################################




mgcv.boxlet.func <- function(whkmNpars, cluster, Nsurv=1000, v=1, sx=max(1, round(w*100)),
		sy=max(1, round(h*100)), min.width=1e-8, est.df=T, dfval=missing(), use.musum=F,
                                   plotit=F, covariates=jeff.covariates, Nmax=NULL)
{
	# mgcv.boxlet.func 6/5/09
	#
	# NOTE: (March 2010) - I can't see any evidence that this has actually been
	# coded for rectangle regions (1 x v) as opposed to square regions (1 x 1):
	# don't apply to rectangle regions 1 x v without checking this first!!
	#
	# Like GAM/expt.boxlet.func and GAM/many.expt.boxlet.func, using mgcv instead of gam.
	# Performs Nsurv surveys with the "try" function to catch errors.
	#
	# If est.df=T, uses Simon's GCV to estimate the df.  Otherwise uses the supplied dfval.
	# Suggest could try df = k*m / 10, (one tenth the number of observations), but this is untested.
	#
	# use.musum : flag to decide whether to use musum (gained from the GAM) as the total number of
	# objects for the striplet multinomial distribution, or whether to use Nhat = Dhat*area = Dhat * v for
	# the rectangle area.  Note that Dhat has slightly lower variance, from sims, so default to using Dhat.
	# use.musum = T : use musum.
	# use.musum=F  (the default, UNLIKE for the distance sampling sims): use Dhat * v.

	if(v!=1) stop("I don't think it's properly encoded for v != 1 : better go read the code and check.")

	whkmNpars <- unlist(whkmNpars)
	w <- whkmNpars[1]
	h <- whkmNpars[2]
	k <- whkmNpars[3]
	m <- whkmNpars[4]
	N <- whkmNpars[5]

	print(c(w, h, k, m, N))
	#----------------------------------------------------------------------------------------------------
	# BOXLET SETUP: the setup of the boxlets, up to the calculation of gbmat, is tested and plotted
	# in directory STRIPLETS/BOXLETS with boxlet.setup.func.
	#
	# First set up boxlet quantities that only need to be calculated once globally.

	line.int <- (1 - 2*w) / k
	# Note that row.int needs to involve the height of the survey area, v:
	row.int <- (v - 2*h) / m

	# Old method using Nboxlet.desired to derive sx and sy: now superceded.
	# Use the value Nboxlet.desired to derive sx and sy:
	# sx = the desired number of striplets per half-width (an integer >= 1);
	# sy = the desired number of bandlets per half-height (an integer >= 1)
	# sx <- 1 + floor(w * sqrt(Nboxlet.desired/v))
	# sy <- 1 + floor(h * sqrt(Nboxlet.desired/v))

	print(c(sx=sx, sy=sy))

	# Set up bvec.x and bvec.y.
	# bvec.x is the vector of x-direction starting points for the grid of systematic samples.
	# bvec.y is the vector of y-direction starting points for the grid of systematic samples.
	# Each systematic sample has lines at bvec.x[i] + (0:(k-1)) * line.int  for some i,
	# and rows at bvec.y[j] + (0:(m-1)) * row.int  for some j.

	bvec.x <- seq(from = w, to = w + line.int, by = w/sx)
	bvec.y <- seq(from = h, to = h + row.int, by = h/sy)

	# Note that bvec might not finish exactly on w + line.int : it will stop at the nearest place below it.

	# Bx and By are the number of starting values in the grid of systematic samples, in the x and y directions.
	# B is the total number of starting values in the grid of systematic samples.
	Bx <- length(bvec.x)
	By <- length(bvec.y)
	B <- Bx * By

	# Find the boundary x-coordinates of the striplets:
	striplet.bounds.x <- striplet.bounds.func(w=w, k=k, s=sx, min.width=min.width)
	n.striplet.x <- length(striplet.bounds.x) - 1

	# Find the boundary y-coordinates of the striplets:
	striplet.bounds.y <- striplet.bounds.func(w=h, k=m, s=sy, min.width=min.width)
	n.striplet.y <- length(striplet.bounds.y) - 1

	# J is the total number of boxlets:
	J <- n.striplet.x * n.striplet.y

	cat("Number of boxlets tesselating region, J = ", J, "\n")
	cat("Number of grid start-points, B = ", B, "\n")

	# midvec is a data frame with J rows such that:
	# midvec$x gives the x-coordinates of the mid-points of the J boxlets;
	# midvec$y gives the y-coordinates of the mid-points of the J boxlets.
	#
	# For example, if striplet.bounds.x = c(0, 0.5, 1) and striplet.bounds.y = c(0, 0.2, 0.4, 0.6)
	# then midvec would be:
	#       	  x   	 y
	#  1 	0.25 	0.1
	#  2 	0.75 	0.1
	#  3 	0.25 	0.3
	#  4 	0.75 	0.3
	#  5 	0.25 	0.5
	#  6 	0.75 	0.5

	midvec.xpoints <- (striplet.bounds.x[-1] + striplet.bounds.x[-(n.striplet.x+1)])/2
	midvec.ypoints <- (striplet.bounds.y[-1] + striplet.bounds.y[-(n.striplet.y+1)])/2
	midvec <- expand.grid(x=midvec.xpoints, y=midvec.ypoints)


	# Next line for rectangle regions only:
	boxlet.area.vec <- rectangle.area2d.func(xlow=striplet.bounds.x[-(n.striplet.x+1)], xhi=striplet.bounds.x[-1],
			ylow=striplet.bounds.y[-(n.striplet.y+1)], yhi=striplet.bounds.y[-1])

	# Now note the following things that do not depend upon the individual surveys:
	# gbmat : matrix with J=Nboxlet rows and B columns, such that column b, gbmat[,b] is a vector
	#	of length J giving the mean detection probability for objects in this boxlet.
	# 	For quadrat counts, this will equate to the area of the boxlet that is covered by the grid at b:
	#	if boxlet j is fully covered (as it will be for rectangular boxlets in a rectangular region)
	#	then gbmat[j, b] will be either 0 (if j is not included in the grid at b) or 1 (if j is included).
	#	For irregularly shaped scenarios, for example circular plots, gbmat[j,b] will be the
	#	proportion of boxlet j's area that is covered - or might also include some detection probability.
	#	If boxlet j is not included in the grid under start-point b, then
	#	gbmat[j, b] = 0 and thus this boxlet contributes 0 to the sum of A(b) later on.
	# 	The matrix is set up to be J x B instead of the reverse for ease of coding:
	#	saves some matrix multiplications because we can then just multiply muvec * gbmat
	#	to get muvec * gbmat[,b] for each column b.
	#
	# For the rectangular region with rectangular quadrats, the output of striplet.bounds.func ensures that
	# every boxlet is wholly in or wholly out of the grid based on b.  Therefore it's sufficient that a boxlet
	# midvec$x is < w from a grid line, and that midvec$y is < h from a grid row, in order for this
	# boxlet to be included in the grid.  All other boxlets get gbmat[b, j] = 0.

	gbmat <- matrix(0, nrow=J, ncol=B)

	# b.ind is the overall index going from 1 to B, corresponding to b.indx and b.indy.
	b.ind <- 0
	for(b.indy in (1:By)){
		b.y <- bvec.y[b.indy]
		rows.grid <- seq(b.y, length=m, by=row.int)

		for(b.indx in (1:Bx)){
			b.x <- bvec.x[b.indx]
			lines.grid <- seq(b.x, length=k, by = line.int)
			b.ind <- b.ind + 1

			# For each boxlet mid-point, find its minimum distance from the grid lines
			# and its minimum distance from the grid rows:
			min.dist.x <- sapply(midvec$x, function(xval) min(abs(xval - lines.grid)))
			min.dist.y <- sapply(midvec$y, function(yval) min(abs(yval - rows.grid)))

			# A boxlet is included in the sample if the distance of the midpoint from the grid lines
			# is less than w, AND the distance of the midpoint from the grid rows is less than h.
			# This follows because each striplet is either wholly included or wholly
			# excluded.  gbmat[,b] for the other striplets is already initialised to 0.
			boxlets.included.b <- (1:J)[min.dist.x < w & min.dist.y < h]
			gbmat[boxlets.included.b, b.ind] <- 1
		}
	}


	#-------------------------------------------------------------------------------------------
	# Setup is complete.  Now run surveys.


	if(is.null(Nmax)) Nmax <- round(qpois(1-1e-10, N) + 0.25*N)
	survest.func <- function(surv){
                 		print(surv)

		jeffpops <- one.2d.ipp.func(Npoints=N, cluster=cluster, covariates=covariates, Nmax=Nmax)
		xvals <- jeffpops$x
                		yvals <- jeffpops$y
                                   ## Spit out an error if Nmax underdid it and there aren't N points here:
                                   if(any(is.na(xvals))) stop("Wrong number of points: quitting this survey")

		## The survey b-values are where the grid starts: "b" for "beginning".
		## b.x is the position of the centerline for the leftmost strip in the grid.
		## b.y is the position of the centerline for the bottom row in the grid.
		## b.x ~ Uniform(w, w+line.int).
		## b.y ~ Uniform(h, h+row.int).
		## uvals.x and uvals.y give the centrelines for the single realised survey.

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


		#---------------------------------------------------------------------------
		# Survey is complete.
		#
		#----------------------------------------------------------------------
		# VARIANCE ESTIMATION
		#----------------------------------------------------------------------
		# Striplet estimator:
		#----------------------------------------------------------------------

		# Area of rectangle region full-boxes with x-centrelines at uvals.x and half-width w,
		# and y-centrelines at uvals.y and half-height h: this was found earlier as area.surv.vec.
		# Note that this has the following ordering:
		#  	area.surv.vec elt	xline	yrow
		# 	  	1	1	1
		# 	  	2	2	1
		# 	  	3	3	1
		# 	  	:	:	1
		# 	  	k	k	1
		# 	  	k+1	1	2
		# 	  	k+2	2	2
		# 	  	:	:	:
		# whereas boxcount.mat has the following arrangement:
		# 	(  1     2     3     4     ...   k )
		#	(k+1 k+2 k+3  k+4 ... 2k )
		#		.....
		# 	(		 mk)
		# So we need to unpack boxcount.mat as as.vector(t(boxcount.mat)) to
		# get the correspondence correct:

		a.offset <- area.surv.vec
		nspotted <- as.vector(t(boxcount.mat))

		if(Nsurv==1 | plotit==T){
			par(mfrow=c(3,2))
			jeff.covariates.plot(covar=covariates)
			## draw.boxes.func(uvals.x=uvals.x, uvals.y=uvals.y, w=w, h=h)
			##
			plot(xvals, yvals, pch=16)
			draw.boxes.func(uvals.x=uvals.x, uvals.y=uvals.y, w=w, h=h, col=2)
		}

		# GAM to find fitted number of objects per unit area density, lambda:
		box.dat <- expand.grid(x=uvals.x, y=uvals.y)
		box.dat$n <-  nspotted
		box.dat$loga.offset <- log(a.offset)

		if(est.df==T) box.gam <- gam(n~offset(loga.offset) + s(x, y), family=poisson(link=log), data=box.dat)
		else box.gam <- gam(n~offset(loga.offset) + s(x, y, fx=TRUE, k=dfval+1),
			family=poisson(link=log), data=box.dat)

                                   if(Nsurv==1 | plotit==T) print(box.gam)
                                   # fx=TRUE, k=dfval+1 is stolen from the BTO code:

		# The 2-d predict.gam functions seem to be terrible when predicting for a
		# fine grid, both for gam() and mgcv().  So instead of using predict.gam with
		# newdata, we'll use predict.gam on the original data frame, then simply
		# interpolate the fitted surface.

		# box.lambda is the fitted object-available density: expected number of
		# objects available per unit area, predicted at the centrelines of the
		# boxes in the survey:

		box.terms <- predict.gam(box.gam, type="terms")
		box.lambda <- exp(box.terms + attributes(box.terms)$constant)

		# Now we need to interpolate this surface to the boxlet midpoints in midvec.
		# Do not try to extrapolate beyond the original data: interp seems to
		# crash when extrap=T anyway!  It also seems to crash when using spline interp,
		# via linear=F.

		# Create a surrogate data frame, surdat, for the coordinates of the interp.
		# When the boxlet midpoints are beyond the range of observed data, this
		# surrogate replaces the point for prediction by the closest within-range point.
		sur.x <- midvec.xpoints
		sur.y <- midvec.ypoints

		sur.x[sur.x < min(uvals.x)] <- min(uvals.x)
		sur.x[sur.x > max(uvals.x)] <- max(uvals.x)

		sur.y[sur.y < min(uvals.y)] <- min(uvals.y)
		sur.y[sur.y > max(uvals.y)] <- max(uvals.y)

		boxlet.lambda.mat <- interp(box.dat$x, box.dat$y, box.lambda, xo=sur.x, yo=sur.y)$z
		# boxlet.lambda.mat is presented as a matrix, such that as.vector(boxlet.lambda.mat)
		# unpacks it in the correct way for the boxlet grid in boxlet.dat.
		# Note that using interp on sur.x and sur.y, and then treating the answer as if it were
		# created on the genuine boxlet midpoints,
		# attributes the prediction using the surrogate x and y values to the required boxlet x and y values:

		boxlet.lambda <- as.vector(boxlet.lambda.mat)

		# boxlet.lambda is now a vector in the same ordering as the boxlet data frame midvec:
		# so we could say midvec$lambda <- boxlet.lambda and each row would then have a
		# (x, y, lambda) triple all correctly corresponding.

		# Now calculate muvec, which is the expected number of objects available per boxlet:
		muvec <- boxlet.lambda * boxlet.area.vec
		musum <- sum(muvec)

		# Finally calculate the boxlet quantities we're looking for:
		# Abvec = A(b) for each b in 1:B, where A(b) = sum_{j} mu_j g_j(b),
		# Here, g_j(b) = 0 for any boxlets j not included in the grid based on b,
		# so the contribution of such boxlets is 0.  Thus A(b) is equally well a sum over all j
		# or a sum over j \in S(b) : they are equivalent because the extra terms in the first are all 0.
		#

		mu.times.g.mat <- muvec * gbmat
		Abvec <- apply(mu.times.g.mat, 2, sum)

		# Put everything together into the striplet estimator:

		if(use.musum) var.n.striplet <- 1/B * sum(Abvec + Abvec^2*(1-1/musum))- (1/B * sum(Abvec))^2
		else{
			Nhat <- v * Dhat
			var.n.striplet <- 1/B * sum(Abvec + Abvec^2*(1-1/Nhat))- (1/B * sum(Abvec))^2
		}
		# Then Dhat = n / area.surv, so varD.striplet = var.n.striplet/area.surv^2.
		varD.striplet <- var.n.striplet/area.surv^2

		box.dat$lambda <- box.lambda
		boxlet.dat <- midvec
		boxlet.dat$lambda <- boxlet.lambda
		boxlet.dat$muvec <- muvec

		if(Nsurv==1 | plotit==T){
			##
			cat("muvec:\n")
			print(summary(muvec))
			image.grid.func(box.dat, zname="lambda")
			persp.grid.func(box.dat, zname="lambda")
			##
			image.grid.func(boxlet.dat, zname="lambda")
			persp.grid.func(boxlet.dat, zname="lambda")
			## Use the plots below to check whether the bilinear interp is tracking the
			## black-point box predictions in x and y directions:
			## zbig <- max(c(box.dat$lambda, boxlet.dat$lambda))
			## plot(boxlet.dat$x, boxlet.dat$lambda, cex=0.5, pch=16, col=2, ylim=c(0, zbig))
                                                    ## points(box.dat$x, box.dat$lambda, pch=16)
                                                    ## title(main="Marginal x : boxlets red, boxes black")
                                                    ##
                                                    ## plot(boxlet.dat$y, boxlet.dat$lambda, pch=16, cex=0.5, 
			                                                      # col=2, ylim=c(0, zbig))
                                                    ## points(box.dat$y, box.dat$lambda, pch=16)
                                                    ## title(main="Marginal y : boxlets red, boxes black")
			##
		}

		#------------------------------------------------------------------------------------------

		# ------------------------------------------------------------------------------------------
		# OTHER VARIANCE ESTIMATORS:
		# P1 : random-sample-based estimator for var(Dhat).  Equivalent to P1 or R1 in the ERVAR
		# paper.
		# O1 : Russell Millar's estimator with overlapping post-stratification. Each stratum is a 2x2 block of
		# points.
		# P4 : Marcello D'Orazio's estimator, P4 in the ERVAR paper. P4.hvd and P4.hv are respectively where
		# the "neighbours" of a point do (hvd) and don't (hv) include the diagonal neighbours.
		# In both cases they include the horizontal (h) and vertical (v) neighbours.
		#

		# ------------------------------------------------------------------------------------------
		# ESTIMATOR P1:
		# Dhat = sum_{i=1}^{km}  n_i / area.surv
		# var(Dhat) = 1/area.surv^2  * sum_{i=1}^{km} var(n_i)
		# 	   = 1/area.surv^2  * 1/(km-1) * km * var(n_1)
		# 	   = 1/area.surv^2  * 1/(km-1) * km * sum_{i=1}^{km} (n_i - nbar)^2
		# 	   = 1/area.surv^2  * km/(km-1) * sum_{i=1}^{km} (n_i - nbar)^2

		nbar <- mean(boxcount.mat)
		nvec <- as.vector(boxcount.mat)
		varD.P1 <- (1/area.surv^2) * k*m/(k*m-1) * sum((nvec-nbar)^2)

		# ------------------------------------------------------------------------------------------
		# ESTIMATOR O1:
		# First construct the strata and find stratum-specific variances.
		# For the rectangular region, this is very easy:
		strat.sum <- 0
		# strat.sum = sum_{strata from 1 to k*m) strat.var
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

		# ------------------------------------------------------------------------------------------
		# ESTIMATOR P4:  D'Orazio's estimator.
		# Try it two ways: one where "adjacent" points do not include diagonally adjacent (P1.hv)
		# and one where they do (P4.hvd).

		# Introduce a matrix boxcount.extend:
		# boxcount.extend = 	NA NA NA
		#			NA mat NA
		#			NA NA NA
		# where "mat" is the original boxcount.mat.

		boxcount.extend <- cbind(NA, boxcount.mat, NA)
		boxcount.extend <- rbind(NA, boxcount.extend, NA)

		delta.sum.hvd <- 0
		delta.nsq.sum.hvd <- 0
		delta.sum.hv <- 0
		delta.nsq.sum.hv <- 0

		for(xline in 1:k){
			for(yrow in 1:m){
				n.centre <- boxcount.mat[yrow, xline]
				# Need to add 1 to xs and ys because there's a column and row "0"
				# full of NAs taking up index 1.
				xs <- seq(xline-1, xline+1)+1
				ys <- seq(yrow-1, yrow+1)+1
				# mini.mat is a 3 x 3 matrix of neighbours of n.centre:
				mini.mat <- boxcount.extend[ys, xs]
				#----------------------------------------------
				# First the hvd calculations:
				mini.vec.hvd <- as.vector(mini.mat)
				mini.vec.hvd <- mini.vec.hvd[!is.na(mini.vec.hvd)]
				# Remember that mini.vec is one element too long,
				# because it includes the central element.
				# This extra point must be removed from deltasum.hvd,
				# but contributes nothing to delta.nsq.sum.hvd (because it's 0).
				delta.sum.hvd <- delta.sum.hvd + length(mini.vec.hvd) - 1
				delta.nsq.sum.hvd <- delta.nsq.sum.hvd + sum((mini.vec.hvd - n.centre)^2)
				#---------------------------------------------
				# Now the hv calculations:
				mini.mat.hv <- mini.mat
				mini.mat.hv[1, 1] <- NA
				mini.mat.hv[3, 1] <- NA
				mini.mat.hv[1, 3] <- NA
				mini.mat.hv[3, 3] <- NA
				mini.vec.hv <- as.vector(mini.mat.hv)
				mini.vec.hv <- mini.vec.hv[!is.na(mini.vec.hv)]
				# Remember that mini.vec.hv is one element too long,
				# because it includes the central element.
				# This extra point must be removed from deltasum.hv,
				# but contributes nothing to delta.nsq.sum.hv (because it's 0).
				delta.sum.hv <- delta.sum.hv + length(mini.vec.hv) - 1
				delta.nsq.sum.hv <- delta.nsq.sum.hv + sum((mini.vec.hv - n.centre)^2)

			}
		}
		varD.P4.hvd <- (1/area.surv^2) * k * m * delta.nsq.sum.hvd / (2 * delta.sum.hvd)
		varD.P4.hv <- (1/area.surv^2) * k * m * delta.nsq.sum.hv / (2 * delta.sum.hv)

		# From experiment, O1 and P4.hvd seem to be VERY similar (try a pairs plot):
		# more similar than P4.hvd and P4.hv, for example.  Wonder if that's always true.
		# The way they are derived is probably pretty similar.


		#-----------------------------------------------------------------------------------
		# Finally compile all estimators into the results from this survey:

		surv.result <- c(n = ntot, musum=musum, Dhat = Dhat,
				varD.striplet=varD.striplet, varD.P1 = varD.P1, varD.O1=varD.O1,
				varD.P4.hvd=varD.P4.hvd, varD.P4.hv=varD.P4.hv,
                                                                      b.x=surv.bval.x, b.y=surv.bval.y)
                                   # print(surv.result)
                                   if(plotit==T)readline()

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



#####################################################################################

oneres.boxplot <- function(res, sd.plot=T)
{
	# oneres.boxplot 4/5/09
	# Input the name directly,
	# EXAMPLE:
	# del <- mgcv.expt.boxlet.func(...)
	# oneres.boxplot(del, sd.plot=T)
	#

	trueval <- var(res$Dhat)

	comps <- list(res$striplet, res$P1, res$O1, res$P4.hvd, res$P4.hv)
	y.hi <- max(c(unlist(comps), trueval))
	y.lo <- min(c(unlist(comps), trueval)) - 1

	if(sd.plot == T) {
		comps <- lapply(comps, sqrt)
		y.hi <- sqrt(y.hi)
		y.lo <- sqrt(y.lo + 1) - 0.4
		trueval <- sqrt(trueval)
	}
	# Use help(bxp) to find some of the extra arguments available to enter
	# in pars=list() argument to boxplot.

	bp <- boxplot(comps, names = c("striplet", "P1", "O1", "P4.hvd", "P4.hv"), ylim = c(y.lo, y.hi),
		pars=list(medlty="blank", boxfill="light blue", outwex = 0.3, boxwex = 0.6,
			outpch=NA, outlty="solid", staplewex=0.6, staplelwd=2), cex.axis=1.8)
	abline(h = trueval, col = 2, lwd = 4)
	for(i in 1:length(comps))
		lines(c(i - 0.4, i + 0.4), rep(mean(comps[[i]]), 2), lwd = 2)

	print(c(true=trueval, means=lapply(comps, mean)))
}

################################################################################

boxlet.simgrid.func <- function(plotit=F){
	# boxlet.simgrid.func 6/5/09
	# Two-dimensional (non-exact) analogue of the settings that Jeff used in 1-dimension
	# for his simulations.  Using cover=0.5 seems a bit ridiculous when doing quadrat counts,
	# so change this to 0.25.
	# Fixes 	k=10 / 20;
	#	cover=0.04 / 0.25
	#	En = 75 / 400
	# and also m = 5 / 16
	# Then if cover = 0.04, sets w = 0.1/k and h=0.1/m
	# If cover = 0.25, sets
	#		w <- 0.5*1.25/(2*k)
	#		h <- 0.5/(1.25*2*m)
	# (These are selected to give rational w and h, and to multiply to the desired value for cover = 4whkm.)
	#
	# Then N = En/cover.
	# Thus (w, h, k, m, N) are recovered.
	#
	# USAGE:
	# boxlet.simgrid <- boxlet.simgrid.func()

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
		# Make N rounded to the nearest integer:
		N <- round(En / P)
		rachel.grid <- rbind(rachel.grid, c(w, h, k, m, N))
	}
	rachel.grid <- as.data.frame(rachel.grid)
	names(rachel.grid) <- c("w", "h", "k", "m", "N")
	rachel.grid <- rachel.grid[order(rachel.grid$N),]
	rownames(rachel.grid) <- 1:ngrid

	if(plotit){
		par(mfrow=c(2,2))
		for(i in 1:ngrid){
			k <- rachel.grid$k[i]
			m <- rachel.grid$m[i]
			N <- rachel.grid$N[i]
			w <- rachel.grid$w[i]
			h <- rachel.grid$h[i]
			pop <- one.2d.ipp.func(Npoints=N, cluster=T, covariates=jeff.covariates)
			plot(pop$x, pop$y, pch=16, main=i, xlim=c(0, 1), ylim=c(0, 1))
			draw.boxes.func(uvals.x=seq(w, by=(1-2*w)/k, length=k),
					uvals.y=seq(h, by=(1-2*h)/m, length=m), w=w, h=h, col=1)
			if(i%%4==0 & i!=16) readline()
		}
	}


	rachel.grid

}


#########################################################################

mgcv.boxlet.wrap <- function()
{
	# mgcv.boxlet.wrap 6/5/09
	# Wrapper for mgcv.boxlet.func : like a night.func()
	#
	for(sim in 1:16){
		name.sim <- paste("mgcv.noclus.sim", sim, sep="")
		print(name.sim)
		ass.only(name.sim, mgcv.boxlet.func(whkmNpars=boxlet.simgrid[sim,], cluster=F, Nsurv=3000))
		savetemp("mgcv.noclus.sim")
	}
	for(sim in 1:16){
		name.sim <- paste("mgcv.clus.sim", sim, sep="")
		print(name.sim)
		ass.only(name.sim, mgcv.boxlet.func(whkmNpars=boxlet.simgrid[sim,], cluster=T, Nsurv=3000))
		savetemp("mgcv.clus.sim")
	}

}

#########################################################################

mgcv.boxplot <- function(clus, sd.plot=T, Dtrue=F)
{
	# mgcv.boxplot 7/5/09
	# Plots results from simulations in directory BOXLETS/MGCV:
	#
	# EVENTUALLY, when coded ....
	# The mean estimated CV is printed over each boxplot, and the CI coverage is printed under
	# each boxplot.
	#
	# clus = T : the populations were generated from Jeff's functions with overdispersion (clustering);
	# clus = F : no overdispersion, so a genuine inhomogeneous Poisson process.
	#
	# Dtrue is an argument to jeff.Dhat.confintcover.func, to say whether the "true" value of D
	# is to be taken as the actual density of objects in the region (Dtrue = T), or the mean estimated
	# density of objects in the region (Dtrue=F).  The latter could be seen as the "true" value of D
	# under the repeated survey framework.
	#
	#
	# EXAMPLE:
	# load("SAVE.mgcv.noclus.sim")
	# mgcv.boxplot(clus=F)
	#

	par(mfrow=c(2,2))
	if(clus==T) namestem <- "mgcv.clus.sim"
	if(clus==F) namestem <- "mgcv.noclus.sim"

	object.names <- paste(namestem, 1:16, sep="")

	Npic <- 16
	picno <- 0
	resplot.func <- function(pic)
	{
		res <- eval(parse(text = object.names[pic]))
		parvec <- attributes(res)$inputs
		whkmNpars <- boxlet.simgrid[pic,]
		w <- whkmNpars[1]
		h <- whkmNpars[2]
		k <- whkmNpars[3]
		m <- whkmNpars[4]
		N <- whkmNpars[5]
		cover.percent <- 4 * w * h * k * m *100
		En <- round(cover.percent/100 * N)

		parname <- paste(pic, ": k=", k, " m=", m, " cov=", cover.percent, "% En=", En,
				" clus=", clus, sep="")

		trueval <- var(res$Dhat)

		comps <- list(res$striplet, res$P1, res$O1, res$P4.hvd, res$P4.hv)

		y.hi <- max(c(unlist(comps), trueval))
		y.lo <- min(c(unlist(comps), trueval)) - 1

		if(sd.plot == T) {
			comps <- lapply(comps, sqrt)
			y.top <- max(c(unlist(comps), sqrt(var(res$Dhat))))
			y.bottom <- min(c(unlist(comps), sqrt(var(res$Dhat))))
			y.hi <- y.top + (y.top - y.bottom)/4
			y.lo <- y.bottom - (y.top - y.bottom)/5
			trueval <- sqrt(trueval)
		}
		# Use help(bxp) to find some of the extra arguments available to enter
		# in pars=list() argument to boxplot.

		bp <- boxplot(comps, names = c("striplet", "P1", "O1", "P4.hvd", "P4.hv"), ylim = c(y.lo, y.hi),
			pars=list(medlty="blank", boxfill="light blue", outwex = 0.3, boxwex = 0.6,
				outpch=NA, outlty="solid", staplewex=0.6, staplelwd=2), cex.axis=1.5)
		abline(h = trueval, col = 2, lwd = 4)
		for(i in 1:length(comps))
			lines(c(i - 0.4, i + 0.4), rep(mean(comps[[i]]), 2), lwd = 2)
		title(main=parname, cex.main=1.5)

		#-----------------------------------------------------
		if(F){
		# Add CI coverage for nominal 95% intervals, and mean %CV, to the plot:
		confres <- jeff.Dhat.confintcover.func(simrow=pic, cluster=clus, conf=0.95, Dtrue=Dtrue,
				satter.striplet=satter.striplet)
		cic.vec <- confres["ci.cover",]
		if(any(names(cic.vec)!=c("striplet", "S1", "S2", "O1", "O2", "R2")))
			stop("Make sure cic.vec is ordered in the same order as the boxplots!!")
		cic.vec <- cic.vec[c("striplet", "S1", "O1", "R2")]

		# Make into integer percentage:
		cic.pc <- round(100 * cic.vec)
		for(i in 1:length(comps))
			text(i, y.lo + (y.hi - y.lo)/15, cic.pc[i], cex = 1.5)

		# Similarly for CV:
		cv.vec <- confres["mean.cv.D",]
		cv.vec <- cv.vec[c("striplet", "S1", "O1", "R2")]
		cv.pc <- round(100 * cv.vec)
		for(i in 1:length(comps))
			text(i, y.hi - (y.hi - y.lo)/15, cv.pc[i], cex = 1.5)
		#-----------------------------------------------------

		pic.results <- c(k=k, cover=cover.percent, En=En, gbar=gbar, true=trueval,
				striplet.mean=mean(comps[[1]]), striplet.sd=sqrt(var(comps[[1]])),
				S1.mean=mean(comps[[2]]), S1.sd=sqrt(var(comps[[2]])),
				O1.mean=mean(comps[[3]]), O1.sd=sqrt(var(comps[[3]])),
				R2.mean=mean(comps[[4]]), R2.sd=sqrt(var(comps[[4]])))
		}
		if(pic%%4 == 0 & pic!=Npic) readline("Press enter for next screen...\n")
		# pic.results
	}
	# resmat <- matrix(0, nrow=Npic, ncol=13)
	# for(i in 1:Npic) resmat[i,] <- resplot.func(i)
	# resmat <- as.data.frame(resmat)
	# names(resmat) <- c("k", "cover", "En", "gbar", "true", "striplet.mean", "striplet.sd",
	# 		"S1.mean", "S1.sd", "O1.mean", "O1.sd", "R2.mean", "R2.sd")
	# resmat
	for(i in 1:Npic) resplot.func(i)

}


###################################################################
mgcv.property.func <- function(cmdfunc, cmdarg, clus, ...)
{
	# mgcv.property.func 7/5/09
	# Quick way to scroll through all objects of type mgcv.clus.sim* or mgcv.noclus.sim*
	# and rattle off the same property for each of them.
	# EXAMPLES:
	# 1. To plot a histogram of Dhat for each noclus object:
	# par(mfrow=c(4,4))
	# mgcv.property.func(hist, "Dhat", clus=F, col=5)
	#
	# 2. To print the percentage error, (mean(Dhat) - N)/N * 100 for each noclus object:
	# mgcv.property.func(print, "(mean(Dhat)-N)/N*100", clus=F)
	#
	# 3. To print the dimensions of the objects themselves, e.g. to check number of failures via "try":
	# mgcv.property.func(print, "dim(obj)", clus=F)

	if(clus==T) namesvec <- mgcv.clus.names
	if(clus==F) namesvec <- mgcv.noclus.names
	for(i in 1:16){
		obj <- eval(parse(text=namesvec[i]))
		for(j in 1:length(names(obj)))
			assign(names(obj)[j], eval(parse(text=paste(namesvec[i], "$", names(obj)[j], sep=""))))
		assign("w", as.numeric(attributes(obj)$inputs["w"]))
		assign("h", as.numeric(attributes(obj)$inputs["h"]))
		assign("k", as.numeric(attributes(obj)$inputs["k"]))
		assign("m", as.numeric(attributes(obj)$inputs["m"]))
		assign("N", as.numeric(attributes(obj)$inputs["N"]))

		arg.in <- eval(parse(text=cmdarg))
		names(arg.in) <- NULL
		cmdfunc(arg.in, ...)

	}
}


####################################################################
