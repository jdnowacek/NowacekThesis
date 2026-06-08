## Conducts line transect distance sampling surveys with systematic lines on triangular shaped survey region.
## For each survey, estimates detection parameter theta and Dhat.  Then calculates the estimated var(er) and
## var(Dhat) using striplets as well as all the other estimators from Fewster et al (Biometrics, 2009).

## All code by Rachel Fewster, r.fewster@auckland.ac.nz, please don't distribute without permission.

## DEMO:
## For a quick run (1000 surveys):
##           tryme.res <- triangle.Dhat.func(c(w=0.02, theta=0.01, k=20, N=1000, superpop=2), Nsurv=1000)
##           triangle.oneres.boxplot(tryme.res)
## Should look good.  This is the run shown in Figure 3, second row, third panel ("Cover=80%")
## in the striplets paper.

## Have a look at the distribution of objects implied by "superpop=2" in the simulations above:
##           generate.superpop.func(superpop=2)

## The striplet work is done in the functions striplet.bounds.func and triangle.Dhat.func.
## striplet.bounds.func is called from within triangle.Dhat.func.

## In triangle.Dhat.func, the work is split into two parts: first there is a striplet setup procedure which
## only needs to be done once.  Then there is the striplet estimation procedure inside the subfunction
## survest.func: this involves fitting the GAM on the particular survey data generated, and calculating
## the striplet variance estimate for each survey.

## Outputs from triangle.Dhat.func are below: hopefully self-explanatory. Some notes: "er" stands for
## "encounter rate" and is n/L, varER.estimator is the estimated var(n/L) from that sim.
## Dhat=estimated density, that=estimated theta ("that"=theta-hat).
##
##  names(tryme.res)
##  [1] "n"              "L"              "er"             "musum"
##  [5] "dispersion"     "Dhat"           "that"           "f0hat"
##  [9] "var.f0hat.H2"   "varD.striplet"  "varD.R2"        "varD.R3"
## [13] "varD.S1"        "varD.S2"        "varD.O1"        "varD.O2"
## [17] "varER.striplet" "varER.R2"       "varER.R3"       "varER.S1"
## [21] "varER.S2"       "varER.O1"       "varER.O2"



library(gam)


halfnorm.func <- function(distvec, theta){
	exp( - distvec^2/(2 * theta^2))
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
	# Called from triangle.Dhat.func.
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
	# print(remove.elts)
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


###################################################################

gtheta.func <- function(w, theta)
{
	# calculates gtheta for the half-normal detection function
	# gtheta = 1/w * int_0^w halfnorm.func(x, theta) dx.
	#
	(sqrt(2 * pi) * theta * (pnorm(w/theta) - 0.5))/w
}


########################################################################

generate.pop.func <- function(x.beta1, x.beta2, y.beta1, y.beta2, plotit = T, Ntot = 1000)
{
	# generate.pop.func  8/12/04
	# EXAMPLES:
	# pop1 <- generate.pop.func(1, 1, 1, 1)
	# pop2 <- generate.pop.func(1, 4, 1, 1)
	# pop3 <- generate.pop.func(4, 1, 1, 1)
	# pop4 <- generate.pop.func(1, 1, 1, 4)
	# pop5 <- generate.pop.func(1, 1, 4, 1)
	# pop6 <- generate.pop.func(4, 1, 1, 4)
	#
	x <- rbeta(3 * Ntot, x.beta1, x.beta2)
	y <- rbeta(3 * Ntot, y.beta1, y.beta2)
	pop <- data.frame(x = x[y < x], y = y[y < x])
	while(nrow(pop)<Ntot){
		x <- c(x, rbeta(Ntot, x.beta1, x.beta2))
		y <- c(y, rbeta(Ntot, y.beta1, y.beta2))
		pop <- data.frame(x = x[y < x], y = y[y < x])
	}
	pop <- pop[1:Ntot,  ]
	# if(is.na(mean(pop$x)))
	#	stop("increase sample size")
	if(plotit == T) {
		plot(pop$x, pop$y, xlim = c(0, 1), ylim = c(0, 1), pch = 16,
			xlab = "", ylab = "", cex = 0.5, col=1)
		# points(pop$x, pop$y, col = 2, cex = 0.5)
		# print(pop)
	}
	pop
}


###################################################################


generate.superpop.func <- function(superpop, plotit = T, Ntot = 1000)
{
	# generate.superpop.func 8/4/09
	# Wrapper for generate.pop.func so only the superpop needs to be entered.

	switch(superpop,
		"1" = generate.pop.func(1, 1, 1, 1, plotit = plotit, Ntot = Ntot),
		"2" = generate.pop.func(1, 4, 1, 1, plotit = plotit, Ntot = Ntot),
		"3" = generate.pop.func(4, 1, 1, 1, plotit = plotit, Ntot = Ntot),
		"4" = generate.pop.func(1, 1, 1, 4, plotit = plotit, Ntot = Ntot),
		"5" = generate.pop.func(1, 1, 4, 1, plotit = plotit, Ntot = Ntot),
		"6" = generate.pop.func(4, 1, 1, 4, plotit = plotit, Ntot = Ntot),
		default = stop("superpop must be 1, 2, 3, 4, 5, or 6."))

}



#########################################################################


triangle.area.func <- function(boundary.low, boundary.high){
	# triangle.area.func 4/4/09
	# One of potentially many functions to calculate the area of a set of strips or
	# striplets, for a particular shape of survey region (triangle in this case),
	# given the boundary x-locations of the set.
	# boundary.low is a vector of the lower boundaries of the strips or striplets.
	# boundary.low is a vector of the upper boundaries of the strips or striplets.

	# EXAMPLE:
	# Take the output from striplet.boundary.func, striplet.pos: this is a single
	# vector where the upper boundary of the first striplet coincides with the lower
	# boundary of the second striplet.   There are Nstriplets altogether.  So apply this
	# function using:
	# 	triangle.area.func(striplet.pos[-(Nstriplets+1)], striplet.pos[-1])

	# Alternatively, if we have a set of full strips starting at positions uvec = (u_1, ..., u_k),
	# each of half-width w, then apply this function using
	# 	triangle.area.func(uvec - w, uvec + w)

	# In the triangle case, calculating the area is very easy, but this code is produced
	# so that the same function structure can be used to find the area for much more
	# complicated regions.
	#----------------------------------------------------------------------------------

	# The area of a striplet with lower bound a and upper bound b is
	# 0.5 * (b-a) * (b+a) = 0.5 * (b^2 - a^2)

	0.5 * (boundary.high^2 - boundary.low^2)

}


###################################################################

triangle.length.func <- function(line.pos){
	# triangle.length.func 4/4/09
	# One of potentially many functions to calculate the line lengths for a
	# set of transect lines, for a particular shape of survey region (triangle in this case),
	# given the line x-locations of the set.
	# line.pos is the x-coordinate (between 0 and 1) of the transect lines in the set.
	#
	# NOTE: triangle.length.func gives lvec, not L (in contrast to earlier versions).
	#
	# EXAMPLE:
	# 	lvec <- triangle.length.func(uvec)

	# In the triangle case, calculating the length is very easy, but this code is produced
	# so that the same function structure can be used to find the length for much more
	# complicated regions.
	#----------------------------------------------------------------------------------

	line.pos

}



#####################################################################################


triangle.Dhat.func <- function(wthetakNpoppars, Nsurv=1000, s=round(2000*w), min.width=1e-8, dfval=4){
	# triangle.Dhat.func 13/7/09
	#
	# Estimates theta and computes Dhat and varhat(Dhat) using the different varhat(n/L) methods.
	#
	# Takes the argument wthetakNpoppars = c(w, theta, k, N, superpop) for ease from BATCH/wrap running
	# or consistency with the Jeff sims.
	#
	# NOTE: uses triangle.length.func and triangle.area.func.
	#
	# For triangle survey area, takes Nsurv surveys and constructs:
	#	- true var(n/L) from Nsurv surveys
	#	- estimated var(n/L) using striplet estimator
	#	- estimated var(n/L) from estimators R2, R3, S1, S2, O1, and O2
	#	- Dhat, f0hat, and varhat(f0hat) = var.f0hat.H2 where H2 uses the 2nd derivative
	#		rather than the 1st derivative squared to approximate the Hessian for
	#		theta-hat.
	#	- varhat(Dhat) using estimators striplet, R2, R3, S1, S2, O1, O2
	#
	# w = full strip half-width
	# theta = true detection parameter for half-normal detection function (but to be estimated here)
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
	# dfval is the df for the smooth term in the gam() call.
	#
	# --------------------------------------------------------------------------------------------------
	# Set up parameters:

	wthetakNpoppars <- unlist(wthetakNpoppars)
	w <- wthetakNpoppars[1]
	theta <- wthetakNpoppars[2]
	k <- wthetakNpoppars[3]
	N <- wthetakNpoppars[4]
	superpop <- wthetakNpoppars[5]

	print(c(w, theta, k, N))
	#----------------------------------------------------------------------------------------------------
	# First set up striplet quantities that only need to be calculated once globally.

	line.int <- (1 - 2*w) / k

	bvec <- seq(from = w, to = w + line.int, by = w/s)
	# bvec is the vector of starting points for the grid of systematic samples.
	# Each systematic sample is bvec[i] + (0:(k-1)) * line.int  for some i.

	# Note that bvec might not finish exactly on w + line.int : it will stop at the nearest place below it.

	B <- length(bvec)
	# B is the number of starting values in the grid of systematic samples.

	# Find the boundary x-coordinates of the striplets:
	striplet.bounds <- striplet.bounds.func(w=w, k=k, s=s, min.width=min.width)
	J <- length(striplet.bounds) - 1
	midvec <- (striplet.bounds[-1] + striplet.bounds[-(J+1)])/2
	# midvec is the vector of striplet mid-points.  J is the number of striplets.

	# Next line for triangle region only:
	striplet.area.vec <- triangle.area.func(striplet.bounds[-(J+1)], striplet.bounds[-1])

	# Now note the following things that do not depend upon the individual surveys:
	# Lbvec : vector of length B such that Lbvec[b] = total line length for survey beginning at b
	# dist.jb.mat : matrix with J rows and B columns, such that column b, dist.jb.mat[,b] is a vector
	#	of length J=Nstriplets giving the distance of (the midpoint of) each striplet from the grid.
	#	Once theta is estimated, this matrix can then be used to generate gbmat quickly
	#	for each survey (see within-survey code below).
	#	For ease of coding, we dispense with summing over j \in S(b) and instead use
	#	dist.jb.mat (followed by gbmat) to do all the work.  If striplet j is not included in the grid
	#	under start-point b, then we set dist.jb.mat[j, b] = -1 and later gbmat[j, b] = 0,
	#	and thus this striplet contributes 0 to the sum of A(b).
	# 	The matrix is set up to be J x B instead of the reverse for ease of coding:
	#	saves some matrix multiplications because we can later just multiply muvec * gbmat
	#	to get muvec * gbmat[,b] for each column b.
	#
	# Note that striplets are set up by striplet.bounds.func so that every striplet should be either
	# wholly in or wholly out of the grid based on b.  Therefore it's sufficient that a striplet
	# mid-point is < w from a grid line in order for this striplet to be included in the grid.  All other
	# striplets get dist.jb.mat[j, b] = -1, and later gbmat[j, b] = 0.

	# Initialise dist.jb.mat to -1 everywhere:
	dist.jb.mat <- (-matrix(1, nrow=J, ncol=B))
	Lbvec <- rep(0, B)

	for(b.ind in (1:B)){
		b <- bvec[b.ind]
		lines.grid <- seq(b, length=k, by = line.int)

		# Calculate L(b):
		Lbvec[b.ind] <- sum(triangle.length.func(lines.grid))

		# For each striplet mid-point, find its minimum distance from the grid lines:
		min.dist.midpoints <- sapply(midvec, function(x) min(abs(x - lines.grid)))

		# A striplet is included in the grid if the distance of the midpoint from the grid lines
		# is less than w.  This follows because each striplet is either wholly included or wholly
		# excluded.  The included striplets will get dist.jb.mat equal to the distance of the striplet
		# from the grid.  dist.jb.mat for the other striplets is already initialised to -1.
		striplets.included.b <- (1:J)[min.dist.midpoints < w]
		dist.jb.mat[striplets.included.b, b.ind] <- min.dist.midpoints[striplets.included.b]

	}


	#-------------------------------------------------------------------------------------------
	# Setup is complete.  Now run the surveys.

	survest.func <- function(surv){
		#print(surv)

		xvals <- switch(superpop,
			"1" = generate.pop.func(1, 1, 1, 1, plotit = F, Ntot = N),
			"2" = generate.pop.func(1, 4, 1, 1, plotit = F, Ntot = N),
			"3" = generate.pop.func(4, 1, 1, 1, plotit = F, Ntot = N),
			"4" = generate.pop.func(1, 1, 1, 4, plotit = F, Ntot = N),
			"5" = generate.pop.func(1, 1, 4, 1, plotit = F, Ntot = N),
			"6" = generate.pop.func(4, 1, 1, 4, plotit = F, Ntot = N),
			default = stop("population must be 1, 2, 3, 4, 5, or 6."))$x

		## The survey b-value is where the grid starts: "b" for "beginning".
		## b is the position of the centerline for the leftmost strip in the grid.
		## b ~ Uniform(w, w+line.int).
		## uvals give the transect positions for the single realised survey.

		surv.bval <- runif(1, w, w + line.int)
		uvals <- surv.bval + seq(from = 0, by = line.int, length = k)

		# line length for the triangle region for this survey:
		lvec <- triangle.length.func(uvals)
		L <- sum(lvec)

		# There are probably improvements that could be made in the survey coding.
		# Leave it for now...
		spotted.func <- function(u)
		{
			x.avail <- xvals[abs(xvals - u) <= w]
			# Hard-coding the half-normal detection function here, for (possible) efficiency:
			seen <- rbinom(length(x.avail), 1, exp( - (x.avail -u)^2/(2 * theta^2)))
			x.seen <- x.avail[seen == 1]
			x.seen
		}
		nspotted <- rep(0, k)
		all.dists <- numeric(0)
		for(i in 1:k) {
			# need to waltz with nspotted[i] in order to avoid
			# complications due to the twin possibilities of spotted.i
			# being numeric(0) (length 0) and being NA (length 1):
			spotted.i <- spotted.func(uvals[i])
			all.dists <- c(all.dists, abs(spotted.i - uvals[i]))
			nspotted[i] <- length(spotted.i)
			if(nspotted[i] > 0)
				nspotted[i] <- length(na.omit(spotted.i))
		}
		ntot <- sum(nspotted)

		#---------------------------------------------------------------------------
		# Survey is complete.  Now carry out estimation of theta:

		## Estimate theta NOT by maximizing the likelihood numerically,
		## but by solving the score equation for dL/dtheta = 0 using uniroot:

		# hist(all.dists)
		# hist.max <- hist(all.dists, plot = F)$counts[1]
		# hist.scale <- hist.max/halfnorm.func(x = 0, theta = theta)
		# lines4(seq(0, w, length = 50), halfnorm.func(x = seq(0, w, length = 50),
		# 	theta = theta) * hist.scale, lwd = 2)

		score.func <- function(tc)
		{
			# tc is the theta argument over which maximization is to occur
			sum(all.dists^2)/tc^3 - ntot/tc + (ntot * w * dnorm(w/tc, 0,
				1))/(tc^2 * (pnorm(w/tc, 0, 1) - 0.5))
		}

		tc.upper <- 100 * w
		tc.lower <- 0.01 * w
		## Problem with solving the score function: if estimated theta
		## is perfect detection, it's not identifiable and the score
		## never reaches 0.  Sort out this problem by arbitrarily setting
		## theta-hat=100: not ideal but gets us there.
		##
		if((score.func(tc.upper) > 0) & (score.func(tc.lower) > 0)) {
			tc.res <- list(root = 100, message = "Setting theta-hat=100")
		}
		else tc.res <- uniroot(score.func, interval=c(tc.lower, tc.upper), tol=1e-10)

		that <- tc.res$root

		#----------------------------------------------------------------------------------
		## Now find the Hessian H(theta) using
		## H.2nd.deriv = - E( d^2 f(X, theta) / d theta^2 )
		## as in the ERVAR paper.  Don't bother with the H.1st.deriv.sq method for this work.
		## Here, X is a single observed distance;
		## f(X, theta) = g(X, theta) / int_0^w g(r, theta) dr   (the distance pdf)
		## and the Hessian (which is a scalar here) is estimated by
		## H.2nd.deriv = 1/n sum_{i=1}^n d^2 f(x_i, theta) / d theta^2
		##
		## See ERVAR notes "var(Dhat)" especially the final summary on p.8.
		##

		dnorm.that <- dnorm(w/that, 0, 1)
		pnorm.that <- pnorm(w/that, 0, 1)
		pnorm.m.half <- pnorm.that - 0.5
		H.2nd.deriv.brackets <- -2 + w^2/that^2 + (w * dnorm.that)/(that * pnorm.m.half)
		H.2nd.deriv.thetaterm <- 1/that^2 + (w * dnorm.that * H.2nd.deriv.brackets)/(that^3 * pnorm.m.half)
		H.2nd.deriv <- -1/ntot * sum((-3 * all.dists^2)/that^4 + H.2nd.deriv.thetaterm)


		## Now find f0hat = 1/(sqrt(2*pi) * that * pnorm.m.half)
		## and f0hat.1st.deriv = d f0hat /d theta,
		## and Dhat = (1/2)* (n/L) * f0hat: see p. 53-54 of Intro to Dist Sampling.  This expression
		## for Dhat does NOT depend on the triangle region: it's general.  The 2 comes from
		## the covered area being 2wL and the encounter rate being just n/L.

		f0hat <- 1/(sqrt(2 * pi) * that * pnorm.m.half)
		f0hat.1st.deriv <- ((w * dnorm.that)/(that * pnorm.m.half) - 1)/(sqrt(2 * pi) * that^2 * pnorm.m.half)
		erhat <- ntot/L
		Dhat <- 0.5 * erhat * f0hat


		## Now var(f0hat) = f0hat.1st.deriv^2 / (ntot * H)
		var.f0hat.H2ndderiv <- f0hat.1st.deriv^2/(ntot * H.2nd.deriv)

		##
		## The variance of Dhat will be enumerated after all the different
		## var(n/L)'s have been found, at the end of the function.
		##


		#----------------------------------------------------------------------
		# ENCOUNTER RATE VARIANCE ESTIMATION
		#----------------------------------------------------------------------
		# Striplet estimator:
		#----------------------------------------------------------------------

		# Estimate gtheta using theta = that:
		gthet.hat <- gtheta.func(w=w, theta=that)
		g.offset <- rep(gthet.hat, k)

		# Area of triangle region full-strips with centrelines at uvals and half-width w:
		a.offset <- triangle.area.func(uvals - w, uvals + w)


		# GAM to find fitted number of objects per unit area density, lambda:
		dat.gam <- data.frame(n = nspotted, u = uvals)

		formula.name <- paste("n ~ offset(log(g.offset) + log(a.offset)) + s(u, df = ", dfval, ")", sep="")
		form.expr <- eval(parse(text=formula.name))

		# Fit the GAM using quasi to get the dispersion parameter out of it:
		# have checked (on Serengeti sims) that the results are identical to family=poisson.
		fit.gam <- gam(form.expr, family=quasipoisson(link=log), data=dat.gam)
		dispersion.surv <- summary(fit.gam)$dispersion

		# plot(uvals, nspotted, pch=16, col=2)
		# lines(uvals, fitted(fit.gam))

		# If we wanted to calculate lambda.strip, where
		# lambda.strip is the fitted object-available density: expected number of
		# objects available per unit area, predicted at the centrelines of the
		# strips in the survey, then we would use the following:
		#
		# terms.gam <- predict.Gam(fit.gam, type="terms")
		# lambda.strip <- exp(terms.gam + attributes(terms.gam)$constant)
                                   # plot(uvals, lambda.strip, type="l")
                                   #
		# However, we actually want to calculate lambda at the locations of the
		# striplet midpoints, not the strip midpoints.  So we need to predict.Gam
		# with newdata given by the striplet midpoints.

		newdat.gam <- data.frame(u = midvec, g.offset=rep(gthet.hat, J), a.offset = striplet.area.vec)
		striplet.pred.gam <- predict.Gam(fit.gam, type="terms", newdata=newdat.gam)
		lambda.striplet <- exp(striplet.pred.gam + attributes(striplet.pred.gam)$constant)

		# Now calculate muvec, which is the expected number of objects available per striplet:
		muvec <- as.vector(lambda.striplet * striplet.area.vec)
		musum <- sum(muvec)

		# Finally calculate the striplet quantities we're looking for:
		# Abvec = A(b) for each b in bvec, where A(b) = sum_{j} mu_j g_j(b),
		# Here, g_j(b) = 0 for any striplets j not included in the grid based on b,
		# so the contribution of such striplets is 0.  Thus A(b) is equally well a sum over all j
		# or a sum over j \in S(b) : they are equivalent because the extra terms in the first are all 0.
		#
		# gbmat.surv is the estimated gbmat for this survey, using the estimated theta "that"
		# as the parameter of the detection function.
		# gbmat.surv is a matrix with J=Nstriplet rows and B columns, such that column b,
		# gbmat.surv[,b] is a vector of length J giving the midpoint detection probability for objects
		# in this striplet assuming that the detection parameter is equal to that.

		gbmat.surv <- matrix(0, nrow=J, ncol=B)
		gbmat.surv[dist.jb.mat >= 0] <- halfnorm.func(distvec=dist.jb.mat[dist.jb.mat >= 0], theta=that)

		mu.times.g.mat <- muvec * gbmat.surv
		Abvec <- apply(mu.times.g.mat, 2, sum)

		# Put everything together into the striplet estimator:

		var.striplet <- 1/B * sum((Abvec + Abvec^2*(1-1/musum) )/Lbvec^2)- (1/B * sum(Abvec/Lbvec))^2

		#-----------------------------------------------------------------------------------
		# Other estimators from the first ERVAR paper:

		#-----------------------------------------------------------------------------------
		# Estimator R2, based on the assumption of a random sample of lines:

		var.R2 <- (k * sum(lvec^2 * (nspotted/lvec - ntot/L)^2))/(L^2 * (k -1))

		#-----------------------------------------------------------------------------------
		# Estimator R3: included because this is still used in Distance5:

		var.R3 <- 1/(L * (k - 1)) * sum(lvec * (nspotted/lvec - ntot/L)^2)

		#-----------------------------------------------------------------------------------
		# Stratified estimators with non-overlapping strata: S1 and S2:
		#
		# First group the lines into strata, so that all strata have two lines,
		# but the last stratum has three if there is an odd number of lines:
		#
		H <- floor(k/2)
		k.h <- rep(2, H)
		if(k %% 2 > 0)
			k.h[H] <- 3
		end.strat <- cumsum(k.h)
		begin.strat <- cumsum(k.h) - k.h + 1

		sum.S1 <- 0
		sum.S2 <- 0
		for(h in 1:H) {
			nvec.strat <- nspotted[begin.strat[h]:end.strat[h]]
			lvec.strat <- lvec[begin.strat[h]:end.strat[h]]
			nbar.strat <- mean(nvec.strat)
			lbar.strat <- mean(lvec.strat)
			##########################
			## S1 calculations:
			inner.strat.S1 <- sum((nvec.strat - nbar.strat - (ntot/L) *
				(lvec.strat - lbar.strat))^2)
			sum.S1 <- sum.S1 + k.h[h]/(k.h[h] - 1) * inner.strat.S1
			##########################
			## S2 calculations: note that we use estimator R2 within
			## each stratum:
			L.strat <- sum(lvec.strat)
			var.strat.S2 <- k.h[h]/(L.strat^2 * (k.h[h] - 1)) * sum(
				lvec.strat^2 * (nvec.strat/lvec.strat - nbar.strat/lbar.strat)^2)
			sum.S2 <- sum.S2 + L.strat^2 * var.strat.S2
		}
		var.S1 <- sum.S1/L^2
		var.S2 <- sum.S2/L^2


		#-----------------------------------------------------------------------------------
		# Stratified estimators with overlapping strata: O1, O2:

		lvec.1 <- lvec[-k]
		lvec.2 <- lvec[-1]
		nvec.1 <- nspotted[-k]
		nvec.2 <- nspotted[-1]
		ervec.1 <- nvec.1/lvec.1
		ervec.2 <- nvec.2/lvec.2

		#########################
		## O1 calculations:
		overlap.varterm <- (nvec.1 - nvec.2 - ntot/L * (lvec.1 - lvec.2))^2
		var.O1 <- k/(2 * L^2 * (k - 1)) * sum(overlap.varterm)
		#########################
		## O2 calculations:
		V.overlap.R2 <- ((lvec.1 * lvec.2)/(lvec.1 + lvec.2))^2 * (ervec.1 -ervec.2)^2
		var.O2 <- (2 * k)/(L^2 * (k - 1)) * sum(V.overlap.R2)

		# --------------------------------------------------------------------------------------
		##
		## Now find the variance of Dhat.  Each calculation uses the H2ndderiv
		## method of approximating the Hessian.  Don't bother with the H1sq method
		## for this work.

		varD.striplet <- Dhat^2 * (var.striplet/erhat^2 + var.f0hat.H2ndderiv/f0hat^2)
		##
		varD.R2 <- Dhat^2 * (var.R2/erhat^2 + var.f0hat.H2ndderiv/f0hat^2)
		##
		varD.R3 <- Dhat^2 * (var.R3/erhat^2 + var.f0hat.H2ndderiv/f0hat^2)
		##
		varD.S1 <- Dhat^2 * (var.S1/erhat^2 + var.f0hat.H2ndderiv/f0hat^2)
		##
		varD.S2 <- Dhat^2 * (var.S2/erhat^2 + var.f0hat.H2ndderiv/f0hat^2)
		##
		varD.O1 <- Dhat^2 * (var.O1/erhat^2 + var.f0hat.H2ndderiv/f0hat^2)
		##
		varD.O2 <- Dhat^2 * (var.O2/erhat^2 + var.f0hat.H2ndderiv/f0hat^2)
		##

		#-----------------------------------------------------------------------------------
		# Finally compile all estimators into the results from this survey:

		surv.result <- c(n = ntot, L = L, er = ntot/L, musum=musum, dispersion=dispersion.surv,
				Dhat = Dhat, that = that,
				f0hat = f0hat, var.f0hat.H2 = var.f0hat.H2ndderiv,
				varD.striplet=varD.striplet, varD.R2 = varD.R2, varD.R3=varD.R3,
				varD.S1=varD.S1, varD.S2=varD.S2, varD.O1=varD.O1, varD.O2=varD.O2,
				varER.striplet=var.striplet, varER.R2 = var.R2, varER.R3 = var.R3,
				varER.S1=var.S1, varER.S2=var.S2, varER.O1=var.O1, varER.O2=var.O2)
		surv.result
	}
	all.results <- t(sapply(1:Nsurv, survest.func))
	row.names(all.results) <- 1:Nsurv
	all.results <- data.frame(all.results)
	names(all.results) <- c("n", "L", "er", "musum", "dispersion", "Dhat", "that", "f0hat", "var.f0hat.H2",
				"varD.striplet", "varD.R2", "varD.R3", "varD.S1", "varD.S2",
				"varD.O1", "varD.O2",
				"varER.striplet", "varER.R2", "varER.R3", "varER.S1", "varER.S2",
				"varER.O1", "varER.O2")
	attributes(all.results)$inputs <- c(w, theta, k, N, superpop)
	names(attributes(all.results)$inputs) <- c("w", "theta", "k", "N", "superpop")
	all.results


}


##############################################################################

triangle.Dhat.confintcover.func  <- function(res, conf = 0.95, Dtrue = F, plotit=F, satter.striplet=F)
{
	# triangle.Dhat.confintcover.func 13/7/09
	#
	# Takes a single simulation result "res" (output from triangle.Dhat.func) and finds the
	# CI coverage for density, D.
	# The triangle aspect of things only enters when computing Dtrue using Dtrue=2N (see below).
	# Search for all instances of "true.D.triangle" to modify.  Can generalise this function by tidying this up.
	#
	# The true density is equal to N/0.5=2N objects per unit area, because we have a triangle region.
	# HOWEVER, in some samples, true density is estimated with bias.
	# So it could be considered unfair on our variance estimators to use true D: makes it
	# look as if the coverage is bad, when in fact it's the erhat that's got the problem.
	# So we default to use Dtrue=F, and use the mean of the estimated D's instead,
	# as being the "true" thing we try to estimate in THIS inferential framework
	# of repeated surveys.  If Dtrue = T, then Dtrue is given the value = N/2 (=D for triangle region).
	#
	# [Comment from the ERVAR paper about the same issue:
	# Using Dtrue=F does at least make the CI coverage consistent
	# with the direction of the variance estimation, i.e. overestimated variance leads
	# to CIC > 0.95.  With Dtrue=T, get overestimated variance and CIC as low as 0.89,
	# which doesn't seem to make sense.]  For the striplet paper - using Dtrue = F has the similar
	# advantage that when the striplet estimator is bang-on the true variance, and has very
	# small variance, then the coverage reported is 95% rather than 93% or 97%!
	#
	#
	# Calculates the confidence interval coverage in the repeated survey framework
	# according to the log-Normal approximation for Dhat.
	# If plotit==T, also plots the histograms showing log-Normal distributions for Dhat:
	# 	- in red lines, show the log-Normal distribution with meanlog = mean(log(estimated Dhats))
	# 	- in dashed blue lines, show the density of the log-Normal distribution with mean
	# 	  equal to the true D.
	# The difference equates to the difference between using Dtrue=F and Dtrue=T: are we trying to
	# cover the (possibly imperfect) mean of Dhat in the repeated survey framework, or trying to cover
	# the true known value of D?  Obviously the latter in real life, BUT, to decide upon the performance
	# of our variance estimators, the former is more suitable, so that imperfections in a different estimation
	# process (that of Dhat) don't make our variance estimator look inappropriately bad.
	#
	# satter.striplet : if satter.striplet = T, the Satterthwaite approximation to the df is used for ALL
	# estimators, INCLUDING striplets.
	# If satter.striplet=F (the default), the Satterthwaite approximation is used for all other estimators,
	# R2, R3, O1, O2, S1, S2; but NOT for striplets.  For striplets, the df is taken to be infinity,
	# so that the distribution used to construct the confidence intervals is genuinely log-Normal rather
	# than log-Normal with a funny student-t-bit thrown in.
	# Could also dispense with the Satterthwaite approximation for the other estimators too - it doesn't
	# make much difference.  The worry is that the Satterthwaite approx is the only way to distinguish
	# between the higher df of O1, O2 relative to S1, S2.  So dispensing with it for these estimators seems
	# to be a worry because it means ignoring the whole df problem altogether.  The precision of the striplet
	# estimator doesn't seem to be particularly affected by the number of lines, so worrying about the df
	# seems to be more of a hindrance than a help.

	res.inputs <- attributes(res)$inputs
	w <- res.inputs["w"]
	theta <- res.inputs["theta"]
	k <- res.inputs["k"]
	N <- res.inputs["N"]
	superpop <- res.inputs["superpop"]

	Dhat <- res$Dhat

	# Density equals N/0.5 = 2N in the triangle region.
	true.D.triangle <- 2*N

	if(Dtrue == F) Dtrue <- mean(Dhat)
	if(Dtrue == T) Dtrue <- true.D.triangle

	erhat <- res$er
	f0hat <- res$f0hat
	t.point <- 1 - (1 - conf)/2

	# ----------------------------------------------------------------------------
	# Striplet estimator: df for the ENCOUNTER RATE component
	# ----------------------------------------------------------------------------
	# Well, this is anybody's guess!  From the simulations, using Dtrue=F,
	# the best results are by making er.df.striplet large (e.g. 1000) so that
	# effectively the Satterthwaite approximation is not being used at all.
	# So run this function using satter.striplet=F.
	# But if we DO want to use Satterthwaite then try this one, with 4 being
	# the df of the GAM used to fit the striplet method:

	er.df.striplet <- k - 4

	# ----------------------------------------------------------------------------
	# Estimators R2 and R3: df for the ENCOUNTER RATE component
	# ----------------------------------------------------------------------------
	# H = number of "strata"
	H.R <- 1
	er.df.R <- k - 1

	# ----------------------------------------------------------------------------
	# Estimators S1, S2: df for the ENCOUNTER RATE component
	# ----------------------------------------------------------------------------
	H.S <- floor(k/2)
	k.h.S <- rep(2, H.S)
	if(k %% 2 > 0)
		k.h.S[H.S] <- 3
	er.df.S <- sum(k.h.S - 1)

	# ----------------------------------------------------------------------------
	# Estimators O1, O2: df for the ENCOUNTER RATE component
	# ----------------------------------------------------------------------------
	H.O <- k - 1
	k.h.O <- rep(2, H.O)
	er.df.O <- sum(k.h.O - 1)
	# ----------------------------------------------------------------------------

	nvec <- res$n
	Ntrials <- length(Dhat)

	if(plotit){
		# Plot the lognormal distribution to see whether it fits Dhat well:
		hist(Dhat, probability=T, col=5)
		x.lo <- min(Dhat)
		x.hi <- max(Dhat)
		x.vals <- seq(x.lo-10, x.hi+10, length=100)

		# In red lines, show the density of the log-Normal distribution with mean-log equal to
		# the mean of the log of the estimated Dhats:
		lines(x.vals, dlnorm(x.vals, meanlog=mean(log(Dhat)), sdlog=sqrt(var(log(Dhat)))), col=2, lwd=2)

		# In dashed blue lines, show the density of the log-Normal distribution with mean-log
		# equal to log of the true D, rather than the mean of the estimated Dhats.
		lines(x.vals, dlnorm(x.vals, meanlog=log(true.D.triangle), sdlog=sqrt(var(log(Dhat)))),
			col=4, lwd=2, lty=2)

	}

	estimator.coverage.func <- function(est, er.t.df)
	{
		# Establish all the required quantities for this
		# estimator:
		# varD is the vector of estimated variances with this estimator.
		varD <- unlist(res[paste("varD.", est, sep = "")])
		cv.D.sq <- varD/Dhat^2

		#----------------------------------------------------------------
		# Satterthwaite method for all estimators:
		# if satter.striplet = F, deal with this case after all the rest.

		cv.f0hat.sq <- res$var.f0hat.H2/f0hat^2

		## Results for the er variance component might (or might not) have been
		## returned directly with "res".  But for convenience, calculate cv(er)
		## from the other components:
		cv.er.sq <- cv.D.sq - cv.f0hat.sq

		## Satterthwaite df on p.78 of Distance Sampling book:
		## Use cv.*.sq^2 to get cv.* ^ 4.
		satterthwaite.df.est <- cv.D.sq^2/(cv.er.sq^2/er.t.df + cv.f0hat.sq^2/(nvec - 1))
		# print(c(est=est, er.df=er.t.df, mean.sat.df=mean(satterthwaite.df.est)))
		t.df <- qt(t.point, unlist(satterthwaite.df.est))
		# print(summary(t.df))

		## Calculate the CI multiplier for this estimator (if satter.striplet=F this will be
		## overwritten in a moment):
		C.est <- exp(t.df * sqrt(log(1 + cv.D.sq)))

		## NOW, if it is the striplet estimator and satter.striplet=F, then replace C.est
		## with the striplet estimator:
		if(est=="striplet" & satter.striplet==F){
			z.alpha <- qnorm(t.point, 0, 1)
			C.est <- exp(z.alpha * sqrt(log(1 + cv.D.sq)))
		}

		lower.conf <- Dhat/C.est
		upper.conf <- Dhat * C.est
		coverage <- length((1:Ntrials)[(lower.conf <= Dtrue) & (upper.conf >= Dtrue)])

		# Uncomment these for a breakdown of incorrect CI coverage:
		# miss.from.below <- length((1:Ntrials)[upper.conf < Dtrue])
		# miss.from.above <- length((1:Ntrials)[lower.conf > Dtrue])
		# if(est=="striplet")catline(est, "cover: ", coverage, ";  miss from below: ",
		#	miss.from.below, "; above: ", miss.from.above)

		# Create summaries:
		mean.cv.D <- mean(sqrt(cv.D.sq))
		mean.cv.er <- mean(sqrt(cv.er.sq))
		mean.cv.f0hat <- mean(sqrt(cv.f0hat.sq))
		mean.ci.width <- mean(upper.conf - lower.conf)
		c(cover=coverage/Ntrials, mean.cv.er=mean.cv.er, mean.cv.f0hat,
                                                    mean.cv.D = mean.cv.D, mean.ci.width=mean.ci.width,
			mean.D=mean(Dhat), true.D = true.D.triangle)
	}
	estnames.striplet <- c("striplet")
	estnames.R <- c("R2", "R3")
	estnames.S <- c("S1", "S2")
	estnames.O <- c("O1", "O2")
	coverage.striplet <- sapply(estnames.striplet, estimator.coverage.func, er.t.df = er.df.striplet)
	coverage.R <- sapply(estnames.R, estimator.coverage.func, er.t.df = er.df.R)
	coverage.S <- sapply(estnames.S, estimator.coverage.func, er.t.df = er.df.S)
	coverage.O <- sapply(estnames.O, estimator.coverage.func, er.t.df = er.df.O)
	##
	results <- cbind(coverage.striplet, coverage.S, coverage.O, coverage.R)
	results <- as.data.frame(results)
	names(results) <- c(estnames.striplet, estnames.S, estnames.O, estnames.R)
	rownames(results) <- c("ci.cover", "mean.cv.er", "mean.cv.f0hat", "mean.cv.D",
                               "mean.ci.width", "mean.D", "true.D")
	results
}


#########################################################################
triangle.oneres.boxplot <- function(res, sd.plot=T, Dtrue=F, satter.striplet=F, paperplot=F)
{
	# triangle.oneres.boxplot 13/7/09
	# Boxplot for a single simulation result, "res", output from triangle.Dhat.func.
	# The mean estimated CV is printed over each boxplot, and the CI coverage is printed under
	# each boxplot.
	#
                 # If paperplot=T, customizes to be plotted in the paper via paper.triangle.plot().  In particular,
                 # it only plots 1000 outcomes in the boxplots, although means, CI covers, CVs, shown on the plot
                 # are all computed on the basis of the full data.
                 #
                 # Dtrue is an argument to triangle.Dhat.confintcover.func, to say whether the "true" value of D
	# is to be taken as the actual density of objects in the region (Dtrue = T), or the mean estimated
	# density of objects in the region (Dtrue=F).  The latter could be seen as the "true" value of D
	# under the repeated survey framework.
	#
	# The triangle aspect of things is just in calling triangle.Dhat.confintcover.func
	# and triangle.mean.cover: otherwise this function is general.
	#
	# EXAMPLE:
	# tryme.res <- triangle.Dhat.func(c(w=0.02, theta=0.01, k=20, N=1000, superpop=2), Nsurv=100)
	# triangle.oneres.boxplot(tryme.res)
	#

	parvec <- attributes(res)$inputs
	w <- parvec["w"]
	theta <- parvec["theta"]
	k <- parvec["k"]
	N <- parvec["N"]
	superpop <- parvec["superpop"]
	nbar <- round(mean(res$n))
	triangle.mean.cover <- round(mean(2*w*res$L/0.5)*100)

	parname <- paste("cover=", triangle.mean.cover, "%, nbar=", nbar,
			", k=", k, ", N=", N, ", superpop=", superpop, sep="")

	trueval <- var(res$Dhat)
	comps <- list(res$varD.striplet, res$varD.S1, res$varD.S2, res$varD.O1, res$varD.O2)
                 if(paperplot==T) compsplot <- lapply(comps, function(cpt)cpt[1:1000])
                 else compsplot <- comps

	if(sd.plot == T) {
		comps <- lapply(comps, sqrt)
		compsplot <- lapply(compsplot, sqrt)
		trueval <- sqrt(trueval)
	}

	y.hi <- max(c(unlist(compsplot), trueval)) + 0.2*diff(range(c(unlist(compsplot), trueval)))
	y.lo <- min(c(unlist(compsplot), trueval)) - 0.2*diff(range(c(unlist(compsplot), trueval)))
	yrange <- y.hi - y.lo

	# Use help(bxp) to find some of the extra arguments available to enter
	# in pars=list() argument to boxplot.

                 bp <- boxplot(compsplot, names = c("striplet", "S1", "S2", "O1", "O2"), ylim = c(y.lo, y.hi),
		pars=list(medlty="blank", boxfill="gray82", outwex = 0.3, boxwex = 0.6,
		outpch=NA, outlty="solid", staplewex=0.6, staplelwd=1, whisklty=1, whisklwd=0.2,
		outlwd=0.2), cex.axis=1.2)

	abline(h = trueval, col = 1, lwd = 3)
	for(i in 1:length(comps))
		lines(c(i - 0.4, i + 0.4), rep(mean(comps[[i]]), 2), lwd = 2)

	#-----------------------------------------------------
	# Add CI coverage for nominal 95% intervals, and mean %CV, to the plot:
	confres <- triangle.Dhat.confintcover.func(res=res, conf=0.95, Dtrue=Dtrue, satter.striplet=satter.striplet)
	print(confres)

	cic.vec <- confres["ci.cover",]
	cic.vec <- cic.vec[c("striplet", "S1", "S2", "O1", "O2")]

	print(cic.vec)
	if(any(names(cic.vec)!=c("striplet", "S1", "S2", "O1", "O2")))
		stop("Make sure cic.vec is ordered in the same order as the boxplots!!")


	# Make into integer percentage:
	cic.pc <- round(100 * cic.vec)
	for(i in 1:length(comps)) text(i, y.lo + (y.hi - y.lo)/15, cic.pc[i], cex = 1.3)

	# Similarly for mean CI width:
	ciwidth.vec <- confres["mean.ci.width",]
	ciwidth.vec <- ciwidth.vec[c("striplet", "S1", "S2", "O1", "O2")]
	ciwidth.round <- round(ciwidth.vec, -1)
	for(i in 1:length(comps)) text(i, y.hi - (y.hi - y.lo)/15, ciwidth.round[i], cex = 1.3)
	#-----------------------------------------------------

	if(!paperplot)title(main=parname, cex.main=1.5)

}

