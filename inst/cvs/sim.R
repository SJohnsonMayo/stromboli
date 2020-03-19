sim.data <- function(n=50, p=30, nsim=100, fam="ln", rho=0.2, seed=21) {
	if (fam=="ln") require(mvtnorm)
	if (fam=="dir") require(gtools)
	set.seed(seed)
	bet <- c(1, -0.8, 0.6, 0, 0, -1.5, -0.5, 1.2, rep(0, p - 8))
	alp <- c(rep(p/2, 5), rep(1, p - 5))
	sig <- matrix(, p, p)
	sig <- rho^abs(row(sig) - col(sig))
	data <- vector("list", nsim)
	for (i in 1:nsim) {
		x <- switch(fam, ln={x <- exp(rmvnorm(n, log(alp), sig)); x/rowSums(x)}, dir=rdirichlet(n, alp))
		z <- log(x)
		y <- rnorm(n, z %*% bet, 0.5)
		data[[i]] <- list(y=y, z=z, bet=bet)
		cat("Sim.", i, "done.\n")
	}
	save(data, file=paste("data/dat", n, "_", p, "_", rho, fam, ".rda", sep=""))
	data <- vector("list", nsim)
	for (i in 1:nsim) {
		x <- switch(fam, ln={x <- exp(rmvnorm(n, log(alp), sig)); x/rowSums(x)}, dir=rdirichlet(n, alp))
		z <- log(x)
		y <- rnorm(n, z %*% bet, 0.5)
		data[[i]] <- list(y=y, z=z, bet=bet)
		cat("Sim.", i, "done.\n")
	}
	save(data, file=paste("data/new", n, "_", p, "_", rho, fam, ".rda", sep=""))
}

sim.comp <- function(n=50, p=30, nsim=100, fam="ln", rho=0.2, seed=42) {
	set.seed(seed)
	load(paste("data/dat", n, "_", p, "_", rho, fam, ".rda", sep=""))
	result <- vector("list", nsim)
	for (i in 1:nsim) {
		# Proposed, GIC
		ans <- gic.cdmm(data[[i]]$y, data[[i]]$z, type="ft")
		result[[i]]$bet.gic <- ans$bet; result[[i]]$lam.gic <- ans$lam
		result[[i]]$int.gic <- ans$int
		# Proposed + refit
		isel <- result[[i]]$bet.gic != 0
		result[[i]]$bet.rft <- result[[i]]$bet.gic
		result[[i]]$int.rft <- result[[i]]$int.gic
		if (any(isel)) {
			ans <- cdmm(data[[i]]$y, as.matrix(data[[i]]$z[, isel]), 0)
			result[[i]]$bet.rft[isel] <- as.numeric(ans$sol)
			result[[i]]$int.rft <- ans$int
		}
		# Lasso (i)
		iref <- sample(1:p, 1)
		zd <- data[[i]]$z[, -iref] - data[[i]]$z[, iref]
		ans <- gic.cdmm(data[[i]]$y, zd, type="ft", constr=FALSE)
		result[[i]]$bet1.gic <- -sum(ans$bet)
		if (iref > 1) result[[i]]$bet1.gic <- c(ans$bet[1:(iref - 1)], result[[i]]$bet1.gic)
		if (iref < p) result[[i]]$bet1.gic <- c(result[[i]]$bet1.gic, ans$bet[iref:(p - 1)])
		result[[i]]$lam1.gic <- ans$lam; result[[i]]$int1.gic <- ans$int
		result[[i]]$iref <- iref
		# Lasso (ii)
		ans <- gic.cdmm(data[[i]]$y, data[[i]]$z, type="ft", constr=FALSE)
		result[[i]]$bet2.gic <- ans$bet; result[[i]]$lam2.gic <- ans$lam
		result[[i]]$int2.gic <- ans$int
		# Lasso (ii) + refit
		isel <- result[[i]]$bet2.gic != 0
		result[[i]]$bet2.rft <- result[[i]]$bet2.gic
		result[[i]]$int2.rft <- result[[i]]$int2.gic
		if (any(isel)) {
			ans <- cdmm(data[[i]]$y, as.matrix(data[[i]]$z[, isel]), 0)
			result[[i]]$bet2.rft[isel] <- as.numeric(ans$sol)
			result[[i]]$int2.rft <- ans$int
		}
		## Proposed, CV
		#ans <- cv.cdmm(data[[i]]$y, data[[i]]$z, type="1se")
		#result[[i]]$bet.cv <- ans$bet; result[[i]]$lam.cv <- ans$lam
		#result[[i]]$int.cv <- ans$int; result[[i]]$foldid <- ans$foldid
		## Lasso (i)
		#iref <- sample(1:p, 1)
		#zd <- data[[i]]$z[, -iref] - data[[i]]$z[, iref]
		#ans <- cv.cdmm(data[[i]]$y, zd, foldid=result[[i]]$foldid, type="1se", constr=FALSE)
		#result[[i]]$bet1.cv <- -sum(ans$bet)
		#if (iref > 1) result[[i]]$bet1.cv <- c(ans$bet[1:(iref - 1)], result[[i]]$bet1.cv)
		#if (iref < p) result[[i]]$bet1.cv <- c(result[[i]]$bet1.cv, ans$bet[iref:(p - 1)])
		#result[[i]]$lam1.cv <- ans$lam; result[[i]]$int1.cv <- ans$int
		## Lasso (ii)
		#ans <- cv.cdmm(data[[i]]$y, data[[i]]$z, foldid=result[[i]]$foldid, type="1se", constr=FALSE)
		#result[[i]]$bet2.cv <- ans$bet; result[[i]]$lam2.cv <- ans$lam
		#result[[i]]$int2.cv <- ans$int
		cat("Sim.", i, "done.\n")
	}
	save(result, file=paste("result/res", n, "_", p, "_", rho, fam, ".rda", sep=""))
}
