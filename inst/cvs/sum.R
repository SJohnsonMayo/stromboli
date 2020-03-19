#n <- 50; p <- 30; rho <- 0.2
fam <- "ln"
load(paste("data/new", n, "_", p, "_", rho, fam, ".rda", sep=""))
load(paste("result/res", n, "_", p, "_", rho, fam, ".rda", sep=""))
nsim <- length(result)

pe.gic <- numeric(nsim)
l1.gic <- numeric(nsim); l2.gic <- numeric(nsim); linf.gic <- numeric(nsim)
fp.gic <- numeric(nsim); fn.gic <- numeric(nsim)
pe.rft <- numeric(nsim)
l1.rft <- numeric(nsim); l2.rft <- numeric(nsim); linf.rft <- numeric(nsim)
pe1.gic <- numeric(nsim)
l11.gic <- numeric(nsim); l21.gic <- numeric(nsim); linf1.gic <- numeric(nsim)
fp1.gic <- numeric(nsim); fn1.gic <- numeric(nsim)
pe2.gic <- numeric(nsim)
l12.gic <- numeric(nsim); l22.gic <- numeric(nsim); linf2.gic <- numeric(nsim)
fp2.gic <- numeric(nsim); fn2.gic <- numeric(nsim)
pe2.rft <- numeric(nsim)
l12.rft <- numeric(nsim); l22.rft <- numeric(nsim); linf2.rft <- numeric(nsim)

#pe.cv <- numeric(nsim)
#l1.cv <- numeric(nsim); l2.cv <- numeric(nsim); linf.cv <- numeric(nsim)
#fp.cv <- numeric(nsim); fn.cv <- numeric(nsim)
#pe1.cv <- numeric(nsim)
#l11.cv <- numeric(nsim); l21.cv <- numeric(nsim); linf1.cv <- numeric(nsim)
#fp1.cv <- numeric(nsim); fn1.cv <- numeric(nsim)
#pe2.cv <- numeric(nsim)
#l12.cv <- numeric(nsim); l22.cv <- numeric(nsim); linf2.cv <- numeric(nsim)
#fp2.cv <- numeric(nsim); fn2.cv <- numeric(nsim)

for (i in 1:nsim) {
	# Proposed, GIC
	pe.gic[i] <- sum((data[[i]]$y - result[[i]]$int.gic - data[[i]]$z %*% result[[i]]$bet.gic)^2)/n
	del <- result[[i]]$bet.gic - data[[i]]$bet
	l1.gic[i] <- sum(abs(del)); l2.gic[i] <- sum(del^2); linf.gic[i] <- max(abs(del))
	fp.gic[i] <- sum(result[[i]]$bet.gic != 0 & data[[i]]$bet == 0)
	fn.gic[i] <- sum(result[[i]]$bet.gic == 0 & data[[i]]$bet != 0)
	# Proposed + refit
	pe.rft[i] <- sum((data[[i]]$y - result[[i]]$int.rft - data[[i]]$z %*% result[[i]]$bet.rft)^2)/n
	del <- result[[i]]$bet.rft - data[[i]]$bet
	l1.rft[i] <- sum(abs(del)); l2.rft[i] <- sum(del^2); linf.rft[i] <- max(abs(del))
	# Lasso (i)
	pe1.gic[i] <- sum((data[[i]]$y - result[[i]]$int1.gic - data[[i]]$z %*% result[[i]]$bet1.gic)^2)/n
	del <- result[[i]]$bet1.gic - data[[i]]$bet
	l11.gic[i] <- sum(abs(del)); l21.gic[i] <- sum(del^2); linf1.gic[i] <- max(abs(del))
	fp1.gic[i] <- sum(result[[i]]$bet1.gic != 0 & data[[i]]$bet == 0)
	fn1.gic[i] <- sum(result[[i]]$bet1.gic == 0 & data[[i]]$bet != 0)
	# Lasso (ii)
	pe2.gic[i] <- sum((data[[i]]$y - result[[i]]$int2.gic - data[[i]]$z %*% result[[i]]$bet2.gic)^2)/n
	del <- result[[i]]$bet2.gic - data[[i]]$bet
	l12.gic[i] <- sum(abs(del)); l22.gic[i] <- sum(del^2); linf2.gic[i] <- max(abs(del))
	fp2.gic[i] <- sum(result[[i]]$bet2.gic != 0 & data[[i]]$bet == 0)
	fn2.gic[i] <- sum(result[[i]]$bet2.gic == 0 & data[[i]]$bet != 0)
	# Lasso (ii) + refit
	pe2.rft[i] <- sum((data[[i]]$y - result[[i]]$int2.rft - data[[i]]$z %*% result[[i]]$bet2.rft)^2)/n
	del <- result[[i]]$bet2.rft - data[[i]]$bet
	l12.rft[i] <- sum(abs(del)); l22.rft[i] <- sum(del^2); linf2.rft[i] <- max(abs(del))
	## Proposed, CV
	#pe.cv[i] <- sum((data[[i]]$y - result[[i]]$int.cv - data[[i]]$z %*% result[[i]]$bet.cv)^2)/n
	#del <- result[[i]]$bet.cv - data[[i]]$bet
	#l1.cv[i] <- sum(abs(del)); l2.cv[i] <- sum(del^2); linf.cv[i] <- max(abs(del))
	#fp.cv[i] <- sum(result[[i]]$bet.cv != 0 & data[[i]]$bet == 0)
	#fn.cv[i] <- sum(result[[i]]$bet.cv == 0 & data[[i]]$bet != 0)
	## Lasso (i)
	#pe1.cv[i] <- sum((data[[i]]$y - result[[i]]$int1.cv - data[[i]]$z %*% result[[i]]$bet1.cv)^2)/n
	#del <- result[[i]]$bet1.cv - data[[i]]$bet
	#l11.cv[i] <- sum(abs(del)); l21.cv[i] <- sum(del^2); linf1.cv[i] <- max(abs(del))
	#fp1.cv[i] <- sum(result[[i]]$bet1.cv != 0 & data[[i]]$bet == 0)
	#fn1.cv[i] <- sum(result[[i]]$bet1.cv == 0 & data[[i]]$bet != 0)
	## Lasso (ii)
	#pe2.cv[i] <- sum((data[[i]]$y - result[[i]]$int2.cv - data[[i]]$z %*% result[[i]]$bet2.cv)^2)/n
	#del <- result[[i]]$bet2.cv - data[[i]]$bet
	#l12.cv[i] <- sum(abs(del)); l22.cv[i] <- sum(del^2); linf2.cv[i] <- max(abs(del))
	#fp2.cv[i] <- sum(result[[i]]$bet2.cv != 0 & data[[i]]$bet == 0)
	#fn2.cv[i] <- sum(result[[i]]$bet2.cv == 0 & data[[i]]$bet != 0)
}

cat("Lasso (i), GIC\n")
cat("PE   ", mean(pe1.gic), " (", sd(pe1.gic)/sqrt(nsim), ")\n", sep="")
cat("L1   ", mean(l11.gic), " (", sd(l11.gic)/sqrt(nsim), ")\n", sep="")
cat("L2   ", mean(l21.gic), " (", sd(l21.gic)/sqrt(nsim), ")\n", sep="")
cat("Linf ", mean(linf1.gic), " (", sd(linf1.gic)/sqrt(nsim), ")\n", sep="")
cat("FP   ", mean(fp1.gic), " (", sd(fp1.gic)/sqrt(nsim), ")\n", sep="")
cat("FN   ", mean(fn1.gic), " (", sd(fn1.gic)/sqrt(nsim), ")\n", sep="")
cat("Lasso (ii)\n")
cat("PE   ", mean(pe2.gic), " (", sd(pe2.gic)/sqrt(nsim), ")\n", sep="")
cat("L1   ", mean(l12.gic), " (", sd(l12.gic)/sqrt(nsim), ")\n", sep="")
cat("L2   ", mean(l22.gic), " (", sd(l22.gic)/sqrt(nsim), ")\n", sep="")
cat("Linf ", mean(linf2.gic), " (", sd(linf2.gic)/sqrt(nsim), ")\n", sep="")
cat("FP   ", mean(fp2.gic), " (", sd(fp2.gic)/sqrt(nsim), ")\n", sep="")
cat("FN   ", mean(fn2.gic), " (", sd(fn2.gic)/sqrt(nsim), ")\n", sep="")
cat("Proposed\n")
cat("PE   ", mean(pe.gic), " (", sd(pe.gic)/sqrt(nsim), ")\n", sep="")
cat("L1   ", mean(l1.gic), " (", sd(l1.gic)/sqrt(nsim), ")\n", sep="")
cat("L2   ", mean(l2.gic), " (", sd(l2.gic)/sqrt(nsim), ")\n", sep="")
cat("Linf ", mean(linf.gic), " (", sd(linf.gic)/sqrt(nsim), ")\n", sep="")
cat("FP   ", mean(fp.gic), " (", sd(fp.gic)/sqrt(nsim), ")\n", sep="")
cat("FN   ", mean(fn.gic), " (", sd(fn.gic)/sqrt(nsim), ")\n", sep="")
cat("Lasso (ii) + refit\n")
cat("PE   ", mean(pe2.rft), " (", sd(pe2.rft)/sqrt(nsim), ")\n", sep="")
cat("L1   ", mean(l12.rft), " (", sd(l12.rft)/sqrt(nsim), ")\n", sep="")
cat("L2   ", mean(l22.rft), " (", sd(l22.rft)/sqrt(nsim), ")\n", sep="")
cat("Linf ", mean(linf2.rft), " (", sd(linf2.rft)/sqrt(nsim), ")\n", sep="")
cat("Proposed + refit\n")
cat("PE   ", mean(pe.rft), " (", sd(pe.rft)/sqrt(nsim), ")\n", sep="")
cat("L1   ", mean(l1.rft), " (", sd(l1.rft)/sqrt(nsim), ")\n", sep="")
cat("L2   ", mean(l2.rft), " (", sd(l2.rft)/sqrt(nsim), ")\n", sep="")
cat("Linf ", mean(linf.rft), " (", sd(linf.rft)/sqrt(nsim), ")\n", sep="")

#cat("Lasso (i), CV\n")
#cat("PE   ", mean(pe1.cv), " (", sd(pe1.cv)/sqrt(nsim), ")\n", sep="")
#cat("L1   ", mean(l11.cv), " (", sd(l11.cv)/sqrt(nsim), ")\n", sep="")
#cat("L2   ", mean(l21.cv), " (", sd(l21.cv)/sqrt(nsim), ")\n", sep="")
#cat("Linf ", mean(linf1.cv), " (", sd(linf1.cv)/sqrt(nsim), ")\n", sep="")
#cat("FP   ", mean(fp1.cv), " (", sd(fp1.cv)/sqrt(nsim), ")\n", sep="")
#cat("FN   ", mean(fn1.cv), " (", sd(fn1.cv)/sqrt(nsim), ")\n", sep="")
#cat("Lasso (ii)\n")
#cat("PE   ", mean(pe2.cv), " (", sd(pe2.cv)/sqrt(nsim), ")\n", sep="")
#cat("L1   ", mean(l12.cv), " (", sd(l12.cv)/sqrt(nsim), ")\n", sep="")
#cat("L2   ", mean(l22.cv), " (", sd(l22.cv)/sqrt(nsim), ")\n", sep="")
#cat("Linf ", mean(linf2.cv), " (", sd(linf2.cv)/sqrt(nsim), ")\n", sep="")
#cat("FP   ", mean(fp2.cv), " (", sd(fp2.cv)/sqrt(nsim), ")\n", sep="")
#cat("FN   ", mean(fn2.cv), " (", sd(fn2.cv)/sqrt(nsim), ")\n", sep="")
#cat("Proposed\n")
#cat("PE   ", mean(pe.cv), " (", sd(pe.cv)/sqrt(nsim), ")\n", sep="")
#cat("L1   ", mean(l1.cv), " (", sd(l1.cv)/sqrt(nsim), ")\n", sep="")
#cat("L2   ", mean(l2.cv), " (", sd(l2.cv)/sqrt(nsim), ")\n", sep="")
#cat("Linf ", mean(linf.cv), " (", sd(linf.cv)/sqrt(nsim), ")\n", sep="")
#cat("FP   ", mean(fp.cv), " (", sd(fp.cv)/sqrt(nsim), ")\n", sep="")
#cat("FN   ", mean(fn.cv), " (", sd(fn.cv)/sqrt(nsim), ")\n", sep="")
