
library(ape)
library(phylowise)
library(nlme)


args = commandArgs(trailingOnly = TRUE)

out.file=args[1]

if (length(args) > 1){
	SKIP.PGLS=as.logical(args[2])
}else{
	SKIP.PGLS=FALSE
}

if (length(args) > 2){
	SKIP.PPC=as.logical(args[3])
}else{
	SKIP.PPC=FALSE
}





cat(paste("saving to", out.file, "\n"))
cat(paste("SKIP.PGLS =", SKIP.PGLS, "\n"))
cat(paste("SKIP.PPC =", SKIP.PPC, "\n"))



# Load a random tree out of the folder
sample.newick = function(dir){
	files = list.files(dir, pattern=".newick")
	f = sample(files, 1)
	tree = read.tree(paste0(dir, "/", f))
	tree$tip.label = paste0("taxon", tree$tip.label)
	tree
}



ntrials=100000



# Expected number of substitutions per unit of time (tree has height of 1)
clock.models = c("TD")
clock.rate.mean=2000
clock.rate.sd=2




max.wall.time = 10


# False discovery rate
trait.methods = c("BM")
tree.dirs = c("../newick/")



youngest=FALSE
nested=TRUE
time.tree.cond=TRUE
n.ppc.trials = 100




setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits
out.df = data.frame(sampleNr = numeric(0), tree.dir = character(0), clock.sigma=numeric(0), ntaxa=numeric(0), clock.model = character(0), beta=numeric(0), trait.method = character(0), nested=logical(0), youngest=logical(0),
	maximise=logical(0), logY=logical(0), standardise=numeric(0), number.of.subst=numeric(0), dmin=numeric(0), dspan=numeric(0), time.tree.cond=logical(0), pgls.ou.p=numeric(0), pgls.rho.p=numeric(0), pgls.blomberg.p=numeric(0), pic.p.st=numeric(0), pic.p.tt=numeric(0), pgls.lambda.p=numeric(0), ls.p=numeric(0), pgls.p.st=numeric(0), pgls.p.tt=numeric(0), 
	ppc.p1=numeric(0), bma.prob1=numeric(0), bma.prob.all=numeric(0), rho = numeric(0), n=numeric(0), t=numeric(0), trait.noise=character(0), rate.noise=character(0))

for (trial in 1:ntrials){


	if (trial %% 1 == 0){
		cat(paste(trial, "/", ntrials, "\n"))
		write.table(out.df, out.file, sep="\t", quote=F, row.names=F)
	}

	for (tree.dir in tree.dirs){


		# Simulate a time tree with height up to 1. Extinct taxa removed
		time.tree = sample.newick(dir=tree.dir)

		
		

		# Drop some at random so we can see the power when n<30 (the smallest tree size)
		if (runif(1,0,1) < 0.01){
			keep = rpois(1, 7) + 3 # Every few trials can keep just 10 or so tips to help with the regression
			nn = length(time.tree$tip.label)
			if (keep > nn){
				keep = nn
			}

			s = sample(1:nn, size=keep)
			drop = setdiff(1:nn, s)

			cat(paste("dropping", length(drop), "to get", keep, "\n"))

			if (length(drop) > 0){
				time.tree = drop.tip(time.tree, drop)
			}
		}

		dmat.time = phylowise::getDistanceMatrix(time.tree)


		# Simulate a subst tree
		for (clock.model in clock.models){


			# These parameters give a clock sd of 0.49 on average
			theta = rexp(1, 0.1)
			V = 0.1
			beta = rexp(1, 2)
			clock.sigma = sqrt(V*2*theta) 


			if (clock.model == "AC"){
				beta = 0
			}

			

			for (trait.method in trait.methods){

				# Simulate traits under this method
				traits = phylowise::simulateTrait(time.tree, sigma=1)
				traits = scale(traits)[,1]


					number.of.subst = rlnorm(1, log(clock.rate.mean) - 0.5*clock.rate.sd^2, clock.rate.sd)
					result = phylowise::simulateSubstitutions(time.tree, traits=traits, beta=beta, sigma=clock.sigma, theta=theta, number.of.subst=number.of.subst, method="TD")
					subst.tree.true = result$subst.tree.true
					subst.tree = result$subst.tree.est
					node.rates.est = result$node.rates.est
					node.rates.true = result$node.rates.true
				 	number.of.subst.total = sum(subst.tree$edge.length)



				for (rate.noise in c("poisson")){


					#for (noise in c("none", "rounding", "gaussian")){
					for (noise in c("none")){						
					#for (noise in c("none")){

						traits.noise = traits

						if (noise == "rounding"){
							traits.noise = round(exp(traits.noise), 1)
							traits.noise = ifelse(traits.noise == 0, 0.01, traits.noise)
							traits.noise = log(traits.noise)
						}

						if (noise == "gaussian"){
							traits.noise = traits.noise + rnorm(length(traits.noise), 0, 0.2)
						}



						nn = length(time.tree$tip.label)
						leaf.rates.est = node.rates.est[1:nn]
						leaf.rates.true = node.rates.true[1:nn]
						leaf.traits = traits.noise[1:nn]

						leaf.rates.est = log(leaf.rates.est) 
						leaf.rates.true = log(leaf.rates.true)



						for (regression.rates.est in c(TRUE)){


							if (regression.rates.est){
								leaf.rates = leaf.rates.est
							}else{
								leaf.rates = leaf.rates.true
							}




							drop = which(leaf.rates == -Inf | is.na(leaf.rates))
							subst.tree.nonzero.rates = subst.tree
							time.tree.nonzero.rates = time.tree
							leaf.traits.nonzero = leaf.traits



							
							

							if (length(drop) > 0){
								subst.tree.nonzero.rates = drop.tip(subst.tree, drop)
								time.tree.nonzero.rates = drop.tip(time.tree, drop)
								leaf.rates = leaf.rates[-drop]
								leaf.traits.nonzero = leaf.traits[-drop]
							}


							


							ntaxa = length(leaf.rates)




							if (!SKIP.PGLS){

								if (length(leaf.rates) < 3){
									next
								}
								

								ls = cor.test(leaf.traits.nonzero, leaf.rates)
								ls.p = ls$p.value
								#print(paste("ls.p=", signif(ls.p, 2), " for R2=", signif(R2.real, 2), " ntaxa = ", ntaxa))



								# PGLS on subst tree 
								pgls.p.st = tryCatch({
									pgls = gls(leaf.rates~leaf.traits.nonzero, correlation=corBrownian(1,subst.tree.nonzero.rates), method="ML")
									pgls.result = summary(pgls)
									
									# Return a p-value of 1 if the slope is negative
									slope = pgls.result$coefficients[2]
									ifelse(slope <= -Inf, 1, pgls.result$tTable[2,4])
								}, error = function(e) {
								  	NA
								})




						

								# PGLS on time tree 
								pgls.p.tt = tryCatch({
									pgls = gls(leaf.rates~leaf.traits.nonzero, correlation=corBrownian(1,time.tree.nonzero.rates), method="ML")
									pgls.result = summary(pgls)
									
									# Return a p-value of 1 if the slope is negative
									slope = pgls.result$coefficients[2]
									ifelse(slope <= -Inf, 1, pgls.result$tTable[2,4])
								}, error = function(e) {
								  	NA
								})



								# PGLS on time tree lambda with multiple starting points
								pgls.lambda.p=NA
								for(initVal in seq(from=1, to=0.1,by=-0.1)){

									setTimeLimit(elapsed = max.wall.time)   # Stop after a few seconds of wall-clock time
									pgls.lambda.p = tryCatch({
										pgls = gls(leaf.rates~leaf.traits.nonzero, correlation=corPagel(initVal, time.tree.nonzero.rates,fixed = FALSE), method="ML")
										pgls.result = summary(pgls)
										
										# Return a p-value of 1 if the slope is negative
										slope = pgls.result$coefficients[2]
										ifelse(slope <= -Inf, 1, pgls.result$tTable[2,4])
									}, error = function(e) {
									  	NA
									})
									setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits

									if (!is.na(pgls.lambda.p)){
										break
									}


								}

								if (is.na(pgls.lambda.p)){
									print("lambda did not converge")
								}
								




								# OU process
								pgls.ou.p=NA
								for(initVal in c(0.01, 0.1, 0.5, 1.0, 2, 5, 10, 50, 100, 1000, 10000)){

									setTimeLimit(elapsed = max.wall.time)   # Stop after a few seconds of wall-clock time
									pgls.ou.p = tryCatch({
										pgls = gls(leaf.rates~leaf.traits.nonzero, correlation=corMartins(initVal, time.tree.nonzero.rates,fixed = FALSE), method="ML")
										pgls.result = summary(pgls)
										
										# Return a p-value of 1 if the slope is negative
										slope = pgls.result$coefficients[2]
										ifelse(slope <= -Inf, 1, pgls.result$tTable[2,4])
									}, error = function(e) {
									  	NA
									})
									setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits

									if (!is.na(pgls.ou.p)){
										break
									}


								}

								if (is.na(pgls.ou.p)){
									print("ou did not converge")
								}





								# Grafen rho process
								pgls.rho.p=NA
								for(initVal in c(0.01, 0.1, 0.5, 1.0, 2, 5, 10, 50, 100)){

									setTimeLimit(elapsed = max.wall.time)   # Stop after a few seconds of wall-clock time
									pgls.rho.p = tryCatch({
										pgls = gls(leaf.rates~leaf.traits.nonzero, correlation=corGrafen(initVal, time.tree.nonzero.rates,fixed = FALSE), method="ML")
										pgls.result = summary(pgls)
										
										# Return a p-value of 1 if the slope is negative
										slope = pgls.result$coefficients[2]
										ifelse(slope <= -Inf, 1, pgls.result$tTable[2,4])
									}, error = function(e) {
									  	NA
									})
									setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits

									if (!is.na(pgls.rho.p)){
										break
									}


								}

								if (is.na(pgls.rho.p)){
									print("rho did not converge")
								}



								# Blomberg rho process
								pgls.blomberg.p=NA
								for(initVal in c(0.0001, 0.01, 0.1, 0.5, 1.0, 2, 5, 10, 50, 100, 1000, 10000)){

									setTimeLimit(elapsed = max.wall.time)   # Stop after a few seconds of wall-clock time
									pgls.blomberg.p = tryCatch({
										pgls = gls(leaf.rates~leaf.traits.nonzero, correlation=corBlomberg(initVal, time.tree.nonzero.rates,fixed = FALSE), method="ML")
										pgls.result = summary(pgls)

										# Return a p-value of 1 if the slope is negative
										slope = pgls.result$coefficients[2]
										ifelse(slope <= -Inf, 1, pgls.result$tTable[2,4])
									}, error = function(e) {
									  	NA
									})
									setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits

									if (!is.na(pgls.blomberg.p)){
										break
									}


								}

								if (is.na(pgls.blomberg.p)){
									print("blomberg did not converge")
								}



								# Subst tree PIC
								# pic1 = pic(leaf.rates, subst.tree.nonzero.rates)
								# pic2 = pic(leaf.traits.nonzero, subst.tree.nonzero.rates)
								# pic_model = lm(pic1 ~ pic2 - 1)
								# pic.res = summary(pic_model)
								# pic.p.st = as.numeric(pic.res$coefficients[1,4])
								pic.p.st = NA


								# Time tree PIC
								pic1 = pic(leaf.rates, time.tree.nonzero.rates)
								pic2 = pic(leaf.traits.nonzero, time.tree.nonzero.rates)
								pic_model = lm(pic1 ~ pic2 - 1)
								pic.res = summary(pic_model)
								pic.p.tt = as.numeric(pic.res$coefficients[1,4])


								
								#pgls.lambda = gls(leaf.rates~leaf.traits, correlation=corPagel(1,subst.tree,fixed = FALSE), method="ML")


								out.df2 = data.frame(sampleNr = trial, tree.dir = tree.dir, clock.sigma=clock.sigma, ntaxa=ntaxa, youngest=NA, number.of.subst=number.of.subst.total, clock.model = clock.model, 
									trait.method = trait.method, nested=NA, maximise=NA, 
										logY=TRUE, standardise=NA, dmin=NA, dspan=NA, beta=beta, time.tree.cond=NA, 
										pic.p.st=pic.p.st, pic.p.tt=pic.p.tt, pgls.lambda.p=pgls.lambda.p, ls.p=ls.p, pgls.ou.p=pgls.ou.p, pgls.rho.p=pgls.rho.p, pgls.blomberg.p=pgls.blomberg.p,
										pgls.p.tt=pgls.p.tt, pgls.p.st=pgls.p.st, ppc.p1=NA,  bma.prob1=NA, bma.prob.all=NA, rho=NA, n=length(leaf.rates), t=NA, rate.noise=rate.noise, trait.noise=noise)
								

								out.df = rbind(out.df, out.df2)


								sub.df = out.df[out.df$R2.real < 0.1 & !is.na(out.df$ls.p),]
								if (noise == "none" && nrow(sub.df) > 1){
									#print(paste("OLS P(p<0.01)=", signif(sum(sub.df$ls.p < 0.01) / nrow(sub.df), 3), "for n=", nrow(sub.df)))
									#hist(sub.df$ls.p)
								}

							}

						}



						if (SKIP.PPC){
							next
						}

						for (youngest in c(FALSE)){
							for (maximise in c(FALSE, TRUE)){
								
								dmin = 0.05# runif(1, 0, 0.5)
								dspan = 10
								dmax = dmin * dspan


								pvals = numeric(0)
								nvals = numeric(0)
								rhovals = numeric(0)

								bma.probs = numeric(0)

								ntaxa = length(subst.tree$tip.label)

								for (t1 in 1:n.ppc.trials){


									# Sample taxa - time tree
									setTimeLimit(elapsed = max.wall.time)   # Stop after a few seconds of wall-clock time
									pairs.df = tryCatch({
										phylowise::sampleTaxonPairs(subst.tree=subst.tree, window.tree=time.tree, covariate=traits.noise, youngest=youngest, dist.min=dmin, dist.max=dmax, nested=nested, maximise=maximise, distance.matrix=dmat.time, verbose=F)
									}, error = function(e) {
										print(e)
										setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits
										print(paste("Timeout 3", dmin, dmax, nested, maximise, time.tree.cond))
									  	data.frame(x=1)
									})
									setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits
									
									#plot(pairs.df$trait1 - pairs.df$trait2, pairs.df$d1 - pairs.df$d2)


									d.all = c(pairs.df$d1, pairs.df$d2)
									if (any(is.na(d.all))){
										cat(paste0("there are some NA values\n"))
										print(pairs.df)
										next
									}



									if (any(abs(d.all) == Inf)) {
										cat(paste0("there are some Inf values\n"))
										print(pairs.df)
										next
									}

									if (nrow(pairs.df) < 5){
										next
									}


									# PPC test
									test.results = phylowise::PPC.test(pairs.df, standardise=2, logY=TRUE, extreme.value.threshold=0)

									if (is.na(test.results$p)){
										next
									}

									pvals = c(pvals, test.results$p)
									nvals = c(nvals, test.results$n)
									rhovals = c(rhovals, test.results$rho)

									dat = test.results$data.df

									ppc.prob = as.numeric(test.results$bma.probs["trait"])
									ppc.bf = as.numeric(test.results$bma.bf["trait"])
									
									bma.probs = c(bma.probs, ppc.prob)


								}





								if (length(pvals) == 0){
									next
								}

								pvals = pvals[!is.na(pvals)]
								if (length(pvals) == 0){
									next
								}




								ppc.p1 = pvals[1] # Just take the first p value, no averaging


								bma.probs = bma.probs[!is.na(bma.probs)]
								if (length(bma.probs) == 0){
									bma.prob1 = NA
									bma.prob.all = NA
								}else{
									bma.prob1 = bma.probs[1]

									ymin=1e-3
									bma.probs[bma.probs>1-ymin] = 1-ymin
									bma.probs[bma.probs<ymin] = ymin

									bma.prob.all = plogis(mean(qlogis(bma.probs)))
								}
								
								

								out.df2 = data.frame(sampleNr = trial, tree.dir = tree.dir, clock.sigma=clock.sigma, ntaxa=ntaxa, clock.model = clock.model, trait.method = trait.method, youngest=youngest, nested=nested, maximise=maximise, 
															time.tree.cond=time.tree.cond, beta=beta, 
															 number.of.subst=number.of.subst.total, logY=TRUE, standardise=2,
															  dmin=dmin, dspan=dspan, ls.p=NA, pic.p.st=NA, pic.p.tt=NA, pgls.lambda.p=NA, pgls.p.tt=NA, pgls.ou.p=NA, pgls.rho.p=NA,  pgls.blomberg.p=NA,
															  pgls.p.st=NA, ppc.p1=ppc.p1, bma.prob1=bma.prob1, bma.prob.all=bma.prob.all, rho=median(rhovals), n=median(nvals), t=NA, rate.noise=rate.noise, trait.noise=noise)





								out.df = rbind(out.df, out.df2)
						
				
							}
						}

					}
				}
			}

		}

	}


}


