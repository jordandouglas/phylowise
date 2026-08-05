


library(ape)
library(phylowise)
library(nlme)


args = commandArgs(trailingOnly = TRUE)
#args=c("false.tsv", FALSE, FALSE)

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
clock.models = c("AC") 
clock.rate.mean=2000
clock.rate.sd=2


# Window size
dmin.subst=10
dmax.subst=200
dmin.time=0.05
dmax.time=0.15




dmin.all.time = c(0.05)
dmin.all = c(50)
dspan.all = c(5)

nwindows1 = length(dmin.all.time)
nwindows2 = length(dspan.all)

# Statistical test options
extreme.value.threshold=0


max.wall.time = 10


# False discovery rate
trait.methods = c("BM") 
tree.dirs = c("../newick/")

setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits
out.df = data.frame(sampleNr = numeric(0), tree.dir = character(0), clock.sigma=numeric(0), R2.true=numeric(0), ntaxa=numeric(0), clock.model = character(0), beta=numeric(0), trait.method = character(0), nested=logical(0), youngest=logical(0),
	maximise=logical(0), logY=logical(0), standardise=numeric(0), number.of.subst=numeric(0), dmin=numeric(0), dspan=numeric(0), time.tree.cond=logical(0), pgls.ou.p=numeric(0), pgls.rho.p=numeric(0), pgls.blomberg.p=numeric(0), pic.p.st=numeric(0), pic.p.tt=numeric(0), pgls.lambda.p=numeric(0), ls.p=numeric(0), pgls.p.st=numeric(0), pgls.p.tt=numeric(0), 
	ppc.p=numeric(0), rho = numeric(0), n=numeric(0), t=numeric(0), trait.noise=character(0), rate.noise=character(0))


time.tree.cond=TRUE
for (trial in 1:ntrials){


	if (trial %% 1 == 0){
		cat(paste(trial, "/", ntrials, "\n"))
		write.table(out.df, out.file, sep="\t", quote=F, row.names=F)
	}

	for (tree.dir in tree.dirs){



		# Simulate a time tree with height up to 1. Extinct taxa removed
		time.tree = sample.newick(dir=tree.dir)
		dmat.time = phylowise::getDistanceMatrix(time.tree)
		ntaxa = length(time.tree$tip.label)


		# Simulate a subst tree
		for (clock.model in clock.models){

			clock.model="TD"
			beta = 0
			if (clock.model != "TD"){
				clock.sigma.mean = 0.3
				clock.sigma.sd = 0.4
				clock.sigma = rlnorm(1, log(clock.sigma.mean) - 0.5*clock.sigma.sd^2, clock.sigma.sd)
			}else{
				theta = rexp(1, 0.1)
				V = 0.2
				clock.sigma = sqrt(V*2*theta) 
			}



			for (trait.method in trait.methods){

				# Simulate traits under this method
				traits = phylowise::simulateTrait(time.tree, sigma=1)
				traits = scale(traits)


				

				number.of.subst = rlnorm(1, log(clock.rate.mean) - 0.5*clock.rate.sd^2, clock.rate.sd)
				result = phylowise::simulateSubstitutions(time.tree, traits=traits, beta=0, theta=theta, sigma=clock.sigma, number.of.subst=number.of.subst, method=clock.model)
				subst.tree.true = result$subst.tree.true
				subst.tree.est = result$subst.tree.est
				node.rates.est = result$node.rates.est
				node.rates.true = result$node.rates.true


				for (rate.noise in c("none", "poisson")){

					if (rate.noise == "none"){
						subst.tree = subst.tree.true
					}

					if (rate.noise == "poisson"){
						subst.tree = subst.tree.est
					}

					number.of.subst.total = sum(subst.tree$edge.length)



					# Distance matrix now (to save time)
					dmat.subst = phylowise::getDistanceMatrix(subst.tree)


					for (noise in c("none", "rounding", "gaussian")){

						traits.noise = traits

						if (noise == "rounding"){
							traits.noise = round(exp(traits.noise), 1)
							traits.noise = ifelse(traits.noise == 0, 0.01, traits.noise)
							traits.noise = log(traits.noise)
						}

						if (noise == "gaussian"){
							traits.noise = traits.noise + rnorm(length(traits.noise), 0, 0.2)
						}


						for (logY in c(TRUE)){

							# PGLS test from the nlme package
							leaf.rates.est = node.rates.est[1:length(subst.tree$tip.label)]
							leaf.rates.true = node.rates.true[1:length(subst.tree$tip.label)]
							leaf.traits = traits.noise[1:length(subst.tree$tip.label)]

							if (logY){
								leaf.rates.est = log(leaf.rates.est) 
								leaf.rates.true = log(leaf.rates.true)
							}




							if (rate.noise == "none"){
								leaf.rates = leaf.rates.true
							}

							if (rate.noise == "poisson"){
								leaf.rates = leaf.rates.est
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
							

							if (!SKIP.PGLS){


								ls = cor.test(leaf.traits.nonzero, leaf.rates)
								ls.p = ls$p.value


								if (time.tree.cond){
									basis.tree = time.tree.nonzero.rates
								}else{
									basis.tree = subst.tree.nonzero.rates
								}


								# PGLM on subst tree 
								setTimeLimit(elapsed = max.wall.time)   # Stop after a few seconds of wall-clock time
								pgls.p.st = tryCatch({
									pgls = gls(leaf.rates~leaf.traits.nonzero, correlation=corBrownian(1,subst.tree.nonzero.rates), method="ML")
									pgls.result = summary(pgls)
									pgls.result$tTable[2,4]
								}, error = function(e) {
									print("Timeout 1")
								  	NA
								})
								setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits


								# PGLM on subst tree 
								setTimeLimit(elapsed = max.wall.time)   # Stop after a few seconds of wall-clock time
								pgls.p.tt = tryCatch({
									pgls = gls(leaf.rates~leaf.traits.nonzero, correlation=corBrownian(1,basis.tree), method="ML")
									pgls.result = summary(pgls)
									pgls.result$tTable[2,4]
								}, error = function(e) {
									print("Timeout 2")
								  	NA
								})
								setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits



								# PGLS on time tree lambda with multiple starting points
								setTimeLimit(elapsed = max.wall.time)   # Stop after a few seconds of wall-clock time
								pgls.lambda.p=NA
								for(initVal in seq(from=1, to=0.1,by=-0.1)){

									pgls.lambda.p = tryCatch({
										pgls = gls(leaf.rates~leaf.traits.nonzero, correlation=corPagel(initVal, basis.tree,fixed = FALSE), method="ML")
										pgls.result = summary(pgls)
										pgls.result$tTable[2,4]
									}, error = function(e) {
									  	NA
									})

									if (!is.na(pgls.lambda.p)){
										break
									}


								}

								if (is.na(pgls.lambda.p)){
									print("lambda did not converge")
								}
								
								setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits



								# OU process
								setTimeLimit(elapsed = max.wall.time)   # Stop after a few seconds of wall-clock time
								pgls.ou.p=NA
								for(initVal in c(0.01, 0.1, 0.5, 1.0, 2, 5, 10, 50, 100, 1000, 10000)){

									pgls.ou.p = tryCatch({
										pgls = gls(leaf.rates~leaf.traits.nonzero, correlation=corMartins(initVal, basis.tree,fixed = FALSE), method="ML")
										pgls.result = summary(pgls)
										pgls.result$tTable[2,4]
									}, error = function(e) {
									  	NA
									})

									if (!is.na(pgls.ou.p)){
										break
									}


								}

								if (is.na(pgls.ou.p)){
									print("ou did not converge")
								}
								setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits





								# Grafen rho process
								setTimeLimit(elapsed = max.wall.time)   # Stop after a few seconds of wall-clock time
								pgls.rho.p=NA
								for(initVal in c(0.01, 0.1, 0.5, 1.0, 2, 5, 10, 50, 100)){

									pgls.rho.p = tryCatch({
										pgls = gls(leaf.rates~leaf.traits.nonzero, correlation=corGrafen(initVal, basis.tree,fixed = FALSE), method="ML")
										pgls.result = summary(pgls)
										pgls.result$tTable[2,4]
									}, error = function(e) {
									  	NA
									})

									if (!is.na(pgls.rho.p)){
										break
									}


								}

								if (is.na(pgls.rho.p)){
									print("rho did not converge")
								}
								setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits



								# Blomberg rho process
								setTimeLimit(elapsed = max.wall.time)   # Stop after a few seconds of wall-clock time
								pgls.blomberg.p=NA
								for(initVal in c(1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.1, 0.5, 1.0, 2, 5, 10, 50, 100, 1000, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e10)){

									pgls.blomberg.p = tryCatch({
										pgls = gls(leaf.rates~leaf.traits.nonzero, correlation=corBlomberg(initVal, basis.tree,fixed = FALSE), method="ML")
										pgls.result = summary(pgls)
										pgls.result$tTable[2,4]
									}, error = function(e) {
									  	NA
									})

									if (!is.na(pgls.blomberg.p)){
										break
									}


								}

								if (is.na(pgls.blomberg.p)){
									print("blomberg did not converge")
								}
								setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits

								


								pic.p.st = NA


								# Time tree PIC
								pic1 = pic(leaf.rates, time.tree.nonzero.rates)
								pic2 = pic(leaf.traits.nonzero, time.tree.nonzero.rates)
								pic_model = lm(pic1 ~ pic2 - 1)
								pic.res = summary(pic_model)
								pic.p.tt = as.numeric(pic.res$coefficients[1,4])




								out.df2 = data.frame(sampleNr = trial, tree.dir = tree.dir, clock.sigma=clock.sigma, R2.true=NA, ntaxa=ntaxa, youngest=NA, number.of.subst=number.of.subst.total, clock.model = clock.model, trait.method = trait.method, nested=NA, maximise=NA, 
										logY=logY, standardise=NA, dmin=NA, dspan=NA, beta=beta, time.tree.cond=NA, 
										pic.p.st=pic.p.st, pic.p.tt=pic.p.tt, pgls.lambda.p=pgls.lambda.p, ls.p=ls.p, pgls.ou.p=pgls.ou.p, pgls.rho.p=pgls.rho.p, pgls.blomberg.p=pgls.blomberg.p,
										pgls.p.tt=pgls.p.tt, pgls.p.st=pgls.p.st, ppc.p=NA, rho=NA, n=length(leaf.rates), t=NA, rate.noise=rate.noise, trait.noise=noise)
								out.df = rbind(out.df, out.df2)

							}

							
							if (SKIP.PPC){
								next
							}

							for (nested in c(TRUE, FALSE)){
								for (maximise in c(FALSE, TRUE)){
									for (youngest in c(FALSE)){

										for (d1 in 1:nwindows1){
											for (d2 in 1:nwindows2){


												standardise = ifelse(time.tree.cond, 2, 1)


												if (time.tree.cond){
													dmin = dmin.all.time[d1]
												}else{
													dmin = dmin.all[d1]
												}

												dspan = dspan.all[d2]
												dmax = dmin * dspan
												#print(paste(dmin, dmax))


												pvals = numeric(0)
												nvals = numeric(0)
												rhovals = numeric(0)



												# Sample taxa - time tree
												setTimeLimit(elapsed = max.wall.time)   # Stop after a few seconds of wall-clock time
												pairs.df = tryCatch({
													if (time.tree.cond){
														phylowise::sampleTaxonPairs(subst.tree=subst.tree, window.tree=time.tree, covariate=traits.noise, dist.min=dmin, dist.max=dmax, nested=nested, youngest=youngest, maximise=maximise, distance.matrix=dmat.time, verbose=F)
													}else{
														phylowise::sampleTaxonPairs(subst.tree=subst.tree, window.tree=subst.tree, covariate=traits.noise, dist.min=dmin, dist.max=dmax, nested=nested, youngest=youngest, maximise=maximise, distance.matrix=dmat.subst, verbose=F)
													}
												}, error = function(e) {
													print(e)
													setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits
													print(paste("Timeout 3", dmin, dmax, nested, maximise, time.tree.cond))
												  	data.frame(x=1)
												})
												setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)  # Reset limits
												
												#plot(pairs.df$trait1 - pairs.df$trait2, pairs.df$d1 - pairs.df$d2)

												if (is.na(pairs.df[1]) || nrow(pairs.df) < 3){
													next
												}


												# PPC test
												test.results = phylowise::PPC.test(pairs.df, standardise=standardise, logY=logY, extreme.value.threshold=extreme.value.threshold)


												pvals = c(pvals, test.results$p)
												nvals = c(nvals, test.results$n)
												rhovals = c(rhovals, test.results$rho)


												pvals = pvals[!is.na(pvals)]
												if (length(pvals) == 0){
													next
												}


												if (length(pvals) == 1){
													p_combined = pvals
												}else{
													p_combined = p.hmp(pvals, L=length(pvals))
												}


												out.df2 = data.frame(sampleNr = trial, tree.dir = tree.dir, clock.sigma=clock.sigma, R2.true=NA, ntaxa=ntaxa, clock.model = clock.model, trait.method = trait.method, youngest=youngest, nested=nested, maximise=maximise, 
														time.tree.cond=time.tree.cond, beta=beta, 
														 number.of.subst=number.of.subst.total, logY=logY, standardise=standardise,
														  dmin=dmin, dspan=dspan, ls.p=NA, pic.p.st=NA, pic.p.tt=NA, pgls.lambda.p=NA, pgls.p.tt=NA, pgls.ou.p=NA, pgls.rho.p=NA,  pgls.blomberg.p=NA,
														  pgls.p.st=NA, ppc.p=p_combined, rho=median(rhovals), n=median(nvals), t=NA, rate.noise=rate.noise, trait.noise=noise)
												out.df = rbind(out.df, out.df2)
										
											}
										}

									}
								}
							}
						}
					}
				}
			}

		}

	}


}

