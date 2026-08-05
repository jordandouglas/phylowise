


library(ape)
library(phylowise)





plot.all = function(min.obs=10, ks.threshold=0.01) {



	ff = list.files(".", pattern="false.+tsv")
	cat(paste("Loading", ff[1], "\n"))
	disco.df = read.table(ff[1], sep="\t", header=T)
	if (length(ff) > 1){
		for (f in ff[2:length(ff)]){
			cat(paste("Loading", f, "\n"))
			disco1.df = read.table(f, sep="\t", header=T)
			disco.df = rbind(disco.df, disco1.df)
		}
	}
	

	

	disco.df = disco.df[disco.df$n >= min.obs,]
	pgls.df = disco.df[!is.na(disco.df$ls.p),]



	write.label = function(txt, srt=0, cex=1.0){
		plot(0,0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), axes=F, xaxs="i", yaxs="i")

		adj=0.5
		x = 0.5
		if (srt == 90){
			adj = c(0.5, 1)
			x = 0.1
		}
		text(x, 0.5, txt, cex=cex, srt=srt, adj=adj)
	}

	plot.hist = function(p, ylim=1.2, sig.threshold=0.01, rate.noise="none", max.nrow = 10000){


		m1=0.5
		m2=1.1
		m3 = (m1+m2)/2

		# Left aligned
		if (rate.noise == "none"){
			par(mar=c(m3,m2,m3,m1))
		}else {
			par(mar=c(m3,m1,m3,m2))
		}

		p = p[!is.na(p)]


		if (length(p) > max.nrow){
			keep = sample(1:length(p), max.nrow)
			p = p[keep]
		}

		print(length(p))

		if (length(p) < 2){
			plot(0,0,type="n")
		}

		if (length(p) > 5){
			ks1 = ks.test(p, alternative="less", "punif", 0, 1)
			ks2 = ks.test(p, alternative="greater", "punif", 0, 1)
			p1 = signif(ks1$p.value, 3)
			p2 = signif(ks2$p.value, 3)
			if (is.na(p1)){
				p1=1
			}
			if (is.na(p2)){
				p2=1
			}
			col1 = ifelse(p1 < ks.threshold, low.col, ifelse(p2 < ks.threshold, high.col, good.col))
			
		}else{
			col1 = good.col
		}
		col2 = paste0(col1, "99")
		hist(p, breaks=breaks, probability=T,
		 xaxs="i", yaxs="i", axes=F, xlab="", ylab="", main="", col=col2, border=col1)

		if (col1 != good.col){
			psig = signif(100*sum(p <= sig.threshold) / length(p), 2)
			mtext(paste0("P(p<0.01)=", psig, "%"), cex=0.5)
		}
		

	}


	################
	# PRETTY PLOTS #
	################

	rate.noise.all = c("none", "poisson")
	trait.noise.all = c("none", "rounding", "gaussian")
	

	good.col="#d3d3d3"
	low.col="#56B4E9"
	high.col="#E69F00"

	
	breaks=10

	demo.cex=0.7
	data.cex=0.85

	

	png(paste0("FP.png"), width=3000, height=3000, res=435)
	par(xpd=T)
	layout(matrix(c(1,  2,  2,  3,  3,  4,  4,0,
					0,  0,  0,  0,  0,  0,  0,0,
					0,  5,  5,  6,  6,  7, 7,84,
					0,  8,  9,  10,  11,  12,  13,85,
					14:20,86,
					21:27,87,
					28:34,88,
					35:41,89,
					rep(0,8),
					42:48,90,
					49:55,91,
					56:62,92,
					63:69,93,
					70:76,94,
					77:83,95), nrow=15, ncol=8, byrow=T), width=c(0.5,2,2,2,2,2,2,0.8), height=c(4, 0.2, 0.5,0.5, 2,2,2,2,0.5,2,2,2,2,2,2))

	# layout(matrix(c(61,  1,  1,  2,  2,  3,  3,0,
	# 				0,  0,  0,  0,  0,  0,  0,0,
	# 				0,  4,  4,  5,  5,  6, 6,0,
	# 				7,  9, 10, 13, 14, 17, 18,0,
	# 				7, 11, 12, 15, 16, 19, 20,0,
	# 				8, 21, 22, 25, 26, 29, 30,0,
	# 				8, 23, 24, 27, 28, 31, 32,0,
	# 				0,  62,  63,  64,  65,  66,  67,0,
	# 				57, 33, 45, 37,  49,  41,  53,0,
	# 				58, 34, 46, 38,  50,  42,  54,0,
	# 				59, 35, 47, 39,  51,  43,  55,0,
	# 				60, 36, 48, 40,  52,  44,  56,0), nrow=12, ncol=8, byrow=T), width=c(0.5,2,2,2,2,2,2,0.5), height=c(3.5,0.5,0.5,2,2,2,2,0.5,1.5,1.5,1.5,1.5))


	


	par(mar=c(0,0,0,0))
	write.label("Legend", srt=90)


	# Demonstration plots
	par(mar=c(2,3,1.5,1))
	plot(0,0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,2), axes=F, xaxs="i", yaxs="i")
	mtext("Inflation of type I error", cex=demo.cex)
	lines(c(0, 1), c(2, 0), lwd=7, col=paste0(high.col, "55"))
	lines(c(0, 1), c(2, 0), lwd=5, col=paste0(high.col, "aa"))
	lines(c(0, 1), c(2, 0), lwd=3, col=paste0(high.col, ""))
	axis(1, at=c(0, 1))
	axis(2, las=2, at=c(0, 1, 2))
	text(0.5, -0.5, "p-value", cex=1)
	text(-0.17, 1, "Probability density", cex=1, srt=90)


	par(mar=c(2,2,1.5,2))
	plot(0,0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,2), axes=F, xaxs="i", yaxs="i")
	mtext("Well-calibrated", cex=demo.cex)
	lines(c(0, 1), c(1, 1), lwd=7, col=paste0(good.col, "55"))
	lines(c(0, 1), c(1, 1), lwd=5, col=paste0(good.col, "aa"))
	lines(c(0, 1), c(1, 1), lwd=3, col=paste0(good.col, ""))
	axis(1, at=c(0, 1))
	axis(2, las=2, at=c(0, 1, 2))
	text(0.5, -0.5, "p-value", cex=1)

	plot(0,0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,2), axes=F, xaxs="i", yaxs="i")
	mtext("Deflation of type I error", cex=demo.cex)
	lines(c(0, 1), c(0, 2), lwd=7, col=paste0(low.col, "55"))
	lines(c(0, 1), c(0, 2), lwd=5, col=paste0(low.col, "aa"))
	lines(c(0, 1), c(0, 2), lwd=3, col=paste0(low.col, ""))
	axis(1, at=c(0, 1))
	axis(2, las=2, at=c(0, 1, 2))
	text(0.5, -0.5, "p-value", cex=1)



	# Column names
	par(mar=c(0,0,0,0))
	write.label("No trait error", cex=1.3)
	write.label("Rounding error", cex=1.3)
	write.label("Gaussian error", cex=1.3)



	write.label("No molecular error", cex=1)
	write.label("Poisson error", cex=1)

	write.label("No molecular error", cex=1)
	write.label("Poisson error", cex=1)

	write.label("No molecular error", cex=1)
	write.label("Poisson error", cex=1)




	
	# PPC
	par(mar=c(0,0,0,0))
	write.label("PPC", srt=90)
	for (noise in trait.noise.all){
		for (rate.noise in rate.noise.all){
			sub.df = disco.df[disco.df$rate.noise == rate.noise & disco.df$trait.noise == noise,]
			sub.df = sub.df[sub.df$nested == FALSE & sub.df$maximise == FALSE,]
			plot.hist(sub.df$ppc.p, ylim=3, rate.noise=rate.noise)
		}
	}


	# PPC nest
	par(mar=c(0,0,0,0))
	write.label("PPC\n+nst", srt=90)
	for (noise in trait.noise.all){
		for (rate.noise in rate.noise.all){
			sub.df = disco.df[disco.df$rate.noise == rate.noise & disco.df$trait.noise == noise,]
			sub.df = sub.df[sub.df$nested == TRUE & sub.df$maximise == FALSE,]
			plot.hist(sub.df$ppc.p, ylim=3, rate.noise=rate.noise)
		}
	}

	# PPC max
	par(mar=c(0,0,0,0))
	write.label("PPC\n+max", srt=90)
	for (noise in trait.noise.all){
		for (rate.noise in rate.noise.all){
			sub.df = disco.df[disco.df$rate.noise == rate.noise & disco.df$trait.noise == noise,]
			sub.df = sub.df[sub.df$nested == FALSE & sub.df$maximise == TRUE,]
			plot.hist(sub.df$ppc.p, ylim=3, rate.noise=rate.noise)
		}
	}

	# PPC nestmax
	par(mar=c(0,0,0,0))
	write.label("PPC\n+nst +max", srt=90)
	for (noise in trait.noise.all){
		for (rate.noise in rate.noise.all){
			sub.df = disco.df[disco.df$rate.noise == rate.noise & disco.df$trait.noise == noise,]
			sub.df = sub.df[sub.df$nested == TRUE & sub.df$maximise == TRUE,]
			plot.hist(sub.df$ppc.p, ylim=3, rate.noise=rate.noise)
		}
	}



	# OLS
	par(mar=c(0,0,0,0))
	write.label("OLS", srt=90)
	for (noise in trait.noise.all){
		for (rate.noise in rate.noise.all){
			sub.df = disco.df[disco.df$rate.noise == rate.noise & disco.df$trait.noise == noise,]
			plot.hist(sub.df$ls.p, ylim=3, rate.noise=rate.noise)
		}
	}


	# Brownian
	par(mar=c(0,0,0,0))
	write.label("PGLS\nBrownian", srt=90)
	for (noise in trait.noise.all){
		for (rate.noise in rate.noise.all){
			sub.df = disco.df[disco.df$rate.noise == rate.noise & disco.df$trait.noise == noise,]
			plot.hist(sub.df$pgls.p.tt, ylim=3, rate.noise=rate.noise)
		}
	}


	# OU
	par(mar=c(0,0,0,0))
	write.label("PGLS\nMartins", srt=90)
	for (noise in trait.noise.all){
		for (rate.noise in rate.noise.all){
			sub.df = disco.df[disco.df$rate.noise == rate.noise & disco.df$trait.noise == noise,]
			plot.hist(sub.df$pgls.ou.p, ylim=3, rate.noise=rate.noise)
		}
	}


	# Lambda λ
	par(mar=c(0,0,0,0))
	write.label("PGLS\nPagel", srt=90)
	for (noise in trait.noise.all){
		for (rate.noise in rate.noise.all){
			sub.df = disco.df[disco.df$rate.noise == rate.noise & disco.df$trait.noise == noise,]
			plot.hist(sub.df$pgls.lambda.p, ylim=3, rate.noise=rate.noise)
		}
	}



	

	# Rho
	par(mar=c(0,0,0,0))
	write.label("PGLS\nGrafen", srt=90)
	for (noise in trait.noise.all){
		for (rate.noise in rate.noise.all){
			sub.df = disco.df[disco.df$rate.noise == rate.noise & disco.df$trait.noise == noise,]
			plot.hist(sub.df$pgls.rho.p, ylim=3, rate.noise=rate.noise)
		}
	}



	# Blomberg
	par(mar=c(0,0,0,0))
	write.label("PGLS\nBlomberg", srt=90)
	for (noise in trait.noise.all){
		for (rate.noise in rate.noise.all){
			sub.df = disco.df[disco.df$rate.noise == rate.noise & disco.df$trait.noise == noise,]
			plot.hist(sub.df$pgls.blomberg.p, ylim=3, rate.noise=rate.noise)
		}
	}




	par(mar=c(0,0,0,0))
	plot(0,0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), axes=F, xaxs="i", yaxs="i")



	dev.off()




}



plot.all(min.obs=20, ks.threshold=0.01)




