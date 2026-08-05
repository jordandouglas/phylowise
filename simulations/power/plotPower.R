
library(IDPmisc)


bone="#696969"
model.cols = c(ppcNst="#0072B2", ppcNstMax="#E69F00", ols="#111111", 
				pgls="#009E73", lambda="#56B4E9", grafen="#D55E00", blomberg="#CC79A7", martins="#F0E442")


show.bf=TRUE
bfmin=-2
bfmax=3

ymin=1e-3

youngest=FALSE
poly=FALSE
sig.threshold=0.001
prior.prob=0.5

pvars = c("bma.prob1", "bma.prob.all")
beta.var = "beta"
pvar = pvars[2]

nvar = "n"

hist.col1="white"
hist.col2="#E3DAC999"
hist.col2="black"


nbreaks=20

min.obs=40
min.subst=10000
min.r2=0.4

beta.max=2
n.max=400

max.nrow = 10000


print.line.pval = function(x, p, col, lty, xlim=c(0,1), logging=F, binomial=T){



	keep = !is.na(p) & !is.na(x)
	x = x[keep]
	p = p[keep]

	if (length(p) > max.nrow){
		keep = sample(1:length(p), max.nrow)
		x = x[keep]
		p = p[keep]
	}

	print(length(p))

	if (length(p) < 20){
		return("")
	}


	y = ifelse(p < sig.threshold, 1, 0)
	

	vals = seq(from=xlim[1], to=xlim[2], length=500)

	minx = min(x)
	vals = vals[vals > minx]

	if (logging){
		vals.log = log(vals)
		vals.log = vals.log[vals.log != -Inf]
		vals = exp(vals.log)
		xt = log(x)
	}else{
		vals.log=vals
		xt=x
	}


	if (binomial){
		m = glm(y~xt, family="binomial")
	}else{
		m = lm(y~xt + I(xt^2) + I(xt^3) + I(xt^4)  + I(xt^5))
	}

	
	pred = predict(m, data.frame(xt=vals.log), interval = "confidence", level=0.95)


	if (binomial){

		yy = plogis(pred)
		lines(vals, yy, col=paste0(col, "44"), lty=lty, lwd=4.5)
		lines(vals, yy, col=paste0(col, "88"), lty=lty, lwd=3.5)
		lines(vals, yy, col=col, lty=lty, lwd=3)
	
	}else{

		if (poly){
			lwr = as.numeric(pred[,"lwr"])
			upr = as.numeric(pred[,"upr"])
			polygon(c(vals[1], vals, rev(vals), vals[1]), c(lwr[1], lwr, rev(upr), upr[1]), col=paste0(col, "66"), border=NA)
		}
		lines(vals, pred[,"fit"], col=col, lty=lty, lwd=3)

	}




}



print.line = function(x, p, col, lty, xlim=c(0,1), logging=F){



	keep = !is.na(p) & !is.na(x) 
	x = x[keep]
	p = p[keep]


	if (length(p) > max.nrow){
		keep = sample(1:length(p), max.nrow)
		x = x[keep]
		p = p[keep]
	}

	print(length(p))

	if (length(p) < 20){
		return("")
	}


	y = p# ifelse(p < sig.threshold, 1, 0)
	

	vals = seq(from=xlim[1], to=xlim[2], length=500)

	minx = min(x)
	vals = vals[vals > minx]

	if (logging){

		

		vals.log = log(vals)
		vals.log = vals.log[vals.log != -Inf]
		vals = exp(vals.log)
		xt = log(x)
	}else{
		vals.log=vals
		xt=x
	}



	# Logit transform for probabilities, but not bayes factors
	y[y>1-ymin] = 1-ymin
	y[y<ymin] = ymin
	if (!show.bf){
		z = qlogis(y)
	}else{

		z = (y/(1-y)) / (prior.prob / (1-prior.prob)) 
		z = log(z, 10)

	}



	m = lm(z~xt + I(xt^2) + I(xt^3))# + I(xt^4) + I(xt^5))
	
	pred = predict(m, data.frame(xt=vals.log), interval = "confidence", level=0.95)

	if (!show.bf){
		pred[,"fit"] = plogis(pred[,"fit"])
		pred[,"lwr"] = plogis(pred[,"lwr"])
		pred[,"upr"] = plogis(pred[,"upr"])
	}

	par(xpd=T)
	if (poly){
		lwr = as.numeric(pred[,"lwr"])
		upr = as.numeric(pred[,"upr"])
		polygon(c(vals[1], vals, rev(vals), vals[1]), c(lwr[1], lwr, rev(upr), upr[1]), col=paste0(col, "66"), border=NA)
		lines(vals, pred[,"fit"], col=col, lty=lty, lwd=3)
	}else{
		lines(vals, pred[,"fit"], col=paste0(col, "44"), lty=lty, lwd=4.5)
		lines(vals, pred[,"fit"], col=paste0(col, "88"), lty=lty, lwd=3.5)
		lines(vals, pred[,"fit"], col=col, lty=lty, lwd=3)
	}
	
	par(xpd=F)




}



plot.bg = function(xlim, show.legend=F){


	cc = "#E3DAC955" # bone

	if (show.bf){
		
		rect(xlim[1], -1, xlim[2]+1, 1, col="#eeeeee44")
		lines(xlim, c(0, 0), lty="dashed", lwd=1)

		if (show.legend){
			text(xlim[1] + (xlim[2]-xlim[1])*0.02,  1.2, "Strong support for", adj=c(0, 0), cex=0.75)
			text(xlim[1] + (xlim[2]-xlim[1])*0.02, -1.2, "Strong support against", adj=c(0, 1), cex=0.75)
		}

	}else{

		bf = 10/11
		rect(xlim[1], 1-bf, xlim[2]+1, bf, col=cc)
		text(xlim[1] + (xlim[2]-xlim[1])*0.02,  bf+0.01, "BF > +1", adj=c(0, 0), cex=0.8)
		text(xlim[1] + (xlim[2]-xlim[1])*0.02, 1-bf-0.01, "BF < −1", adj=c(0, 1), cex=0.8)

		lines(xlim, c(prior.prob, prior.prob), lty="dashed")
		text(xlim[1] + (xlim[2]-xlim[1])*0.02, prior.prob+0.05, "Prior", adj=0, cex=0.8)


		


	}

}



plot.beta.pval = function(sub.df,  pvar="ppc.p1",  xlim=c(0,beta.max)){

	plot(0,0, type="n", axes=F, xaxs="i", yaxs="i", xlab="", ylab="", main="", xlim=xlim, ylim=c(0,1.05))


	s.df = sub.df[!is.na(sub.df$pgls.p.tt),]

	
	
	# PGLS blomberg
	print.line.pval(s.df$beta, s.df$pgls.blomberg.p, col=model.cols["blomberg"], lty="solid", xlim=xlim)

	# PGLS grafen
	print.line.pval(s.df$beta, s.df$pgls.rho.p, col=model.cols["grafen"], lty="solid", xlim=xlim)

	# PGLS pagel
	print.line.pval(s.df$beta, s.df$pgls.lambda.p, col=model.cols["lambda"], lty="solid", xlim=xlim)


	# PGLS martins
	print.line.pval(s.df$beta, s.df$pgls.ou.p, col=model.cols["martins"], lty="solid", xlim=xlim)

	# OLS
	print.line.pval(s.df$beta, s.df$ls.p, col=model.cols["ols"], lty="solid", xlim=xlim)

		

	# PGLS brownian
	print.line.pval(s.df$beta, s.df$pgls.p.tt, col=model.cols["pgls"], lty="solid", xlim=xlim)
	
	# PPC nomax
	s.df = sub.df[!is.na(sub.df$maximise) & sub.df$maximise==FALSE,]
	print.line.pval(s.df$beta, s.df[,pvar], col=model.cols["ppcNst"], lty="solid", xlim=xlim)

	# PPC max
	s.df = sub.df[!is.na(sub.df$maximise) & sub.df$maximise==TRUE,]
	print.line.pval(s.df$beta, s.df[,pvar], col=model.cols["ppcNstMax"], lty="solid", xlim=xlim)

	axis(1)
	axis(2, las=2)


	# Arrow
	arrowY = 0.2
	Arrows(0.45, arrowY, 0.1, arrowY, open=FALSE, sh.adj=1, sh.lwd=2, size=0.5, sh.col="#696969", h.lwd=0.8)
	text(0.26, 0.28, "False\npositives", adj=c(0,0.15), cex=0.75)

}



plot.n.pval = function(sub.df,  pvar="ppc.p1",  xlim=c(0,n.max)){

	plot(0,0, type="n", axes=F, xaxs="i", yaxs="i", xlab="", ylab="", main="", xlim=xlim, ylim=c(0,1.05))


	s.df = sub.df[!is.na(sub.df$pgls.p.tt),]

	

	# PGLS blomberg
	print.line.pval(s.df[,nvar], s.df$pgls.blomberg.p, col=model.cols["blomberg"], lty="solid", xlim=xlim)

	# PGLS grafen
	print.line.pval(s.df[,nvar], s.df$pgls.rho.p, col=model.cols["grafen"], lty="solid", xlim=xlim)

	# PGLS pagel
	print.line.pval(s.df[,nvar], s.df$pgls.lambda.p, col=model.cols["lambda"], lty="solid", xlim=xlim)

	# PGLS martins
	print.line.pval(s.df[,nvar], s.df$pgls.ou.p, col=model.cols["martins"], lty="solid", xlim=xlim)

	# OLS
	print.line.pval(s.df[,nvar], s.df$ls.p, col=model.cols["ols"], lty="solid", xlim=xlim)
		

	xlim2 = c(0,250)


	# PGLS brownian
	print.line.pval(s.df[,nvar], s.df$pgls.p.tt, col=model.cols["pgls"], lty="solid", xlim=xlim)
	
	# PPC nomax
	s.df = sub.df[!is.na(sub.df$maximise) & sub.df$maximise==FALSE,]
	print.line.pval(s.df[,nvar], s.df[,pvar], col=model.cols["ppcNst"], lty="solid", xlim=xlim2)

	# PPC max
	s.df = sub.df[!is.na(sub.df$maximise) & sub.df$maximise==TRUE,]
	print.line.pval(s.df[,nvar], s.df[,pvar], col=model.cols["ppcNstMax"], lty="solid", xlim=xlim2)

	axis(1, at=c(0, 100, 200, 300, 400))
	axis(2, las=2)

}




plot.subst.pval = function(sub.df,  pvar="ppc.p1",  xlim=c(2,7)){

	plot(0,0, type="n", axes=F, xaxs="i", yaxs="i", xlab="", ylab="", main="", xlim=xlim, ylim=c(0,1.05))


	s.df = sub.df[!is.na(sub.df$pgls.p.tt),]

	
	

	# PGLS blomberg
	print.line.pval(log(s.df$number.of.subst, 10), s.df$pgls.blomberg.p, col=model.cols["blomberg"], lty="solid", xlim=xlim)

	# PGLS grafen
	print.line.pval(log(s.df$number.of.subst, 10), s.df$pgls.rho.p, col=model.cols["grafen"], lty="solid", xlim=xlim)

	# PGLS pagel
	print.line.pval(log(s.df$number.of.subst, 10), s.df$pgls.lambda.p, col=model.cols["lambda"], lty="solid", xlim=xlim)

	# PGLS martins
	print.line.pval(log(s.df$number.of.subst, 10), s.df$pgls.ou.p, col=model.cols["martins"], lty="solid", xlim=xlim)
		

	# OLS
	print.line.pval(log(s.df$number.of.subst, 10), s.df$ls.p, col=model.cols["ols"], lty="solid", xlim=xlim)


	# PGLS brownian
	print.line.pval(log(s.df$number.of.subst, 10), s.df$pgls.p.tt, col=model.cols["pgls"], lty="solid", xlim=xlim)
	
	# PPC nomax
	s.df = sub.df[!is.na(sub.df$maximise) & sub.df$maximise==FALSE,]
	print.line.pval(log(s.df$number.of.subst, 10), s.df[,pvar], col=model.cols["ppcNst"], lty="solid", xlim=xlim)

	# PPC max
	s.df = sub.df[!is.na(sub.df$maximise) & sub.df$maximise==TRUE,]
	print.line.pval(log(s.df$number.of.subst, 10), s.df[,pvar], col=model.cols["ppcNstMax"], lty="solid", xlim=xlim)

	axis(1, at=c(2:7), labels=c("10²", "10³", "10⁴", "10⁵", "10⁶", "10⁷"))
	axis(2, las=2)

}





plot.beta = function(sub.df, main="", sig.threshold, pvar, beta.var, show.legend=F, xlim=c(0,1), show.ylab=FALSE){


	 
	if (show.bf){
		ylim = c(bfmin,bfmax)
	}else{
		ylim=c(0,1)
	}

	
	plot(0,0, type="n", axes=F, xaxs="i", yaxs="i", xlab="", ylab="", main=main, xlim=xlim, ylim=ylim)

	xlim2 = c(0, xlim[2])

	plot.bg(xlim, show.legend=T)

	s.df = sub.df[!is.na(sub.df$pgls.p.tt),]


	# PPC nst
	s.df = sub.df[!is.na(sub.df$maximise) & sub.df$maximise==FALSE & sub.df$youngest == youngest & sub.df$nested==TRUE,]
	print.line(s.df[,beta.var], s.df[,pvar], col=model.cols["ppcNst"], lty="solid", xlim=xlim2)


	# PPC nst max
	s.df = sub.df[!is.na(sub.df$maximise) & sub.df$maximise==TRUE  & sub.df$youngest == youngest & sub.df$nested==TRUE,]
	print.line(s.df[,beta.var], s.df[,pvar], col=model.cols["ppcNstMax"], lty="solid", xlim=xlim2)


	axis(1)
	axis(2, las=2)

}


plot.nobs = function(sub.df, main="", pvar, beta.var, sig.threshold, xlim=c(0,n.max), show.ylab=FALSE){

	if (show.bf){
		ylim = c(bfmin,bfmax)
	}else{
		ylim = c(0,1)
	}

	plot(0,0, type="n", axes=F, xaxs="i", yaxs="i", xlab="", ylab="", main=main, xlim=xlim, ylim=ylim)

	plot.bg(xlim)

	xlim2=c(0,250)

	s.df = sub.df[!is.na(sub.df$pgls.p.tt),]


	# PPC nst
	s.df = sub.df[!is.na(sub.df$maximise) & sub.df$maximise==FALSE & sub.df$youngest == youngest & sub.df$nested==TRUE,]
	print.line(s.df[,nvar], s.df[,pvar], col=model.cols["ppcNst"], lty="solid", xlim=xlim2)


	# PPC nst max
	s.df = sub.df[!is.na(sub.df$maximise) & sub.df$maximise==TRUE & sub.df$youngest == youngest & sub.df$nested==TRUE,]
	print.line(s.df[,nvar], s.df[,pvar], col=model.cols["ppcNstMax"], lty="solid", xlim=xlim2)


	

	axis(1, at=c(0, 100, 200, 300, 400))
	axis(2, las=2)

}


plot.subst = function(sub.df, main="", sig.threshold, pvar, beta.var, xlim=c(2,7), show.ylab=FALSE){


	if (show.bf){
		ylim = c(bfmin,bfmax)
	}else{
		ylim = c(0,1)
	}

	plot(0,0, type="n", axes=F, xaxs="i", yaxs="i", xlab="", ylab="", main=main, xlim=xlim, ylim=ylim)

	plot.bg(xlim)


	s.df = sub.df[!is.na(sub.df$pgls.p.tt),]



	# PPC nst
	s.df = sub.df[!is.na(sub.df$maximise) & sub.df$maximise==FALSE & sub.df$youngest == youngest & sub.df$nested==TRUE,]
	print.line(log(s.df$number.of.subst, 10), s.df[,pvar], col=model.cols["ppcNst"], lty="solid", xlim=xlim, logging=F)


	# PPC nst max
	s.df = sub.df[!is.na(sub.df$maximise) & sub.df$maximise==TRUE  & sub.df$youngest == youngest & sub.df$nested==TRUE,]
	print.line(log(s.df$number.of.subst, 10), s.df[,pvar], col=model.cols["ppcNstMax"], lty="solid", xlim=xlim, logging=F)



	axis(1, at=c(2:7), labels=c("10²", "10³", "10⁴", "10⁵", "10⁶", "10⁷"))
	axis(2, las=2)

}


add.xlab = function(xlab){

	par(xpd=T)

	usr = par("usr")

	# temporarily set coordinates to 0-1
	par(usr = c(0, 1, 0, 1))

	text(0.5, -0.35, xlab, adj=0.5, cex=1)

	# restore original coordinates
	par(usr = usr)
	par(xpd=F)

}



add.ylab = function(xlab, hist=FALSE){

	par(xpd=T)

	usr = par("usr")

	# temporarily set coordinates to 0-1
	par(usr = c(0, 1, 0, 1))

	x=ifelse(hist, -0.02, -0.27) #-0.21 for bf
	text(x, 0.5, xlab, adj=0.5, cex=1, srt=90)

	# restore original coordinates
	par(usr = usr)
	par(xpd=F)

}



add.label = function(letter, is.hist=FALSE){

	#is.hist=F

	par(xpd=T)

	usr = par("usr")
	letter = toupper(letter)

	# temporarily set coordinates to 0-1
	par(usr = c(0, 1, 0, 1))

	if (is.hist){
		text(-0.25, 1.05, letter, adj=0, cex=1.3)
	}else{
		text(-0.25, 1.07, letter, adj=0, cex=1.3)
	}
	
	# restore original coordinates
	par(usr = usr)
	par(xpd=F)

}

write.label = function(txt, srt=0, cex=1.0){
	plot(0,0, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), axes=F, xaxs="i", yaxs="i")
	text(0.5, 0.5, txt, cex=cex, srt=srt)
}


plot.all = function(sig.threshold=0.01, min.obs=20, min.subst=1000, min.r2=0.2, poly=T,  beta.var=c("R2.log", "R2.real"), pvar = c("ppc.p1", "ppc.p.wilson", "ppc.p.harmonic")) {

	


	beta.var = beta.var[1]
	pvar = pvar[1]
	


	ff = list.files(".", pattern="power.+tsv")
	disco.df = read.table(ff[1], sep="\t", header=T)
	cat(paste("Loading", ff[1], "\n"))
	if (length(ff) > 1){
		for (f in ff[2:length(ff)]){
			cat(paste("Loading", f, "\n"))
			disco1.df = read.table(f, sep="\t", header=T)
			disco.df = rbind(disco.df, disco1.df)
		}
	}

	

	true.df = disco.df[disco.df$clock.model == "TD",]


	r2.df = true.df[true.df[,nvar] >= min.obs & true.df$number.of.subst >= min.subst,]
	nobs.df = true.df[true.df[,beta.var] >= min.r2 & true.df$number.of.subst >= min.subst,]
	subst.df = true.df[true.df[,beta.var] >= min.r2 & true.df[,nvar] >= min.obs,]





	png(paste0("power.png"), width=3000, height=2200, res=400)


	layout(matrix(c(14,11,12,13,
					15,1,2,3,
					16,4,5,6,
					17,7,8,9,
					17,7,18,9,
					10,10,10,10), nrow=6, ncol=4, byrow=T), height=c(0.3,2,2,1,1,0.5), width=c(0.3,5,5,5))


	

	par(mar=c(4, 4, 1, 1))

	# Bayes factor R2
	sub.df = r2.df[r2.df$trait.noise=="none",]
	plot.beta(sub.df, main="", show.ylab=T, sig.threshold=sig.threshold, pvar=pvar, beta.var=beta.var, xlim=c(0, beta.max))
	add.label("A")


	# Bayes factor nobs
	sub.df = nobs.df[nobs.df$trait.noise=="none",]
	plot.nobs(sub.df, main="", show.ylab=T, sig.threshold=sig.threshold,  pvar=pvar, beta.var=beta.var)
	add.label("B")


	# Bayes factor nsubst
	sub.df = subst.df[subst.df$trait.noise=="none",]
	plot.subst(sub.df, main="", show.ylab=T, sig.threshold=sig.threshold,  pvar=pvar, beta.var=beta.var)
	add.label("C")




	# p-value R2
	sub.df = r2.df[r2.df$trait.noise=="none",]
	plot.beta.pval(sub.df)
	add.label("D")

	# p-value nobs
	sub.df = nobs.df[nobs.df$trait.noise=="none",]
	plot.n.pval(sub.df)
	add.label("E")

	# p-value subst
	sub.df = subst.df[subst.df$trait.noise=="none",]
	plot.subst.pval(sub.df)
	add.label("F")


	# Histograms
	par(mar=c(2, 4, 1, 1))


	# R2 hist
	bb = r2.df[,beta.var]
	bb = bb[bb < beta.max]
	h = hist(bb, breaks=nbreaks, plot=F)
	cols = ifelse(h$mids < min.r2, hist.col1, hist.col2)
	hist(bb, main="", xlab="", probability=T, axes=F, yaxs="i", ylab="", breaks=nbreaks, xlim=c(0,beta.max), col=cols)
	x1 = h$mids[h$mids >= min.r2][2]
	y1 = h$density[h$mids >= min.r2][1]
	text(x1, y1, "Cutoff — other columns are\nconditioned on this range", adj=c(0,0.5), cex=0.75)
	axis(1)
	add.label("G", is.hist=TRUE)



	# N hist
	par(mar=c(2, 4, 1, 1))
	nobs.df1 = nobs.df[!is.na(nobs.df$ls.p),]
	nn = nobs.df1[,nvar]
	nn = nn[nn <= n.max]
	h = hist(nn, breaks=nbreaks,  plot=F)
	cols = ifelse(h$mids < min.obs, hist.col1, hist.col2)
	hist(nn,  xlab="", probability=T, axes=F, yaxs="i", ylab="", breaks=nbreaks, main="", xlim=c(0,n.max), col=cols)
	axis(1, at=c(0, 100, 200, 300, 400))
	mtext("Taxa (OLS and PGLS)", cex=0.6)
	add.label("H", is.hist=TRUE)



	par(mar=c(2, 4, 1, 1))


	# Subst hist
	ss = log(subst.df$number.of.subst, 10)
	ss = ss[ss >= 2 & ss <= 7]
	cutoff = log(min.subst, 10)
	h = hist(ss, breaks=nbreaks, plot=F)
	cols = ifelse(h$mids < cutoff, hist.col1, hist.col2)
	hist(ss, main="", xlab="", probability=T, axes=F, yaxs="i", ylab="", breaks=nbreaks, xlim=c(2,7), col=cols)
	#axis(1, at=c(1,2,3,4,5,6,7), labels=c("10¹", "10²", "10³", "10⁴", "10⁵", "10⁶", "10⁷"))
	axis(1, at=c(2:7), labels=c("10²", "10³", "10⁴", "10⁵", "10⁶", "10⁷"))
	add.label("J", is.hist=TRUE)
	


	# Legend
	par(mar=c(0,0,0,0))
	plot(0,0, type="n", axes=F, xaxs="i", yaxs="i", xlab="", ylab="", main="", xlim=c(0,1), ylim=c(0,1))

	x1=0.05
	x2=0.3
	x3=0.55
	x4=0.8
	y1=0.5

	dy=-0.3
	marx1 = 0.05
	marx2 = 0.06


	text(x1, 0.8, "Lines of best fit", adj=0, cex=0.8)


	# PPC
	y = y1
	x = x1
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["ppcNst"], "44"), lwd=4.5)
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["ppcNst"], "88"), lwd=3.5)
	lines(c(x, x+marx1), c(y,y), col=model.cols["ppcNst"], lwd=3)
	text(x+marx2, y, "PPC +nst", adj=0)



	# PPC max
	y = y + dy
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["ppcNstMax"], "44"), lwd=4.5)
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["ppcNstMax"], "88"), lwd=3.5)
	lines(c(x, x+marx1), c(y,y), col=model.cols["ppcNstMax"], lwd=3)
	text(x+marx2, y, "PPC +nst +max", adj=0)



	# OLS# PGLS brownian
	y = y1
	x = x2
	
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["ols"], "44"), lwd=4.5)
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["ols"], "88"), lwd=3.5)
	lines(c(x, x+marx1), c(y,y), col=model.cols["ols"], lwd=3)
	text(x+marx2, y, "OLS", adj=0)



	
	y = y + dy
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["pgls"], "44"), lwd=4.5)
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["pgls"], "88"), lwd=3.5)
	lines(c(x, x+marx1), c(y,y), col=model.cols["pgls"], lwd=3)
	text(x+marx2, y, "PGLS Brownian", adj=0)


	y = y1
	x = x3

	# PGLS pagel
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["lambda"], "44"), lwd=4.5)
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["lambda"], "88"), lwd=3.5)
	lines(c(x, x+marx1), c(y,y), col=model.cols["lambda"], lwd=3)
	text(x+marx2, y, "PGLS Pagel", adj=0)


	# PGLS grafen
	y = y + dy
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["grafen"], "44"), lwd=4.5)
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["grafen"], "88"), lwd=3.5)
	lines(c(x, x+marx1), c(y,y), col=model.cols["grafen"], lwd=3)
	text(x+marx2, y, "PGLS Grafen", adj=0)


	x = x4
	y = y1

	# PGLS martins
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["martins"], "44"), lwd=4.5)
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["martins"], "88"), lwd=3.5)
	lines(c(x, x+marx1), c(y,y), col=model.cols["martins"], lwd=3)
	text(x+marx2, y, "PGLS Martins", adj=0)


	# PGLS blomberg
	y = y + dy
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["blomberg"], "44"), lwd=4.5)
	lines(c(x, x+marx1), c(y,y), col=paste0(model.cols["blomberg"], "88"), lwd=3.5)
	lines(c(x, x+marx1), c(y,y), col=model.cols["blomberg"], lwd=3)
	text(x+marx2, y, "PGLS Blomberg", adj=0)




	print("X")

	par(mar=c(0,4,0,1))
	write.label("Effect size β", cex=1.2)
	write.label("Number of taxa or taxon-pairs", cex=1.2)
	write.label("Number of substitutions", cex=1.2)

	par(mar=c(0,0,0,0))
	write.label("", cex=1.0)


	print("Y")
	par(mar=c(2,0,0.5,0))
	if (show.bf){
		ylab1 = "Bayes factor log₁₀ K"
	}else{
		ylab1 = "Posterior probability"
	}


	write.label(ylab1, cex=1.2, srt=90)
	write.label(paste0("Power (p<", sig.threshold, ")"), cex=1.2, srt=90)
	write.label("Probability density", cex=1.2, srt=90)


	# N hist
	par(mar=c(2, 4, 1, 1))
	nobs.df1 = nobs.df[is.na(nobs.df$ls.p),]
	nn = nobs.df1[,nvar]
	nn = nn[nn <= n.max]
	h = hist(nn, breaks=nbreaks,  plot=F)
	cols = ifelse(h$mids < min.obs, hist.col1, hist.col2)
	hist(nn,  xlab="", probability=T, axes=F, yaxs="i", ylab="", breaks=nbreaks, main="", xlim=c(0,n.max), col=cols)
	axis(1, at=c(0, 100, 200, 300, 400))
	mtext("Taxon pairs (PPC)", cex=0.6)
	add.label("I", is.hist=TRUE)


	dev.off()



}


plot.all(min.obs=min.obs, min.subst=min.subst, min.r2=min.r2, sig.threshold=sig.threshold, poly=T, beta.var=beta.var, pvar=pvar)




