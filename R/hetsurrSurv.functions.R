#estimates P(T>t|W=w)
pred.smooth.surv.w <- function(x.ref, delta.ref, w.ref, w.apply, myt, extrapolate = T, h.use = NULL, warn.support = FALSE)
  { 
    if(is.null(h.use)) {
    		bwini = bw.nrd(w.ref)
    		n.s = length(w.ref)
    		bw <- bwini/(n.s^0.1)

    	}
    	else {bw=h.use}
    kerni.ss = Kern.FUN(zz=w.ref,zi=w.apply,bw)           
    tmpind = (x.ref<=myt)&(delta.ref==1); tj = x.ref[tmpind]; 
    kerni.1 = t(t(kerni.ss))
    pihamyt0.tj.ss = helper.si(tj, "<=", x.ref, Vi=t(kerni.1)) ## n.tj x n.ss matrix ##   
    dLamhat.tj.ss = t((kerni.1[,tmpind]))/pihamyt0.tj.ss; 
    ret = apply(dLamhat.tj.ss,2,sum)
    Phat.ss  =exp(-ret)
    if(sum(is.na(Phat.ss))>0 & extrapolate){
    	if(!warn.support) {print(paste("Note: ", sum(is.na(Phat.ss)), " values extrapolated."))}
    	c.mat = cbind(w.apply, Phat.ss)
    	for(o in 1:length(Phat.ss)) {
    		if(is.na(Phat.ss[o])){
    			distance = abs(w.apply - w.apply[o])
    			c.temp = cbind(c.mat, distance)
    			c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where predication is not na
    			new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    			Phat.ss[o] = new.est[1]   #in case there are multiple matches
    	}
    }
	}
    return(Phat.ss)
    }


#estimates P(T>t|T>t_0, S=s, W=w)
#only send reference data that is already conditioned on X>t_0
pred.smooth.surv.w.s <- function(x.ref, delta.ref, w.ref, s.ref, w.apply, s.apply, h.s, h.w, myt, extrapolate = T, kerni.ss.s, tmpind, tj)
  { 
        
    kerni.ss.w = Kern.FUN(zz=w.ref,zi=w.apply,h.w)
    kerni.ss = t(as.vector(kerni.ss.w)*t(kerni.ss.s))           
    
    kerni.1 = t(t(kerni.ss))
    pihamyt0.tj.ss = helper.si(tj, "<=", x.ref, Vi=t(kerni.1)) ## n.tj x n.ss matrix ##   
    dLamhat.tj.ss = t((kerni.1[,tmpind]))/pihamyt0.tj.ss; 
    ret = apply(dLamhat.tj.ss,2,sum)
    Phat.ss  =exp(-ret)
  
    
    if(sum(is.na(Phat.ss))>0 & extrapolate){
    	print(paste("Note: ", sum(is.na(Phat.ss)), " values extrapolated."))
    	c.mat = cbind(s.apply, Phat.ss)
    	for(o in 1:length(Phat.ss)) {
    		if(is.na(Phat.ss[o])){
    			distance = abs(s.apply - s.apply[o])
    			c.temp = cbind(c.mat, distance)
    			c.temp = c.temp[!is.na(c.temp[,2]),]  #all rows where prediction is not na
    			new.est = c.temp[c.temp[,3] == min(c.temp[,3]), 2]
    			Phat.ss[o] = new.est[1]   #in case there are multiple matches
    	}
    }
	}
    return(Phat.ss)
    }

R.surv.s.w <- function(xone, xzero, deltaone, deltazero, sone, szero, wone, wzero, w.grd, myt, landmark, extrapolate = T,  h.0=NULL,h.1=NULL, h.w=NULL, h.s=NULL,h.w.1=NULL){
	
	#delta_w
	tau.1.w = pred.smooth.surv.w(x.ref=xone, delta.ref = deltaone, w.ref=wone, w.apply = w.grd, myt=myt,  h.use = h.1)
	tau.0.w = pred.smooth.surv.w(x.ref=xzero, delta.ref = deltazero, w.ref=wzero, w.apply = w.grd, myt=myt, h.use = h.0)
	delta.w = tau.1.w - tau.0.w
	S.1= tau.1.w
	S.0 = tau.0.w

	#delta.s.w
	
	if(is.null(h.w)) {h.w = 3*bw.nrd(wone[xone>landmark])/((length(wone[xone>landmark]))^(0.1))}
	if(is.null(h.s)) {h.s = 3*bw.nrd(sone[xone>landmark])/((length(sone[xone>landmark]))^(0.1))}
	if(is.null(h.w.1)) {h.w.1 = 3*bw.nrd(wzero[xzero>landmark])/((length(wzero[xzero>landmark]))^(0.1))}

	tau.10.w = pred.smooth.surv.w(x.ref=xzero, delta.ref = deltazero, w.ref=wzero, w.apply = w.grd, myt=landmark, h.use = h.1)

	x.ref=xone[xone>landmark]; delta.ref = deltaone[xone>landmark]; w.ref=wone[xone>landmark];s.ref=sone[xone>landmark]
	
	id=order(x.ref)
    x.ref=x.ref[id]
    delta.ref=delta.ref[id]
    w.ref=w.ref[id]
    s.ref=s.ref[id]
    
    tmpind = (x.ref<=myt & delta.ref==1) 
    tj = x.ref[tmpind] 
    kerni.ss.s = Kern.FUN(zz=s.ref,zi=szero[xzero>landmark],h.s)
    
	tau.0.tsw  = sapply(w.grd,pred.smooth.surv.w.s,x.ref=x.ref, delta.ref = delta.ref, w.ref=w.ref, s.ref=s.ref, s.apply = szero[xzero>landmark], myt=myt, h.s = h.s, h.w = h.w, extrapolate=T, kerni.ss.s = kerni.ss.s, tmpind=tmpind, tj=tj)
	      
	k.weight=t(Kern.FUN(wzero, w.grd, h.w.1))

    int.v.s = apply(tau.0.tsw*k.weight[xzero>landmark,],2,sum) /apply(k.weight[xzero>landmark,],2,sum)*tau.10.w 

	
	delta.s.w = int.v.s - tau.0.w
	R.s.w = 1-delta.s.w/delta.w
	return(list("R.s.w" = R.s.w, "delta.s.w" = delta.s.w, "delta.w"=delta.w,"S.1" = S.1, "S.0"=S.0))
	
}

R.main.estimate <- function(xone, xzero, deltaone, deltazero, sone, szero, wone, wzero, w.grd, myt, landmark, type = "cont", var= FALSE, test = FALSE, extrapolate = T,  h.0=NULL,h.1=NULL, h.w=NULL, h.s=NULL,h.w.1=NULL){
	
	if(type == "cont") {
	R.f = R.surv.s.w(xone=xone, xzero=xzero, deltaone=deltaone, deltazero=deltazero, sone=sone, szero=szero, wone=wone, wzero=wzero, w.grd=w.grd, myt=myt, landmark=landmark)
	
	results.list = list("w.values" = w.grd, "R.s.w" = R.f$R.s.w, "delta.w"=R.f$delta.w, "delta.s.w" = R.f$delta.s.w)
	
	if(var | test){
		my.delta = function(w.want, thegrid = w.grd, thedelta){
		rrr = approx(thegrid, thedelta, xout = w.want, method = "linear")
		return(rrr$y)
		}	
		r.term = 1-integrate(my.delta, thedelta = R.f$delta.s.w, min(w.grd), max(w.grd))$value/integrate(my.delta, thedelta = R.f$delta.w, min(w.grd), max(w.grd))$value
		r.d = R.f$R.s.w - r.term

		#BOOTSTRAP
		b.num = 300
		R.mat.boot = matrix(nrow = b.num, ncol=length(w.grd))
		delta.mat.boot = matrix(nrow = b.num, ncol=length(w.grd))
		delta.s.mat.boot = matrix(nrow = b.num, ncol=length(w.grd))
		r.d.boot = matrix(nrow = b.num, ncol=length(w.grd))
		n1 = length(xone)
		n0 = length(xzero)
		for(uuu in 1:b.num) {
			ind.boot.1 = sample(1:n1, n1, replace = T)
			ind.boot.0 = sample(1:n0, n0, replace = T)
			R.b = R.surv.s.w(xone=xone[ind.boot.1], xzero=xzero[ind.boot.0], deltaone=deltaone[ind.boot.1], deltazero=deltazero[ind.boot.0], sone=sone[ind.boot.1], szero=szero[ind.boot.0], wone=wone[ind.boot.1], wzero=wzero[ind.boot.0], w.grd=w.grd, myt=myt, landmark=landmark, extrapolate = T)
			R.mat.boot[uuu,] = R.b$R.s.w
			delta.mat.boot[uuu,] = R.b$delta.w
			delta.s.mat.boot[uuu,] = R.b$delta.s.w
		
			r.d.boot[uuu,] = rep(NA, length(r.d.boot[uuu,]))
			tryCatch({
			r.term.b = 1-integrate(my.delta, thedelta = R.b$delta.s.w, min(w.grd), max(w.grd))$value/integrate(my.delta, thedelta = R.b$delta.w, min(w.grd), max(w.grd))$value
			r.d.b = R.b$R.s.w - r.term.b
			r.d.boot[uuu,] = r.d.b
			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
		}	
		r.d.boot = r.d.boot[!is.na(r.d.boot[,1]),]
		b.num = sum(!is.na(r.d.boot[,1]))	
		sd.R.p = apply(R.mat.boot, 2, mad) 
		sd.r.d = apply(r.d.boot, 2, mad)
	
		delta.s.sd = apply(delta.s.mat.boot, 2, mad)
		delta.sd = apply(delta.mat.boot, 2, mad)
		
		results.list = list("w.values" = w.grd, "R.s.w" = R.f$R.s.w, "delta.w"=R.f$delta.w, "delta.s.w" = R.f$delta.s.w, "sd.R" = as.vector(sd.R.p), "sd.delta" = as.vector(delta.sd), "sd.delta.s" = as.vector(delta.s.sd))
	}
	
	if(test){
		tstar = abs(r.d/sd.r.d)
		t.est = max(tstar)
	
		aa = covRob(r.d.boot, corr=TRUE,estim="mcd")$cov
		ww=mvrnorm(1000, mu = rep(0,length(w.grd)),aa)
		max.this = apply(abs(ww), 1, max)
		pval = mean(max.this >= t.est)

		ww=mvrnorm(1000, mu = rep(0,length(w.grd)),diag(1, nrow=length(w.grd)))
		max.this = apply(abs(ww), 1, max)
		pval.con = mean(max.this >= t.est)
		results.list = c(results.list, "pval.omnibus" = pval, "pval.con.omnibus" = pval.con)
	}
	}
	if(type == "discrete"){
	if(is.null(w.grd)){
			w.grd = sort(unique(c(wone, wzero)))
		}
		delta.diff.w = vector(length = length(w.grd))
		delta.s.diff.w = vector(length = length(w.grd))
    	for(jj in 1:length(w.grd)) {
		delta.diff.w[jj] = delta.surv.estimate(xone[wone==w.grd[jj]], xzero[wzero==w.grd[jj]], deltaone[wone==w.grd[jj]], deltazero[wzero==w.grd[jj]], t = myt)$delta
		delta.s.diff.w[jj] = delta.s.surv.estimate(xone[wone==w.grd[jj]], xzero[wzero==w.grd[jj]], deltaone[wone==w.grd[jj]], deltazero[wzero==w.grd[jj]], sone[wone==w.grd[jj]], szero[wzero==w.grd[jj]], landmark=landmark, t = myt)
    	}
		
		R.s.w = 1-delta.s.diff.w/delta.diff.w 	
		if(var | test) {
			ll.1 = length(xone)
			ll.0 = length(xzero)
			boot.num = 500
			num.w = length(w.grd)
    		boot.mat = matrix(nrow = boot.num,ncol = length(w.grd)*3)
				for(kk in 1:boot.num) {
					index.boot.1 = sample(1:ll.1, ll.1, replace = T) 
					index.boot.0 = sample(1:ll.0, ll.0, replace = T) 
					s1.b = sone[index.boot.1]
					x1.b = xone[index.boot.1]
					delta1.b = deltaone[index.boot.1]
					u1.b = wone[index.boot.1]
					
					s0.b = szero[index.boot.0]
					x0.b = xzero[index.boot.0]
					delta0.b = deltazero[index.boot.0]
					u0.b = wzero[index.boot.0]

					#			delta as a function of w
					delta.diff.w.boot = vector(length = length(w.grd))
					delta.s.diff.w.boot = vector(length = length(w.grd))
    				for(jj in 1:length(w.grd)) {
						delta.diff.w.boot[jj] = delta.surv.estimate(x1.b[u1.b==w.grd[jj]], x0.b[u0.b==w.grd[jj]], delta1.b[u1.b==w.grd[jj]], delta0.b[u0.b==w.grd[jj]], t = myt)$delta
						delta.s.diff.w.boot[jj] = delta.s.surv.estimate(x1.b[u1.b==w.grd[jj]], x0.b[u0.b==w.grd[jj]], delta1.b[u1.b==w.grd[jj]], delta0.b[u0.b==w.grd[jj]], s1.b[u1.b==w.grd[jj]], s0.b[u0.b==w.grd[jj]], landmark=landmark, t = myt)
    				}
					boot.mat[kk,1:num.w] = delta.diff.w.boot
 					boot.mat[kk,(num.w+1):(num.w+num.w)] = delta.s.diff.w.boot
					boot.mat[kk,(num.w*2+1):(num.w*3)] = 1-delta.s.diff.w.boot/delta.diff.w.boot	
		
				} #closes for loop
			var.results.mat = diag(var(boot.mat))
			se.delta.w = sqrt(var.results.mat[1:num.w])
			se.delta.w.s = sqrt(var.results.mat[(num.w+1):(num.w+num.w)])
			se.R.w.s = sqrt(var.results.mat[(num.w*2+1):(num.w*3)])
			#testing
			#make contrast matrix
			cont = diag(1,nrow = length(w.grd)-1, ncol = length(w.grd)) + cbind(rep(0,length(w.grd)-1), diag(-1,nrow = length(w.grd)-1, ncol = length(w.grd)-1))
			vec.R = R.s.w
			delta.test = cont %*% as.matrix(vec.R)
			sand = solve(cont %*% var(boot.mat[,(length(w.grd)*2+1):(length(w.grd)*3)]) %*% t(cont))
			G = t(delta.test) %*% sand %*% delta.test
			test.stat = as.numeric(G)
			discrete.p.value = as.numeric(1-pchisq(G, num.w-1))
		}
		results.list = list("w.values" = w.grd, "R.w.s" = R.s.w, "delta.w" = delta.diff.w, "delta.s.w" = delta.s.diff.w)
		if(var) {results.list = list("w.values" = w.grd, "R.w.s" = R.s.w, "delta.w" = delta.diff.w, "delta.s.w" = delta.s.diff.w,"sd.R" = se.R.w.s, "sd.delta" = se.delta.w, "sd.delta.s" = se.delta.w.s)}
		if(test) {results.list = c(results.list, "pval.discrete" = discrete.p.value)}
	} #close discrete
		
	return(results.list)
}

test.multiplet <- function(t.mult, xone, xzero, deltaone, deltazero, sone, szero, wone, wzero, w.grd, landmark, extrapolate = T,  h.0=NULL,h.1=NULL, h.w=NULL, h.s=NULL,h.w.1=NULL){

	my.delta = function(w.want, thegrid = w.grd, thedelta){
		rrr = approx(thegrid, thedelta, xout = w.want, method = "linear")
		return(rrr$y)
		}	
	big.delta = rep(0, length(w.grd))	
	big.delta.s= rep(0, length(w.grd))
	for(kk in 1:length(t.mult)){
		R.f = R.surv.s.w(xone=xone, xzero=xzero, deltaone=deltaone, deltazero=deltazero, sone=sone, szero=szero, wone=wone, wzero=wzero, w.grd=w.grd, myt=t.mult[kk], landmark=landmark)
		big.delta = big.delta + R.f$delta.w
		big.delta.s = big.delta.s + R.f$delta.s.w
	}
	R.bigt = 1-(big.delta.s)/(big.delta)
	r.term = 1-integrate(my.delta, thedelta = big.delta.s, min(w.grd), max(w.grd))$value/integrate(my.delta, thedelta = big.delta, min(w.grd), max(w.grd))$value
	r.d.bigt = R.bigt  - r.term
	
	#BOOTSTRAP
	n1 = length(xone)
	n0=length(xzero)
	b.num = 300
	delta.bigt.boot = matrix(0,nrow = b.num, ncol=length(w.grd))
	delta.s.bigt.boot = matrix(0,nrow = b.num, ncol=length(w.grd))
	d.bigt.boot = matrix(nrow = b.num, ncol=length(w.grd))
	for(uuu in 1:b.num) {
		ind.boot.1 = sample(1:n1, n1, replace = T)
		ind.boot.0 = sample(1:n0, n0, replace = T)
		for(qq in 1:length(t.mult)) {
			t.use=t.mult[qq]
			R.b = R.surv.s.w(xone=xone[ind.boot.1], xzero=xzero[ind.boot.0], deltaone=deltaone[ind.boot.1], deltazero=deltazero[ind.boot.0], sone=sone[ind.boot.1], szero=szero[ind.boot.0], wone=wone[ind.boot.1], wzero=wzero[ind.boot.0], w.grd=w.grd, myt=t.use, landmark=landmark, extrapolate = T)
			delta.bigt.boot[uuu,] = delta.bigt.boot[uuu,] + R.b$delta.w
			delta.s.bigt.boot[uuu,] = delta.s.bigt.boot[uuu,] + R.b$delta.s.w
		}
		R.bigt.b = 1-delta.s.bigt.boot[uuu,]/delta.bigt.boot[uuu,]
		tryCatch({
		r.bigt.term.b = 1-integrate(my.delta, thedelta = delta.s.bigt.boot[uuu,], min(w.grd), max(w.grd))$value/integrate(my.delta, thedelta = delta.bigt.boot[uuu,], min(w.grd), max(w.grd))$value
		r.d.b = R.bigt.b - r.bigt.term.b
		d.bigt.boot[uuu,] = r.d.b
		}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
	}	
	d.bigt.boot = d.bigt.boot[!is.na(d.bigt.boot[,1]),]
	b.num = sum(!is.na(d.bigt.boot[,1]))	
	sd.d.bigt.b = apply(d.bigt.boot, 2, mad)

	tstar = abs(r.d.bigt/sd.d.bigt.b)
	t.est = max(tstar)
	
	aa = covRob(d.bigt.boot, corr=TRUE,estim="mcd")$cov
	ww=mvrnorm(1000, mu = rep(0,length(w.grd)),aa)
	max.this = apply(abs(ww), 1, max)
	pval.mult = mean(max.this >= t.est)
	
	ww=mvrnorm(1000, mu = rep(0,length(w.grd)),diag(1, nrow=length(w.grd)))
	max.this = apply(abs(ww), 1, max)
	pval.mult.con = mean(max.this >= t.est)

return(list("pval.multi" = pval.mult, "pval.con.multi" = pval.mult.con))
}
