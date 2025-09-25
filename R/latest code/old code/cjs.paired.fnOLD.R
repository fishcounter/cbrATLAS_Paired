# estimate paired CJS w/o adjusting for tag life.  Censoring may or may not be included in history
# INPUT
# detects: [list of 2]
# detect.in: [0|1|2] detection history in wide ATlAS format, w/o rel.group or bin columns
# L.in: vector of mean probabilities tag is working, given detection at each site.
#    if [L.in = NULL]: function will estimate the CJS parameters w/o tag failure (Li = 1)
# seeds.in: optional starting estimates for CJS model, includes additional p(censoring) per site at end 
#    parameters are Si, pi, li where i is 1:(total number of sites-1)
#    di (p(censoring)), is not estimated in this version
# se.out: (T/F) calculate se(cjs var) from variance matrix of parameters using Hessian?
# common.start.in: [0 to 2*(num.period-1)] indicates location common parameters start.0=all unique
# OUTPUT
# cjs.params: 2-column matrix with CJS estimates and standard errors (if estimated).

cjs.paired.fn=function(detects.in,L.in=NULL,seeds.in=NULL,se.out=F,d=NULL,bootstrap=F,common.start.in=0){
	print("intta cjs.paired.fn")

	if(is.null(seeds.in)){seeds.in=list(NULL,NULL)}else{
		seeds.tmp=list(NULL,NULL)
		seeds.tmp[[1]]=seeds.in[1:(length(seeds.in)/2)]
		seeds.tmp[[2]]=seeds.in[((length(seeds.in)/2)+1):length(seeds.in)]
		seeds.in=seeds.tmp}
	counts.in=list(NULL,NULL)
	obs.hist=list(NULL,NULL)
	d=list(NULL,NULL)
	if(is.null(L.in)){L.in=list(NULL,NULL)}

	for(g in 1:2){
		# call each release separately (1=upper release; 2=lower release)
		# common lower release likelihood inputs will be set to upper release corresponding upper release likelihood inputs 

		detect.in=detects.in[[g]]
		num.period=dim(detect.in)[2]
		if(is.null(L.in[[g]])){L.in[[g]]=rep(1,num.period)} # if no tag failure, L = 100% for all sites
  
		# summary of counts.in/history and unique detection matrix in same order
		counts.in[[g]] = thist0(detect.in)$count
		obs.hist[[g]] = thist0(detect.in)$hist.matrix

		if(is.null(seeds.in[[g]])){
			# seed estimates do not estimate survival/detection with censoring present.
			seeds.in[[g]] = rep(0.5, num.period * 2 - 1)
			names(seeds.in[[g]])=c(paste0(rep(c("S","p"),rep(num.period-1,2)),rep(1:(num.period-1),2),g),paste0("l",g))
		
    
		if (num.period > 2) {
			seeds.in[[g]][num.period] = sum(detect.in[apply(detect.in[,-1]==1, 1, sum) > 0, 1] == 1)/sum(apply(detect.in[, -1] == 1, 1, sum) > 0)
			for (i in (2:(num.period - 2))) {
				seeds.in[[g]][num.period + i - 1] = sum(detect.in[apply(detect.in[,-c(1:i)] == 1, 1, sum) > 0, i] == 1)/
				sum(apply(detect.in[,-c(1:i)] == 1, 1, sum) > 0)
				}
			seeds.in[[g]][num.period * 2 - 2] = sum(detect.in[detect.in[,num.period] == 1, num.period - 1] == 1)/sum(detect.in[,num.period] == 1)
			}else{seeds.in[[g]][num.period] = sum(detect.in[detect.in[, 2] == 1,1] == 1)/sum(detect.in[, 2] == 1)}

		seeds.in[[g]][(2 * num.period) - 1] = sum(detect.in[, num.period] == 1)/(sum(detect.in[, num.period - 1] == 1)/seeds.in[[g]][2 * (num.period - 1)])
		seeds.in[[g]][1] = (sum(detect.in[, 1] == 1)/dim(detect.in)[1])/seeds.in[[g]][num.period]
		if (num.period > 2) {
			for (i in 2:(num.period - 1)) {
				seeds.in[[g]][i] = (sum(detect.in[, i] == 1)/seeds.in[[g]][i + num.period - 1])/(sum(detect.in[, i - 1] == 1)/seeds.in[[g]][i + num.period - 2])                
				}
			}
		

		# print("p(censored)")
		d[[g]]=rep(0,num.period-1)
		for(i in (1:(num.period-1))){
			d[[g]][i]= sum(detect.in[,i]==2)/sum(detect.in[,i]>0)
			}
#			print(d[[g]])
		names(d[[g]])=c(paste0(rep(paste0(rep("d",rep(num.period-1)))),rep(g,rep(num.period-1)),c(1:(num.period-1))))
			
			
		}

		# adjust S,lambda by expected tag-life.  L.in
		seeds.in[[g]][c(1:(num.period - 1), 2 * num.period - 1)] = seeds.in[[g]][c(1:(num.period - 1), 2 * num.period - 1)]/L.in[[g]]
		# print("seeds done ")
		}
	## maximum likelihood estimate of CJS parameters
	parscale.in=1e-04
	reltol.in=1e-30
	maxit.in=2e+06

### print inputs 	
#print("entering paired.cjs.lik 1st time");print(paste0(c("common.start.in= ",as.numeric(common.start.in))),sep="")
# if common.start.in is > 0, set lower release to equal upper release starting at first common site to end
	paired.cjs.out = optim(par=unlist(seeds.in), fn=paired.cjs.lik, counts.in = counts.in, num.period=num.period, use.hist=obs.hist,L=L.in,d.in=d,
		common.start=as.numeric(common.start.in), method = "BFGS", hessian = F, control = list(parscale = rep(parscale.in, 2*(num.period * 2 - 1)), 
		reltol = reltol.in,maxit=maxit.in,ndeps=rep(parscale.in, 2*(num.period * 2 - 1))))

	params=paired.cjs.out$par
	if(common.start.in>0){
		common.param=c(seq(1,(length(params)-1),2),seq(2,(length(params)-1),2),length(params))

#	print("back into cjs.paired.fn")
#	print(paired.cjs.out$par)		
#	print("set common estimates")
#	print(paired.cjs.out$par[c(1:(length(paired.cjs.out$par)/2))][!(common.param<(common.start.in+1))])
#
#	paired.cjs.out$par[length(paired.cjs.out$par)/2+c(1:(length(paired.cjs.out$par)/2))][!(common.param<(common.start.in+1))]=
#			paired.cjs.out$par[c(1:(length(paired.cjs.out$par)/2))][!(common.param<(common.start.in+1))] # set 2nd release common params  = 1st release common params
		}  
#print("Done with paired cjs w equal params")		
	paired.cjs.out$par=correct.fn(paired.cjs.out$par)

#print("made it to correction")
	##################

	#print(cjs.out$par)
	if(se.out){
	browser()
		# rerun to get hessian
		paired.cjs.out = optim(paired.cjs.out$par, fn=paired.cjs.lik, counts.in = counts.in, num.period = num.period, use.hist=obs.hist,L=L.in,d.in=d,common.start=common.start.in,
			method = "BFGS", hessian = T, control = list(parscale = rep(parscale.in, 2*(num.period * 2 - 1)), reltol = reltol.in,maxit=maxit.in,ndeps=rep(parscale.in, 2*(num.period * 2 - 1))))
		paired.cjs.out$par=correct.fn(paired.cjs.out$par)

		# calculate standard errors on CJS estimates, based on Hessian
		temp = paired.cjs.out$hessian
		singular = colSums(temp) == 0
		var.vec = rep(0, 2*(num.period * 2 - 1))
		var.vec[!singular] = diag(solve(temp[!singular, !singular]))

		# format estimates into a matrix, 2 columns. 1 row per parameter
		paired.cjs.param = matrix(c(paired.cjs.out$par, sqrt(var.vec)), ncol = 2)

		if(common.start.in>0){  
#			common.param=c(seq(1,(length(params[[1]])-1),2),seq(2,(length(params[[1]])-1),2),length(params[[1]]))
#			print(common.param)
#		paired.cjs.param[dim(paired.cjs.param)[1]/2+c(1:(dim(paired.cjs.param)[1]/2)),][!(common.param<(common.start.in+1)),]=
#				paired.cjs.param[c(1:(dim(paired.cjs.param)[1]/2)),][!(common.param<(common.start.in+1)),] # set 2nd release common params  = 1st release common params
			}  

	# note: after paired.cjs.param	(line 112), estimates match ATLAS, but missing last site SE (rel2 lambda SE)
		row.names(paired.cjs.param)=names(seeds.in)
		dimnames(paired.cjs.param)[[2]] = c("Estimate", "s.e.")
		}else{
			# format estimates into a matrix, 1 columns. 1 row per parameter
			paired.cjs.param = matrix(c(paired.cjs.out$par), ncol = 1)
			row.names(paired.cjs.param)=names(unlist(seeds.in))
#			print("last paired.cjs.param est in cjs.paired.fn")
#			print(paired.cjs.param)
			}
	print("outta cjs.paired.fn")
	print(list(paired.cjs.param=paired.cjs.param, d=d))
	return(list(paired.cjs.param=paired.cjs.param, d=d))
	}


