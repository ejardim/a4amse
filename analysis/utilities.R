#====================================================================
# 20260429EJ
# msemodules
#====================================================================

#--------------------------------------------------------------------
# oem: generates catches and index in the data year to be used as
# input in the assessment
#--------------------------------------------------------------------

sca.oem <- function(stk, deviances, observations, args, tracking) {

  # TODO needs more work to remove the index OE, for now index OE is mandatory

	#dataYears <- 1:(args$ay-args$y0-args$data_lag+1)
	dataYears <- args$y0:args$dy
	mxy <- ac(max(dataYears))
	assessmentYear <- ac(args$ay)

	# catch.n
	# note it's adding 1 individual to avoid sca from crashing
	catch.n(observations$stk)[,mxy] <- catch.n(stk)[,mxy] * deviances$stk$catch.n[,mxy] + 1
	catch(observations$stk)[,mxy] <- computeCatch(observations$stk[,mxy])
	stk0 <- observations$stk[,ac(dataYears)]

	# indices
	idx0 <- observations$idx
	for (idx_count in 1:length(observations$idx)){
		TS <- mean(range(observations$idx[[idx_count]])[c("startf", "endf")])
		ages <- dimnames(observations$idx[[idx_count]])$age
		i0 <- (stock.n(stk)[,mxy] * exp((-m(stk)[,mxy] - harvest(stk)[,mxy]) * TS))[ages]
		i0 <- i0 * deviances$idx[[idx_count]][,mxy]
		if(any(i0==0)) i0[i0==0] <- min(i0[i0>0])/2
		index(observations$idx[[idx_count]])[,mxy] <- i0
		idx0[[idx_count]] <- observations$idx[[idx_count]][,ac(range(observations$idx[[idx_count]])['minyear']:mxy)]
	}

	# return
	list(stk=stk0, idx=idx0, observations=observations, tracking=tracking)
}
