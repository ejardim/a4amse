#====================================================================
# 20260429EJ
# msemodules
#====================================================================

#--------------------------------------------------------------------
# oem: generates catches and index in the data year to be used as
# input in the assessment
#--------------------------------------------------------------------

sca.oem <- function(stk, deviances, observations, args, tracking, log=TRUE) {

  # TODO needs more work to remove the index OE, for now index OE is mandatory

	#dataYears <- 1:(args$ay-args$y0-args$data_lag+1)
	dataYears <- args$y0:args$dy
	mxy <- ac(max(dataYears))
	assessmentYear <- ac(args$ay)

	# catch.n
	# note it's adding 1 individual to avoid sca from crashing
	if(isTRUE(log)){
		catch.n(observations$stk)[,mxy] <- catch.n(stk)[,mxy] * exp(deviances$stk$catch.n[,mxy]) + 1
	} else {
		catch.n(observations$stk)[,mxy] <- catch.n(stk)[,mxy] * deviances$stk$catch.n[,mxy] + 1
	}
	catch(observations$stk)[,mxy] <- computeCatch(observations$stk[,mxy])
	stk0 <- observations$stk[,ac(dataYears)]

	# indices
	idx0 <- observations$idx
	for (idx_count in 1:length(observations$idx)){
		TS <- mean(range(observations$idx[[idx_count]])[c("startf", "endf")])
		ages <- dimnames(observations$idx[[idx_count]])$age
		i0 <- (stock.n(stk)[,mxy] * exp((-m(stk)[,mxy] - harvest(stk)[,mxy]) * TS))[ages]
		if(isTRUE(log)){
			i0 <- i0 * index.q(observations$idx[[idx_count]])[,mxy] * exp(deviances$idx[[idx_count]][,mxy])
		} else {
			i0 <- i0 * index.q(observations$idx[[idx_count]])[,mxy] * deviances$idx[[idx_count]][,mxy]
		}
		index(observations$idx[[idx_count]])[,mxy] <- i0
		observations$idx[[idx_count]] <- replaceZeros(observations$idx[[idx_count]])
		idx0[[idx_count]] <- observations$idx[[idx_count]][,ac(range(observations$idx[[idx_count]])['minyear']:mxy)]
	}

	# return
	list(stk=stk0, idx=idx0, observations=observations, tracking=tracking)
}

effort.is <- function(stk, ctrl, Fref=yearMeans(fbar(stk)[, ac(seq(dy - nyears, dy))]),
  nyears=args$nsqy, constraints=missing, fbyfleet=missing, args, tracking) {

  # CHECK ctrl sets F or effort
  if(!ctrl$quant %in% c("fbar", "f", "effort"))
    stop("'effort.is' can only accept ctrl set on 'f', 'fbar' or 'effort'")

  # EXTRACT args
  spread(args)

	# target to reach defined by HCR
	trgt <- ctrl$value
  # multiplier
	mult <- trgt / c(unitMeans(Fref))
	# deals with 0/0
	mult[is.na(mult)] <- 0
    # catch constraints
    if(!missing(constraints)){
      mult[mult<c(1-constraints["low"])] <- c(1-constraints["low"])
      mult[mult>c(1+constraints["upp"])] <- c(1+constraints["upp"])
    }
    ctrl$value <- mult*unitMeans(fbar(stk)[,ac(dy)])
  # new control file, in relative terms
  #ctrl <- fwdControl(year = ctrl$year, quant = ctrl$quant, value = mult, relYear = ctrl$year - data_lag)

  # effort reduction by fleet
    if(!missing(fbyfleet)){
      # something needs to happen
    }

    return(list(ctrl = ctrl, tracking = tracking))

}

tac_sca.is <- function (stk, ctrl, args, output = "catch", recyrs = -2, dtaclow = NA, dtacupp = NA, fmin = 0, reuse = TRUE, initac = metrics(stk, output)[, ac(iy - data_lag)], tracking)
{
    spread(args)
    fut <- fwdWindow(stk, end = mys[length(mys)], nsq = nsqy)
    id <- dimnames(stk)$year
    if (!is.list(recyrs)) {
        recyrs <- list(recyrs)
    }
    for (i in recyrs) {
        if (is(i, "character")) {
            id <- id[!id %in% i]
        }
        else if (all(i < 0)) {
            if (length(i) == 1)
                id <- rev(rev(id)[-seq(abs(i))])
            else id <- rev(rev(id)[i])
        }
        else if (all(i > 0)) {
            id <- rev(rev(id)[seq(abs(i))])
        }
    }
    recyrs <- id
    if (!all(recyrs %in% dimnames(stk)$year)) {
        stop("'recyrs' cannot be found in input stk")
    }
    gmnrec <- exp(yearMeans(log(unitSums(rec(stk))[, recyrs])))
    srr <- predictModel(model = rec ~ a, params = FLPar(a = gmnrec))
    track(tracking, "gmrec.isys", ay) <- gmnrec
    ftar <- c(ctrl$value)
    track(tracking, "fbar.isys", ay) <- ftar
    if (management_lag > 0) {
        fsq <- yearMeans(unitMeans(fbar(stk)[, ac((dy-2):dy)]))
        if (data_lag == 0) {
            fctrl <- fwdControl(list(year = mys, quant = "fbar",
                value = c(ftar)))
        }
        else {
            fctrl <- fwdControl(c(lapply(seq(dy + 1, mys[1] -
                1), function(y) list(year = y, quant = "fbar",
                value = c(fsq))), lapply(mys, function(y) list(year = y,
                quant = "fbar", value = c(ftar)))))
        }
    }
    else {
        fctrl <- fwdControl(list(year = ay, quant = "fbar", value = ftar))
    }
    fut <- fwd(fut, sr = srr, control = fctrl)
    id <- tracking[metric == "rule.hcr" & year == ay, data >
        2] & c(fbar(fut)[, ac(ay + management_lag)] > fmin)
    if (isTRUE(reuse) | toupper(reuse) == "C") {
        TAC <- areaSums(unitSums(expand(catch(fut)[, ac(mys)[1]],
            year = seq(length(mys)))))
    }
    else {
        TAC <- areaSums(unitSums(catch(fut)[, ac(mys)]))
    }
    if (ay == iy)
        prev_tac <- rep(c(initac), length = args$it)
    else prev_tac <- c(tracking[metric == "isys" & year == ay])
    if (!is.na(dtacupp)) {
        iter(TAC, id) <- pmin(c(iter(TAC, id)), prev_tac[id] *
            dtacupp)
    }
    if (!is.na(dtaclow)) {
        iter(TAC, id) <- pmax(c(iter(TAC, id)), prev_tac[id] *
            dtaclow)
    }
    ctrl <- fwdControl(lapply(seq(length(mys)), function(x) list(year = mys[x],
        quant = output, value = TAC[, x])))
    if (management_lag == 0) {
        ctrl <- merge(ctrl, fwdControl(list(year = mys[length(mys)] +
            1, quant = "catch", value = 0)))
    }
    return(list(ctrl = ctrl, tracking = tracking))
}

f.phcr <- function(stk, frp="f0.1", model="missing", interval, args, tracking) {

    # args
    spread(args)

  # RUN brp() with or without SR fit
	if(ay == iy | (ay - iy) %% interval == 0){
		if(!missing(model)){
			sr0 <- fmle(as.FLSR(stk, model=model))
			hcrpars <- refpts(brp(FLBRP(stk, sr0)))[tolower(frp),"harvest"]
		} else {
			hcrpars <- refpts(brp(FLBRP(stk)))[tolower(frp),"harvest"]
		}
	} else {
		hcrpars <- FLPar(tracking[metric=="phcr" & year==ay-1, "data"])
	}
	tracking[metric=="phcr" & year==ay, "data"] <- c(hcrpars)
	list(hcrpars=hcrpars, tracking=tracking)
} # }}}

f.hcr <- function(stk, hcrpars, args, tracking){
    spread(args)
    # rule
	if(!is(hcrpars, "FLQuant"))
    hcrpars <- FLQuant(c(hcrpars), dimnames=list(iter=dimnames(stk@catch)$iter))

  # create control file
    ctrl <- fwdControl(year=seq(ay + management_lag, ay + frq), quant="fbar", value=c(hcrpars))

    tracking[metric=="hcr" & year==ay, "data"] <- c(hcrpars)
  # return
	list(ctrl=ctrl, tracking=tracking)
} # }}}
