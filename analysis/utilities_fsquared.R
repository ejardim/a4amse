# utilities.R - Extra functions
# /utilities.R

# Distributed under the terms of the EUPL-1.2

sampling.oem <- function(stk, deviances, observations, args, tracking,
  oe=c("both","index","catch")) {

  # TODO needs more work to remove the index OE, for now index OE is mandatory

	#dataYears <- 1:(args$ay-args$y0-args$data_lag+1)
	dataYears <- args$y0:args$dy
	mxy <- ac(max(dataYears))
	assessmentYear <- ac(args$ay)

	# carry on stock information in the observations for "short-cut" approach
	stock.n(observations$stk)[,assessmentYear] <- stock.n(stk)[,assessmentYear]
	harvest(observations$stk)[,assessmentYear] <- harvest(stk)[,assessmentYear]

	# catch.n
	# note it's adding 1 individual to avoid sca from crashing
	if(any(oe %in% c("both","catch"))){
		catch.n(observations$stk)[,mxy] <- catch.n(stk)[,mxy] *
      deviances$stk$catch.n[,mxy] + 1
		catch(observations$stk)[,mxy] <- computeCatch(observations$stk[,mxy])
		stk0 <- observations$stk[,ac(dataYears)]
	}

	# indices
	if(any(oe %in% c("both","index"))){
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
	}

	# return
	list(stk=stk0, idx=idx0, observations=observations, tracking=tracking)
} # }}}



# ICES performance statistics {{{ ----

icestats <- list(

  # C
  C=list(~C, name="Catch (t)",
    desc="Catch inn tonnes"),

  # L
  L=list(~L, name="Landings (t)",
    desc="Landings in tonnes"),

  # D
  D=list(~L, name="Discards (t)",
    desc="Discards in tonnes"),

  # SB
  SB=list(~SB, name="SB (t)",
    desc="Spawning stock biomass"),

  # R
  R=list(~R, name="Recruits (1e3)",
    desc="Recruitment in numbers"),
  
  # F
  F=list(~F, name="bar(F)",
    desc="Fishing mortality"),

  # cv(C)
  cvC=list(~sqrt(iterVars(C)) / iterMeans(C), name="cv(C)",
    desc="CV of catch per year"),

  # AVVC
  AAVC=list(~abs(C[, -1] - C[, -dim(C)[2]]) / C[, -dim(C)[2]],
    name="AAV(C)", desc="Average annual variability in catch"),

  # IACC
  IACC=list(~100 * (C[, -1] - C[, -dim(C)[2]]) / C[, -dim(C)[2]],
    name="IAC(C)",
    desc="Percentage inter-annual change in catch"),

  # IACB
  IACB=list(~100 * (SB[, -1] - SB[, -dim(SB)[2]]) / SB[, -dim(SB)[2]], name="IAC(B)", desc="Percentage inter-annual change in biomass"),

  # P(SB<SBlim)
  PBlim=list(~iterMeans((SB/Blim) < 1)  , name="P(SB<SB[lim])",
    desc="Probability that spawner biomass is below Blim"),

  # P(SB>SBtrigger)
  PBtrigger=list(~iterMeans((SB/Btrigger) > 1), name="P(SB>B[trigger])",
    desc="Probability that spawner biomass is above Btrigger"),

  # P(SB < SBlim) at least once
  risk2 = list(~iterMeans((SB / Blim) < 1) > 0,
    name="once(P(SB<B[limit]))",
    desc="ICES Risk 2, probability that spawner biomass is above Blim once")
)

# }}}
