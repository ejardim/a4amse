#====================================================================
# 20260429EJ
#====================================================================

# load libraries and functions
library(FLa4a)
library(FLasher)
library(FLBRP)
library(ggplotFL)
library(TAF)
library(mse)
library(msemodules)
source("utilities_fsquared.R")

# load data
# stock
load("../data/HKE_1_5_6_7_STK.Rdata")
stk <- hke.stk
# index
load("../data/HKE_1_5_6_7_IDX.Rdata")
idx <- FLIndices(idx=hke.idx)

#====================================================================
# SETUP
#====================================================================

# Name of the stock
stkname <- name(stk)
# TODO: Recruitment models to be used in the OM conditioning
srmodels <- c("bevholt") # segreg, bevholt, ricker
# Initial year of projections
iy <- dims(stk)$maxyear
# First year of data
y0 <- dims(stk)$minyear
# Years to be used to compute SPR0 for stock-recruitment model
spryrs <- seq(dims(stk)$minyear, iy)
# Data lag
dl <- 1
# Management lag
ml <- 1 
# Assessment frequency
af <- 1
# Data year
dy <- iy - dl
# Final year
fy <- iy + 10
# Years to compute probability metrics
pys <- seq(fy - 5, fy)
# How many years from the past to condition the future
conditioning_ny <- 5
# CV for SSB to add uncertainty in the shortcut estimator
bcv_sa <- 0.5
# CV for F to add uncertainty in the shortcut estimator
fcv_sa <- 0.5
# Years for geometric mean in short term forecast
recyrs_mp <- -2
# TODO: Blim and Btrigger
Blim <- 34788
Btrigger <- 48340
refpts <- FLPar(c(Blim = Blim, Btrigger = Btrigger))
# TODO: no. of cores to use in parallel, defauls to 2/3 of those in machine
cores <- 2#floor(availableCores() * 0.5)
# TODO: F search grid
fg_mp <- seq(0, 3, length=9)
# Number of iterations (minimum of 25 for testing, 500 for final)
it <- 2
# Random seed
set.seed(987)

# PARALLEL setup via doFuture
if(os.linux()) {
  plan(multicore, workers=cores)
} else {
  plan(multisession, workers=cores)
}

options(doFuture.rng.onMisuse="ignore")

#====================================================================
# OM conditioning
#====================================================================

#--------------------------------------------------------------------
# assessment
#--------------------------------------------------------------------

# survey catchability submodel
qmod <- list(~ I(1/(1+exp(-age))))

# fishing mortality submodel
fmod <- ~s(age, k = 4) + s(year, k = 8) + te(age, year, k = c(3, 10))

# fit
fit <- sca(hke.stk, hke.idx, fmodel=fmod, qmodel=qmod)

# om with uncertainty
om <- hke.stk + simulate(fit, it)

# Stock-recruitment relationship(s)
om.sr <- fmle(as.FLSR(om, model="bevholt"), control = list(trace = 0))
om.srdevs <- rlnormar1(it, sdlog=sd(residuals(om.sr)), years=seq(dy, fy))

# BRP
om.brp <- brp(FLBRP(om, sr=om.sr))

#--------------------------------------------------------------------
# BUILD FLom, OM FLR object
#--------------------------------------------------------------------
om <- FLom(stock=om, refpts=refpts(om.brp), model="bevholt",
  params=params(om.sr), deviances=om.srdevs, name=stkname)

# SETUP om future: average of most recent years set by conditioning_ny
om <- fwdWindow(om, end=fy, nsq=conditioning_ny)

#====================================================================
# OEM
# Needs to create objects with the same dimensions as the OM
#====================================================================
#--------------------------------------------------------------------
# deviances for indices using q estimated by the model
#--------------------------------------------------------------------
idcs <- FLIndices()
for (i in 1:length(idx)){
	i.q0 <- predict(fit)$qmodel[[i]]
	i.q <- window(i.q0, end=fy)
	i.q[,ac((iy):fy)] <- i.q[,ac(dy)]
	i.fit <- window(index(fit)[[i]], end=fy)
	idx_temp <- FLIndex(index=i.fit, index.q=i.q)
	range(idx_temp)[c("startf", "endf")] <- range(idx[[i]])[c("startf", "endf")]
	idcs[[i]] <- idx_temp
}
names(idcs) <- names(idx)

#--------------------------------------------------------------------
# deviances for catches
#--------------------------------------------------------------------

# create object from OM structure
catch.dev <- catch.n(stock(om))
# get observation error variance from fit assuming fixed over time
catch.dev[] <- yearMeans(sqrt(predict(fit)$vmodel$catch))
# randomize
catch.dev <- rnorm(it, 0, sd=catch.dev)

#--------------------------------------------------------------------
# build OEM object
#--------------------------------------------------------------------

idxDev <- lapply(idcs, index.q)
names(idxDev) <- c("index.q", "index.q")
stkDev <- FLQuants(catch.n=catch.dev)

# deviances
dev <- list(idx=idxDev, stk=stkDev)
# observations
# WARNING: note we're selecting one index only
obs <- list(idx=idcs, stk=stk)

# OEM
oem <- FLoem(method=mse::sampling.oem, observations=obs, deviances=dev)

#====================================================================
# MP
#====================================================================

# SET intermediate year + start of stks, lags and frequency
mseargs <- list(iy=iy, fy=fy, data_lag=dl, management_lag=ml, frq=af)

# SETUP standard ICES advice rule
arule <- mpCtrl(

# (est)imation method: shortcut.sa + SSB deviances
  est = mseCtrl(method=sca.sa, args=list(fmodel=fmod, qmodel=qmod, update=FALSE)),

  # hcr: hockeystick (fbar ~ ssb | lim, trigger, target, min)
  hcr = mseCtrl(method=hockeystick.hcr,
    args=list(lim=0, trigger=0, target=0.5, min=0,
    metric="ssb", output="fbar")),

  # (i)mplementation (sys)tem: tac.is (C ~ F)
  isys = mseCtrl(method=tac.is, args=list(recyrs=recyrs_mp))
)

#====================================================================
# Run simulations
#====================================================================

tes0 <- mp(om, oem, ctrl=arule, args=mseargs)

# RUN over Ftarget grid
fgrid <- mps(om, ctrl=arule, args=mseargs, hcr=list(target=fg_mp),
  names=paste0("F", fg_mp))

# PLOT
plot(om, fgrid)
