#====================================================================
# 20260429EJ
# MSE template for a full feedback model with a4a sca
#====================================================================

# load libraries and functions
library(FLa4a)
library(FLBRP)
library(ggplotFL)
library(mse)
library(msemodules)
source("utilities.R")

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
# Recruitment model to be used in the OM conditioning
srmodel <- "segreg" # segreg, bevholt, ricker
# Initial year of projections
iy <- dims(stk)$maxyear
# First year of data
y0 <- dims(stk)$minyear
# Data lag
dl <- 1
# Management lag
ml <- 1 
# Assessment frequency
af <- 1
# Data year
dy <- iy - dl
# NUmber of years to projections
npy <- 10
# Final year
fy <- iy + npy
# Years to compute probability metrics
pys <- seq(fy - 5, fy-1)
# How many years from the past to condition the future
conditioning_ny <- 5
# Years for geometric mean in short term forecast
recyrs_mp <- -2
# Number of iterations (minimum of 25 for testing, 500 for final)
it <- 1
# Random seed
set.seed(987)

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
fit <- simulate(fit, it)
om <- hke.stk + fit

# Stock-recruitment relationship(s)
om.sr <- fmle(as.FLSR(qapply(om, iterMedians), model=srmodel), control = list(trace = 0))
om.srdevs <- rlnormar1(it, sdlog=sd(residuals(om.sr)), years=seq(dy, fy))

# BRP
om.brp <- brp(FLBRP(om, sr=om.sr))
rp <- remap(refpts(om.brp))

#--------------------------------------------------------------------
# BUILD FLom, OM FLR object
#--------------------------------------------------------------------
om <- FLom(stock=om, refpts=rp, model=srmodel,
  params=params(om.sr), deviances=om.srdevs, name=stkname)

# SETUP om future: average of most recent years set by conditioning_ny
om <- fwdWindow(om, end=fy, nsq=conditioning_ny)

#====================================================================
# OEM
# Needs to create objects with the same dimensions as the OM
#====================================================================
#--------------------------------------------------------------------
# deviances for indices using q vmodel(s)
#--------------------------------------------------------------------
idcs <- FLIndices()
idxDev <- FLQuants()

for (i in 1:length(idx)){
  # catchability predicted by the model
  i.q0 <- predict(fit)$qmodel[[i]]
  # use median, no estimation error
  i.q0[] <- iterMedians(i.q0)
  i.q <- window(i.q0, end=fy)
  i.q[,ac((iy):fy)] <- i.q[,ac(dy)]
  i.fit <- window(index(fit)[[i]], end=fy)
  idx_temp <- FLIndex(index=i.fit, index.q=i.q)
  range(idx_temp)[c("startf", "endf")] <- range(idx[[i]])[c("startf", "endf")]
  idcs[[i]] <- idx_temp
  idxDev[[i]] <- index(idx_temp)
  idxDev[[i]][] <- yearMeans(sqrt(predict(fit)$vmodel[[i+1]]))
  # randomize
  idxDev[[i]] <- rnorm(it, 0, sd=idxDev[[i]])
}
names(idcs) <- names(idxDev) <- names(idx)

#--------------------------------------------------------------------
# deviances for catches using catch vmodel
#--------------------------------------------------------------------
# create object from OM structure
catch.dev <- catch.n(stock(om))
# get observation error variance from fit assuming fixed over time
catch.dev[] <- yearMeans(sqrt(predict(fit)$vmodel$catch))
# randomize
catch.dev <- rnorm(it, 0, sd=catch.dev)
stkDev <- FLQuants(catch.n=catch.dev)

#--------------------------------------------------------------------
# build OEM object
#--------------------------------------------------------------------
# deviances
dev <- list(idx=idxDev, stk=stkDev)
# observations
obs <- list(idx=idcs, stk=stock(om))
# FLoem object
oem <- FLoem(method=sca.oem, observations=obs, deviances=dev)

#====================================================================
# IEM
# not much to say, just assuming our own ignorance ...
#====================================================================
noise <- catch(stock(om))
noise[] <- 0.1
noise <- rlnorm(it, 0, sd=noise)
iem <- FLiem(method=noise.iem, args=list(noise=noise))

#====================================================================
# MP
#====================================================================

# SET intermediate year + start of stks, lags and frequency
mseargs <- list(iy=iy, fy=fy, data_lag=dl, management_lag=ml, frq=af)

# SETUP standard ICES advice rule
arule <- mpCtrl(

  # estimation method: full feedback (a4a) sca
  est = mseCtrl(method=sca.sa, args=list(fmodel=fmod, qmodel=qmod, update=FALSE)),

  # parametrizing the HCR
  phcr = mseCtrl(method=f.phcr, args=list(interval=3)),

  # hcr: hockeystick (fbar ~ ssb | lim, trigger, target, min)
#   hcr = mseCtrl(method=hockeystick.hcr,
#         #args=list(lim=0.25*refpts(om.brp)["fmax","ssb"], trigger=0.5*refpts(om.brp)["fmax","ssb"],
#         args=list(lim=0, trigger=0, target=1, min=0, metric="ssb", output="fbar")),

  # hcr: fixed F
  hcr = mseCtrl(method=f.hcr),

  # (i)mplementation (sys)tem: tac.is (C ~ F)
  #isys = mseCtrl(method=tac_sca.is, args=list(recyrs=recyrs_mp))
  isys = mseCtrl(method=effort.is, args=list(nyrs=3))
)

#====================================================================
# Run simulations
#====================================================================

tes0 <- mp(om, oem, iem, ctrl=arule, args=mseargs)

# # RUN over Ftarget grid
# TODO: F search grid
# fg_mp <- seq(0, 3, length=9)
# fgrid <- mps(om, ctrl=arule, args=mseargs, hcr=list(target=fg_mp),
#   names=paste0("F", fg_mp))
#
# # PLOT
# plot(om, fgrid)
