#====================================================================
# EJ (20251118)
# MSE for hake in 1-7 benchmark tests
#====================================================================

#====================================================================
# load libraries and data and setup run
#====================================================================
library(FLa4a)
library(FLasher)
library(FLBRP)
library(ggplotFL)

# stock
load("../data/HKE_1_5_6_7_STK.Rdata")
# index
load("../data/HKE_1_5_6_7_IDX.Rdata")
hke.idx <- FLIndices(ibts=hke.idx)
tos <- mean(range(hke.idx[[1]])[c("startf","endf")]) # time of survey
iar <- ac(range(hke.idx[[1]])["min"]:range(hke.idx[[1]])["max"]) # index age range

# set up
rng <- range(hke.stk)
ny <- 5
yrs <- rng["maxyear"]+0:(ny-1)
its <- 25
seed <- 123
af <- 1 # advice frequency
frefpt <- "f0.1"
estimator.lst <- split(yrs, yrs)
perceivedstock.lst <- split(yrs, yrs)
tracking <- FLQuant(dimnames=list(quant=c("om.f", "mp.f", paste("mp",frefpt, sep="."), "mp.catch_advice"), year=yrs, iter=1:its))

#====================================================================
# OM conditioning based on stock assessment
# Uncertainty is introduced by MCMC fitting process
# Uncertainty is propagated through the SR model into reference points
#====================================================================
# survey catchability submodel
qmod <- list(~ I(1/(1+exp(-age))))

# fishing mortality submodel
fmod <- ~s(age, k = 4) + s(year, k = 8) + te(age, year, k = c(3, 10))

# fit
fit <- sca(hke.stk, hke.idx, fmodel=fmod, qmodel=qmod)
fit <- simulate(fit, its, seed=seed)
pred <- predict(fit)
om <- hke.stk + fit
om.sr <- fmle(as.FLSR(om, model="bevholt"), control = list(trace = 0))
om.brp <- brp(FLBRP(om, sr=om.sr))
plot(om)

#====================================================================
# Observation error model based on stock assessment estimates
# Uncertainty is introduced by adding observation error to catches
# and survey
# This process creates a stock object with estimation and observation
# errors up to maxyear and the observation error estimates for the
# future. All objects have range up to the end of the analyis.
#====================================================================

fit.oe <- simulate(fit, its, seed=seed, obserror=TRUE)
pred.oe <- predict(fit.oe)

# index
idx.oe <- hke.idx + fit.oe
idx.oe <- window(idx.oe, end=rng["maxyear"]+ny)

flq0 <- idx.oe$ibts@index
flq0[] <- pred$qmodel[[1]][,1] # q constant over time, no oe
idx.oe$ibts@index.q <- flq0

flq0[] <- pred$vmodel$ibts[,1] # q oe constant over time
q.oe <- flq0

# catch
stk.oe <- hke.stk + fit.oe
stk.oe <- fwdWindow(stk.oe, end=rng["maxyear"]+ny)
c.oe <- catch.n(stk.oe)
c.oe[] <- pred$vmodel$catch[,1] # catch oe constant over time

#====================================================================
# Management procedure
#====================================================================

# estimator setup
fmod <- ~s(age, k=3) + s(year, k=20)
fmod <- ~s(age, k = 4) + s(year, k = 8) + te(age, year, k = c(3, 10))
qmod <- defaultQmod(idx.oe)
srmod <- defaultSRmod(stk.oe)
n1mod <- defaultN1mod(stk.oe)
vmod <- defaultVmod(stk.oe, idx.oe)

#--------------------------------------------------------------------
# Year 1, providing advice for year 2, data up to previous year (1-1)
#--------------------------------------------------------------------

ay <- rng["maxyear"]
dy <- ac(ay - af)
py <- ay + af
# store om value for later analysis
tracking["om.f",ac(ay)] <- fbar(om)[,dy]

# OEM
# Observation model generates data up to year 1-1
# Already sorted out in conditioning
idx <- window(idx.oe, end=dy)
stk <- window(stk.oe, end=dy)

# MP
# Estimator is sca with FLa4a
estimator <- sca(stk, idx, fmodel=fmod, qmodel=qmod, srmodel=srmod, n1model=n1mod, vmodel=vmod, fit="MP")
stk0 <- stk + estimator
# beverton and holt stock recruitment is fitted every year
sr0 <- fmle(as.FLSR(window(stk0, start=2015), model="geomean"), control = list(trace = 0))
# HCR is fmax as target, also estimated every year
ftrg <- refpts(brp(FLBRP(stk0, sr0)))[frefpt, "harvest"]
# Projection
yrs <- c(ay:py)
ctrl <- fwdControl(list(year=yrs, value=ftrg, quant="fbar"))
stk_fwd <- fwdWindow(stk0, end=py)
stk_fwd <- fwd(stk_fwd, control=ctrl, sr=sr0)
# catch advice is based on the projected catch if the target is applied
catch_advice <- catch(stk_fwd)[,ac(py)]

# IEM
# No implementation model

# Update tracking matrix and lists
tracking["mp.f",ac(ay)] <- fbar(stk0)[,dy]
tracking[paste("mp",frefpt, sep="."),ac(ay)] <- ftrg
tracking["mp.catch_advice",ac(ay)] <- catch_advice
estimator.lst[[ac(ay)]] <- fit
perceivedstock.lst[[ac(ay)]] <- stk0

#--------------------------------------------------------------------
# Year 2, providing advice for year 3, data up to previous year (2-1)
#--------------------------------------------------------------------

ay <- ay + 1
dy <- ac(ay - af)
py <- ay + af
# store om value for later analysis
tracking["om.f",ac(ay)] <- fbar(om)[,dy]

# Update OM with previous decision
ctrl <- fwdControl(list(year=ay, value=catch_advice, quant="catch"))
om <- fwdWindow(om, end=ay)
om <- fwd(om, control=ctrl, sr=om.sr)

# OEM
# Already sorted out in conditioning
idx <- window(idx.oe, end=ay-1)
stk <- window(stk.oe, end=ay-1)

# MP
estimator <- sca(stk, idx, fmodel=fmod, qmodel=qmod, srmodel=srmod, n1model=n1mod, vmodel=vmod, fit="MP")
stk0 <- stk2 <- stk + estimator
sr0 <- fmle(as.FLSR(window(stk0, start=2015), model="geomean"), control = list(trace = 0))
ftrg <- refpts(brp(FLBRP(stk0, sr0)))[frefpt, "harvest"]
yrs <- c(ay:py)
ctrl <- fwdControl(list(year=yrs, value=ftrg, quant="fbar"))
stk_fwd <- fwdWindow(stk0, end=py)
stk_fwd <- fwd(stk_fwd, control=ctrl, sr=sr0)
catch_advice <- catch(stk_fwd)[,ac(py)]

# IEM

# Update tracking matrix
tracking["mp.f",ac(ay)] <- fbar(stk0)[,dy]
tracking[paste("mp",frefpt, sep="."),ac(ay)] <- ftrg
tracking["mp.catch_advice",ac(ay)] <- catch_advice
estimator.lst[[ac(ay)]] <- fit
perceivedstock.lst[[ac(ay)]] <- stk0

#--------------------------------------------------------------------
# Year 3, providing advice for year 4, data up to previous year (3-1)
#--------------------------------------------------------------------

ay <- ay + 1
dy <- ac(ay - af)
py <- ay + af
# store om value for later analysis
tracking["om.f",ac(ay)] <- fbar(om)[,dy]

# Update OM with previous decision
ctrl <- fwdControl(list(year=ay, value=catch_advice, quant="catch"))
om <- fwdWindow(om, end=ay)
om <- fwd(om, control=ctrl, sr=om.sr)

# OEM
idx0 <- (stock.n(om)[,dy] * exp(-(harvest(om)[,dy]+m(om)[,dy])*tos))[iar] * index.q(idx.oe[[1]])[,dy]
index(idx.oe[[1]])[,dy] <- exp(log(idx0) + rnorm(its, 0, sqrt(q.oe[,dy])))
idx <- window(idx.oe, end=ay-1)

stk <- window(om, end=dy)
stk[,dy]@catch.n <- exp(log(stk[,dy]@catch.n) + rnorm(its, 0, sqrt(c.oe[,dy])))

# MP
estimator <- sca(stk, idx, fmodel=fmod, qmodel=qmod, srmodel=srmod, n1model=n1mod, vmodel=vmod, fit="MP")
stk0 <- stk3 <- stk + estimator
sr0 <- fmle(as.FLSR(window(stk0, start=2015), model="geomean"), control = list(trace = 0))
ftrg <- refpts(brp(FLBRP(stk0, sr0)))[frefpt, "harvest"]
yrs <- c(ay:py)
ctrl <- fwdControl(list(year=yrs, value=ftrg, quant="fbar"))
stk_fwd <- fwdWindow(stk0, end=py)
stk_fwd <- fwd(stk_fwd, control=ctrl, sr=sr0)
catch_advice <- catch(stk_fwd)[,ac(py)]

# IEM

# Update tracking matrix
tracking["mp.f",ac(ay)] <- fbar(stk0)[,dy]
tracking[paste("mp",frefpt, sep="."),ac(ay)] <- ftrg
tracking["mp.catch_advice",ac(ay)] <- catch_advice
estimator.lst[[ac(ay)]] <- fit
perceivedstock.lst[[ac(ay)]] <- stk0

#--------------------------------------------------------------------
# Year 4, providing advice for year 5, data up to previous year (4-1)
#--------------------------------------------------------------------

ay <- ay + 1
dy <- ac(ay - af)
py <- ay + af
# store om value for later analysis
tracking["om.f",ac(ay)] <- fbar(om)[,dy]

# Update OM with previous decision
ctrl <- fwdControl(list(year=ay, value=catch_advice, quant="catch"))
om <- fwdWindow(om, end=ay)
om <- fwd(om, control=ctrl, sr=om.sr)

# OEM
idx0 <- (stock.n(om)[,dy] * exp(-(harvest(om)[,dy]+m(om)[,dy])*tos))[iar] * index.q(idx.oe[[1]])[,dy]
index(idx.oe[[1]])[,dy] <- exp(log(idx0) + rnorm(its, 0, sqrt(q.oe[,dy])))
idx <- window(idx.oe, end=ay-1)

stk <- window(om, end=dy)
stk[,dy]@catch.n <- exp(log(stk[,dy]@catch.n) + rnorm(its, 0, sqrt(c.oe[,dy])))

# MP
estimator <- sca(stk, idx, fmodel=fmod, qmodel=qmod, srmodel=srmod, n1model=n1mod, vmodel=vmod, fit="MP")
stk0 <- stk4 <- stk + estimator
sr0 <- fmle(as.FLSR(window(stk0, start=2015), model="geomean"), control = list(trace = 0))
ftrg <- refpts(brp(FLBRP(stk0, sr0)))[frefpt, "harvest"]
yrs <- c(ay:py)
ctrl <- fwdControl(list(year=yrs, value=ftrg, quant="fbar"))
stk_fwd <- fwdWindow(stk0, end=py)
stk_fwd <- fwd(stk_fwd, control=ctrl, sr=sr0)
catch_advice <- catch(stk_fwd)[,ac(py)]

# IEM

# Update tracking matrix
tracking["mp.f",ac(ay)] <- fbar(stk0)[,dy]
tracking[paste("mp",frefpt, sep="."),ac(ay)] <- ftrg
tracking["mp.catch_advice",ac(ay)] <- catch_advice
estimator.lst[[ac(ay)]] <- fit
perceivedstock.lst[[ac(ay)]] <- stk0

#--------------------------------------------------------------------
# Year 5
#--------------------------------------------------------------------

ay <- ay + 1
dy <- ac(ay - af)
py <- ay + af
# store om value for later analysis
tracking["om.f",ac(ay)] <- fbar(om)[,dy]

# Update OM with previous decision
ctrl <- fwdControl(list(year=ay, value=catch_advice, quant="catch"))
om <- fwdWindow(om, end=ay)
om <- fwd(om, control=ctrl, sr=om.sr)

#====================================================================
# Analysis of results
#====================================================================

#--------------------------------------------------------------------
# kobe is your friend
#--------------------------------------------------------------------
# ssb and F relative to MSY reference points
ssb_ssbmsy <- ssb(om)/refpts(om.brp)["msy", "ssb"]
f_fmsy <- fbar(om)/refpts(om.brp)["msy", "harvest"]

ssb_ssbmsy <- window(iterMedians(ssb(om)/refpts(om.brp)["msy", "ssb"]), start=2020)
f_fmsy <- window(iterMedians(fbar(om)/refpts(om.brp)["msy", "harvest"]), start=2020)

# code to compute kobe plot
ggplot(mapping=aes(y=c(f_fmsy), x=c(ssb_ssbmsy))) +
  # Add quadrants
  geom_rect(aes(xmin = 1, xmax = Inf, ymin = 0, ymax = 1), fill = "green", alpha = 0.5) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), fill = "yellow", alpha = 0.5) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = 1, ymax = Inf), fill = "red", alpha = 0.5) +
  geom_rect(aes(xmin = 1, xmax = Inf, ymin = 1, ymax = Inf), fill = "yellow", alpha = 0.5) +
  # Reference lines
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  # Points
  geom_point(size = 2) +
  # lines
  geom_path(arrow = arrow(type = "open", length = unit(0.15, "inches")), linewidth = 0.5) +
  # Labels and theme
  labs(
    x = expression(B / B[MSY]),
    y = expression(F / F[MSY]),
  ) +
  theme_minimal()

#--------------------------------------------------------------------
# Did the assessment track the OM?
#--------------------------------------------------------------------

plot(log(tracking["mp.f"])~log(tracking["om.f"]), ylab="f estimated in the mp", xlab="F in the om", pch=19)

perceivedstock.lst[[5]] <- om
names(perceivedstock.lst) <- c(1:4, "om")
plot(FLStocks(perceivedstock.lst))

#--------------------------------------------------------------------
# Were the reference points stable?
#--------------------------------------------------------------------

bwplot(data~year, data=tracking[paste("mp",frefpt, sep=".")], ylab=frefpt, xlab="Year")

#--------------------------------------------------------------------
# How was catch advice in relation to the status of the stock?
# (note there's no evaluation of status in the HCR)
#--------------------------------------------------------------------

plot(tracking["mp.catch_advice"]~I(tracking["mp.f"]/tracking[paste("mp",frefpt, sep=".")]), ylab="catch advice", xlab="F status in the mp", pch=19)

#====================================================================
# Extensions
#====================================================================

- Use a different stock assessment model
- Use a different assessment frequency
- Don't reestimate reference points every year
- Use a different catchability for the survey
- Introduce uncertainty in catchability, or recruitment
- ...

