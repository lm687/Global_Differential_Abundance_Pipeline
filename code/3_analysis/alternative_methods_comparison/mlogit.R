## https://cran.r-project.org/web/packages/mlogit/vignettes/c5.mxl.html
##' in mlogit it is the correlations of the independent variables, and not of the dependent variables,
##' which are included in the random effects
##' in the example below we pnly have two choices, A and B, so it's binomial?


library("mlogit")
data("Train", package = "mlogit")
Train$choiceid <- 1:nrow(Train)
Tr <- dfidx(Train, choice = "choice", varying = 4:11, sep = "_",
            opposite = c("price", "comfort", "time", "change"),
            idx = list(c("choiceid", "id")), idnames = c("chid", "alt"))
Tr$price <- Tr$price / 100 * 2.20371
Tr$time <- Tr$time / 60
Train.ml <- mlogit(choice ~ price + time + change + comfort | - 1, Tr)
coef(summary(Train.ml))

table(Train$choice) ## binomial?

View(Train)
coef(Train.ml)[- 1] / coef(Train.ml)[1]

## uncorrelated mixed effects model
Train.mxlu <- mlogit(choice ~ price + time + change + comfort | - 1, Tr,
                     panel = TRUE, rpar = c(time = "n", change = "n", comfort = "n"), R = 100,
                     correlation = FALSE, halton = NA, method = "bhhh")
names(coef(Train.mxlu))

## correlated model
Train.mxlc <- update(Train.mxlu, correlation = TRUE)
names(coef(Train.mxlc))
marg.ut.time <- rpar(Train.mxlc, "time")
summary(marg.ut.time)
vcov(Train.mxlc, what = "rpar")

summary(vcov(Train.mxlc, what = "rpar", type = "cor"))

### four choices
data("RiskyTransport", package = "mlogit")
RT <- dfidx(RiskyTransport, choice = "choice", idx = list(c("chid", "id"), "mode"),
            idnames = c("chid", "alt"))
ml.rt <- mlogit(choice ~ cost + risk  + seats + noise + crowdness +
                  convloc + clientele | 0, data = RT, weights = weight)

coef(ml.rt)[c("risk", "cost")]

mx.rt <- mlogit(choice ~ cost + risk  + seats + noise + crowdness +
                  convloc + clientele | 0, data = RT, weights = weight,
                rpar = c(cost = 'zbt', risk = 'zbt'), R = 100, halton = NA, panel = TRUE)

library("texreg")
htmlreg(list('Multinomial logit' = ml.rt, 'Mixed logit' = mx.rt),
        digits = 3, float.pos = "hbt", label = "tab:risktr", single.row = TRUE,
        caption = "Transportation choices.")
