require(asreml)
require(tidyverse)
#data()
data(nin89)
nin89 



#ASReml-R Reference Manual Version 4 P.45
rcb.asr1 <- asreml(yield ~ Variety + Rep, residual = ~idv(units), data = nin89,
                   na.action = na.method(x = "omit"))
summary(rcb.asr1)$aic

sp.asr2 <- asreml(yield ~ Variety, residual = ~idv(Column):ar1(Row), data = nin89)
summary(sp.asr2)$aic
plot.varioGram(varioGram(sp.asr2))

sp.asr2a <- asreml(yield ~ Variety, residual = ~ar1v(Column):ar1(Row), data = nin89)
summary(sp.asr2a)$aic
plot.varioGram(varioGram(sp.asr2a))

sp.asr2b <- asreml(yield ~ Variety, random = ~idv(units), 
                   residual = ~ar1v(Column):ar1(Row), 
                   data = nin89,
                   maxiter = 100)
summary(sp.asr2b)$aic
plot.varioGram(varioGram(sp.asr2b))

sp.asr3 <- asreml(yield ~ Variety, random = ~ar1v(Column):ar1(Row), residual = ~idv(units), data = nin89)
summary(sp.asr3)$aic
plot.varioGram(varioGram(sp.asr3))

nin89$nloc <- nin89$nloc-3
nin89
