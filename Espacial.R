require(asreml)
#data()
data(nin89)
nin89
#ASReml-R Reference Manual Version 4 P.45
rcb.asr1 <- asreml(yield ~ Variety + Rep, residual = ~idv(units), data = nin89, na.method.X="omit")
summary(rcb.asr1)
rcb.asr1a <- asreml(yield ~ Variety, random = ~idv(Rep), residual = ~idv(units), data = nin89, na.method.X="omit")
summary(rcb.asr1a)
sp.asr2 <- asreml(yield ~ Variety, residual = ~idv(Column):ar1(Row), data = nin89)
summary(sp.asr2)
sp.asr2a <- asreml(yield ~ Variety, residual = ~ar1v(Column):ar1(Row), data = nin89)
summary(sp.asr2a)
sp.asr2b <- asreml(yield ~ Variety, random = ~idv(units), residual = ~ar1v(Column):ar1(Row), data = nin89)
summary(sp.asr2b)
sp.asr3 <- asreml(yield ~ Variety, random = ~ar1v(Column):ar1(Row), residual = ~idv(units), data = nin89)
summary(sp.asr3)

nin89$nloc <- nin89$nloc-3
nin89