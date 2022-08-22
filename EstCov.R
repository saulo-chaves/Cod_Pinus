library(asreml)

data1 = read.table("EstCov.txt", header = T)

data1$amb = as.factor(data1$amb)
data1$gen = as.factor(data1$gen)
data1$rept = as.factor(data1$rept)


# Modelo de simetria composta: Resíduos homogêneos
m1 = asreml(fixed = y ~ rept:amb,
            random = ~gen + gen:amb,
            data = data1)
summary(m1)$varcomp
summary(m1)$aic
summary(m1)$bic

# Modelo de simetria composta: Resíduos heterogêneos
m2 = asreml(fixed = y ~ rept:amb,
            random = ~gen + gen:amb,
            residual = ~dsum(~id(units)|amb),
            data = data1)
summary(m2)$varcomp
summary(m2)$aic
summary(m2)$bic

# Modelo de simetria composta heterogênea
m3 = asreml(fixed = y ~ rept:amb,
            random = ~corh(amb):gen,
            residual = ~dsum(~id(units)|amb),
            data = data1)
summary(m3)$varcomp
summary(m3)$aic
summary(m3)$bic


# Modelo não-estruturado
m4 = asreml(fixed = y ~ rept:amb,
            random = ~us(amb):gen,
            residual = ~dsum(~id(units)|amb),
            data = data1)
summary(m4)$varcomp
summary(m4)$aic

