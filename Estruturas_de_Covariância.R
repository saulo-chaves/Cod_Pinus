library(asreml)

## Carregando o conjunto de dados
data = read.table("EstCov.txt", header = T)
data$amb = as.factor(data$amb)
data$gen = as.factor(data$gen)
data$rept = as.factor(data$rept)

## Modelagem das estruturas de covariância

### Modelo 1: Simetria composta para genótipo e identidade para resíduo

m1 = asreml(fixed = y ~ rept:amb,
            random = ~gen + gen:amb,
            data = data)
summary(m1)$varcomp

### Modelo 2: Simetria composta para genótipo e bloco diagonal para resíduo

m2 = asreml(fixed = y ~ rept:amb,
            random = ~gen + gen:amb,
            residual = ~dsum(~id(units)|amb),
            data = data)
summary(m2)$varcomp

### Modelo 3: Simetria composta heterogênea para genótipo e bloco diagonal para resíduo

m3 = asreml(fixed = y ~ rept:amb,
            random = ~corh(amb):gen,
            residual = ~dsum(~id(units)|amb),
            data = data)
summary(m3)$varcomp


### Modelo 4: Estrutura multivariada para genótipo e bloco diagonal para resíduo

m4 = asreml(fixed = y ~ rept:amb,
            random = ~us(amb):gen,
            residual = ~dsum(~id(units)|amb),
            data = data)
summary(m4)$varcomp

### Seleção do modelo

data.frame("Modelo" = seq(1,4),
           "AIC" = c(summary(m1)$aic,summary(m2)$aic,summary(m3)$aic,
                     summary(m4)$aic))
