# 3. Análise espacial

require(asreml)

## Carregando o conjunto de dados

data_esp = read.csv("espacial.csv", header = T,sep=";")

data_esp$Variety = as.factor(data_esp$Variety)
data_esp$Block = as.factor(data_esp$Block)
data_esp$Row = as.factor(data_esp$Row)
data_esp$Col= as.factor(data_esp$Col)

## Carregando a função para calcular o AIC ajustado ($AIC_c$)

icREML <- function(fm, scale=1) {
  if(!is.list(fm)) stop(" Models need to be in a list\n")
  if(is.null(names(fm))) namesfm <- paste("fm", 1:length(fm))
  else namesfm <- names(fm)
  require(asreml)
  asreml.options(Cfixed = TRUE, gammaPar=FALSE)
  fm <- lapply(fm, function(el) {
    if(is.null(el$Cfixed)) {
      out <- update(el, maxit=1) }
    else out <- el
    out})
  logl <- lapply(fm, function(el) el$loglik)
  summ <- lapply(fm, function(el) summary(el, coef=TRUE)$coef.fixed)
  which.X0 <- lapply(summ, function(el) !is.na(el[, "z.ratio"]))
  p.0 <- lapply(which.X0, function(el) sum(el))
  Cfixed <- lapply(fm, function(el) el$Cfixed)
  logdet <- lapply(1:length(fm), function(el, Cfixed, which.X0, scale) {
    log(prod(svd(as.matrix(scale*Cfixed[[el]][which.X0[[el]], 
                                              which.X0[[el]]]))$d))
  }, Cfixed, which.X0, scale)
  vparam <- lapply(fm, function(el) summary(el)$varcomp)
  q.0 <- lapply(vparam, function(el) sum(!(el$bound == "F" | el$bound == "B")))
  b.0 <- lapply(vparam, function(el) sum(el$bound == "F" | el$bound == "B"))
  logl <- lapply(1:length(fm), function(el, logl, logdet, p.0) {
    logl[[el]] - logdet[[el]]/2}, logl, logdet,p.0)
  aic <- unlist(lapply(1:length(fm), function(el, logl, p.0, q.0) {
    -2*logl[[el]] + 2*(p.0[[el]] + q.0[[el]])}, logl, p.0, q.0))
  bic <- unlist(lapply(1:length(fm), function(el, logl, p.0, q.0, fm) {
    -2*logl[[el]] + log(fm[[el]]$nedf+p.0[[el]])*(p.0[[el]] + q.0[[el]])},
    logl, p.0, q.0, fm))
  results <- data.frame(model=namesfm, loglik = unlist(logl), p=unlist(p.0),
                        q=unlist(q.0), b = unlist(b.0), AIC = aic, BIC = bic, logdet=unlist(logdet))
  row.names(results) <- 1:dim(results)[1]
  invisible(results)
}

## Modelo 1: Sem ajuste espacial

m1 <- asreml(yield ~ 1,
             random = ~ Variety + Block,
             data = data_esp)

plot(varioGram(m1))

## Modelo 2: Auto-regressivo

m2 <- asreml(yield ~ 1,
             random = ~ Variety + Block,
             residual = ~ ar1(Col):ar1(Row),
             data = data_esp)

plot(varioGram(m2))

## Modelo 3: Auto-regressivo + Pepita

m3 <- asreml(yield ~ 1,
             random = ~ Variety + Block + units,
             residual = ~ ar1(Col):ar1(Row),
             data = data_esp)

plot(varioGram(m3))

## Modelo 4: Anteriores + Aleatório de linha e coluna


m4 <- asreml(yield ~ 1,
             random = ~ Variety + Block + Row + Col + units,
             residual = ~ ar1(Col):ar1(Row),
             data = data_esp)
m4 <- update(m4)

plot(varioGram(m4))

## Modelo 5: Anteriores + Efeito linear fixo de linha e coluna


m5 <- asreml(yield ~ lin(Row) + lin(Col),
             random = ~ Variety + Block + Row + Col + units,
             residual = ~ ar1(Col):ar1(Row),
             data = data_esp)

plot(varioGram(m5))

## Modelo 6: Anteriores + Spline no sentido da coluna

m6 <- asreml(yield ~ lin(Row) + lin(Col),
             random = ~ Variety + Block + Row + Col + spl(col) + units,
             residual = ~ ar1(Col):ar1(Row),
             data = data_esp)

plot(varioGram(m6))

## AIC_c

mods = list(m1=m1, m2=m2, m3=m3, m4=m4, m5=m5, m6=m6)

selection <- icREML(mods)
selection
