## Modelo Fator Analítico Multiplicativo Misto

### Carregando o conjunto de dados

dataset = read.csv("https://raw.githubusercontent.com/Kaio-Olimpio/Probability-for-GEI/master/maize_dataset.csv", 
                   sep = ",")

dataset <- transform(dataset, env = factor(Location), gen = factor(Hybrid), rept = factor(Rep), 
                     block = factor(Block))

num.env = nlevels(factor(dataset$env))
num.gen = nlevels(factor(dataset$gen))
name.env = levels(factor(dataset$env))
name.gen = levels(factor(dataset$gen))

### Seleção do modelo

m1 = asreml(fixed=GY~1 + rept:env + env,
            random=~fa(env,1):gen  + diag(env):block,
            residual= ~dsum(~id(units)|env),
            data=dataset,
            maxiter=100)
m2 = asreml(fixed=GY~1 + rept:env + env,
            random=~fa(env,2):gen  + diag(env):block,
            residual= ~dsum(~id(units)|env),
            data=dataset,
            maxiter=100)
m3 = asreml(fixed=GY~1 + rept:env + env,
            random=~fa(env,3):gen  + diag(env):block,
            residual= ~dsum(~id(units)|env),
            data=dataset,
            maxiter=100)
m4 = asreml(fixed=GY~1 + rept:env + env,
            random=~fa(env,4):gen  + diag(env):block,
            residual= ~dsum(~id(units)|env),
            data=dataset,
            maxiter=100)
m5 = asreml(fixed=GY~1 + rept:env + env,
            random=~fa(env,5):gen  + diag(env):block,
            residual= ~dsum(~id(units)|env),
            data=dataset,
            maxiter=100)

aicm1 = summary(m1)$aic
sum1 = summary(m1)$varcomp
fa1.loadings = sum1[grep('fa1', rownames(sum1)),1]
spvar = as.matrix(diag(sum1[grep("var",row.names(sum1)),1]))
lamblamb = fa1.loadings %*% t(fa1.loadings)
expvar1 = (sum(diag(lamblamb))/sum(diag(lamblamb + spvar)))*100

aicm2 = summary(m2)$aic
sum2 = summary(m2)$varcomp
fa1.loadings = sum2[grep('fa1', rownames(sum2)),1]
fa2.loadings = sum2[grep('fa2', rownames(sum2)),1]
mat.loadings2 = as.matrix(cbind(fa1.loadings, fa2.loadings))
spvar = as.matrix(diag(sum2[grep("var",row.names(sum2)),1]))
lamblamb = mat.loadings2 %*% t(mat.loadings2)
expvar2 = (sum(diag(lamblamb))/sum(diag(lamblamb + spvar)))*100
svdmat = svd(mat.loadings2)
mat.loadings2.star = -1*mat.loadings2 %*% svdmat$v
lamblamb.star = mat.loadings2.star %*% t(mat.loadings2.star)
expvar2.star = (sum(diag(lamblamb.star))/sum(diag(lamblamb.star + spvar)))*100

aicm3 = summary(m3)$aic
bicm3 = summary(m3)$bic
sum3 = summary(m3)$varcomp
fa1.loadings = sum3[grep('fa1', rownames(sum3)),1]
fa2.loadings = sum3[grep('fa2', rownames(sum3)),1]
fa3.loadings = sum3[grep('fa3', rownames(sum3)),1]
mat.loadings3 = as.matrix(cbind(fa1.loadings, fa2.loadings,fa3.loadings))
spvar = as.matrix(diag(sum3[grep("var",row.names(sum3)),1]))
lamblamb = mat.loadings3 %*% t(mat.loadings3)
expvar3 = (sum(diag(lamblamb))/sum(diag(lamblamb + spvar)))*100
svdmat = svd(mat.loadings3)
mat.loadings3.star = -1*mat.loadings3 %*% svdmat$v
lamblamb.star = mat.loadings3.star %*% t(mat.loadings3.star)
expvar3.star = (sum(diag(lamblamb.star))/sum(diag(lamblamb.star + spvar)))*100

aicm4 = summary(m4)$aic
bicm4 = summary(m4)$bic
sum4 = summary(m4)$varcomp
fa1.loadings = sum4[grep('fa1', rownames(sum4)),1]
fa2.loadings = sum4[grep('fa2', rownames(sum4)),1]
fa3.loadings = sum4[grep('fa3', rownames(sum4)),1]
fa4.loadings = sum4[grep('fa4', rownames(sum4)),1]
mat.loadings4 = as.matrix(cbind(fa1.loadings, fa2.loadings,fa3.loadings,
                                fa4.loadings))
spvar = as.matrix(diag(sum4[grep("var",row.names(sum4)),1]))
lamblamb = mat.loadings4 %*% t(mat.loadings4)
expvar4 = (sum(diag(lamblamb))/sum(diag(lamblamb + spvar)))*100
svdmat = svd(mat.loadings4)
mat.loadings4.star = -1*mat.loadings4 %*% svdmat$v
lamblamb.star = mat.loadings4.star %*% t(mat.loadings4.star)
expvar4.star = (sum(diag(lamblamb.star))/sum(diag(lamblamb.star + spvar)))*100

aicm5 = summary(m5)$aic
sum5 = summary(m5)$varcomp
fa1.loadings = sum5[grep('fa1', rownames(sum5)),1]
fa2.loadings = sum5[grep('fa2', rownames(sum5)),1]
fa3.loadings = sum5[grep('fa3', rownames(sum5)),1]
fa4.loadings = sum5[grep('fa4', rownames(sum5)),1]
fa5.loadings = sum5[grep('fa5', rownames(sum5)),1]
mat.loadings5 = as.matrix(cbind(fa1.loadings, fa2.loadings,fa3.loadings,
                                fa4.loadings,fa5.loadings))
spvar = as.matrix(diag(sum5[grep("var",row.names(sum5)),1]))
lamblamb = mat.loadings5 %*% t(mat.loadings5)
expvar5 = (sum(diag(lamblamb))/sum(diag(lamblamb + spvar)))*100
svdmat = svd(mat.loadings5)
mat.loadings5.star = -1*mat.loadings5 %*% svdmat$v
lamblamb.star = mat.loadings5.star %*% t(mat.loadings5.star)
expvar5.star = (sum(diag(lamblamb.star))/sum(diag(lamblamb.star + spvar)))*100


fit = data.frame("Model" = c('FA1', 'FA2','FA3','FA4','FA5'),
                 'AIC' = c(aicm1, aicm2, aicm3, aicm4, aicm5),
                 'ExpVar.Star' = c(expvar1, expvar2.star,expvar3.star,expvar4.star,expvar5.star))
fit

### Obtenção da matriz de covariâncias genéticas

#### Matriz de cargas fatoriais 

sum5 = summary(m5)$varcomp
fa1.loadings = sum5[grep('fa1', rownames(sum5)),1]
fa2.loadings = sum5[grep('fa2', rownames(sum5)),1]
fa3.loadings = sum5[grep('fa3', rownames(sum5)),1]
fa4.loadings = sum5[grep('fa4', rownames(sum5)),1]
fa5.loadings = sum5[grep('fa5', rownames(sum5)),1]
mat.loadings5 = as.matrix(cbind(fa1.loadings, fa2.loadings,fa3.loadings,
                                fa4.loadings,fa5.loadings))

#### Rotação da matriz de cargas fatoriais

svdmat = svd(mat.loadings5)
mat.loadings5.star = -1*mat.loadings5 %*% svdmat$v

#### Matriz de variâncias específicas

spvar = as.matrix(diag(sum5[grep("var",row.names(sum5)),1]))

#### Matriz de covariâncias genéticas
  
gencov = mat.loadings5.star %*% t(mat.loadings5.star) + spvar

### Correlações genéticas entre ambientes

rownames(gencov) = colnames(gencov) = name.env

corr = cov2cor(gencov)

Heatmap(corr,col=colorRampPalette(brewer.pal(8, "YlOrBr"))(25),
        column_dend_height = unit(2, "cm"), 
        clustering_method_rows = "complete",
        row_dend_width = unit(2, "cm"), 
        clustering_method_columns = "complete",
        heatmap_legend_param = list(title="Correlation",
                                    at=c(-1,0,1),labels=c("-1","0","1")),
        row_names_gp = gpar(fontsize = 11),
        column_names_gp = gpar(fontsize = 11))

### Obtenção dos EBLUPs

#### Vetor de escores fatoriais

coef5 = summary(m5,coef = T)$coef.random

fa1.scores = coef5[grep("Comp1",row.names(coef5)),1];
names(fa1.scores) = sub("fa(env, 4)_Comp1:gen_","",names(fa1.scores),fixed=T)
fa2.scores = coef5[grep("Comp2",row.names(coef5)),1];
names(fa2.scores) = sub("fa(env, 4)_Comp2:gen_","",names(fa2.scores),fixed=T)
fa3.scores = coef5[grep("Comp3",row.names(coef5)),1];
names(fa3.scores) = sub("fa(env, 4)_Comp3:gen_","",names(fa3.scores),fixed=T)
fa4.scores = coef5[grep("Comp4",row.names(coef5)),1];
names(fa4.scores) = sub("fa(env, 4)_Comp4:gen_","",names(fa4.scores),fixed=T)
fa5.scores = coef5[grep("Comp4",row.names(coef5)),1];
names(fa5.scores) = sub("fa(env, 4)_Comp4:gen_","",names(fa5.scores),fixed=T)

fa.scores = rbind(as.matrix(fa1.scores),as.matrix(fa2.scores),
                  as.matrix(fa3.scores),as.matrix(fa4.scores),
                  as.matrix(fa5.scores))

#### Rotação do vetor de escores fatoriais 
  
fa.scores.star = -kronecker(t(svdmat$v), diag(num.gen))%*%fa.scores 
rownames(fa.scores.star) = rep(name.gen,5)

fa1.scores.star = fa.scores.star[1:num.gen,1]
fa2.scores.star = fa.scores.star[(num.gen+1):(num.gen*2),1]
fa3.scores.star = fa.scores.star[(num.gen*2+1):(num.gen*3),1]
fa4.scores.star = fa.scores.star[(num.gen*3+1):(num.gen*4),1]
fa5.scores.star = fa.scores.star[(num.gen*4+1):(num.gen*5),1]

#### EBLUPs marginais

EBLUPs_marg = (kronecker(mat.loadings5.star,diag(num.gen))) %*% fa.scores.star 
EBLUPs_marg = data.frame("Environment" = rep(name.env,each = num.gen),
                         "Gen" = rep(name.gen,num.env),
                         "EBLUP_marg" = EBLUPs_marg)


### Regressões latentes

EBLUPs_marg$FL_1 = rep(mat.loadings5.star[,1], each = num.gen)
ggplot(subset(EBLUPs_marg, Gen %in% c("G9","G32","G20","G23")), 
       aes(x=FL_1,y=EBLUP_marg,colour=Gen))+
  geom_smooth(aes(x=FL_1,y=EBLUP_marg),method=lm, fill = NA)+
  facet_wrap(~Gen, ncol = 4, as.table = T, dir = 'v', scales = 'free_x')+
  geom_vline(xintercept = mean(mat.loadings5.star[,1]),linetype="dashed") + 
  geom_point()+
  xlab("Cargas do primeiro fator")+ylab("EBLUP")+
  labs(color = "Genótipo", title = "EBLUP x Cargas do primeiro fator")

EBLUPs_marg$EBLUPs_marg_f2 = EBLUPs_marg$EBLUP_marg - 
  (kronecker(mat.loadings5.star[,1],diag(num.gen))) %*%
  as.matrix(fa1.scores.star)
EBLUPs_marg$FL_2 = rep(mat.loadings5.star[,2], each = num.gen)

ggplot(subset(EBLUPs_marg, Gen %in% c("G9","G32","G20","G23")), 
       aes(x=FL_2,y=EBLUPs_marg_f2,colour=Gen))+
  geom_smooth(aes(x=FL_2,y=EBLUPs_marg_f2),method=lm, fill = NA)+
  facet_wrap(~Gen, ncol = 4, as.table = T, dir = 'v', scales = 'free_x')+
  geom_vline(xintercept = mean(mat.loadings5.star[,1]),linetype="dashed") + 
  geom_point()+
  xlab("Cargas do segundo fator")+
  ylab(expression(EBLUP[m] - lambda[1]^'*' * f[1]^'*'))+
  labs(color = "Genótipo", title = "EBLUP x Cargas do segundo fator")

EBLUPs_marg$EBLUPs_marg_f3 = EBLUPs_marg$EBLUP_marg - 
  ((kronecker(mat.loadings5.star[,1],diag(num.gen)))
   %*% as.matrix(fa1.scores.star)+
     (kronecker(mat.loadings5.star[,2],diag(num.gen))) 
   %*% as.matrix(fa2.scores.star))
EBLUPs_marg$FL_3 = rep(mat.loadings5.star[,3], each = num.gen)

ggplot(subset(EBLUPs_marg, Gen %in% c("G9","G32","G20","G23")), 
       aes(x=FL_3,y=EBLUPs_marg_f3,colour=Gen))+
  geom_smooth(aes(x=FL_3,y=EBLUPs_marg_f3),method=lm, fill = NA)+
  facet_wrap(~Gen, ncol = 4, as.table = T, dir = 'v', scales = 'free_x')+
  geom_vline(xintercept = mean(mat.loadings5.star[,3]),linetype="dashed") + 
  geom_point()+
  xlab("Cargas do terceiro fator")+
  ylab(expression(EBLUP[m] - (lambda[1]^'*' * f[1]^'*' + 
                                lambda[2]^'*' * f[2]^'*')))+
  labs(color = "Genótipo",title = "EBLUP x Cargas do terceiro fator")


EBLUPs_marg$EBLUPs_marg_f4 = EBLUPs_marg$EBLUP_marg - 
  ((kronecker(mat.loadings5.star[,1],diag(num.gen))) 
   %*% as.matrix(fa1.scores.star)+
     (kronecker(mat.loadings5.star[,2],diag(num.gen))) 
   %*% as.matrix(fa2.scores.star) + 
     (kronecker(mat.loadings5.star[,3],diag(num.gen))) 
   %*% as.matrix(fa3.scores.star))
EBLUPs_marg$FL_4 = rep(mat.loadings5.star[,4], each = num.gen)

ggplot(subset(EBLUPs_marg, Gen %in% c("G9","G32","G20","G23")), 
       aes(x=FL_4,y=EBLUPs_marg_f4,colour=Gen))+
  geom_smooth(aes(x=FL_4,y=EBLUPs_marg_f4),method=lm, fill = NA)+
  facet_wrap(~Gen, ncol = 4, as.table = T, dir = 'v', scales = 'free_x')+
  geom_vline(xintercept = mean(mat.loadings5.star[,3]),linetype="dashed") + 
  geom_point()+
  xlab("Cargas do quarto fator")+
  ylab(expression(EBLUP[m] - (lambda[1]^'*' * f[1]^'*' +
                                lambda[2]^'*' * f[2]^'*' + lambda[3]^'*' * f[3]^'*')))+ 
  labs(color = "Genótipo",title = "EBLUP x Cargas do quarto fator")


EBLUPs_marg$EBLUPs_marg_f5 = EBLUPs_marg$EBLUP_marg - 
  ((kronecker(mat.loadings5.star[,1],diag(num.gen))) %*%
     as.matrix(fa1.scores.star)+
     (kronecker(mat.loadings5.star[,2],diag(num.gen))) %*%
     as.matrix(fa2.scores.star) + 
     (kronecker(mat.loadings5.star[,3],diag(num.gen))) %*% 
     as.matrix(fa3.scores.star) +
     (kronecker(mat.loadings5.star[,4],diag(num.gen))) %*% 
     as.matrix(fa4.scores.star))
EBLUPs_marg$FL_5 = rep(mat.loadings5.star[,5], each = num.gen)

ggplot(subset(EBLUPs_marg, Gen %in% c("G9","G32","G20","G23")), 
       aes(x=FL_5,y=EBLUPs_marg_f5,colour=Gen))+
  geom_smooth(aes(x=FL_5,y=EBLUPs_marg_f5),method=lm, fill = NA)+
  facet_wrap(~Gen, ncol = 4, as.table = T, dir = 'v', scales = 'free_x')+
  geom_vline(xintercept = mean(mat.loadings5.star[,3]),linetype="dashed") + 
  geom_point()+
  xlab("Cargas do quinto fator")+
  ylab(expression(EBLUP[m] - (lambda[1]^'*' * f[1]^'*' +
                                lambda[2]^'*' * f[2]^'*' + lambda[3]^'*' * f[3]^'*' +
                                lambda[4]^'*' * f[4]^'*')))+ 
  labs(color = "Genótipo",title = "EBLUP x Cargas do quinto fator")

### Ferramentas de seleção

#### Performance geral x Estabilidade

PG = mean(mat.loadings5.star[,1]) * as.matrix(fa1.scores.star)
rest = (EBLUPs_marg$EBLUP_marg - 
          (kronecker(mat.loadings5.star[,1],diag(num.gen))) %*%
          fa1.scores.star)^2
STA = data.frame("gen" = rep(name.gen, num.env),
                 "ST" = rest)
STA = STA %>% group_by(gen) %>% summarise(ST = sqrt(mean(ST)))
plot1 = data.frame("gen"=name.gen,
                   "PG" = PG,
                   "ST" = STA$ST)

ggplot(data=plot1)+
  geom_point(aes(x = ST, y = PG, color=gen),size=2)+
  geom_vline(xintercept = 0,linetype="dashed", colour = "black")+
  geom_hline(yintercept = 0,linetype="dashed")+
  scale_x_continuous(name = "Estabilidade", limits = c(0,1))+
  ylab("Performance Geral") +
  gghighlight(max(ST)<=0.50,max(PG)>=0.25,use_direct_label = F)+
  geom_label_repel(aes(x = ST, y = PG, label = gen), size = 4, box.padding = 1)+
  labs(color = "Genótipo")


#### Performance geral x Responsividade aos fatores

mat.loadings5.star.df = as.data.frame(mat.loadings5.star)

colnames(mat.loadings5.star.df) = c("fa1.loadings","fa2.loadings","fa3.loadings",
                                    "fa4.loadings","fa5.loadings")

fa2.mednegload = mean(mat.loadings5.star.df[mat.loadings5.star.df$fa2.loadings <0,]$fa2.loadings)
fa3.mednegload = mean(mat.loadings5.star.df[mat.loadings5.star.df$fa3.loadings <0,]$fa3.loadings)
fa4.mednegload = mean(mat.loadings5.star.df[mat.loadings5.star.df$fa4.loadings <0,]$fa4.loadings)
fa5.mednegload = mean(mat.loadings5.star.df[mat.loadings5.star.df$fa5.loadings <0,]$fa5.loadings)

fa2.medposload = mean(mat.loadings5.star.df[mat.loadings5.star.df$fa2.loadings >0,]$fa2.loadings)
fa3.medposload = mean(mat.loadings5.star.df[mat.loadings5.star.df$fa3.loadings >0,]$fa3.loadings)
fa4.medposload = mean(mat.loadings5.star.df[mat.loadings5.star.df$fa4.loadings >0,]$fa4.loadings)
fa5.medposload = mean(mat.loadings5.star.df[mat.loadings5.star.df$fa5.loadings >0,]$fa5.loadings)


respfa2 = (fa2.medposload-fa2.mednegload)*fa2.scores.star
respfa3 = (fa3.medposload-fa3.mednegload)*fa3.scores.star
respfa4 = (fa4.medposload-fa4.mednegload)*fa4.scores.star
respfa5 = (fa5.medposload-fa5.mednegload)*fa5.scores.star


plot2 = data.frame("gen" = name.gen,
                   "PG" = PG,
                   "respfa2" = respfa2,
                   "respfa3" = respfa3,
                   "respfa4" = respfa4,
                   "respfa5" = respfa5)

ggplot(data=plot2)+
  geom_point(aes(x = respfa2, y = PG, color=gen),size=2)+
  geom_vline(xintercept = 0,linetype="dashed", colour = "black")+
  geom_hline(yintercept = 0,linetype="dashed")+
  ylab("Performance geral") + 
  scale_x_continuous(name = "Responsividade ao segundo fator",
                     limits = c(-1.5,1.5))+
  labs(color = "Genótipo")

ggplot(data=plot2)+
  geom_point(aes(x = respfa3, y = PG, color=gen),size=2)+
  geom_vline(xintercept = 0,linetype="dashed", colour = "black")+
  geom_hline(yintercept = 0,linetype="dashed")+
  ylab("Performance geral") + 
  scale_x_continuous(name = "Responsividade ao terceiro fator", 
                     limits = c(-1.5,1.5))+
  labs(color = "Genótipo")


ggplot(data=plot2)+
  geom_point(aes(x = respfa4, y = PG, color=gen),size=2)+
  geom_vline(xintercept = 0,linetype="dashed", colour = "black")+
  geom_hline(yintercept = 0,linetype="dashed")+
  ylab("Performance Geral") + 
  scale_x_continuous(name = "Responsividade ao quarto fator",
                     limits = c(-1.5,1.5))+
  labs(color = "Genótipo")


ggplot(data=plot2)+
  geom_point(aes(x = respfa5, y = PG, color=gen),size=2)+
  geom_vline(xintercept = 0,linetype="dashed", colour = "black")+
  geom_hline(yintercept = 0,linetype="dashed")+
  ylab("Performance Geral") + 
  scale_x_continuous(name = "Responsividade ao quinto fator",
                     limits = c(-1.5,1.5))+
  labs(color = "Genótipo")
