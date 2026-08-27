rm(list=ls())

# Library -----------------------------------------------------------------
library(tidyverse)
library(BGLR)
library(AGHmatrix)
library(ComplexHeatmap)
library(plotly)

# Phenotypic and genomic data ---------------------------------------------------------
dat = readRDS("Data/G2F.RDS")
str(dat)
geno.code = readRDS(file = "Data/mrk_filtered.RDS")

# ggplot(data = as.data.frame(table(dat$gen, dat$env)),
#        aes(x = Var2, y = Var1, fill = as.factor(Freq))) + 
#   geom_tile(aes(alpha = as.factor(Freq))) + 
#   theme_minimal() + 
#   theme(axis.text = element_blank(), legend.title = element_blank(), 
#         legend.position = 'top') +
#   scale_fill_manual(values = c('#b2bc63','#10342d'), labels = c("Abscent", "Present")) +
#   scale_alpha_manual(values = c(.6, 1), labels = c("Abscent", "Present")) + 
#   labs(x = "Environment (Year-Location combination)", y = "Hybrid")

dat = droplevels(dat[which(dat$year > 2021 & dat$loc %in% c("GAH1","IAH1","IAH2","ILH1", "INH1", "MNH1", "MOH1", "MOH2")),])
dat = dat[which(dat$gen %in% rownames(geno.code)),]
set.seed(7)
dat = dat[which(dat$gen %in% sample(unique(dat$gen),150)),]
geno.code = geno.code[which(rownames(geno.code) %in% dat$gen),]
all(unique(dat$gen) %in% rownames(geno.code))

# aux = sort(rowSums(table(dat$gen, dat$env)), TRUE)
# aux = names(aux[which(aux >= 15)])
# set.seed(4588)
# dat = dat[which(dat$gen %in% sample(aux,150)),]

ggplot(data = dat, aes(
  x = reorder(env, blue),
  y = reorder(gen, blue),
  fill = blue
)) +
  geom_tile() +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 90),
        legend.position = "top", legend.title = element_blank()) +
  scale_fill_viridis_c(option = "inferno") +
  labs(x = "Environment", y = "Hybrid")

# G matrix ---------------------------------------------------------
Gmat = Gmatrix(SNPmatrix = geno.code, method = "VanRaden")

Heatmap(
  Gmat,
  col = viridis::inferno(20),
  show_row_names = F,
  show_column_names = F,
  heatmap_legend_param = list(title = "IBS")
)

pc = prcomp(Gmat) 
pc.df = rbind(
  as.data.frame(pc$rotation) |> 
    rownames_to_column("ID") |> 
    mutate(comp = "load")
)
# ggplot(dat = pc.df, aes(x = PC1, y = PC2)) + 
#   geom_point(size = 3, alpha = .5, colour = "#b2bc63") +
#   geom_vline(xintercept = 0, linetype = "dashed") + 
#   geom_hline(yintercept = 0, linetype = 'dashed') +
#   coord_fixed() + 
#   theme_minimal() + 
#   theme(text= element_text(size = 18)) + 
#   labs(x = paste0("PC1 (", round(pc$sdev[1]^2/sum(pc$sdev^2)*100,2), "%)"),
#        y= paste0("PC2 (", round(pc$sdev[2]^2/sum(pc$sdev^2)*100,2), "%)"))

plot_ly(
  pc.df,
  x = ~ PC1,
  y = ~ PC2,
  z = ~ PC3
)  |> add_markers(
  text = ~ paste0("Hybrid: ", ID),
  hoverinfo = "text",
  marker = list(color = '#b2bc63', opacity = .62)
) |>  layout(scene = list(
  xaxis = list(
    title = paste0("PC1 (", round(pc$sdev[1]^2/sum(pc$sdev^2)*100,2), "%)"),
    showgrid = TRUE,
    range = c(min(c(
      range(pc.df$pc1), range(pc.df$pc2), range(pc.df$pc3)
    )), max(c(
      range(pc.df$pc1), range(pc.df$pc2), range(pc.df$pc3)
    )))
  ),
  yaxis = list(
    title = paste0("PC2 (", round(pc$sdev[2]^2/sum(pc$sdev^2)*100,2), "%)"),
    showgrid = TRUE,
    range = c(min(c(
      range(pc.df$pc1), range(pc.df$pc2), range(pc.df$pc3)
    )), max(c(
      range(pc.df$pc1), range(pc.df$pc2), range(pc.df$pc3)
    )))
  ),
  zaxis = list(
    title = paste0("PC3 (", round(pc$sdev[3]^2/sum(pc$sdev^2)*100,2), "%)"),
    showgrid = TRUE,
    range = c(min(c(
      range(pc.df$pc1), range(pc.df$pc2), range(pc.df$pc3)
    )), max(c(
      range(pc.df$pc1), range(pc.df$pc2), range(pc.df$pc3)
    )))
  )
))

# Model -------------------------------------------------------------------

## Adjusting the dataset to consider missing values

nrow(dat)
grid = expand.grid(
  gen = unique(dat$gen),
  env = unique(dat$env),
  stringsAsFactors = FALSE
)
dat = merge(grid, dat, by = c("gen", "env"), all.x = TRUE)
nrow(dat)

## Building the covariance matrices
dat$gen = factor(dat$gen, levels = rownames(Gmat), ordered = TRUE)
dat = dat[order(dat$env, dat$gen),]
ZG = model.matrix(~-1 + gen, data = dat)
rownames(ZG) = dat$gen
colnames(ZG) = gsub("gen", "", colnames(ZG))
ZG[1:5, 1:5]

ZE = model.matrix(~-1 + env, data = dat)
rownames(ZE) = dat$env
colnames(ZE) = gsub("env", "", colnames(ZE))
ZE[1:5,1:5]

ZGZ = ZG %*% Gmat %*% t(ZG)
rownames(ZGZ) -> aux1
dim(ZGZ)
ZEZ = tcrossprod(ZE)
dim(ZEZ)
GEI = ZGZ * ZEZ
rownames(GEI) -> aux3
dim(GEI)
GEI = eigen(GEI)
rownames(GEI$vectors) = aux3
save(GEI, file = "ETA_GEI.RDA")
ZGZ = eigen(ZGZ)
rownames(ZGZ$vectors) = aux1
save(ZGZ, file = "ETA_ZGZ.RDA")

ETA = list(
  env = list(X = ZE, model = 'FIXED'),
  hyb = list(V = ZGZ$vectors, d = ZGZ$values, model = 'RKHS'),
  gei = list(V = GEI$vectors, d = GEI$values, model = 'RKHS')
)

model1 = BGLR(y = dat$blue, ETA = ETA, nIter = 12000, burnIn = 2000, 
              thin = 10, verbose = TRUE, saveAt = "complete_met1_")
save(model1, file = "met_gs1.RDA")


## Cholesly decomposition
L_G = t(chol(Gmat + diag(x=0.0001, nrow = nrow(Gmat))))
X_g = ZG %*% L_G
X_ge = model.matrix(~ -1 + X_g:dat$env)

ETA2 = list(
  env = list(X = ZE, model = 'FIXED'),
  hyb = list(X = X_g, model = 'BRR'),
  gei = list(X = X_ge, model = 'BRR')
)

model2 = BGLR(y = dat$blue, ETA = ETA2, nIter = 12000, burnIn = 2000, 
             thin = 10, verbose = TRUE, saveAt = "complete_met2_")
save(model2, file = "met_gs2.RDA")


#### Comparing both models
model1$varE
model2$varE

var_g1  = model1$ETA$hyb$varU
var_g2  = model2$ETA$hyb$varU

var_ge1 <- model1$ETA$gei$varU
var_ge2 <- model2$ETA$gei$varU

var_e1  <- model1$varE
var_e2  <- model2$varE

var_total1 <- var_g1 + var_ge1 + var_e1
var_total2 <- var_g2 + var_ge2 + var_e2


res = data.frame(
  gen = dat$gen,
  env = dat$env,
  blue = dat$blue,
  pred1 = model1$yHat,
  g1 = model1$ETA$hyb$u,
  ge1 = model1$ETA$gei$u,
  pred2 = model2$yHat
  # g2 = model2$ETA$hyb$b,
  # ge2 = model2$ETA$gei$b
)

cor(res$pred1, res$pred2)
# cor(res$g1, res$g2)
# cor(res$ge1, res$ge2)

gemat1 <- tapply(res$pred1, INDEX = list(res$gen, res$env), FUN = mean)
gemat2 <- tapply(res$pred2, INDEX = list(res$gen, res$env), FUN = mean)

model1$fit$DIC
model2$fit$DIC

cor(res$blue, res$pred1, use = "complete.obs")
cor(res$blue, res$pred2, use = "complete.obs")


ggplot(data = res, aes(x = gen, y = g1))+
  geom_point(alpha = .75, colour = "#10342d") + 
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())+
  labs(x = "Hybrid", y = "GEBV")

ggplot(data = res, aes(x = gen, y = env, fill = pred1))+ 
  geom_tile(alpha = .75) + 
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank(),
        legend.position = "top")+
  labs(x = "Hybrid", y = "Environment", fill = "GEBV") +
  scale_fill_viridis_c(option = "inferno") 


# Checking the convergence
library(coda)

# Lendo os arquivos de traço da MCMC salvos no disco
trace_varE   <- read.table("complete_met1_varE.dat")[,1]
trace_varG   <- read.table("complete_met1_ETA_hyb_varU.dat")[,1]
trace_varGEI <- read.table("complete_met1_ETA_gei_varU.dat")[,1]

# Gráfico da convergência da variância genética
plot(trace_varG, type = "l", main = "Convergência de Var(G)", ylab = "Variância", xlab = "Iteração (thinned)")

varG_chain = mcmc(trace_varG[201:length(trace_varG)])
plot(varG_chain)

effectiveSize(varG_chain)

# 5. Teste de Convergência de Geweke (p-valor deve ser > 0.05)
geweke.diag(varG_chain)


# Cross-validation --------------------------------------------------------
nfolds = 10
nrept = 5 
seed = 8  

# CV2 ==========
cvdata = list()
seed = 7
for (rept in 1:nrept) {
  set.seed(seed * rept)
  cvdata[[rept]] = dat
  cvdata[[rept]]$set = NA
  for (id in unique(dat$gen)) {
    cvdata[[rept]][cvdata[[rept]]$gen == id, 'set'] = sample(1:nfolds,
                                                             size = dim(cvdata[[rept]][cvdata[[rept]]$gen == id, ])[1],
                                                             replace = dim(cvdata[[rept]][cvdata[[rept]]$gen == id, ])[1] > nfolds)
  }
  cvdata[[rept]]$rept = rept
}

cv2 = lapply(cvdata, function(x){
  res.list = list()
  for (i in unique(x$set)) {
    
    yNA = x$blue
    yNA[x$set == i] = NA
    
    mod_cv = BGLR(y = yNA, ETA = ETA2, nIter = 8000, burnIn = 1000, 
                  thin = 10, verbose = FALSE, saveAt = "modCV")
    unlink(list.files(pattern = "modCV"))
    
    res.list[[i]] = data.frame(
      gen = x$gen,
      env = x$env,
      blue = x$blue,
      set = x$set,
      yNA = yNA,
      yhat = mod_cv$yHat
    ) |> filter(set == i)
    
  }
  res.list
})

save(cv2, file = "CV2.RDA")

# CV1 ==========
set.seed(7)
sets = split(
  rep(
    1:nfolds, length(unique(dat$gen)) * nrept
  )[order(runif(length(unique(dat$gen)) * nrept))],
  f = 1:nrept
)
cvdata = lapply(sets, function(x){
  cvdata = dat
  cvdata = merge(cvdata, data.frame(
    gen = unique(dat$gen),
    set = x
  ), by = 'gen')
})
for (i in 1:length(cvdata)) cvdata[[i]]$rept = i


cv1 = lapply(cvdata, function(x){
  res.list = list()
  for (i in unique(x$set)) {
    
    yNA = x$blue
    yNA[x$set == i] = NA
    
    mod_cv = BGLR(y = yNA, ETA = ETA2, nIter = 8000, burnIn = 1000, 
                  thin = 10, verbose = FALSE, saveAt = "modCV")
    unlink(list.files(pattern = "modCV"))
    
    res.list[[i]] = data.frame(
      gen = x$gen,
      env = x$env,
      blue = x$blue,
      set = x$set,
      yNA = yNA,
      yhat = mod_cv$yHat
    ) |> filter(set == i)
    
  }
  res.list
})

save(cv1, file = "CV1.RDA")









