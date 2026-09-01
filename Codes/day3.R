rm(list=ls())

# Library -----------------------------------------------------------------
library(tidyverse)
library(BGLR)
library(AGHmatrix)
library(ComplexHeatmap)
library(plotly)
library(circlize)
library(coda)

# Phenotypic and genomic data ---------------------------------------------------------
dat = readRDS("Data/G2F.RDS")
str(dat)
geno.code = readRDS(file = "Data/mrk_filtered.RDS")

dat = droplevels(dat[which(dat$year > 2021 & dat$loc %in% c("GAH1","IAH1","IAH2","ILH1", "INH1", "MNH1", "MOH1", "MOH2")),])
dat = dat[which(dat$gen %in% rownames(geno.code)),]
set.seed(7)
dat = dat[which(dat$gen %in% sample(unique(dat$gen),150)),]
geno.code = geno.code[which(rownames(geno.code) %in% dat$gen),]
all(unique(dat$gen) %in% rownames(geno.code))

# Environmental data ------------------------------------------------------
weather = read.csv(file = "Data/weather_seasons.csv")
soil = read.csv(file = "Data/Soil.csv")

# Not all environments have soil data, so let's keep only those that have
dat = dat[which(dat$env %in% soil$Env),]
weather = weather[which(weather$Env %in% dat$env),]
soil = soil[which(soil$Env %in% dat$env),]

# Let's take a quick look and the environmental profile of each environment
soil = soil[,-which(apply(soil, 2, function(x) all(is.na(x))))]
soil = soil[,which(!colnames(soil) %in% c("LabID", "Date.Received", "Date.Reported", "E.Depth", "Texture", "Comments"))]

ggplot(data = soil |> pivot_longer(colnames(soil[-c(1,2)])), aes(x = Env, y = value)) +
  facet_wrap(.~name, scales = "free_y") +
  geom_col() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust= .5, hjust = 1))

soil = soil[,which(!colnames(soil) %in% c("X.Na.Sat", "X.H.Sat",
                                          "Texture.No"))]

ggplot(data = soil |> pivot_longer(colnames(soil[-c(1,2)])), aes(x = Env, y = value)) +
  facet_wrap(.~name, scales = "free_y") +
  geom_col() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust= .5, hjust = 1), text = element_text(size = 15))


weather = weather[which(weather$Env %in% dat$env),]
weather$Date = lubridate::ymd(weather$Date)
# Sys.setlocale("LC_TIME", "English") 
ggplot(data = weather |> pivot_longer(colnames(weather[-c(1, 2)])) |>
         filter(Env == "GAH1_2022"),
       aes(x = Date, y = value)) +
  facet_wrap(name ~ ., scales = "free_y") +
  geom_point(color = '#10342d', alpha = .75) +
  geom_line(aes(group = 1), color = '#10342d', alpha = .8) +
  theme_bw() +
  theme(text = element_text(size = 15))

## Checking the range of each variable
apply(weather, 2, function(x) range(x))

## Extracting the environmental markers (summary statistics)
envdata = list()
for (i in colnames(weather)[-c(1:2, 18)]) {
  temp = do.call(rbind, tapply(weather[,i], weather$Env, summary))
  colnames(temp) = c("min", "1qrt", "median", "mean", "3qrt", 'max')
  colnames(temp) = paste(i, colnames(temp), sep = "_")
  envdata[[i]] = temp
  rm(temp)
}
envdata = do.call(cbind, envdata)
envdata = cbind(envdata, 
                PRECTOT_sum = tapply(weather$PRECTOTCORR, weather$Env, sum),
                PRECTOT_mean = tapply(weather$PRECTOTCORR, weather$Env, mean),
                soil[match(rownames(envdata), soil$Env),-c(1,2)])

dim(envdata)

## Overview of the matrix of environmental markers
any(is.na(envdata))
envdata = scale(envdata)

col_fun = colorRamp2(c(min(envdata), 0, max(envdata)), c("blue4", "white", "red4"))

Heatmap(matrix = envdata,
        show_row_dend = FALSE, show_column_dend = FALSE, 
        column_order = colnames(envdata),
        column_split = ifelse(grepl("\\.", colnames(envdata)), "SOIL", gsub("_.*", "", colnames(envdata))),
        show_column_names = FALSE,
        column_title_rot = 45,
        border = TRUE,
        border_gp = gpar(col = "black", lty = 2),
        col = col_fun,
        show_heatmap_legend = FALSE)


## Building the environmental linear kernel
Emat = tcrossprod(envdata)/ncol(envdata)
Emat

col_fun = colorRamp2(c(min(Emat), 0, max(Emat)), c("blue4", "white", "green4"))
Heatmap(matrix = Emat,
        col = col_fun,
        show_heatmap_legend = TRUE,
        heatmap_legend_param = list(title = "Relationship"))

pc = prcomp(Emat) 
pc.df = rbind(
  as.data.frame(pc$rotation) |> 
    rownames_to_column("ID") |> 
    mutate(comp = "load")
)

plot_ly(
  pc.df,
  x = ~ PC1,
  y = ~ PC2,
  z = ~ PC3
)  |> add_markers(
  text = ~ paste0("Environment: ", ID),
  hoverinfo = "text",
  marker = list(color = '#b2bc63', opacity = .62)
) |>  layout(scene = list(
  xaxis = list(
    title = paste0("PC1 (", round(pc$sdev[1]^2/sum(pc$sdev^2)*100,2), "%)"),
    showgrid = TRUE,
    range = c(min(c(
      range(pc.df$PC1), range(pc.df$PC2), range(pc.df$PC3)
    )), max(c(
      range(pc.df$PC1), range(pc.df$PC2), range(pc.df$PC3)
    )))
  ),
  yaxis = list(
    title = paste0("PC2 (", round(pc$sdev[2]^2/sum(pc$sdev^2)*100,2), "%)"),
    showgrid = TRUE,
    range = c(min(c(
      range(pc.df$PC1), range(pc.df$PC2), range(pc.df$PC3)
    )), max(c(
      range(pc.df$PC1), range(pc.df$PC2), range(pc.df$PC3)
    )))
  ),
  zaxis = list(
    title = paste0("PC3 (", round(pc$sdev[3]^2/sum(pc$sdev^2)*100,2), "%)"),
    showgrid = TRUE,
    range = c(min(c(
      range(pc.df$PC1), range(pc.df$PC2), range(pc.df$PC3)
    )), max(c(
      range(pc.df$PC1), range(pc.df$PC2), range(pc.df$PC3)
    )))
  )
))


# Model -------------------------------------------------------------------
Gmat = Gmatrix(SNPmatrix = geno.code, method = "VanRaden")


nrow(dat)
grid = expand.grid(
  gen = unique(dat$gen),
  env = unique(dat$env),
  stringsAsFactors = FALSE
)
dat = merge(grid, dat, by = c("gen", "env"), all.x = TRUE)
nrow(dat)

ggplot(data = dat, aes(
  x = reorder(env, blue),
  y = reorder(gen, blue),
  fill = blue
)) +
  geom_tile() +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(angle = 90),
        legend.position = "top", legend.title = element_blank(), 
        text = element_text(size = 18)) +
  scale_fill_viridis_c(option = "inferno", na.value = "lightgrey") +
  labs(x = paste("Environments: ", length(unique(dat$env))),
       y = paste("Hybrid: ", length(unique(dat$gen))))

## Building the covariance matrices
dat$gen = factor(dat$gen, levels = rownames(Gmat), ordered = TRUE)
dat$env = factor(dat$env, levels = rownames(Emat), ordered = TRUE)
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
aux1 = rownames(ZGZ)
dim(ZGZ)

ZEZ = tcrossprod(ZE)
dim(ZEZ)

ZOZ = ZE %*% Emat %*% t(ZE)
dim(ZEZ)
aux2 = rownames(ZOZ)

GEI = ZGZ * ZEZ
aux3 = rownames(GEI)
dim(GEI)

GEI_EC = ZGZ * ZOZ
aux4 = rownames(GEI_EC)
dim(GEI_EC)

## Eigendecomposition
ZGZ = eigen(ZGZ)
rownames(ZGZ$vectors) = aux1
ZEZ = eigen(ZEZ)
rownames(ZEZ$vectors) = aux2
ZOZ = eigen(ZOZ)
GEI = eigen(GEI)
rownames(GEI$vectors) = aux3
GEI_EC = eigen(GEI_EC)
rownames(GEI_EC$vectors) = aux4

## Cholesky fatoration
L_G = t(chol(Gmat + diag(x=0.0001, nrow = nrow(Gmat))))
L_E = t(chol(Emat))
X_g = ZG %*% L_G
X_e = ZE
X_ec = ZE %*% L_E
X_ge = model.matrix(~ -1 + X_g:dat$env)
X_ge_ec = do.call(cbind, lapply(1:ncol(X_ec), function(j) X_g * X_ec[, j]))


## Model 1: only genotype x environmental marker interaction
ETA1_kernel = list(
  ec = list(V = ZOZ$vectors, d = ZOZ$values, model = 'RKHS'),
  hyb = list(V = ZGZ$vectors, d = ZGZ$values, model = 'RKHS'),
  gei_ec = list(V = GEI_EC$vectors, d = GEI_EC$values, model = 'RKHS')
)

a = Sys.time()
model1 = BGLR(y = dat$blue, ETA = ETA1_kernel, nIter = 12000, burnIn = 2000, 
              thin = 10, verbose = TRUE, saveAt = "complete_met_ec_kernel")
b = Sys.time()
model1$ellapsed = b-a
save(model1, file = "met_ec_kernel.RDA")


ETA1_brr = list(
  ec = list(X = X_ec, model = 'BRR'),
  hyb = list(X = X_g, model = 'BRR'),
  gei_ec = list(X = X_ge_ec, model = 'BRR')
)

a = Sys.time()
model2 = BGLR(y = dat$blue, ETA = ETA1_brr, nIter = 12000, burnIn = 2000, 
              thin = 10, verbose = TRUE, saveAt = "complete_met_ec_brr")
b = Sys.time()
model2$ellapsed = b-a
save(model2, file = "met_ec_brr.RDA")


# Model 2: GxEC and GxE
ETA2_kernel = list(
  ec =  list(V = ZOZ$vectors, d = ZOZ$values, model = 'RKHS'),
  hyb = list(V = ZGZ$vectors, d = ZGZ$values, model = 'RKHS'),
  gei = list(V = GEI$vectors, d = GEI$values, model = 'RKHS'),
  gei_ec = list(V = GEI_EC$vectors, d = GEI_EC$values, model = 'RKHS')
)

a = Sys.time()
model3 = BGLR(y = dat$blue, ETA = ETA2_kernel, nIter = 12000, burnIn = 2000, 
              thin = 10, verbose = TRUE, saveAt = "complete_met_geiec_kernel")
b = Sys.time()
model3$ellapsed = b-a
save(model3, file = "met_geiec_kernel.RDA")


ETA2_brr = list(
  ec = list(X = X_ec, model = 'BRR'),
  hyb = list(X = X_g, model = 'BRR'),
  gei = list(X = X_ge, model = 'BRR'),
  gei_ec = list(X = X_ge_ec, model = 'BRR')
)

a = Sys.time()
model4 = BGLR(y = dat$blue, ETA = ETA2_brr, nIter = 12000, burnIn = 2000, 
              thin = 10, verbose = TRUE, saveAt = "complete_met_geiec_brr")
b = Sys.time()
model4$ellapsed = b-a
save(model4, file = "met_geiec_brr.RDA")


# Ellapsed time
ellapsed = data.frame(
  model = rep(c("GxEC", "GxEC + GxE"), each = 2),
  mechanism = rep(c("RKHS", "CHOL"), times = 2),
  ellapsed = -1*c(model1$ellapsed, model2$ellapsed, model3$ellapsed, model4$ellapsed)
)

ggplot(data = ellapsed, aes(x = model, y = ellapsed, fill = mechanism)) +
  geom_col(position = position_dodge(), color = "black") + 
  theme_bw() + 
  theme(text = element_text(size = 18), legend.position = 'top', 
        legend.title = element_blank()) +
  scale_fill_manual(values = c('#b2bc63','#10342d')) + 
  labs(x = "Model", y = "Ellapsed time (min)")

# Model parameters

modparam = data.frame(
  residual = c(model1$varE, model2$varE, model3$varE, model4$varE),
  ec = c(
    model1$ETA$ec$varU,
    model2$ETA$ec$varB,
    model3$ETA$ec$varU,
    model4$ETA$ec$varB
  ),
  g = c(
    model1$ETA$hyb$varU,
    model2$ETA$hyb$varB,
    model3$ETA$hyb$varU,
    model4$ETA$hyb$varB
  ),
  gxec = c(
    model1$ETA$gei_ec$varU,
    model2$ETA$gei_ec$varB,
    model3$ETA$gei_ec$varU,
    model4$ETA$gei_ec$varB
  ),
  gxe = c(NA, NA, model3$ETA$gei$varU, model4$ETA$gei$varB),
  model = rep(c("GxEC", "GxEC + GxE"), each = 2),
  mechanism = rep(c("RKHS", "CHOL"), times = 2)
)

modparam |> mutate(control = paste(model, mechanism, sep = "_")) %>% 
  group_by(control) |> 
  mutate(varT = ifelse(grepl("GxE_", control),residual + ec + g + gxec + gxe,
                       residual + ec + g + gxec),
         PvR = residual/varT,
         PvEC = ec/varT,
         PvG = g/varT,
         PvGxEC = gxec/varT,
         PvGxE = gxe/varT) |> 
  pivot_longer(PvR:PvGxE) |> 
  ggplot(aes(x = model, y = value)) +
  facet_wrap(.~mechanism) +
  geom_col(aes(fill = name), color = "black") + 
  theme_bw() + 
  labs(x = "Model", y = "Percentage from total variance") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_viridis_d(option="mako", 
                       labels = c('PvR' = "varE", 'PvG' = "varG", 'PvGxE' = "varGEI",
                                  "PvGxEC" = "varE", "PvEC" = "varEC")) +
  theme(text = element_text(size = 18), legend.position = "top", 
        legend.title = element_blank())

# Predictions

res = data.frame(gen = dat$gen, env = dat$env, blue = dat$blue,
                 'GxEC_RKHS' = model1$yHat, 'GxEC_CHOL' = model2$yHat,
                 'GxE_GxEC_RKHS' = model3$yHat, 'GxE_GxEC_CHOL' = model4$yHat)

ggplot(data = res, aes(x = GxEC_RKHS, y = GxEC_CHOL)) + 
  geom_point(alpha = .6, color = '#10342d') +
  labs(x = "RKHS", y = "CHOL (BRR)",
       caption = paste("Correlation = ", round(cor(res$GxEC_RKHS, res$GxEC_CHOL), 4))) + 
  theme_bw() +
  theme(text = element_text(size = 18))

ggplot(data = res, aes(x = GxE_GxEC_RKHS, y = GxE_GxEC_CHOL)) + 
  geom_point(alpha = .6, color = '#10342d') +
  labs(x = "RKHS", y = "CHOL (BRR)",
       caption = paste("Correlation = ", round(cor(res$GxE_GxEC_RKHS, res$GxE_GxEC_CHOL), 4))) + 
  theme_bw() +
  theme(text = element_text(size = 18))

GGally::ggpairs(res[,-c(1:3)]) + theme_bw() + theme(text = element_text(size = 16))

# Convergence

trace_varR = read.table("complete_met_geiec_brrvarE.dat")[,1]
trace_varG = read.table("complete_met_geiec_brrETA_hyb_varB.dat")[,1]
trace_varGEI = read.table("complete_met_geiec_brrETA_gei_varB.dat")[,1]
trace_varEC = read.table("complete_met_geiec_brrETA_ec_varB.dat")[,1]
trace_varGEI_EC = read.table("complete_met_geiec_brrETA_gei_ec_varB.dat")[,1]

varG_chain = mcmc(trace_varG[201:length(trace_varG)])
plot(varG_chain)
effectiveSize(varG_chain)
geweke.diag(varG_chain)

varR_chain = mcmc(trace_varR[201:length(trace_varR)])
plot(varR_chain)
effectiveSize(varR_chain)
geweke.diag(varR_chain)

varEC_chain = mcmc(trace_varEC[201:length(trace_varEC)])
plot(varEC_chain)
effectiveSize(varEC_chain)
geweke.diag(varEC_chain)

varGEI_EC_chain = mcmc(trace_varGEI_EC[201:length(trace_varGEI_EC)])
plot(varGEI_EC_chain)
effectiveSize(varGEI_EC_chain)
geweke.diag(varGEI_EC_chain)

varGEI_chain = mcmc(trace_varGEI[201:length(trace_varGEI)])
plot(varGEI_chain)
effectiveSize(varGEI_chain)
geweke.diag(varGEI_chain)

# Cross-validation --------------------------------------------------------

## CV0: Leave-one-out -----------------------------------------------------
cv0=list()
for (i in levels(dat$env)) {
  
  yNA = dat$blue
  yNA[dat$env == i] = NA
  
  mod1_cv = BGLR(y = yNA, ETA = ETA1_brr, nIter = 8000, burnIn = 1000, 
                thin = 10, verbose = FALSE, saveAt = "modCV")
  unlink(list.files(pattern = "modCV"))
  mod2_cv = BGLR(y = yNA, ETA = ETA2_brr, nIter = 8000, burnIn = 1000, 
                 thin = 10, verbose = FALSE, saveAt = "modCV")
  unlink(list.files(pattern = "modCV"))
  
  cv0[[i]] = data.frame(
    gen = dat$gen,
    env = dat$env,
    blue = dat$blue,
    yNA = yNA,
    yhat_GxEC = mod1_cv$yHat,
    yhat_GxE_GxEC = mod2_cv$yHat
  ) |> filter(env == i)
}
save(cv0, file = "CV0.RDA")


## CV00 ========
nfolds = 5
nrept = 10
seed = 8  

sets = list()
i = 1
repeat{
  set.seed(7 * i)
  sets[[i]] = sample(rep(1:nfolds, length.out = nlevels(dat$gen)))
  i = i + 1
  if(i > nrept) break
}
cvdata = lapply(sets, function(x){
  cvdata = dat
  aux = data.frame(gen = levels(dat$gen), set = x)
  cvdata$set = aux$set[match(cvdata$gen, aux$gen)]
  return(cvdata)
})
for (i in 1:length(cvdata)) cvdata[[i]]$rept = i


cv00 = lapply(cvdata, function(x){
  res.env = list()
  for (j in levels(x$env)) {
    res.list = list()
    for (i in unique(x$set)) {
      yNA = x$blue
      yNA[x$env == j | x$set == i] = NA
      
      mod1_cv = BGLR(y = yNA, ETA = ETA1_brr, nIter = 8000, burnIn = 1000, 
                     thin = 10, verbose = FALSE, saveAt = "modCV")
      unlink(list.files(pattern = "modCV"))
      mod2_cv = BGLR(y = yNA, ETA = ETA2_brr, nIter = 8000, burnIn = 1000, 
                     thin = 10, verbose = FALSE, saveAt = "modCV")
      unlink(list.files(pattern = "modCV"))
      
      res.list[[i]] = data.frame(
        gen = x$gen,
        env = x$env,
        set = x$set,
        blue = x$blue,
        yNA = yNA,
        yhat_GxEC = mod1_cv$yHat,
        yhat_GxE_GxEC = mod2_cv$yHat
      ) |> filter(env == j & set == i)
    }
    res.env[[j]] = do.call(rbind, res.list)
  }
  res.env
})
save(cv00, file = "CV00.RDA")


# CV2 ==========

sets = list()
i = 1
repeat{
  set.seed(987 * i)
  sets[[i]] = sample(rep(1:nfolds, length.out = nrow(dat)))
  i = i + 1
  if(i > nrept) break
}
cvdata = lapply(sets, function(x){
  cvdata = dat
  cvdata$set = x
  return(cvdata)
})
for (i in 1:length(cvdata)) cvdata[[i]]$rept = i


cv2 = lapply(cvdata, function(x){
  res.list = list()
  for (i in unique(x$set)) {
    
    yNA = x$blue
    yNA[x$set == i] = NA
    
    mod1_cv = BGLR(y = yNA, ETA = ETA1_brr, nIter = 8000, burnIn = 1000, 
                   thin = 10, verbose = FALSE, saveAt = "modCV")
    unlink(list.files(pattern = "modCV"))
    mod2_cv = BGLR(y = yNA, ETA = ETA2_brr, nIter = 8000, burnIn = 1000, 
                   thin = 10, verbose = FALSE, saveAt = "modCV")
    unlink(list.files(pattern = "modCV"))
    
    res.list[[i]] = data.frame(
      gen = x$gen,
      env = x$env,
      blue = x$blue,
      set = x$set,
      yNA = yNA,
      yhat_GxEC = mod1_cv$yHat,
      yhat_GxE_GxEC = mod2_cv$yHat
    ) |> filter(set == i)
    
  }
  res.list
})

cv2 = lapply(cv2, function(x) do.call(rbind, x))
save(cv2, file = "CV2_EC.RDA")

# CV1 ==========
sets = list()
i = 1
repeat{
  set.seed(7 * i)
  sets[[i]] = sample(rep(1:nfolds, length.out = nlevels(dat$gen)))
  i = i + 1
  if(i > nrept) break
}
cvdata = lapply(sets, function(x){
  cvdata = dat
  aux = data.frame(gen = levels(dat$gen), set = x)
  cvdata$set = aux$set[match(cvdata$gen, aux$gen)]
  return(cvdata)
})
for (i in 1:length(cvdata)) cvdata[[i]]$rept = i


cv1 = lapply(cvdata, function(x){
  res.list = list()
  for (i in unique(x$set)) {
    
    yNA = x$blue
    yNA[x$set == i] = NA
    
    mod1_cv = BGLR(y = yNA, ETA = ETA1_brr, nIter = 8000, burnIn = 1000, 
                   thin = 10, verbose = FALSE, saveAt = "modCV")
    unlink(list.files(pattern = "modCV"))
    mod2_cv = BGLR(y = yNA, ETA = ETA2_brr, nIter = 8000, burnIn = 1000, 
                   thin = 10, verbose = FALSE, saveAt = "modCV")
    unlink(list.files(pattern = "modCV"))
    
    res.list[[i]] = data.frame(
      gen = x$gen,
      env = x$env,
      blue = x$blue,
      set = x$set,
      yNA = yNA,
      yhat_GxEC = mod1_cv$yHat,
      yhat_GxE_GxEC = mod2_cv$yHat
    ) |> filter(set == i)
    
  }
  res.list
})
cv1 = lapply(cv1, function(x) do.call(rbind, x))
save(cv1, file = "CV1_EC.RDA")







