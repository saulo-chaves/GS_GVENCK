rm(list=ls())

# Library -----------------------------------------------------------------
library(tidyverse)
library(rrBLUP)
library(vcfR)
library(BGLR)

# Phenotypic data ---------------------------------------------------------
dat = readRDS("Data/G2F.RDS")
str(dat)

ggplot(data = as.data.frame(table(dat$gen, dat$env)),
       aes(x = Var2, y = Var1, fill = as.factor(Freq))) + 
  geom_tile(aes(alpha = as.factor(Freq))) + 
  theme_minimal() + 
  theme(axis.text = element_blank(), legend.title = element_blank(), 
        legend.position = 'top') +
  scale_fill_manual(values = c('#b2bc63','#10342d'), labels = c("Abscent", "Present")) +
  scale_alpha_manual(values = c(.6, 1), labels = c("Abscent", "Present")) + 
  labs(x = "Environment (Year-Location combination)", y = "Hybrid")

## Choose only one trial
dat = dat[which(dat$env == "WIH3_2023"),]

# Genomic data ---------------------------------------------------------
vcf = read.vcfR(file = "Data/markers.vcf")
geno = extract.gt(x = vcf)
geno = geno[-grep(',', vcf@fix[,'ALT']),]  # 42 multiallelic marker were removed
vcf@fix = vcf@fix[-grep(',',vcf@fix[,'ALT']),]  # Removing from the fix file
geno = gsub('/', "", geno)
geno[1:5,1:5]

geno = cbind(geno, vcf@fix[,c('REF', 'ALT')])  # Binding reference and alternative allele to marker matrix
geno = apply(geno, 1, function(x){
  y = gsub('0', x['REF'], x)
  z = gsub('1', y['ALT'], y)
  z
})

geno = geno[-which(rownames(geno) %in% c('REF', 'ALT')), ]
geno[grep('CA', geno)] = 'AC'
geno[grep('GA', geno)] = 'AG'
geno[grep('TA', geno)] = 'AT'
geno[grep('GC', geno)] = 'CG'
geno[grep('TC', geno)] = 'CT'
geno[grep('TG', geno)] = 'GT'

# Genotypic and allele frequencies before filtering

Freqsnps = apply(geno, 2, table)
Freqsnps = do.call(rbind, lapply(Freqsnps, data.frame)) |> 
  rownames_to_column("mrk") |> mutate(mrk = gsub('\\..*', '', mrk))

Freqall = apply(geno, 2, function(x){
  aux = as.data.frame(table(do.call(c, strsplit(x = x, split = ''))))
  return(aux)
})
Freqall = do.call(rbind, Freqall) |> 
  rownames_to_column("mrk") |> mutate(mrk = gsub('\\..*', '', mrk))

Freqmis = apply(geno, 2, function(x) mean(is.na(x)))

ggpubr::ggarrange(
  data.frame(Freq = tapply(Freqsnps$Freq, Freqsnps$Var1, sum)/sum(Freqsnps$Freq)) |> 
    rownames_to_column("Genotype") |> ggplot(aes(x = Genotype, y = Freq)) + 
    geom_bar(stat = 'identity', fill = '#b2bc63', color = 'black') +
    theme_bw() +
    scale_y_continuous(labels = scales::percent) + 
    labs(x = 'Genotype', y = 'Frequency'),
  
  data.frame(Freq = tapply(Freqall$Freq, Freqall$Var1, sum)/sum(Freqall$Freq)) |> 
    rownames_to_column("Alleles") |> ggplot(aes(x = Alleles, y = Freq)) + 
    geom_bar(stat = 'identity', fill = '#b2bc63', color = 'black') +
    theme_bw() +
    scale_y_continuous(labels = scales::percent) + 
    labs(x = 'Alleles', y = 'Frequency')
)

# Removing markers with more than 5% of missing data
geno = geno[,-which(Freqmis >= 0.05)]  
vcf@fix = vcf@fix[-which(Freqmis >= 0.05),]  

# Imputing missing data with the mode
for (i in colnames(geno)) {
  if(!any(is.na(geno[,i]))) next else{
    geno[which(is.na(geno[,i])),i] = as.character(Freqsnps[which(Freqsnps$mrk == i), ][which.max(Freqsnps[which(Freqsnps$mrk == i), ][,'Freq']), "Var1"])
  }
  rm(i)
}

# Removing markers with MAF < 0.05
Freqall$Fper = Freqall$Freq / (nrow(geno) * 2)
geno = geno[,-which(Freqall$Fper <= 0.05)]

# Genotypic and allele frequencies after filtering
Freqsnps = apply(geno, 2, table)
Freqsnps = do.call(rbind, lapply(Freqsnps, data.frame)) |> 
  rownames_to_column("mrk") |> mutate(mrk = gsub('\\..*', '', mrk))

Freqall = apply(geno, 2, function(x){
  aux = as.data.frame(table(do.call(c, strsplit(x = x, split = ''))))
  return(aux)
})
Freqall = do.call(rbind, Freqall) |> 
  rownames_to_column("mrk") |> mutate(mrk = gsub('\\..*', '', mrk))

ggpubr::ggarrange(
  data.frame(Freq = tapply(Freqsnps$Freq, Freqsnps$Var1, sum)/sum(Freqsnps$Freq)) |> 
    rownames_to_column("Genotype") |> ggplot(aes(x = Genotype, y = Freq)) + 
    geom_bar(stat = 'identity', fill = '#b2bc63', color = 'black') +
    theme_bw() +
    scale_y_continuous(labels = scales::percent) + 
    labs(x = 'Genotype', y = 'Frequency'),
  
  data.frame(Freq = tapply(Freqall$Freq, Freqall$Var1, sum)/sum(Freqall$Freq)) |> 
    rownames_to_column("Alleles") |> ggplot(aes(x = Alleles, y = Freq)) + 
    geom_bar(stat = 'identity', fill = '#b2bc63', color = 'black') +
    theme_bw() +
    scale_y_continuous(labels = scales::percent) + 
    labs(x = 'Alleles', y = 'Frequency')
)

# Recodification according to the dosage of the reference allele

geno.code = matrix(nrow = nrow(geno), ncol = ncol(geno),
                   dimnames = list(rownames(geno), colnames(geno)))
for (i in colnames(geno)) {
  x = geno[,i]
  geno.code[,i] = nchar(gsub(paste0("[^", vcf@fix[which(vcf@fix[,"ID"] == i), 'REF'], "]"), "", x))
}

# saveRDS(geno.code, file = "Data/mrk_filtered.RDS")

# PCA
pca = prcomp(geno.code) 
pc.df = rbind(
  as.data.frame(pca$rotation) |> 
    rownames_to_column("ID") |> 
    mutate(comp = "load"),
  as.data.frame(pca$x) |> 
    mutate(ID = paste0("ind", 1:nrow(pca$x)),
           comp = "score")
)
ggplot(dat = subset(pc.df, subset = comp == 'score'), aes(x = PC1, y = PC2)) + 
  geom_point(size = 3, alpha = .5, colour = "#b2bc63") +
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = 0, linetype = 'dashed') +
  scale_x_continuous(limits = range(pca$x[,1:2])) +
  scale_y_continuous(limits = range(pca$x[,1:2])) +
  coord_fixed() + 
  theme_minimal() + 
  theme(text= element_text(size = 18)) + 
  labs(x = paste0("PC1 (", round(pca$sdev[1]^2/sum(pca$sdev^2)*100,2), "%)"),
       y= paste0("PC2 (", round(pca$sdev[2]^2/sum(pca$sdev^2)*100,2), "%)"))

# RRBLUP ------------------------------------------------------------------

# Check compatibility pheno -> mrk
all(dat$gen %in% rownames(geno.code))
which(!dat$gen %in% rownames(geno.code))
dat = dat[which(dat$gen %in% rownames(geno.code)),]
all(dat$gen %in% rownames(geno.code))

# Check compatibility mrk -> pheno
all(rownames(geno.code) %in% dat$gen)
which(!rownames(geno.code) %in% dat$gen)
geno.code = geno.code[which(rownames(geno.code) %in% dat$gen),]
all(rownames(geno.code) %in% dat$gen)

# Check the order
all(rownames(geno.code) == dat$gen)

# Recode the mrk matrix
geno_rr = geno.code - 1  # Recoding the marker matrix

# Fit the model
mod = mixed.solve(y = dat$blue, Z = geno_rr) 
blup_rr = mod$u 
data.frame(blup_rr) |>  
  rownames_to_column(var = "Marker") |>  
  rename(BLUP = 'blup_rr') |>  
  ggplot(aes(x = Marker, y = BLUP))+
  geom_point(alpha = .6, colour = "#10342d")+
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank())

GEBV = geno_rr %*% blup_rr

as.data.frame(GEBV) |> rownames_to_column(var = "gen") |> 
  rename(GEBV = V1) |> 
  ggplot(aes(x = gen, y = GEBV))+
  geom_point(alpha = .75, colour = "#10342d") + 
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())+
  labs(x = "Hybrid", y = "GEBV")


# Parameter estimates
Va = mod$Vu * ncol(geno_rr) 
Ve = mod$Ve 
mu = mod$beta 
geno_her = Va/(Va + Ve)
cor(dat$blue, GEBV, use = "complete.obs")

# GBLUP ------------------------------------------------------------------
Gmat = A.mat(geno_rr)

ComplexHeatmap::Heatmap(
  Gmat,
  col = viridis::inferno(20),
  show_row_names = F,
  show_column_names = F,
  heatmap_legend_param = list(title = "IBS")
)

mod2 = mixed.solve(dat$blue, K = Gmat)
GEBV2 = mod2$u
cor(GEBV2, GEBV)

as.data.frame(GEBV2) %>% 
  rownames_to_column(var = "Hybrid") %>% 
  rename(GEBV = 'GEBV2') %>% 
  ggplot(aes(x = Hybrid, y = GEBV))+
  geom_point(alpha = .75, colour = "#10342d") + 
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())+
  labs(x = "Hybrid", y = "GEBV")

# BayesB ------------------------------------------------------------------
# Centre the mrk matrix so the model converge faster
X_centred = scale(geno.code, center = TRUE, scale = FALSE)

ETA = list(list(X = X_centred, model = "BayesB"))

BayesB = BGLR(y = dat$blue, ETA = ETA, nIter = 70000, burnIn = 20000, 
              thin = 10, verbose = FALSE, saveAt = "complete_")
save(BayesB, file = "bayesbmod.RData")


# Checking the convergence
library(coda)

# 1. Ler o arquivo de variância residual salvo pelo BGLR
varE = read.table("complete_varE.dat")[,1]

# 2. Converter para objeto MCMC descartando o burnIn
# Se você rodou nIter=50000, burnIn=10000 e thin=10:
# As primeiras 1000 linhas salvas correspondem ao burn-in (10000 / 10)
varE_chain = mcmc(varE[2001:length(varE)])

# 3. Inspeção Visual (Trace plot deve parecer um "ruído branco")
plot(varE_chain)

# 4. Checar o Tamanho Efetivo da Amostra (ESS)
# O valor de ESS DEVE ser maior que 200 (idealmente > 500)
effectiveSize(varE_chain)

# 5. Teste de Convergência de Geweke (p-valor deve ser > 0.05)
geweke.diag(varE_chain)

# Extracting parameters
BayesB$ETA[[1]]$b  # Shrinkage
BayesB$ETA[[1]]$d # Probability of inclusion: Frequency that the SNP had effect different than 0
BayesB$varE # Residual variance
var(BayesB$yHat)  # Genetic variance (markers)

GEBV_bb = data.frame(Hybrid = dat$gen, GEBV = BayesB$yHat)

ggplot(data = GEBV_bb, aes(x = Hybrid, y = GEBV))+
  geom_point(alpha = .75, colour = "#10342d") + 
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x=element_blank())+
  labs(x = "Hybrid", y = "GEBV")

cor(GEBV, GEBV_bb$GEBV)
cor(GEBV2, GEBV_bb$GEBV)


df_bayesB = data.frame(
  Marker = colnames(X_centred),
  Effect = BayesB$ETA[[1]]$b,  # Efeito da marca
  PIP = BayesB$ETA[[1]]$d,     # Probabilidade de inclusão
  SD = BayesB$ETA[[1]]$SD.b    # Desvio padrão do efeito da marca
)

# Plot de PIP (Identificação de QTLs)
ggplot(df_bayesB, aes(x = 1:nrow(df_bayesB), y = PIP)) +
  geom_point(aes(color = PIP > 0.5), alpha = 0.8) +
  scale_color_manual(values = c("grey", "#10342d")) +
  labs(x = "Marcadores", y = "Probabilidade de Inclusão Posterior (PIP)") +
  theme_minimal()

# Cross-validations ------------------------------------------------------------
nfolds = 5
nrept = 5 
seed = 8  

cvdata = list() 
for (j in 1:nrept) {
  set.seed(seed*j) 
  cvdata[[j]] = dat[,c('gen', 'blue')]
  cvdata[[j]]$set = NA
  for (i in cvdata[[j]]$gen) {
    cvdata[[j]][cvdata[[j]]$gen == i,'set'] = sample(
      1:nfolds, 
      size = 1
    )
  }
}

## RRBLUP ====  
cv_rr = lapply(cvdata, function(x){
  res.list = list()
  for (i in unique(x$set)) {
    
    yNA = x$blue
    yNA[x$set == i] = NA
    
    mod_cv = mixed.solve(y = yNA, Z = geno_rr)
    
    pred = c(mod_cv$beta) + (geno_rr %*% mod_cv$u)
    
    res.list[[i]] = cbind(
      x[x$set == i,],
      yhat = pred[x$set == i,]
    )
  }
  res.list
})

cv_rr = lapply(cv_rr, function(x) do.call(rbind, x))

# Capacidade preditiva (quanto maior, melhor)
pa_rr = lapply(cv_rr, function(x) cor(x$blue, x$yhat, use = 'complete.obs'))

# MSPE: Quanto menor (mais próximo de 0), melhor
mspe_rr = lapply(cv_rr, function(x) sqrt(mean((x$blue - x$yhat)^2, na.rm = T)))

## BayesB ====  
cv_bb = lapply(cvdata, function(x){
  res.list = list()
  for (i in unique(x$set)) {
    
    yNA = x$blue
    yNA[x$set == i] = NA
    
    mod_cv = BGLR(y = yNA, ETA = ETA, nIter = 8000, burnIn = 1000, 
                  thin = 10, verbose = FALSE)
    unlink("*.dat")
    
    res.list[[i]] = data.frame(gen = x$gen, yhat = mod_cv$yHat) |> 
      left_join(x, by = "gen") |> filter(set == i)
  }
  res.list
})

cv_bb = lapply(cv_bb, function(x) do.call(rbind, x))
save(cv_bb, file = "cv_bb.RDA")

# Capacidade preditiva (quanto maior, melhor)
pa_bb = lapply(cv_bb, function(x) cor(x$blue, x$yhat, use = 'complete.obs'))

# MSPE: Quanto menor (mais próximo de 0), melhor
mspe_bb = lapply(cv_bb, function(x) sqrt(mean((x$blue - x$yhat)^2, na.rm = T)))









