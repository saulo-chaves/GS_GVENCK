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

ggplot(data = as.data.frame(table(dat$gen, dat$env)),
       aes(x = Var2, y = Var1, fill = as.factor(Freq))) + 
  geom_tile(aes(alpha = as.factor(Freq))) + 
  theme_minimal() + 
  theme(axis.text = element_blank(), legend.title = element_blank(), 
        legend.position = 'top') +
  scale_fill_manual(values = c('#b2bc63','#10342d'), labels = c("Abscent", "Present")) +
  scale_alpha_manual(values = c(.6, 1), labels = c("Abscent", "Present")) + 
  labs(x = "Environment (Year-Location combination)", y = "Hybrid")

# Filtering the dataset to analyse only the last year

dat = droplevels(dat[which(dat$year %in% c(2023)),])
ggplot(data = as.data.frame(table(dat$gen, dat$env)),
       aes(x = Var2, y = Var1, fill = as.factor(Freq))) + 
  geom_tile(aes(alpha = as.factor(Freq))) + 
  theme_minimal() + 
  theme(axis.text = element_blank(), legend.title = element_blank(), 
        legend.position = 'top') +
  scale_fill_manual(values = c('#b2bc63','#10342d'), labels = c("Abscent", "Present")) +
  scale_alpha_manual(values = c(.6, 1), labels = c("Abscent", "Present")) + 
  labs(x = "Environment (Year-Location combination)", y = "Hybrid")

# Genomic data
geno.code = readRDS(file = "Data/mrk_filtered.RDS")

# Filtering out ungenotyped and unphenotyped individuals
dat = dat[which(dat$gen %in% rownames(geno.code)),]
geno.code = geno.code[which(rownames(geno.code) %in% dat$gen),]

all(unique(dat$gen) %in% rownames(geno.code))

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

## Building the covariance matrices
dat$gen = factor(dat$gen, levels = rownames(Gmat), ordered = TRUE)

ZG = model.matrix(~-1 + gen, data = dat)
ZE = model.matrix(~-1 + env, data = dat)

rownames(ZG) = rownames(ZE) = dat$gen
colnames(ZG) = gsub("gen", "", colnames(ZG))
colnames(ZE) = gsub("env", "", colnames(ZE))

ZGZ = eigen(ZG %*% Gmat %*% t(ZG), )
ZEZ = eigen(tcrossprod(ZE))
GEI = eigen(ZGZ * ZEZ)

ETA = list(
  env = list(V = eigen(ZEZ)$vectors, d = eigen(ZEZ)$values, model = 'RKHS'),
  hyb = list(V = eigen(ZGZ)$vectors, d = eigen(ZGZ)$values, model = 'RKHS'),
  gei = list(V = eigen(GEI)$vectors, d = eigen(GEI)$values, model = 'RKHS')
)

