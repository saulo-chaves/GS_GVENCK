rm(list=ls())

# Library -----------------------------------------------------------------
library(tidyverse)
library(BGLR)
library(AGHmatrix)
library(ComplexHeatmap)
library(plotly)
library(circlize)

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
Sys.setlocale("LC_TIME", "English") 
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
rownames(ZGZ) = aux1
dim(ZGZ)
ZEZ = tcrossprod(ZE)
dim(ZEZ)

ZOZ = ZE %*% Emat %*% t(ZE)
dim(ZEZ)
rownames(ZOZ) -> aux2

GEI = ZGZ * ZEZ
rownames(GEI) -> aux3
dim(GEI)

GEI_EC = ZGZ * ZOZ
rownames(GEI_EC) -> aux4
dim(GEI_EC)



