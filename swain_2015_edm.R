library(ggplot2)
library(tidyverse)
library(rEDM)

## Data for cod, thorny skate, and white hake in the Gulf of St. Lawrence

cod <- read.csv("data/cod.csv")
skate <- read.csv("data/thorny_skate.csv")
hake <- read.csv("data/white_hake.csv")

## Data processing

## Cod
cod_n <- cod %>% group_by(stratum,year) %>% summarize(n())

cod_str <- cod %>%
  group_by(stratum,year) %>% 
  summarize(catch=mean(catch)) %>% 
  mutate(stcatch=(catch-mean(catch))/sd(catch))

segs <- cod_str %>% 
  ungroup()%>%
  mutate(ind = row_number()) %>% 
  group_by(stratum) %>%
  summarise(first=first(ind),last=last(ind)) %>%
  select(-stratum)

cod_simp <- simplex(cod_str$stcatch,lib=as.matrix(segs),E=1:20)

plot(cod_simp$E,cod_simp$mae,type="l")
plot(cod_simp$E,cod_simp$rho,type="l")
cod_tp <- simplex(cod_str$stcatch,tp=1:10,E=14)
plot(cod_tp$tp,cod_tp$rho,type="l")

## Hake
hake_n <- hake %>% group_by(stratum,year) %>% summarize(n())

hake_str <- hake %>%
  group_by(stratum,year) %>% 
  summarize(catch=mean(catch)) %>% 
  mutate(stcatch=(catch-mean(catch))/sd(catch))

hake_segs <- hake_str %>% 
  ungroup()%>%
  mutate(ind = row_number()) %>% 
  group_by(stratum) %>%
  summarise(first=first(ind),last=last(ind)) %>%
  select(-stratum)

hake_simp <- simplex(hake_str$stcatch,lib=as.matrix(hake_segs),E=1:20)

plot(hake_simp$E,hake_simp$mae,type="l")
plot(hake_simp$E,hake_simp$rho,type="l")
hake_tp <- simplex(hake_str$stcatch,tp=1:10,E=15)
plot(hake_tp$tp,hake_tp$rho,type="l")

## Thorny skate
skate_n <- skate %>% group_by(stratum,year) %>% summarize(n())

skate_str <- skate %>%
  group_by(stratum,year) %>% 
  summarize(catch=mean(catch)) %>% 
  mutate(stcatch=(catch-mean(catch))/sd(catch))

skate_segs <- skate_str %>% 
  ungroup()%>%
  mutate(ind = row_number()) %>% 
  group_by(stratum) %>%
  summarise(first=first(ind),last=last(ind)) %>%
  select(-stratum)

skate_simp <- simplex(skate_str$stcatch,lib=as.matrix(skate_segs),E=1:25)

plot(skate_simp$E,skate_simp$mae,type="l")
plot(skate_simp$E,skate_simp$rho,type="l")
skate_tp <- simplex(skate_str$stcatch,tp=1:10,E=19)
plot(skate_tp$tp,skate_tp$rho,type="l")

#### Univariate simplex ####
## With Nndx
cod_Nndx <- cod %>% select(year,Nndx) %>% distinct(Nndx,.keep_all=T)%>% 
  mutate(stcod=(Nndx-mean(Nndx))/sd(Nndx))
cod_simp2 <- simplex(cod_Nndx$stcod)
plot(cod_simp2$E,cod_simp2$mae,type="l")
plot(cod_simp2$E,cod_simp2$rho,type="l")

hake_Nndx <- hake %>% select(year,Nndx) %>% distinct(Nndx,.keep_all=T)%>% 
  mutate(sthake=(Nndx-mean(Nndx))/sd(Nndx))
hake_simp2 <- simplex(hake_Nndx$sthake)
plot(hake_simp2$E,hake_simp2$mae,type="l")
plot(hake_simp2$E,hake_simp2$rho,type="l")

skate_Nndx <- skate %>% select(year,Nndx) %>% distinct(Nndx,.keep_all=T)%>% 
  mutate(stskate=(Nndx-mean(Nndx))/sd(Nndx))
skate_simp2 <- simplex(skate_Nndx$stskate)
plot(skate_simp2$E,skate_simp2$mae,type="l")
plot(skate_simp2$E,skate_simp2$rho,type="l")


## Seals
seal <- cod %>% 
  select(year,kseal) %>% 
  distinct(kseal,.keep_all=T) %>% 
  mutate(stseal=(kseal-mean(kseal))/sd(kseal))
seal_simp2 <- simplex(seal$stseal)
plot(seal_simp2$E,seal_simp2$mae,type="l")
plot(seal_simp2$E,seal_simp2$rho,type="l")

#### Univariate S maps ####
## Cod S-map (best E=2)
cod_smap <- s_map(cod_Nndx$stcod,E=2)
plot(cod_smap$theta,cod_smap$rho,type="l")
# best theta=2

## hake S-map (best E=4)
hake_smap <- s_map(hake_Nndx$sthake,E=4)
plot(hake_smap$theta,hake_smap$rho,type="l")
# best theta=1

## skate S-map
skate_smap <- s_map(skate_Nndx$stskate,E=4)
plot(skate_smap$theta,skate_smap$rho,type="l")
# best theta=0.5

## seal S-map (E=4??) maybe not appropriate
seal_smap <- s_map(seal$stseal,E=4)
plot(seal_smap$theta,seal_smap$rho,type="l")

#### Multivariate analyses ####
## Block
dat <- bind_cols(cod_Nndx,skate_Nndx,hake_Nndx,seal)

# For cod
cod_multi <- block_lnlp(dat,target_column = "stcod",method="s-map",num_neighbors=0,
                        theta=2,columns = c(3,6,9,12),save_smap_coefficients=T)
cod_multi[[1]]$stats
cod_ints <- cod_multi[[1]]$smap_coefficients %>% as.data.frame()
colnames(cod_ints) <- c("cod","skate","hake","seal","constant")

# Boxplot of interactions
# long form
cod_ints_long <- cod_ints %>% 
  gather(key="effect")
cod_ints_plot <- ggplot(cod_ints_long,aes(x=effect,y=value,fill=effect)) +
  geom_boxplot()+
  geom_hline(yintercept=0,linetype=2)+
  ggtitle("Interactive Effects on Cod")+
  theme_minimal()
cod_ints_plot

# For hake
hake_multi <- block_lnlp(dat,target_column = "sthake",method="s-map",num_neighbors=0,
                        theta=1,columns = c(3,6,9,12),save_smap_coefficients=T)
hake_multi[[1]]$stats
hake_ints <- hake_multi[[1]]$smap_coefficients %>% as.data.frame()
colnames(hake_ints) <- c("cod","skate","hake","seal","constant")

# Boxplot of interactions
# long form
hake_ints_long <- hake_ints %>% 
  gather(key="effect")
hake_ints_plot <- ggplot(hake_ints_long,aes(x=effect,y=value,fill=effect)) +
  geom_boxplot()+
  geom_hline(yintercept=0,linetype=2)+
  ggtitle("Interactive Effects on Hake")+
  theme_minimal()
hake_ints_plot

# For skate
skate_multi <- block_lnlp(dat,target_column = "stskate",method="s-map",num_neighbors=0,
                        theta=0.5,columns = c(3,6,9,12),save_smap_coefficients=T)
skate_multi[[1]]$stats
skate_ints <- skate_multi[[1]]$smap_coefficients %>% as.data.frame()
colnames(skate_ints) <- c("cod","skate","hake","seal","constant")

# Boxplot of interactions
# long form
skate_ints_long <- skate_ints %>% 
  gather(key="effect")
skate_ints_plot <- ggplot(skate_ints_long,aes(x=effect,y=value,fill=effect)) +
  geom_boxplot()+
  geom_hline(yintercept=0,linetype=2)+
  ggtitle("Interactive Effects on Skate")+
  theme_minimal()
skate_ints_plot

