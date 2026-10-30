source("NKPC_1.R")
library(dplyr)
library(tidyr)
library(tibble)
library(lubridate)
#############
##Construct output gap
ngdp     = read.csv("china_gdp.csv")
rownames(ngdp) = ngdp[,1]
ngdp = ngdp[,-1]
gdp_q_n  = rep(NA,ncol(ngdp))
ct = 1
for(year in 2023:2006){
  for (Quarter in 4:1){
    gdp_q_n[ct] = paste0(as.character(year),"Q",as.character(as.character(Quarter)))
    ct = ct + 1
  }
}
colnames(ngdp) = gdp_q_n

# Year-to-date nominal GDP -> quarterly flows, before deflation and filtering.
# Columns descend Q4 to Q1; leave Q1 unchanged and difference within each year.
non_q1 = which(!endsWith(colnames(ngdp), "Q1"))
ngdp[,non_q1] = ngdp[,non_q1,drop=FALSE] - ngdp[,non_q1+1L,drop=FALSE]

rgdpt = t(ngdp[,1:48])

for(i in 1:ncol(rgdpt)){
  rgdpt[,i] = rgdpt[,i]/cum_cpi[,i]
}

rgdpt = as.data.frame(rgdpt)

rgdpt <- rgdpt %>%
  rownames_to_column(var = "Date")

rgdpt$Date <- yq(rgdpt$Date)

rgdpt <- rgdpt %>%
  arrange(Date)

for(i in 2:ncol(rgdpt)){
  rgdpt[,i] = nkpc_cycle(rgdpt[,i], quarter(rgdpt$Date))
}

panel_rgdp <- rgdpt %>%
  pivot_longer(cols = -Date, names_to = "region", values_to = "y_gap") %>%
  rename(Date = Date)


uniq = unique(panel_rgdp$region)
for(i in uniq){
  print(i)
}

