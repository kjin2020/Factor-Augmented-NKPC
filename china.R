source("NKPC_1.R")

#############
##Construct inflation
cpi_df   = read.csv("china_cpi.csv")
cpi_mat  = as.matrix(cpi_df[,c(2:ncol(cpi_df))])
cpi_mat_Q= matrix(NA,nrow(cpi_mat),ncol(cpi_mat)/3)
cum_multi = function(x){
  res = 1
  for (i in x){
    res = res * i/100
  }
  return(res*100)
}
for (i in 1:(ncol(cpi_mat)/3)){
  cpi_mat_Q[,i] = apply(cpi_mat[,c(1:3)+3*(i-1)],1,cum_multi)
}
cpi_q_n  = rep(NA,ncol(cpi_mat_Q))
ct = 1
for(year in 2023:2012){
  for (Quarter in 4:1){
    cpi_q_n[ct] = paste0(as.character(year),"Q",as.character(as.character(Quarter)))
    ct = ct + 1
  }
}
cpi_mat_Q1 = cpi_mat_Q-100
for (i in 1:nrow(cpi_mat_Q1)){
  cpi_mat_Q1[i,] = nkpc_quarter_adjust(cpi_mat_Q1[i,], substr(cpi_q_n, 6, 6))
}
cpi_df_Q = as.data.frame(cpi_mat_Q1)

rownames(cpi_df_Q) = cpi_df$Region
colnames(cpi_df_Q) = cpi_q_n

cpi_df_Q_t = t(cpi_df_Q)

cpi_df_Q_t    = cpi_df_Q_t[rev(rownames(cpi_df_Q_t)), ]

start_year    = as.numeric(substr(rownames(cpi_df_Q_t)[1], 1, 4))
start_quarter = as.numeric(substr(rownames(cpi_df_Q_t)[1], 6, 6))

cum_prod = function(x){
  res = rep(1,length(x))
  for(i in 2:length(x)){
    res[length(x)-i+1] = res[length(x)-i+2] * x[length(x)-i+1]
  }
  return(res)
}

cum_cpi = apply(cpi_mat_Q/100,1,cum_prod)

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
#######################
for(i in 1:ncol(rgdpt)){
  rgdpt[,i] = rgdpt[,i]/cum_cpi[,i]
}#####################
my_lag <- function(x, k){
  c(rep(NA, k), x[1:(length(x)-k)])
}
my_lead <-function(x, k){
  c(x[(1+k):length(x)],rep(NA,k))
}
region    = colnames(rgdpt)[1]
panel_rgdp= as.data.frame(cbind(region,rgdpt[,1]))
panel_rgdp[,2]= as.numeric(as.character(panel_rgdp[,2]))
panel_rgdp    = panel_rgdp[rev(rownames(panel_rgdp)), ]
start_year    = as.numeric(substr(rownames(panel_rgdp)[1], 1, 4))
start_quarter = as.numeric(substr(rownames(panel_rgdp)[1], 6, 6))

panel_rgdp[,2]= nkpc_cycle(panel_rgdp[,2], substr(rownames(panel_rgdp), 6, 6))
#max_lag <- 5
#for(j in 1:max_lag){
#  panel_rgdp <- cbind(panel_rgdp, my_lag(panel_rgdp[,2], k=j))
#}
#max_lead <- 1
#for(j in 1:max_lead){
#  panel_rgdp <- cbind(panel_rgdp, my_lead(panel_rgdp[,2], k=j))
#}

#colnames(panel_rgdp) = c("region","rgdp","lag1","lag2","lag3","lag4","lag5","lead1")
colnames(panel_rgdp) = c("region","rgdp")

#panel_rgdp = panel_rgdp[-c(1:5,nrow(panel_rgdp)),]

for(i in 2:ncol(rgdpt)){
  market = colnames(rgdpt)[i]
  tmp    = as.data.frame(cbind(market,rgdpt[,i]))
  tmp[,2]= as.numeric(as.character(tmp[,2]))
  tmp    = tmp[rev(rownames(tmp)), ]
  tmp[,2]= nkpc_cycle(tmp[,2], substr(rownames(tmp), 6, 6))
  
#  max_lag <- 5
#  for(j in 1:max_lag){
#    tmp <- cbind(tmp, my_lag(tmp[,2], k=j))
#  }
#  max_lead <- 1
#  for(j in 1:max_lead){
#    tmp <- cbind(tmp, my_lead(tmp[,2], k=j))
#  }
  colnames(tmp) = c("region","rgdp")
#  tmp = tmp[-c(1:5,nrow(tmp)),]
  
  panel_rgdp = rbind(panel_rgdp,tmp)
}
#############
##Construct real estate investment 
investment_cycle = function(x, quarters){
  # Keep the reported flows; a nonpositive flow leaves this province's cycle missing.
  if (any(!is.finite(x) | x <= 0)) return(rep(NA_real_, length(x)))
  nkpc_cycle(x, quarters)
}
real_inv     = read.csv("real_es.csv", check.names = FALSE)
rownames(real_inv) = rownames(ngdp)
quarter_ends = grep("年(3|6|9|12)月$", names(real_inv))
real_inv_mat_Q = as.matrix(real_inv[,quarter_ends,drop=FALSE])
real_inv_q_n = paste0(sub("年.*$", "", colnames(real_inv_mat_Q)), "Q",
  as.integer(sub(".*年([0-9]+)月$", "\\1", colnames(real_inv_mat_Q))) %/% 3L)
colnames(real_inv_mat_Q) = real_inv_q_n
# Quarter-end observations are year-to-date totals, just as for nominal GDP.
non_q1 = which(!endsWith(real_inv_q_n, "Q1"))
real_inv_mat_Q[,non_q1] = real_inv_mat_Q[,non_q1,drop=FALSE] -
  real_inv_mat_Q[,non_q1+1L,drop=FALSE]
real_inv_mat_Q = t(real_inv_mat_Q[,colnames(ngdp)[1:48],drop=FALSE])

colnames(real_inv_mat_Q) = colnames(rgdpt)
rownames(real_inv_mat_Q) = rownames(rgdpt)

for(i in 1:ncol(real_inv_mat_Q)){
  real_inv_mat_Q[,i] = real_inv_mat_Q[,i]/cum_cpi[,i]
}
###########
region    = colnames(real_inv_mat_Q)[1]
panel_real_estate_inv = as.data.frame(cbind(region,real_inv_mat_Q[,1]))
panel_real_estate_inv[,2]= as.numeric(as.character(panel_real_estate_inv[,2]))
panel_real_estate_inv    = panel_real_estate_inv[rev(rownames(panel_real_estate_inv)), ]

start_year    = as.numeric(substr(rownames(panel_real_estate_inv)[1], 1, 4))
start_quarter = as.numeric(substr(rownames(panel_real_estate_inv)[1], 6, 6))

panel_real_estate_inv[,2]= investment_cycle(panel_real_estate_inv[,2],
  substr(rownames(panel_real_estate_inv), 6, 6))

#max_lag <- 5
#for(j in 1:max_lag){
#  panel_real_estate_inv <- cbind(panel_real_estate_inv, my_lag(panel_real_estate_inv[,2], k=j))
#}
#max_lead <- 1
#for(j in 1:max_lead){
#  panel_real_estate_inv <- cbind(panel_real_estate_inv, my_lead(panel_real_estate_inv[,2], k=j))
#}

colnames(panel_real_estate_inv) = c("region","reinv")
#panel_real_estate_inv = panel_real_estate_inv[-c(1:5,nrow(panel_real_estate_inv)),]

for(i in 2:ncol(real_inv_mat_Q)){
  market = colnames(real_inv_mat_Q)[i]
  tmp    = as.data.frame(cbind(market,real_inv_mat_Q[,i]))
  tmp[,2]= as.numeric(as.character(tmp[,2]))
  tmp    = tmp[rev(rownames(tmp)), ]
  
  tmp[,2]= investment_cycle(tmp[,2], substr(rownames(tmp), 6, 6))
  
#  max_lag <- 5
#  for(j in 1:max_lag){
#    tmp <- cbind(tmp, my_lag(tmp[,2], k=j))
#  }
#  max_lead <- 1
#  for(j in 1:max_lead){
#    tmp <- cbind(tmp, my_lead(tmp[,2], k=j))
#  }
  colnames(tmp) = c("region","reinv")
#  tmp = tmp[-c(1:5,nrow(tmp)),]
  
  panel_real_estate_inv = rbind(panel_real_estate_inv,tmp)
}
panel_real_estate_inv

