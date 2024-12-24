iT = sum(panel_rgdp$region=="Beijing")
N  = length(unique(panel_rgdp$region))

eigen_ratio_test = function(eigenVal, n, k) {
  temp = min(c(n, k))
  len  = ifelse(temp<4,temp-1,bitwShiftR(temp,1))
  print(len)
  if (len == 0) {
    rst = as.vector(eigenVal[k - 1])
    return(rst)
  }
  ratio = rep(NA,len)
  comp = eigenVal[k - 1] / eigenVal[k - 2]
  ratio[1] = comp
  for (i in 2:len) {
    ratio[i] = eigenVal[k - 1 - i] / eigenVal[k - 2 - i]
  }
  return (ratio)
}

uobfMg = function(ztz,N,iT){
  egv      = eigen(ztz)
  eigenVal = egv$values
  mg       = which.max(eigen_ratio_test(eigenVal,iT,iT))
  G_0      = t(sqrt(iT) * t(egv$vectors[,c(1:mg)]))
  Mg       = diag(iT) - G_0 %*%solve(t(G_0)%*%G_0)%*%t(G_0)
  return(list(G_0=G_0,Mg=Mg,mg=mg))
}

create_ztz = function(panel_data,iT,name){
  ztz = matrix(0,iT,iT)
  uniq= unique(panel_data[,name])
  for(i in 1:nrow(uniq)){
    z = panel_data[panel_data[,name] == as.character(uniq[i,]),]
    z = as.matrix(z[, sapply(z, is.numeric)])
    ztz = ztz + crossprod(t(z))
  }
  ztz = ztz/(iT*N)
  return(ztz)
}

#for(i in unique(panel_rgdp$region)){
#  z   = cbind(panel_rgdp[panel_rgdp$region==i,2],panel_real_estate_inv[panel_real_estate_inv$region==i,2])
#  z   = panel_rgdp[panel_rgdp$region==i,2]
#  ztz = ztz + crossprod(t(z))
#}
#ztz = ztz/(iT*N)

Mgmg= uobfMg(create_ztz(panel_rgdp,iT,"region"),N,iT)


first_region = unique(panel_rgdp$region)[1]

#z   = cbind(panel_rgdp[panel_rgdp$region==first_region,2],panel_real_estate_inv[panel_real_estate_inv$region==first_region,2])
z = panel_rgdp[panel_rgdp$region==first_region,3]
z = as.matrix(z[, sapply(z, is.numeric)])
BigZ= Mgmg$Mg%*%z
max_lag <- 5
for(j in 1:max_lag){
#  BigZ <- cbind(BigZ,rbind(matrix(NA,nrow=j,ncol=2),BigZ[1:(nrow(BigZ)-j),c(1,2)]))
  BigZ <- cbind(BigZ,rbind(matrix(NA,nrow=j,ncol=1),as.matrix(BigZ[1:(nrow(BigZ)-j),1])))
}
max_lead <- 1
for(j in 1:max_lead){
#  BigZ <- cbind(BigZ,rbind(BigZ[(j+1):nrow(BigZ),c(1,2)],matrix(NA,nrow=j,ncol=2)))
  BigZ <- cbind(BigZ,rbind(as.matrix(BigZ[(j+1):nrow(BigZ),1]),matrix(NA,nrow=j,ncol=1)))
}

BigZ = BigZ[-c(1:5,nrow(BigZ)),]

region = first_region

BigZ = cbind(region,as.data.frame(BigZ))

newiT = nrow(BigZ)

for (regi in unique(panel_rgdp$region)[2:N]){
#  z   = cbind(panel_rgdp[panel_rgdp$region==regi,2],panel_real_estate_inv[panel_real_estate_inv$region==regi,2])
  z   = panel_rgdp[panel_rgdp$region==regi,3]
  z = as.matrix(z[, sapply(z, is.numeric)])
  tmpBigZ= Mgmg$Mg%*%z
  max_lag <- 5
  for(j in 1:max_lag){
#    tmpBigZ <- cbind(tmpBigZ,rbind(matrix(NA,nrow=j,ncol=2),tmpBigZ[1:(nrow(tmpBigZ)-j),c(1,2)]))
    tmpBigZ <- cbind(tmpBigZ,rbind(matrix(NA,nrow=j,ncol=1),as.matrix(tmpBigZ[1:(nrow(tmpBigZ)-j),1])))
  }
  max_lead <- 1
  for(j in 1:max_lead){
#    tmpBigZ <- cbind(tmpBigZ,rbind(tmpBigZ[(j+1):nrow(tmpBigZ),c(1,2)],matrix(NA,nrow=j,ncol=2)))
    tmpBigZ <- cbind(tmpBigZ,rbind(as.matrix(tmpBigZ[(j+1):nrow(tmpBigZ),1]),matrix(NA,nrow=j,ncol=1)))
  }
  tmpBigZ = tmpBigZ[-c(1:5,nrow(tmpBigZ)),]
  region = regi
  tmpBigZ = cbind(region,as.data.frame(tmpBigZ))
  BigZ = rbind(BigZ,tmpBigZ)
}

cpi_df_Q_t = t(cpi_df_Q)

cpi_df_Q_t    = cpi_df_Q_t[rev(rownames(cpi_df_Q_t)), ]

regions = colnames(cpi_df_Q_t)

cpi_matrix = as.data.frame(cpi_df_Q_t[,1])

max_lag <- 1
for(j in 1:max_lag){
  cpi_matrix <- cbind(cpi_matrix,rbind(matrix(NA,nrow=j,ncol=1),as.matrix(cpi_matrix[1:(nrow(cpi_matrix)-j),1])))
}
max_lead <- 1
for(j in 1:max_lead){
  cpi_matrix <- cbind(cpi_matrix,rbind(as.matrix(cpi_matrix[(j+1):nrow(cpi_matrix),1]),matrix(NA,nrow=j,ncol=1)))
}

cpi_matrix = cpi_matrix[-c(1:5,nrow(cpi_matrix)),]

cpi_matrix = cbind(regions[1],cpi_matrix)

colnames(cpi_matrix) = c("region","pi","lag","lead")

for (regi in regions[2:N]){
  tmpx      = as.data.frame(cpi_df_Q_t[,regi])
  max_lag <- 1
  for(j in 1:max_lag){
    tmpx <- cbind(tmpx,rbind(matrix(NA,nrow=j,ncol=1),as.matrix(tmpx[1:(nrow(tmpx)-j),1])))
  }
  max_lead <- 1
  for(j in 1:max_lead){
    tmpx <- cbind(tmpx,rbind(as.matrix(tmpx[(j+1):nrow(tmpx),1]),matrix(NA,nrow=j,ncol=1)))
  }
  
  tmpx = tmpx[-c(1:5,nrow(tmpx)),]
  
  tmpx = cbind(regi,tmpx)
  
  colnames(tmpx) = c("region","pi","lag","lead")
  
  cpi_matrix = rbind(cpi_matrix,tmpx)
}

X = cbind(cpi_matrix,BigZ[,2])

IV_1st = function(BigZ,X,newiT,N){
  pi   = X[,c(1,2)]
  BigX = X[,c(1,3,4,5)]
  ZZ   = matrix(0,ncol(BigZ)-1,ncol(BigZ)-1)
  ZX   = matrix(0,ncol=ncol(BigX)-1,nrow=ncol(BigZ)-1)
  Zg   = matrix(0,nrow=ncol(BigZ)-1,ncol=1)
  for (regi in unique(BigZ$region)){
    pii  = as.matrix(pi[pi$region==regi,2])
    Zi   = as.matrix(BigZ[BigZ$region==regi,c(2:ncol(BigZ))])
    Xi   = as.matrix(BigX[BigX$region==regi,c(2,3,4)])
    ZZ   = ZZ + t(Zi)%*%Zi
    ZX   = ZX + t(Zi)%*%Xi
    Zg   = Zg + t(Zi)%*%pii
  }
  A = ZX  / (N*iT)
  B = ZZ  / (N*iT)
  g = Zg /  (N*iT)
  theta1 = solve(t(A)%*%solve(B)%*%A) %*% t(A) %*% solve(B) %*% g
  return(theta1)
}

IV_1 = IV_1st(BigZ,X,newiT,N)

u = X[,2] - as.matrix(X[,c(3,4,5)])%*%IV_1

u = cbind(X[,1],as.data.frame(u))

uu= matrix(0,newiT,newiT)
for (i in regions){
  ui = as.matrix(u[u$`X[, 1]`==i,2])
  uu = uu + ui%*%t(ui)
}

uu = uu/(N*newiT)

egv= eigen(uu)
eigenVal = egv$values
mf       = which.max(eigen_ratio_test(eigenVal,newiT,newiT))
F_0      = sqrt(iT) * egv$vectors[,c(1:mf)]

Mf       = diag(newiT)  - F_0 %*% solve(t(F_0)%*%F_0)%*%t(F_0)

IV_2st = function(BigZ,X,newiT,N,Mf){
  pi   = X[,c(1,2)]
  BigX = X[,c(1,3,4,5)]
  ZZ   = matrix(0,ncol(BigZ)-1,ncol(BigZ)-1)
  ZX   = matrix(0,ncol=ncol(BigX)-1,nrow=ncol(BigZ)-1)
  Zg   = matrix(0,nrow=ncol(BigZ)-1,ncol=1)
  for (regi in unique(BigZ$region)){
    pii  = as.matrix(pi[pi$region==regi,2])
    Zi   = as.matrix(BigZ[BigZ$region==regi,c(2:ncol(BigZ))])
    Xi   = as.matrix(BigX[BigX$region==regi,c(2,3,4)])
    ZZ   = ZZ + t(Zi)%*%Mf%*%Zi
    ZX   = ZX + t(Zi)%*%Mf%*%Xi
    Zg   = Zg + t(Zi)%*%Mf%*%pii
  }
  A = ZX  / (N*iT)
  B = ZZ  / (N*iT)
  g = Zg /  (N*iT)
  theta2 = solve(t(A)%*%solve(B)%*%A) %*% t(A) %*% solve(B) %*% g
  return(list(theta2=theta2,A2=A,B2=B))
}

IV_2_obj = IV_2st(BigZ,X,newiT,N,Mf)
IV_2 = IV_2_obj$theta2

rownames(IV_2) = c("pi_lag","pi_lead","y_gap")

sighat = function(A,B,newiT,N,u,BigZ,Mf){
  Delta = matrix(0,ncol(BigZ)-1,ncol(BigZ)-1)
  for (i in regions){
    ui = as.matrix(u[u$`X[, 1]`==i,2])
    sigma2i = t(ui)%*%Mf%*%ui / newiT
    Zi   = as.matrix(BigZ[BigZ$region==i,c(2:ncol(BigZ))])
    Delta  = Delta + sigma2i[1,1]*t(Zi)%*%Mf%*%Zi
  }
  Delta = Delta/(N*newiT)
  print(Delta)
  Sigma = solve(t(A)%*%solve(B)%*%A)%*%t(A)%*%solve(B)%*%Delta%*%solve(B)%*%A%*%solve(t(A)%*%solve(B)%*%A)
  return(Sigma)
}

Sigma1 = sighat(IV_2_obj$A2,IV_2_obj$B2,newiT,N,u,BigZ,Mf)

IV_2/sqrt(diag(Sigma1)) 

Jhat = function(u,Mf,Z,Delta1){
  return(iT^(-1)* (t(u)%*%Mf%*%Z)%*%Delta1%*%(t(Z)%*%Mf%*%u))
}





