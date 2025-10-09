
# Convert lat/lon (degrees) to 3D unit vector
sph_to_cart <- function(latlon) {
  lat <- latlon[,1] * pi/180
  lon <- latlon[,2] * pi/180
  x <- sin(lat) * cos(lon)
  y <- sin(lat) * sin(lon)
  z <- cos(lat)
  cbind(x,y,z)
}

# Covariance matrix with spherical harmonics kernel
cov_matrix <- function(X1, X2, Lmax) {
  n <- nrow(X1)
  m = nrow(X2)
  X1 <- sph_to_cart(X1)
  X2 <- sph_to_cart(X2)        
  
  # cosine of all pairwise angles: n x m
  cosmat <- X1 %*% t(X2)
  
  coeffs <- (2*(1:Lmax)+1)/((1:Lmax)*(2:(Lmax+1)))^2 #coefficients of the polynomials
  
  # Initialize P_0 and P_1
  P_prev <- matrix(1, n, m)              # P0 = 1
  P_curr <- cosmat                       # P1 = cosθ
  
  # Start accumulation
  K <- coeffs[1] * P_curr                # add l=1 term
  
  if (Lmax > 1) {
    for (l in 2:Lmax) {
      P_next <- ((2*l-1)*cosmat*P_curr - (l-1)*P_prev)/l #recurrence formula for the legendre polynomials
      K <- K + coeffs[l] * P_next
      P_prev <- P_curr
      P_curr <- P_next
    }
  }
  return(K)
}


pred_krig = function(df_obs, df_pred, Lmax, noise = 0){
    m = nrow(df_pred)
    n = nrow(df_obs)
    cov_inv = solve(cov_matrix(df_obs[,1:2],df_obs[,1:2], Lmax = Lmax)+ noise^2*diag(nrow = n)) #compute inverse of covariance matrix
    vec_ones = rep(1, n) 
    trend = as.numeric(solve(t(vec_ones)%*%cov_inv%*%vec_ones)%*%t(vec_ones)%*%cov_inv%*%df_obs[,3]) #compute trend (ordinary kriging)
    return(trend + cov_matrix(df_pred, df_obs[,1:2], Lmax = Lmax)%*%cov_inv%*%(df_obs[,3]-rep(trend,n))) #compute prediction
}