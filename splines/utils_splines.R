power_iteration <- function(A = NULL, chol_A = NULL, v0 = NULL, 
                            num_iter = 1000, tol = 0.001, seed = 10) { #Power iteration to estimate the dominant eigenvalue of a matrix A
  # or the minimum eigenvalue using a Cholesky factorization
  # Power iteration to estimate the dominant eigenvalue of a matrix A
  # or an equivalent formulation using a Cholesky factorization of A+I

  if(!is.null(A)) {
    n <- nrow(A)
  } else {
    n <- nrow(chol_A)
  }
  
  # Initialize a random vector of length n
  set.seed(seed)
  b <- runif(n)
  # Normalize the initial vector
  b <- b / sqrt(sum(b^2))
  
  # Iterative power method
  for (i in 1:num_iter) {
    if(!is.null(A)) {
      # Standard power iteration: multiply by A
      b_new <- A %*% b
    } else {
      # Alternative formulation:
      # Solve system using Cholesky factorization and apply orthogonalization
      b_new <- solve(a = chol_A, b = b) - v0 %*% t(v0) %*% b
    }
    
    # Normalize the new vector
    b_new <- b_new / sqrt(sum(b_new^2))
    
    # Check convergence: if the vector has stabilized (within tol)
    if (sqrt(sum((b_new - b)^2)) < tol) {
      break
    }
    # Update the vector for the next iteration
    b <- b_new
  }
  
  # Estimate the corresponding eigenvalue
  if(!is.null(A)) {
    # Rayleigh quotient for matrix A
    lambda <- as.numeric(t(b) %*% A %*% b)
  } else {
    # Equivalent Rayleigh quotient formulation in the minimum eigenvalue case
    lambda <- as.numeric(t(b) %*% (solve(a = chol_A, b = b) - v0 %*% t(v0) %*% b))
    lambda <- 1 / lambda - 1
  }
  
  return(lambda)
}
                              
compute_chol = function(Q){ ## function to compute cholesky of Q + vv^T, when Q is not invertible
  if(is.list(Q)){
    v1 = Matrix(c(rep(0, nrow(Q$Q_1)-1),1), sparse = TRUE) ## To compute cholesky of Q, we introduce a vector v1 and compute cholesky of Q+v1v1^T
    chol1 = Cholesky(Q$Q_1+v1%*%t(v1), LDL = FALSE, perm = FALSE) #compute cholesky of Q+v1v1^T
    chol2 = updown(update = TRUE, C = Matrix(Q$v, sparse = TRUE), L = chol1) #compute cholesky of Q+v1v1^T + vv^T with rank-one update
    chol_Q = updown(update = FALSE, C = v1, L = chol2) #compute cholesky of Q+ vv^T with rank-one downdate
    return(chol_Q)
  }else{return(t(chol(Q)))}
}


augmented_matrices = function(list_tri, theta = rep(0,3), rho = rep(1, 3), var_a0 = NULL, sparse = TRUE){ #Function to compute the matrices sqrt(M)S^2 sqrt(M)^T and the term 1/alpha(Mphi0)(Mphi0)^T
  M_F = compute_M_F(list_tri = list_tri, theta = theta, rho = rho, only_diag = TRUE) #Compute the matrices M and F 
  mat_M = M_F$mat_M
  mat_F = M_F$mat_F
  sqrt_M = sqrt(mat_M) 
  inv_sqrt_M = Diagonal(n=length(mat_M),x = 1/sqrt_M) 
  phi_0_vec = -Matrix(rep(1,length(mat_M))/sqrt(sum(abs(mat_M))), sparse = TRUE) #get the vector phi_0
  N = mat_M*phi_0_vec #get Mphi0
  S = inv_sqrt_M%*%mat_F%*%inv_sqrt_M #compute S
  S2 = S%*%S
  if(is.null(var_a0)){ #if alpha is not specified
      chol_S2 = Cholesky(S2+Diagonal(n = nrow(S2)), LDL = FALSE, perm = FALSE) #cholesky of S2 + Im
      v0 = (sqrt_M*phi_0_vec) #first eigenvector of S
      lmin = power_iteration(chol_A = chol_S2, v0 = v0, seed = 2) #get lowest eigen value
      var_a0 = 0.01/lmin 
}
  sqrt_M =  Diagonal(n=length(mat_M),x = sqrt_M)
  part_1 =  sqrt_M%*%S2%*%sqrt_M 
  Q_tilde = list(Q_1 = part_1, v = 1/sqrt(var_a0)*N) #return a list with sqrt(M)S^2 sqrt(M)^T and the term vector 1/sqrt(alpha)(Mphi0)
  return(list(Q_tilde = Q_tilde, phi_0_vec = phi_0_vec, var_a0 = var_a0)) 
}

get_schur_comp = function(Q_tilde, idx_obs, y, ybis, return_chol = FALSE){ #this function computes (A_n Q_alpha^{-1] A_n^T)^-1 y and  (A_n Q_alpha^{-1] A_n^T)^-1 ybis using Schur Complement
  idx_bar = setdiff(1:ncol(Q_tilde$Q_1), idx_obs) #The indices Ibar 
  chol_Qbar = compute_chol(list(Q_1 = Q_tilde$Q_1[idx_bar, idx_bar], v = Q_tilde$v[idx_bar])) #Compute cholesky of Q_[Ibar,I]
  omega_y = solve(a = chol_Qbar, b =  Q_tilde$Q_1[idx_bar, idx_obs]%*%y+Q_tilde$v[idx_bar]%*%t(Q_tilde$v[idx_obs]%*%y)) #Compute Q_[Ibar,I]^-1 y
  omega_anphi0 = solve(a = chol_Qbar, b =  Q_tilde$Q_1[idx_bar, idx_obs]%*%ybis+Q_tilde$v[idx_bar]%*%t(Q_tilde$v[idx_obs]%*%ybis)) #Compute Q_[Ibar,I]^-1 ybis
  schur_comp_y = Q_tilde$Q_1[idx_obs, idx_obs]%*%y - Q_tilde$Q_1[idx_obs, idx_bar]%*%omega_y+ Q_tilde$v[idx_obs]%*%(t(Q_tilde$v[idx_obs])%*%y)-Q_tilde$v[idx_obs]%*%(t(Q_tilde$v[idx_bar])%*%omega_y) # Compute (A_n Q_alpha^{-1] A_n^T)^-1 y
  schur_comp_anphi0 = Q_tilde$Q_1[idx_obs, idx_obs]%*%ybis - Q_tilde$Q_1[idx_obs, idx_bar]%*%omega_anphi0+ Q_tilde$v[idx_obs]%*%(t(Q_tilde$v[idx_obs])%*%ybis)-Q_tilde$v[idx_obs]%*%(t(Q_tilde$v[idx_bar])%*%omega_anphi0) # Compute (A_n Q_alpha^{-1] A_n^T)^-1 ybis
  res = list(schur_comp_y = schur_comp_y, schur_comp_anphi0 = schur_comp_anphi0)
  if(return_chol){res$chol_Qbar = chol_Qbar}
  return(res)
}
                     
pred_func = function(list_tri, y,AA, idx_obs, var_a0 = NULL, theta = rep(0,3), rho = rep(1,3), noise = 0, list_augmented = NULL, chol_Qbar = NULL, omega_anphi0 = NULL){ #Compute spline prediction at the triangulation points
  if(is.null(list_augmented)){list_augmented = augmented_matrices(list_tri, theta = theta, rho = rho, var_a0 = var_a0)} #Compute matrix Q
  phi_0_vec = list_augmented$phi_0_vec #Get phi0
  var_a0 = list_augmented$var_a0
  Q_tilde = list_augmented$Q_tilde #Get sqrt(M)S^2 sqrt(M)^T and the term 1/alpha(Mphi0)(Mphi0)^T
  ybis = AA%*%phi_0_vec #get Aphi0

  if(noise == 0){ #First scenario
    idx_bar = setdiff(1:ncol(Q_tilde$Q_1), idx_obs) #The indices Ibar 

    if(is.null(chol_Qbar)){
        chol_Qbar = compute_chol(list(Q_1 = Q_tilde$Q_1[idx_bar, idx_bar], v = Q_tilde$v[idx_bar]))
        omega_anphi0 = -solve(a = chol_Qbar, b =  Q_tilde$Q_1[idx_bar, idx_obs]%*%ybis+Q_tilde$v[idx_bar]%*%t(Q_tilde$v[idx_obs]%*%ybis)) #Compute Q_[Ibar,I]^-1 Anphi0 
    }
    omega_y = -solve(a = chol_Qbar, b =  Q_tilde$Q_1[idx_bar, idx_obs]%*%y+Q_tilde$v[idx_bar]%*%t(Q_tilde$v[idx_obs]%*%y)) #Compute Q_[Ibar,I]^-1 y
    u_alp_y = seq_len(nrow(list_tri$points_2d)) #Compute u_alpha(y)
    u_alp_y[idx_bar] = omega_y
    u_alp_y[idx_obs] = y
    u_alp_anphi0 = seq_len(nrow(list_tri$points_2d)) #Compute u_alpha(Anphi0)
    u_alp_anphi0[idx_bar] = omega_anphi0
    u_alp_anphi0[idx_obs] = ybis

      
  }
  else{ #Second scenario
    chol_QA = compute_chol(Q = list(Q_1 = noise^2*Q_tilde$Q_1 + t(AA)%*%AA, v = Q_tilde$v*noise)) #Compute Cholesky decomposition of (tau^2 Q_alph + A^TA)
    u_alp_y = as.matrix(solve(a = chol_QA, b = t(AA)%*%y))#Compute u_alpha(y)
    u_alp_anphi0 = as.matrix(solve(a = chol_QA, b = t(AA)%*%ybis)) #Compute u_alpha(Aphi0)
  }
  a_alp_y = as.numeric(t(Q_tilde$v*sqrt(var_a0))%*%u_alp_y)  
  a_alp_anphi0 = as.numeric(t(Q_tilde$v*sqrt(var_a0))%*%u_alp_anphi0)
  h_u_anphi0 = (u_alp_anphi0 - a_alp_anphi0*phi_0_vec)/(1-a_alp_anphi0)
  res = as.matrix(u_alp_y + a_alp_y*(h_u_anphi0-phi_0_vec+(phi_0_vec-h_u_anphi0)/a_alp_anphi0)) #Compute the spline prediction
  return(res)
}

                              
ll_func = function(param, list_tri, AA, y, var_a0 = NULL, noise = 0, idx_obs = NULL, constant_anisotropies = TRUE){ #Compute the log-likelihood for given paramaters
  if(constant_anisotropies){
      if(length(param)==3){ #Local charts in 2D (theta, phi) or (theta,z) for instance
        theta = param[1] 
        rho = param[2:3]
      }
      else{ #3D charts x,y,z
        theta = param[1:3]
        rho = param[4:6]
      }
  }
  else{
    theta = matrix(param[,1])
    rho = param[,2:3]  
  }
  list_augmented = augmented_matrices(list_tri, theta = theta, rho = rho, var_a0 = var_a0) #Compute matrix Q
  phi_0_vec = list_augmented$phi_0_vec #Get phi0
  var_a0 = list_augmented$var_a0
  Q_tilde = list_augmented$Q_tilde #Get sqrt(M)S^2 sqrt(M)^T and the term 1/alpha(Mphi0)(Mphi0)^T
  ybis = AA%*%phi_0_vec #get Aphi0
  chol_Q = compute_chol(Q_tilde) #Compute Cholesky of sqrt(M)S^2 sqrt(M)^T + 1/alpha(Mphi0)(Mphi0)^T
  if(noise == 0){ #first scenario
    list_schur = get_schur_comp(Q_tilde = Q_tilde, idx_obs = idx_obs, y = y, ybis = ybis, return_chol = TRUE) #Compute  (A_n Q_alpha^{-1] A_n^T)^-1 y and  (A_n Q_alpha^{-1] A_n^T)^-1  Aphi0
    a_alp_y = as.numeric(var_a0*t(phi_0_vec)%*%t(AA)%*%list_schur$schur_comp_y)
    a_alp_anphi0 = as.numeric(var_a0*t(phi_0_vec)%*%t(AA)%*%list_schur$schur_comp_anphi0)
    part_quad = as.numeric(t(y)%*%list_schur$schur_comp_y-1/var_a0*a_alp_y^2/a_alp_anphi0) #Compute l1(param)
    part_det = sum(log(diag(list_schur$chol_Qbar)))-sum(log(diag(chol_Q)))+ log(1-a_alp_anphi0) #Compute l2(param)
  }
  
  else{
    chol_QA = compute_chol(Q = list(Q_1 = noise^2*Q_tilde$Q_1 + t(AA)%*%AA, v = Q_tilde$v*noise)) #Compute Cholesky decomposition of (tau^2 Q_alph + A^TA)
    u_alp_y = as.matrix(solve(a = chol_QA, b = t(AA)%*%y)) #Compute u_alpha(y)
    u_alp_anphi0 = as.matrix(solve(a = chol_QA, b = t(AA)%*%ybis)) #Compute u_alpha(Aphi0)
    a_alp_y = as.numeric(t(Q_tilde$v*sqrt(var_a0))%*%u_alp_y) 
    a_alp_anphi0 = as.numeric(t(Q_tilde$v*sqrt(var_a0))%*%u_alp_anphi0)
    part_quad = as.numeric(1/noise^2*(t(y)%*%y-t(y)%*%AA%*%u_alp_y)-1/var_a0*a_alp_y^2/a_alp_anphi0) #Compute l1(param)
    part_det = (length(y)-length(phi_0_vec))*log(noise^2)-sum(log(diag(chol_Q)))+sum(log(diag(chol_QA)))+log(1-a_alp_anphi0) #Compute l2(param)
  }

  return(-(part_quad+part_det))
}


