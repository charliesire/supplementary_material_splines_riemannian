
get_corners <- function(point, list_tri) {#get the corners of each triangle

  all_triangles <- list_tri$triangles[list_tri$dic_tri[[point]], , drop = FALSE]
  
  corners_list <- lapply(1:nrow(all_triangles), function(i) list_tri$points_3d[all_triangles[i, ], , drop = FALSE])
  
  return(corners_list)
}

############## for anisotropies (see supplementary material)

centroid_from_corners = function(corners){ #get the centroid of a triangle
  return(apply(corners, 2, mean))
}

rot_mat_func = function(theta){ #get a 2D rotation matrix
  return(cbind(c(cos(theta), sin(theta)), c(-sin(theta), cos(theta))))
}

anis_mat_func = function(theta, rho){ #Get the 2D or 3D anisotropy matrix with lengthscales and rotations
  l_rot = list()
  if(length(theta) == 3){
    for(k in 1:3){
      l_rot[[k]] = matrix(0, nrow = 3, ncol=3)
      l_rot[[k]][-k,-k] = rot_mat_func(theta[k])
      l_rot[[k]][k,k] = 1
    }
    mat_rot = l_rot[[1]]%*%l_rot[[2]]%*%l_rot[[3]]
  }
  else{mat_rot = rot_mat_func(theta)}
  if(length(rho)==2){rho = c(rho)}
  mat_diag = diag(1/rho^2)
  return(mat_rot%*%mat_diag%*%t(mat_rot))
}

jacob_sphere = function(x){ #Jacobian matrix of the coordinate change in the sphere case
  rho = sqrt(x[1]^2+x[2]^2)
  r = sqrt(x[1]^2+x[2]^2+x[3]^2)
  res = rbind(c(x[1]*x[3]/r^2/rho, x[2]*x[3]/r^2/rho,-rho/r^2), c(-x[2]/rho^2, x[1]/rho^2,0))
  return(res)
}

jacob_cylindre = function(x){#Jacobian matrix of the coordinate change in the cylinder case
  r2 = x[1]^2+x[2]^2
  res = rbind(c(-x[2]/r2, x[1]/r2,0), c(0,0,1))
  return(res)
}

anis_loc = function(x, anis_mat){ #anisotropy matrice t(J)G^x_T J
  J = jacob(x)
  return(t(J)%*%anis_mat%*%J)
}

########### For the computation of the mass and stiffness matrices

get_MT = function(corners){ #compute the matrix M_T associated to each triangle
  MT = as.matrix(t(rbind(corners[1, ] - corners[3, ],
                         corners[2, ] - corners[3, ])))
  return(MT)
}

get_GT = function(corners, anis_mat){ #Compute the matrix GT associated to each triangle, from the anisotropy matrix t(J)G^x_T J
  MT = get_MT(corners)
  if(ncol(anis_mat)==2){anis_mat = anis_loc(x = centroid_from_corners(corners), anis_mat = anis_mat)}
  GT = t(MT)%*%anis_mat%*%MT
  return(GT)
}


compute_M_F = function(list_tri, theta = rep(0,3), rho = rep(1,3), only_diag = TRUE){ #compute M and F matrices 
  mat_proj = (cbind(c(1,0), c(0,1), c(-1,-1)))
  Ns = nrow(list_tri$points_2d)
  
  if(length(theta)==3){ #if anisotropy in 3D
    anis_mat = anis_mat_func(theta, rho) #compute anisotropy matrix
    det_anis = det(anis_mat) #get determinant 
    inv_anis = solve(anis_mat) #inverse of anisotropy matrix
    fun_tri = function(tri_idx){ #for each triangle
      corners = list_tri$points_3d[list_tri$triangles[tri_idx,],] #get corners
      MT = get_MT(corners) #MT matrix
      PT = solve(t(MT)%*%MT) 
      area2 = 1/sqrt(det(PT)) #squared area of the triangle
      PT = PT%*%t(MT) 
    f_matrix = area2*sqrt(det_anis)*t(mat_proj)%*%PT%*%inv_anis%*%t(PT)%*%mat_proj/2 #contribution of this triangle to the stiffness matrix F
    return(c(area2/6*sqrt(det_anis), as.vector(f_matrix))) #return both contributions of this triangle to M and F
    }
  }
  else{ #if 2D anisotropies
    if(is.null(nrow(theta))){anis_mat = anis_mat_func(theta,rho)} #if theta is constant across all triangles
    else{anis_mat = NULL}
    fun_tri = function(tri_idx){ 
      if(is.null(anis_mat)){anis_mat = anis_mat_func(theta[tri_idx,], rho[tri_idx,])} #if theta and rho are not constant, compute anisotropy matrix for each triangle
      GT = get_GT(list_tri$points_3d[list_tri$triangles[tri_idx,],],anis_mat) #get GT
      det_GT = det(GT) #determinant of GT
      for_M = sqrt(det_GT) #
      f_matrix = sqrt(det_GT)/2*t(mat_proj)%*%solve(GT)%*%mat_proj #contribution of this triangle to the stiffness matrix F
      return(c(for_M/6,as.vector(f_matrix))) #return both contributions of this triangle to M and F
    }
  }


  res = t(sapply(1:nrow(list_tri$triangles),fun_tri)) 
  mat_M <- sapply(1:Ns, function(point){sum(res[,1][list_tri$dic_tri[[point]]])}) #sum all contributions for M

  vec_i = as.vector(t(list_tri$triangles[rep(1:nrow(list_tri$triangles), each = 3),]))
  vec_j = rep(as.vector(t(list_tri$triangles)),each = 3)
  vec_x = as.vector(t(res[,-1]))
  mat_F = sparseMatrix(i=vec_i, j = vec_j, x = vec_x, dims = c(Ns, Ns)) #sum all contributions for F, at the correct locations
  if(!only_diag){mat_M <- Matrix(diag(mat_M), sparse = TRUE)}
  return(list(mat_M = mat_M, mat_F = mat_F))
}

####### To compute the projection matrix


psi <- function(s, i, list_tri) { #compute psi_i(s)
  if (all(list_tri$points_2d[i, ] == s)) { #if s in the node i, then psi_i(s)=1
    return(1)
  } else {
    kk=1
    res = 0
    while(res == 0 & kk <= length(list_tri$dic_tri[[i]])){ #list_tri$dict_tri[[i]] provides the triangles indices on which psi_i has an effect
      res = psi_trunc(s, i, list_tri$dic_tri[[i]][kk], list_tri) #compute the restriction of psi_i on tri_idx, evaluated at s
      kk = kk+1
    }
    return(res)
  }
}


A_func = function(points_mes = NULL, list_tri = NULL, sparse = TRUE, all_in_tri = FALSE, idx_obs = NULL, parallel = FALSE){ #compute the matrix A_n of projection from the observation points to the mesh nodes
  Ns = nrow(list_tri$points_2d)
  if(all_in_tri){ #all_in_tri indicates whether all the observation points are on the triangulation nodes
    A_block = sparseMatrix( #if yes, then A juste contains 0 and 1 at the correct location
      i = 1:length(idx_obs),
      j = idx_obs,
      x = 1,
      dims = c(length(idx_obs), Ns)
    )
  }
  else{ #if all_in_tri is FALSE
    n_points <- nrow(points_mes)
    if(parallel){
        A_block <- foreach(idx = 1:n_points, .combine = 'rbind') %dopar% {
            s <- as.numeric(points_mes[idx, ])         # Coordinates of s
            row <- sapply(1:Ns, function(i) psi(s, i, list_tri))  # compute psi_i(s) for all s
            row   # return each row
        }
    }else{
        A_block <- foreach(idx = 1:n_points, .combine = 'rbind') %do% {
            s <- as.numeric(points_mes[idx, ])         # Coordinates of s
            row <- sapply(1:Ns, function(i) psi(s, i, list_tri))  # compute psi_i(s) for all s
            row   # return each row
        }
      }
  }
  if(!sparse){A_block = as.matrix(A_block)}
  return(A_block)
}





