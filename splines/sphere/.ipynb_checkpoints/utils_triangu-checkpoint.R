build_triangu = function(radius, n_theta = NULL, n_phi = n_theta, theta = NULL, phi = NULL, nodes = NULL){ #build triangulation for sphere
  # radius is the radius of the sphere
  # n_theta is the number of points for theta
  # n_phi is the number of points for phi
  # theta is the vector for theta
  # phi is the vector for phi
  # nodes is the a matrix of all the triangulation nodes
  set.seed(13)
  if(is.null(nodes)){ # if nodes are not provided
      if(is.null(theta)){
        theta <- seq(0.1, 180-0.1, length.out = n_theta)#[-c(1, n_theta+3)]  # Angles en degrés
      }
      if(is.null(phi)){
        phi <- seq(0, 360, length.out = n_phi+1)  # Angles en degrés
      }
      if(!(0 %in% phi)){phi = c(0,phi)} 
      if(!(360%in%phi)){phi = c(phi,360)}
      
      points_2d = expand.grid(theta = theta, phi = phi)
      
      points_2d = points_2d[order(points_2d$phi),]
  }
  else{
      colnames(nodes) = c("theta","phi")
      points_2d = nodes[order(nodes$phi),]
      phi = sort(unique(nodes$phi))
      theta = sort(unique(nodes$theta))
  }
  if(0 %in% (points_2d$theta%%180)){stop("A pole is in the triangulation")} #we do not want theta = 0 or theta = 180

  points_3d = cbind(radius*sin(points_2d[,1]*pi/180)*cos(points_2d[,2]*pi/180), radius*sin(points_2d[,1]*pi/180)*sin(points_2d[,2]*pi/180), radius*cos(points_2d[,1]*pi/180))
  
  triangles <- delaunayn(points_2d)  # Triangulation Delaunay en 2D
  

  fn_unique = function(i){ #for points with phi = 2pi, return the index of the same point with phi = 0. Else just return the index of the point
      if (i %in% (nrow(points_2d)-length(theta)+1):nrow(points_2d)){return(i-(nrow(points_2d)-length(theta)))}
      else{return(i)}
    }
  associate_unique_points = lapply(1:nrow(points_2d), fn_unique) #list that associate to each point the value of fn_unique
  
  idx_unique = unique(unlist(associate_unique_points))
  
  old_points_2d = points_2d #save the dataframe of points including phi = 2pi
  points_2d = points_2d[idx_unique,] #remove points with phi = 2pi 
  points_3d = points_3d[idx_unique,] #remove points with phi = 2pi 
  

  
  old_triangles = triangles #save the triangles including the points with phi = 2pi

  triangles <- apply(triangles, c(1, 2), function(value) associate_unique_points[[value]]) #in the triangles, replace the indices of phi =2pi with the indices of phi = 0
  wrong_triangles = which(apply(triangles,1, function(aa){length(unique(aa))<3})) #If a triangle now has twice the same point, remove it
  if(length(wrong_triangles)>0){
  triangles = triangles[-wrong_triangles,]
  old_triangles = old_triangles[-wrong_triangles,]
  }
  dic_tri = lapply(1:nrow(points_2d),{function(k){c()}}) #dic_tri provides for each point the triangles it belongs to
  for(i in 1:nrow(triangles)){
    for(j in 1:ncol(triangles)){
      dic_tri[[triangles[i,j]]] = c(dic_tri[[triangles[i,j]]], i)
    }
  }
  return(list(triangles = triangles, points_2d = points_2d, points_3d = points_3d, dic_tri = dic_tri, radius = radius, old_triangles = old_triangles, old_points_2d = old_points_2d))
  
}

psi_trunc <- function(s, i, tri_idx, list_tri) { #psi_trunc compute for a given point s, a given function psi_i and a given triangle tri_idx, the value of the restriction of psi_i on tri_idx, evaluated at s
  tri = list_tri$triangles[tri_idx, ]
  if (!(i %in% tri)) { #if i is not in tri, psi_i = 0 on tri
    return(0)
  } else {
    
    if(sum(tri == list_tri$old_triangles[tri_idx,]) != 3){ #if tri contains a point with phi = 360
      tri_2d <- list_tri$old_points_2d[list_tri$old_triangles[tri_idx,], , drop = FALSE] #go to old triangle including phi = 360 and not phi = 0
    }
    else{tri_2d <- list_tri$points_2d[tri, , drop = FALSE]} #else get the triangle
  
    tri_0 <- as.matrix(t(rbind(tri_2d[1, ] - tri_2d[3, ],
                               tri_2d[2, ] - tri_2d[3, ]))) #matrix to change to barycentric coordinates
    s_0 <- as.numeric(s - tri_2d[3, ])
    barycentric_coordinates <- solve(tri_0)%*% s_0 #barycentric coordinates of s in tri
    if ((sum(barycentric_coordinates) > 1) || any(barycentric_coordinates < 0)) {
      return(0)
    } else if (i == tri[[3]]) {
      return(1 - sum(barycentric_coordinates))
    } else if (i == tri[[2]]) {
      return(barycentric_coordinates[2])
    } else {
      return(barycentric_coordinates[1])
    }
  }
}
