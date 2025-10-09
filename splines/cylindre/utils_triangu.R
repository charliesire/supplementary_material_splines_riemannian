
build_triangu = function(min_h, max_h, radius, n_circ, n_height = NULL, nodes = NULL){ #build triangulation for cylinder
  # min_h and max_h are the lower and upper bounds of the longitudinal coordinate
  # radius is the radius of the cylinder
  # n_circ the number of triangulation points on each section
  # Optional n_height is the number of triangulation point for a given theta
  # Optional nodes if the user provides the triangulation points 
  set.seed(13)
  if(is.null(nodes)){
    theta <- seq(0, 360, length.out = n_circ + 1)  # Angles in degrees
    if(is.null(n_height)){ #if n_height not provided, compute it to have almost isosceles triangles
      pas = theta[2]-theta[1]
      pas_z = radius*sqrt(2-2*cos(pas*pi/180))
      reste = pas_z - (max_h- min_h)%%pas_z
      z = seq(min_h-reste/2-pas_z, max_h + reste/2+pas_z, pas_z)
      n_height = length(z)
    }
    else{ # Leave some margin along the longitudinal axis to avoid edge effects
      min_l = min_h
      max_l = max_h
      z = seq(min_l, max_l, l=n_height)
      pas_z = (max_l-min_l)/(n_height-1)
      while((min_h-z[1]) < pas_z){
        min_l = min_l-0.01
        max_l = max_l + 0.01
        pas_z = (max_l-min_l)/(n_height-1)
        z = seq(min_l, max_l, l=n_height)
      }
    }
    points_2d = expand.grid(theta = theta, z = z)
    points_2d = points_2d[order(points_2d$theta),]
    Ns = n_height*n_circ
  }
  else{
    colnames(nodes) = c("theta","z")
    points_2d = nodes[order(nodes$theta),]
    z = sort(unique(nodes$z))
    theta = sort(unique(nodes$theta))
  }
    
  points_3d = cbind(radius*cos(points_2d[,1]*pi/180), radius*sin(points_2d[,1]*pi/180), points_2d[,2])
  
  triangles <- delaunayn(points_2d)  # Triangulation Delaunay en 2D

  fn_unique = function(i){ #for points with theta = 2pi, return the index of the same point with theta = 0. Else just return the index of the point
    if (i %in% (nrow(points_2d)-length(z)+1):nrow(points_2d)){return(i-(nrow(points_2d)-length(z)))}
    else{return(i)}
  }
  
  points_360 = which(points_2d[,1] == 360) # points at theta = 2pi
  warning_360 = which(apply(triangles, 1,function(tri){any(tri %in% points_360)})) #triangles with points at theta = 2pi
  
  associate_unique_points = lapply(1:nrow(points_2d), fn_unique) #list that associate to each point the value of fn_unique
  
  idx_unique = unique(unlist(associate_unique_points))
  
  points_2d = points_2d[setdiff(1:nrow(points_2d),points_360),] #remove points with theta = 2pi 
  points_3d = points_3d[setdiff(1:nrow(points_2d),points_360),]
  
  triangles <- apply(triangles, c(1, 2), function(value) associate_unique_points[[value]]) #in the triangles, replace the indices of theta =2pi with the indices of theta = pi
  
  dic_tri = lapply(1:nrow(points_2d),{function(k){c()}}) #dic_tri provides for each point the triangles it belongs to
  for(i in 1:nrow(triangles)){
    for(j in 1:ncol(triangles)){
      dic_tri[[triangles[i,j]]] = c(dic_tri[[triangles[i,j]]], i)
    }
  }
  return(list(triangles = triangles, points_2d = points_2d, points_3d = points_3d, dic_tri = dic_tri, radius = radius, warning_360 = warning_360))
  
}

psi_trunc <- function(s, i, tri_idx, list_tri) { #psi_trunc compute for a given point s, a given function psi_i and a given triangle tri_idx, the value of the restriction of psi_i on tri_idx, evaluated at s
  tri = list_tri$triangles[tri_idx, ]
  if (!(i %in% tri)) {#if i is not in tri, psi_i = 0 on tri
    return(0)
  } else {
    tri_2d <- list_tri$points_2d[tri, , drop = FALSE] #get triangle
    if(tri_idx %in% list_tri$warning_360){tri_2d[tri_2d[,1] == 0,1] = 360} #if warning_360, put theta = 2pi instead of theta = 0
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
