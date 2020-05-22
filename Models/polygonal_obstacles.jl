""" ---------------------------------
    ToDo: Replace this custom-made signed
          distance function for polytopic
          objects with faster library,
          e.g. FCL.:
          
    github.com/BerkeleyAutomation/python-fcl/
          (ongoing work, please contact
            thomas.lew@stanford.edu 
           for more information.)
"""

mutable struct PolygonalObstacle
    # center
    c # (3,)  

    # half widths
    dx
    dy
    dz

    widths

    # normals
    n_vec # (3,)

    # positions of corners
    p_vec # (3,)
end

function PolygonalObstacle(c::Vector, 
                           dx::Real, dy::Real, dz::Real,
                           n_vec::Vector, p_vec::Vector)
    # full width
    widths = 2. * Vector([dx,dy,dz])

    PolygonalObstacle(c, dx, dy, dz,
                      widths,
                      n_vec, p_vec)
end

function PolygonalObstacle(center::Vector, widths_vec::Vector)
  # center     : (3,)
  # widths_vec : (3,) width along each dimension (2*radius)
  c  = center
  dx = widths_vec[1] / 2.
  dy = widths_vec[2] / 2.
  dz = widths_vec[3] / 2.

  p_vec = compute_corners_positions(c[1], c[2], c[3], dx, dy, dz)
  n_vec = compute_surfaces_normals( c[1], c[2], c[3], dx, dy, dz)
  return PolygonalObstacle(c, dx,dy,dz,
                           n_vec, p_vec)
end

function compute_surfaces_normals(cx::Real, cy::Real, cz::Real, 
                                  dx::Real, dy::Real, dz::Real)
  n1 = Vector([ 0., -1.,  0.])
  n2 = Vector([ 1.,  0.,  0.])
  n3 = Vector([ 0.,  1.,  0.])
  n4 = Vector([-1.,  0.,  0.])
  n5 = Vector([ 0.,  0.,  1.])
  n6 = Vector([ 0.,  0., -1.])
  return [n1, n2, n3, n4, n5, n6]
end

function get_surfaces_normals(obs::PolygonalObstacle)
  # obs::PolygonalObstacle
  n = obs.n_vec
  return n[1], n[2], n[3], n[4], n[5], n[6]
end


function compute_corners_positions(cx::Real, cy::Real, cz::Real, 
                                   dx::Real, dy::Real, dz::Real)
  p125 = Vector([cx + dx, cy - dy, cz + dz])
  p235 = Vector([cx + dx, cy + dy, cz + dz])
  p345 = Vector([cx - dx, cy + dy, cz + dz])
  p145 = Vector([cx - dx, cy - dy, cz + dz])
  p126 = Vector([cx + dx, cy - dy, cz - dz])
  p236 = Vector([cx + dx, cy + dy, cz - dz])
  p346 = Vector([cx - dx, cy + dy, cz - dz])
  p146 = Vector([cx - dx, cy - dy, cz - dz])
  return [p125, p235, p345, p145, p126, p236, p346, p146]
end

function get_corners_positions(obs::PolygonalObstacle)
  # obs::PolygonalObstacle
  p = obs.p_vec
  p125, p235, p345, p145, p126, p236, p346, p146 = (p[1], 
        p[2], p[3], p[4], p[5], p[6], p[7], p[8])
  return p125, p235, p345, p145, p126, p236, p346, p146
end

function B_is_in_zone_i(x::Vector, i::Int, 
                        obs::PolygonalObstacle)
  # x::Vector
  # i::Int
  # obs::PolygonalObstacle
  c = obs.c
  n = obs.n_vec[i]
  dx, dy, dz = obs.dx, obs.dy, obs.dz

  if (i == 1 || i == 3)
    return (dot(x-c,n) >= dy)
  elseif (i == 2 || i == 4)
    return (dot(x-c,n) >= dx)
  elseif (i == 5 || i == 6)
    return (dot(x-c,n) >= dz)
  else
    error("[polygonal_obstacles.jl::is_in_zone_i] i not recognized.")
  end
end

function signed_distance(x::Vector, obs::PolygonalObstacle)
  # x::Vector 
  # obs::PolygonalObstacle

  # extract values
  c                                              = obs.c
  cx, cy, cz                                     = obs.c[1], obs.c[2], obs.c[3]
  dx, dy, dz                                     = obs.dx,   obs.dy,   obs.dz
  p125, p235, p345, p145, p126, p236, p346, p146 = get_corners_positions(obs)
  n1, n2, n3, n4, n5, n6                         = get_surfaces_normals(obs)

  # check in which zones it is
  Bz_vec = np.zeros(6, dtype=bool) # vector of booleans: i-th value is true if x is in i-th zone
  for i = 1:6
    Bz_vec[i] = B_is_in_zone_i(x, i, obs)
  end

  # x is outside a corner
  if (Bz_vec[1] && Bz_vec[2] && Bz_vec[5]) # in 125
    return norm(x-p125)
  elseif (Bz_vec[2] && Bz_vec[3] && Bz_vec[5]) # in 235
    return norm(x-p235)
  elseif (Bz_vec[3] && Bz_vec[4] && Bz_vec[5]) # in 345
    return norm(x-p345)
  elseif (Bz_vec[1] && Bz_vec[4] && Bz_vec[5]) # in 145
    return norm(x-p145)
  elseif (Bz_vec[1] && Bz_vec[2] && Bz_vec[6]) # in 126
    return norm(x-p126)
  elseif (Bz_vec[2] && Bz_vec[3] && Bz_vec[6]) # in 236
    return norm(x-p236)
  elseif (Bz_vec[3] && Bz_vec[4] && Bz_vec[6]) # in 346
    return norm(x-p346)
  elseif (Bz_vec[1] && Bz_vec[4] && Bz_vec[6]) # in 146
    return norm(x-p146)

  # x is outside an edge
  elseif (Bz_vec[1] && Bz_vec[2]) # in 12
    xhat = [x[1], x[2], cz+dz] 
    return norm(xhat-p125)
  elseif (Bz_vec[2] && Bz_vec[3]) # in 23
    xhat = [x[1], x[2], cz+dz] 
    return norm(xhat-p235)
  elseif (Bz_vec[3] && Bz_vec[4]) # in 34
    xhat = [x[1], x[2], cz+dz] 
    return norm(xhat-p345)
  elseif (Bz_vec[1] && Bz_vec[4]) # in 14
    xhat = [x[1], x[2], cz+dz] 
    return norm(xhat-p145)
  elseif (Bz_vec[1] && Bz_vec[5]) # in 15
    xhat = [x[1]+dx, x[3], x[3]] 
    return norm(xhat-p125)
  elseif (Bz_vec[2] && Bz_vec[5]) # in 25
    xhat = [x[1], x[2]-dy, x[3]] 
    return norm(xhat-p125)
  elseif (Bz_vec[3] && Bz_vec[5]) # in 35
    xhat = [x[1]+dx, x[3], x[3]] 
    return norm(xhat-p235)
  elseif (Bz_vec[4] && Bz_vec[5]) # in 45
    xhat = [x[1], x[2]+dy, x[3]] 
    return norm(xhat-p345)
  elseif (Bz_vec[1] && Bz_vec[6]) # in 16
    xhat = [x[1]+dx, x[3], x[3]] 
    return norm(xhat-p126)
  elseif (Bz_vec[2] && Bz_vec[6]) # in 26
    xhat = [x[1], x[2]-dy, x[3]] 
    return norm(xhat-p126)
  elseif (Bz_vec[3] && Bz_vec[6]) # in 36
    xhat = [x[1]+dx, x[3], x[3]] 
    return norm(xhat-p236)
  elseif (Bz_vec[4] && Bz_vec[6]) # in 46
    xhat = [x[1], x[2]+dy, x[3]] 
    return norm(xhat-p346)

  # x is outside a flat edge (not on border)
  elseif Bz_vec[1] # in 1
    return (dot(x-c,n1) - dy)
  elseif Bz_vec[2] # in 2
    return (dot(x-c,n2) - dx)
  elseif Bz_vec[3] # in 3
    return (dot(x-c,n3) - dy)
  elseif Bz_vec[4] # in 4
    return (dot(x-c,n4) - dx)
  elseif Bz_vec[5] # in 5
    return (dot(x-c,n5) - dz)
  elseif Bz_vec[6] # in 6
    return (dot(x-c,n6) - dz)

  # is inside
  else
    return (norm(x-c) - sqrt(dx^2+dy^2+dz^2)) # not quite accurate, but good enough for obstacle avoidance
  end
end

function signed_distance_with_closest_point_on_surface(x::Vector, obs::PolygonalObstacle)
  # x::Vector 
  # obs::PolygonalObstacle

  # extract values
  c                                              = obs.c
  cx, cy, cz                                     = obs.c[1], obs.c[2], obs.c[3]
  dx, dy, dz                                     = obs.dx,   obs.dy,   obs.dz
  p125, p235, p345, p145, p126, p236, p346, p146 = get_corners_positions(obs)
  n1, n2, n3, n4, n5, n6                         = get_surfaces_normals(obs)


  # check in which zones it is
  Bz_vec = falses(6) # vector of booleans: i-th value is true if x is in i-th zone
  for i = 1:6
    Bz_vec[i] = B_is_in_zone_i(x, i, obs)
  end

  # x is outside a corner
  if (Bz_vec[1] && Bz_vec[2] && Bz_vec[5]) # in 125
    return norm(x-p125), p125
  elseif (Bz_vec[2] && Bz_vec[3] && Bz_vec[5]) # in 235
    return norm(x-p235), p235
  elseif (Bz_vec[3] && Bz_vec[4] && Bz_vec[5]) # in 345
    return norm(x-p345), p345
  elseif (Bz_vec[1] && Bz_vec[4] && Bz_vec[5]) # in 145
    return norm(x-p145), p145
  elseif (Bz_vec[1] && Bz_vec[2] && Bz_vec[6]) # in 126
    return norm(x-p126), p126
  elseif (Bz_vec[2] && Bz_vec[3] && Bz_vec[6]) # in 236
    return norm(x-p236), p236
  elseif (Bz_vec[3] && Bz_vec[4] && Bz_vec[6]) # in 346
    return norm(x-p346), p346
  elseif (Bz_vec[1] && Bz_vec[4] && Bz_vec[6]) # in 146
    return norm(x-p146), p146

  # x is outside an edge
  elseif (Bz_vec[1] && Bz_vec[2]) # in 12
    xhat = [x[1], x[2], cz+dz]
    ds = norm(xhat-p125)
    return ds, (x-(xhat-p125))
  elseif (Bz_vec[2] && Bz_vec[3]) # in 23
    xhat = [x[1], x[2], cz+dz]
    ds = norm(xhat-p235)
    return ds, (x-(xhat-p235))
  elseif (Bz_vec[3] && Bz_vec[4]) # in 34
    xhat = [x[1], x[2], cz+dz]
    ds = norm(xhat-p345)
    return ds, (x-(xhat-p345))
  elseif (Bz_vec[1] && Bz_vec[4]) # in 14
    xhat = [x[1], x[2], cz+dz]
    ds = norm(xhat-p145)
    return ds, (x-(xhat-p145))
  elseif (Bz_vec[1] && Bz_vec[5]) # in 15
    xhat = [c[1]+dx, x[2], x[3]]
    ds = norm(xhat-p125)
    return ds, (x-(xhat-p125))
  elseif (Bz_vec[2] && Bz_vec[5]) # in 25
    xhat = [x[1], c[2]-dy, x[3]]
    ds = norm(xhat-p125)
    return ds, (x-(xhat-p125))
  elseif (Bz_vec[3] && Bz_vec[5]) # in 35
    xhat = [c[1]+dx, x[2], x[3]]
    ds = norm(xhat-p235)
    return ds, (x-(xhat-p235))
  elseif (Bz_vec[4] && Bz_vec[5]) # in 45
    xhat = [x[1], c[2]+dy, x[3]]
    ds = norm(xhat-p345)
    return ds, (x-(xhat-p345))
  elseif (Bz_vec[1] && Bz_vec[6]) # in 16
    xhat = [c[1]+dx, x[2], x[3]]
    ds = norm(xhat-p126)
    return ds, (x-(xhat-p126))
  elseif (Bz_vec[2] && Bz_vec[6]) # in 26
    xhat = [x[1], c[2]-dy, x[3]]
    ds = norm(xhat-p126)
    return ds, (x-(xhat-p126))
  elseif (Bz_vec[3] && Bz_vec[6]) # in 36
    xhat = [c[1]+dx, x[2], x[3]]
    ds = norm(xhat-p236)
    return ds, (x-(xhat-p236))
  elseif (Bz_vec[4] && Bz_vec[6]) # in 46
    xhat = [x[1], c[2]+dy, x[3]]
    ds = norm(xhat-p346)
    return ds, (x-(xhat-p346))

  # x is outside a flat edge (not on border)
  elseif Bz_vec[1] # in 1
    ds = (dot(x-c,n1)-dy)
    return ds, (x - n1*ds)
  elseif Bz_vec[2] # in 2
    ds = (dot(x-c,n2)-dx)
    return ds, (x - n2*ds)
  elseif Bz_vec[3] # in 3
    ds = (dot(x-c,n3)-dy)
    return ds, (x - n3*ds)
  elseif Bz_vec[4] # in 4
    ds = (dot(x-c,n4)-dx)
    return ds, (x - n4*ds)
  elseif Bz_vec[5] # in 5
    ds = (dot(x-c,n5)-dz)
    return ds, (x - n5*ds)
  elseif Bz_vec[6] # in 6
    ds = (dot(x-c,n6)-dz)
    return ds, (x - n6*ds)

  # is inside
  else
    # not quite accurate, but good enough for obstacle avoidance
    ds = 1.01*(norm(x-c) - sqrt(dx^2+dy^2+dz^2))
    pt_outside = (x+(x-c)*(-ds))
    _, pt_on_surface = signed_distance_with_closest_point_on_surface(pt_outside, obs)
    ds = -norm(x-pt_on_surface)
    return ds, pt_on_surface
  end
end