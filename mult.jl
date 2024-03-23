#
 # Author:  Lamine Djibo, ldjibo2016@fit.edu
 # Course:  CSE 4250, Section 01, Fall 2023
 #  Project: Proj 01, Multilateration
 #/

# Program to find location of object in 3D

# import libraries
using Printf
# using Statistics

# Functions

# # function to get the Rletters
# function getR_letter(coord::Vector)
#   return sqrt((0-coord[1])^2 + (0-coord[2])^2 + (0-coord[3])^2)
# end

# function to get difference in Rs using the speed of light formula
function getR_between(sat1::Float64, sat2::Float64)
  return (sat1 - sat2)*299792458
end

# function to get coordinates difference
function getCoord_Diff(num1::Float64, num2::Float64)
  return (num1 - num2)
end

# function to return sign
function sgn(num::Float64)
  if num < 0
    return -1
  elseif num == 0
    return 0
  else
    return 1
  end
end



# Initialize an empty array to store the lines
lines = []
count = 1

#object list of times
obj = []


# filename as an argument
if length(ARGS) < 1
  println("No file for input")
  exit(1)
end

# Get the filename from the first argument
file_path = ARGS[1]

#read file into a list of lines
try
  global lines = readlines(file_path)
catch e
  println("Cannot read file: $e")
end

#println(lines)


# sattelite lists with coordinates, index starts at 1
# strip of left and right spaces, split using delimiters and parse to Float64
satI = parse.(Float64, split(strip(lines[1])))
satJ = parse.(Float64, split(strip(lines[2])))
satK = parse.(Float64, split(strip(lines[3])))
satL = parse.(Float64, split(strip(lines[4])))


# convert time to seconds, use index for updating
for o in obj
  for i in 1:length(o)
    o[i] = o[i] / 1e9
  end
end

# if line is not empty append to obj list to be located for
# computing after parsing and convert nanosecond to second
count = 1
for i in 5:(length(lines))
  if lines[i] != " "
    push!(obj, parse.(Float64, split(strip(lines[i]))))
    obj[count] = obj[count] / 1e9
    global count += 1
  end
end

#println(obj)


# Loop for all objects, one by one
for i in 1:length(obj)

  #println(obj[i])
  # get the Rs for each satellite
  R_i = obj[i][1]
  R_j = obj[i][2]
  R_k = obj[i][3]
  R_l = obj[i][4]


  # get the Rs between satellites, TDOA * speed of light
  R_ij = getR_between(R_i, R_j)
  R_ik = getR_between(R_i, R_k)
  R_kj = getR_between(R_k, R_j)
  R_kl = getR_between(R_k, R_l)


  # compute the big difference in x
  x_ki = getCoord_Diff(satK[1], satI[1])
  x_ji = getCoord_Diff(satJ[1], satI[1])
  x_lk = getCoord_Diff(satL[1], satK[1])
  x_jk = getCoord_Diff(satJ[1], satK[1])

  # compute the big difference in y
  y_ki = getCoord_Diff(satK[2], satI[2])
  y_ji = getCoord_Diff(satJ[2], satI[2])
  y_lk = getCoord_Diff(satL[2], satK[2])
  y_jk = getCoord_Diff(satJ[2], satK[2])

  # compute the big difference in z
  z_ki = getCoord_Diff(satK[3], satI[3])
  z_ji = getCoord_Diff(satJ[3], satI[3])
  z_lk = getCoord_Diff(satL[3], satK[3])
  z_jk = getCoord_Diff(satJ[3], satK[3])


  # Compute the big Xs
  X_ijy = R_ij * y_ki - R_ik * y_ji
  X_ikx = R_ik * x_ji - R_ij * x_ki
  X_ikz = R_ik * z_ji - R_ij * z_ki
  X_kjy = R_kj * y_lk - R_kl * y_jk
  X_klx = R_kl * x_jk - R_kj * x_lk
  X_klz = R_kl * z_jk - R_kj * z_lk


  # Compute the big Ss 
  Si2 = satI[1]^2 + satI[2]^2 + satI[3]^2
  Sj2 = satJ[1]^2 + satJ[2]^2 + satJ[3]^2
  Sk2 = satK[1]^2 + satK[2]^2 + satK[3]^2
  Sl2 = satL[1]^2 + satL[2]^2 + satL[3]^2


  # Compute the Rs^2
  Rij2xyz = R_ij^2 + Si2 - Sj2
  Rik2xyz = R_ik^2 + Si2 - Sk2
  Rkj2xyz = R_kj^2 + Sk2 - Sj2
  Rkl2xyz = R_kl^2 + Sk2 - Sl2 
  
  
  # Compute A, B, C, D
  A = X_ikx / X_ijy
  B = X_ikz / X_ijy
  C = X_klx / X_kjy
  D = X_klz / X_kjy


  # Compute E, F
  E = (R_ik * Rij2xyz - R_ij * Rik2xyz) / (2 * X_ijy)
  F = (R_kl * Rkj2xyz - R_kj * Rkl2xyz) / (2 * X_kjy)


  # Compute G, H
  G = (D - B) / (A - C)
  H = (F - E) / (A - C)


  # Compute I, J 
  I = (A * G) + B
  J = (A * H) + E


  # Compute K, L 
  K = Rik2xyz + (2 * x_ki * H) + (2 * y_ki * J)
  L = 2 * ((x_ki * G) + (y_ki * I) + z_ki)


  # compute M, N, O 
  M = (4 * (R_ik^2) * (G^2 + I^2 + 1)) - L^2
  N = (8 * (R_ik^2) * (G*(satI[1]-H) + I*(satI[2]-J) + satI[3])) + (2*L*K)
  O = (4 * (R_ik^2) * ((satI[1]-H)^2 + (satI[2]-J)^2 + satI[3]^2)) - (K^2)

  # println(M)
  # println(N)
  # println(O)
  

  # Compute Q
  signed = sgn(N)
  Q = N + (signed * sqrt(N^2 - (4*M*O)))


  # Compute 2 potential z coordinates, z1 and z2
  z1 = Q / (2*M)
  z2 = (2*O) / Q

  x1 = (G*z1) + H
  x2 = (G*z2) + H

  y1 = (I*z1) + J
  y2 = (I*z2) + J

  # Compute r, distance from origin
  r1 = sqrt(x1^2 + y1^2 + z1^2)
  r2 = sqrt(x2^2 + y2^2 + z2^2)

  @printf("g= %10.2e, h= %10.2e, j= %10.2e, m= %10.2e, o= %10.2e \n", G, H, J, M, O)
  @printf("+) x= %10d, y= %10d, z= %10d; r= %10d \n", round(x1), round(y1), round(z1), round(r1))
  @printf("-) x= %10d, y= %10d, z= %10d; r= %10d \n", round(x2), round(y2), round(z2), round(r2))

  # print("x1:  ")
  # print(x1)
  # print("  y1:  ")
  # print(y1)
  # print("  z1:  ")
  # print(z1)
  # print("  r:  ")
  # println(r1)
  # println("\n")
end





