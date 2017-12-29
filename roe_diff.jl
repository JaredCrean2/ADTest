function RoeSolver_diff{Tmsh, Tsol, Tres}(params::ParamType{2},
                                     q::AbstractArray{Tsol,1},
                                     qg::AbstractArray{Tsol, 1},
                                     aux_vars::AbstractArray{Tres, 1},
                                     nrm::AbstractArray{Tmsh,1},
                                     fluxL_dot::AbstractArray{Tres, 2},
                                     fluxR_dot::AbstractArray{Tres, 2},
                                     use_efix::Int=1)

  # SAT terms are used for ensuring consistency with the physical problem. Its
  # similar to upwinding which adds dissipation to the problem. SATs on the
  # boundary can be thought of as having two overlapping nodes and because of
  # the discontinuous nature of SBP adds some dissipation.

  # Declaring constants
  d1_0 = 1.0
  d0_0 = 0.0
  d0_5 = 0.5
  tau = 1.0
  gamma = params.gamma
  gami = params.gamma_1
  sat_fac = 1  # multiplier for SAT term

  # Begin main executuion
  nx = nrm[1]
  ny = nrm[2]

  # Compute the Roe Averaged states
  # The left state of Roe are the actual solution variables
  # All the _dot variables are wrt q
  fac = d1_0/q[1]
  fac_dot1 = -fac*fac
#  fac_dot1 = -d1_0/(q[1]*q[1])

  uL = q[2]*fac
  uL_dot1 = q[2]*fac_dot1
  uL_dot2 = fac

  vL = q[3]*fac
  vL_dot1 = q[3]*fac_dot1
  vL_dot3 = fac

  phi = d0_5*(uL*uL + vL*vL)
  phi_dot1 = uL*uL_dot1 + vL*vL_dot1
  phi_dot2 = uL*uL_dot2
  phi_dot3 = vL*vL_dot3

  HL = gamma*q[4]*fac - gami*phi # Total enthalpy, H = e + 0.5*(u^2 + v^2) + p/rho,
                                 # where e is the internal energy per unit mass
  HL_dot1 = gamma*q[4]*fac_dot1 - gami*phi_dot1
  HL_dot2 = -gami*phi_dot2
  HL_dot3 = -gami*phi_dot3
  HL_dot4 = gamma*fac  # q[4]_dot = 1

  # The right side of the Roe solver comprises the boundary conditions
  # all the _dot variables are wrt qg now
  fac = d1_0/qg[1]
  fac_dot1 = -fac*fac
#  fac_dot1 = -d1_0/(qg[1]*qg[1])

  uR = qg[2]*fac
  uR_dot1 = qg[2]*fac_dot1
  uR_dot2 = fac

  vR = qg[3]*fac
  vR_dot1 = qg[3]*fac_dot1
  vR_dot3 = fac

  phi = d0_5*(uR*uR + vR*vR)
  phi_dot1 = uR*uR_dot1 + vR*vR_dot1
  phi_dot2 = uR*uR_dot2
  phi_dot3 = vR*vR_dot3


  HR = gamma*qg[4]*fac - gami*phi # Total Enthalpy
  HR_dot1 = gamma*qg[4]*fac_dot1 - gami*phi_dot1
  HR_dot2 = -gami*phi_dot2
  HR_dot3 = -gami*phi_dot3
  HR_dot4 = gamma*fac

  # Averaged states
  sqL = sqrt(q[1])
  sqL_dot1 = 0.5/sqL

  sqR = sqrt(qg[1])
  sqR_dot1 = 0.5/sqR

  fac = d1_0/(sqL + sqR)
  t1 = -1/((sqL + sqR)*(sqL + sqR))
  fac_dotL1 = t1*sqL_dot1
  fac_dotR1 = t1*sqR_dot1

  u = (sqL*uL + sqR*uR)*fac
  t2 = sqR*uR
  t3 = sqL*uL
  t4 = sqL*fac
  t5 = sqR*fac
  u_dotL1 = t3*fac_dotL1 + sqL*fac*uL_dot1 + uL*fac*sqL_dot1 + t2*fac_dotL1
  u_dotR1 = t2*fac_dotR1 + sqR*fac*uR_dot1 + uR*fac*sqR_dot1 + t3*fac_dotR1


  u_dotL2 = t4*uL_dot2
  u_dotR2 = t5*uR_dot2

  v = (sqL*vL + sqR*vR)*fac
  t2 = sqL*vL
  t3 = sqR*vR
  v_dotL1 = t2*fac_dotL1 + sqL*fac*vL_dot1 + vL*fac*sqL_dot1 + t3*fac_dotL1
  v_dotR1 = t3*fac_dotR1 + sqR*fac*vR_dot1 + vR*fac*sqR_dot1 + t2*fac_dotR1

  v_dotL3 = t4*vL_dot3
  v_dotR3 = t5*vR_dot3

  H = (sqL*HL + sqR*HR)*fac
  t2 = sqL*HL
  t3 = sqR*HR
  H_dotL1 = t2*fac_dotL1 + sqL*fac*HL_dot1 + HL*fac*sqL_dot1 + t3*fac_dotL1
  H_dotR1 = t3*fac_dotR1 + sqR*fac*HR_dot1 + HR*fac*sqR_dot1 + t2*fac_dotR1
 
  H_dotL2 = t4*HL_dot2 
  H_dotR2 = t5*HR_dot2

  H_dotL3 = t4*HL_dot3
  H_dotR3 = t5*HR_dot3

  H_dotL4 = t4*HL_dot4
  H_dotR4 = t5*HR_dot4


  dq = params.v_vals2 # zeros(Tsol, 4)
  for i=1:length(dq)
    dq[i] = q[i] - qg[i]
  end

  # dq_dotL* = 1, dq_dotR* = -1, so omit them

  sat = params.sat_vals
#  calcSAT(params, nrm, dq, sat, u, v, H, use_efix)
  
  sat[1] = q[1]*dq[1]
  sat1_dotL1 = dq[1]*1 + q[1]
  sat1_dotR1 = -q[1]

  sat[2] = u*dq[2]
  sat2_dotL1 = dq[2]*u_dotL1
  sat2_dotR1 = dq[2]*u_dotR1


  sat2_dotL2 = u + dq[2]*u_dotL2
  sat2_dotR2 = -u + dq[2]*u_dotR2


  sat[3] = v*dq[3]
  sat3_dotL1 = dq[3]*v_dotL1
  sat3_dotR1 = dq[3]*v_dotR1

  sat3_dotL3 = v + dq[3]*v_dotL3
  sat3_dotR3 = -v + dq[3]*v_dotR3

  sat[4] = H*dq[4]
  sat4_dotL1 = dq[4]*H_dotL1
  sat4_dotR1 = dq[4]*H_dotR1

  sat4_dotL2 = dq[4]*H_dotL2
  sat4_dotR2 = dq[4]*H_dotR2

  sat4_dotL3 = dq[4]*H_dotL3
  sat4_dotR3 = dq[4]*H_dotR3

  sat4_dotL4 = H + dq[4]*H_dotL4
  sat4_dotR4 = -H + dq[4]*H_dotR4

  #euler_flux = params.flux_vals1
  euler_fluxjac = params.euler_fluxjac

  v_vals = params.q_vals
  nrm2 = params.nrm
  nrm2[1] = nx   # why are we assigning to nrm2?
  nrm2[2] = ny

#  convertFromNaturalToWorkingVars(params, q, v_vals)
  calcEulerFlux_diff(params, q, aux_vars, nrm2, euler_fluxjac)

  fluxL_dot[1, 1] = sat_fac*sat1_dotL1 + euler_fluxjac[1, 1]
  fluxL_dot[2, 1] = sat_fac*sat2_dotL1 + euler_fluxjac[2, 1]
  fluxL_dot[3, 1] = sat_fac*sat3_dotL1 + euler_fluxjac[3, 1]
  fluxL_dot[4, 1] = sat_fac*sat4_dotL1 + euler_fluxjac[4, 1]

  fluxL_dot[1, 2] =                      euler_fluxjac[1, 2]
  fluxL_dot[2, 2] = sat_fac*sat2_dotL2 + euler_fluxjac[2, 2]
  fluxL_dot[3, 2] =                    + euler_fluxjac[3, 2]
  fluxL_dot[4, 2] = sat_fac*sat4_dotL2 + euler_fluxjac[4, 2]

  fluxL_dot[1, 3] =                      euler_fluxjac[1, 3]
  fluxL_dot[2, 3] =                      euler_fluxjac[2, 3]
  fluxL_dot[3, 3] = sat_fac*sat3_dotL3 + euler_fluxjac[3, 3]
  fluxL_dot[4, 3] = sat_fac*sat4_dotL3 + euler_fluxjac[4, 3]

  fluxL_dot[1, 4] =                      euler_fluxjac[1, 4]
  fluxL_dot[2, 4] =                      euler_fluxjac[2, 4]
  fluxL_dot[3, 4] =                    + euler_fluxjac[3, 4]
  fluxL_dot[4, 4] = sat_fac*sat4_dotL4 + euler_fluxjac[4, 4]

  fill!(fluxR_dot, 0.0)
  fluxR_dot[1, 1] = sat_fac*sat1_dotR1

  fluxR_dot[2, 1] = sat_fac*sat2_dotR1
  fluxR_dot[3, 1] = sat_fac*sat3_dotR1
  fluxR_dot[4, 1] = sat_fac*sat4_dotR1

  fluxR_dot[2, 2] = sat_fac*sat2_dotR2
  fluxR_dot[4, 2] = sat_fac*sat4_dotR2

  fluxR_dot[3, 3] = sat_fac*sat3_dotR3
  fluxR_dot[4, 3] = sat_fac*sat4_dotR3

  fluxR_dot[4, 4] = sat_fac*sat4_dotR4

  return nothing

end # ends the function eulerRoeSAT

function calcEulerFlux_diff{Tmsh, Tsol, Tres}(params::ParamType{2},
                      q::AbstractArray{Tsol,1},
                      aux_vars::AbstractArray{Tres, 1},
                      dir::AbstractArray{Tmsh},  Fjac::AbstractArray{Tsol,2})
# calculates the Euler flux in a particular direction at a point
# eqn is the equation type
# q is the vector (of length 4), of the conservative variables at the point
# aux_vars is the vector of auxiliary variables at the point
# dir is a vector of length 2 that specifies the direction
# F is populated with the flux Jacobian
# 2D  only


  p_dot = params.p_dot
#  p_dot = zeros(q)  # TODO: replace this with a pre-allocated array
  press = calcPressure_diff(params, q, p_dot)
#  press = getPressure(aux_vars)
#  press = @getPressure(aux_vars)
  fac = 1/q[1]
  U = (q[2]*dir[1] + q[3]*dir[2])*fac
  U_dot1 = -(q[2]*dir[1] + q[3]*dir[2])*fac*fac
  U_dot2 = dir[1]*fac
  U_dot3 = dir[2]*fac

  # F[1] = q[1]*U
  # F[2] = q[2]*U + dir[1]*press
  # F[3] = q[3]*U + dir[2]*press
  # F[4] = (q[4] + press)*U
  Fjac[1, 1] = U + q[1]*U_dot1
  Fjac[2, 1] = q[2]*U_dot1 + dir[1]*p_dot[1]
  Fjac[3, 1] = q[3]*U_dot1 + dir[2]*p_dot[1]
  Fjac[4, 1] = q[4]*U_dot1 + press*U_dot1 + U*p_dot[1]

  Fjac[1, 2] = q[1]*U_dot2
  Fjac[2, 2] = U + q[2]*U_dot2 + dir[1]*p_dot[2]
  Fjac[3, 2] = q[3]*U_dot2 + dir[2]*p_dot[2]
  Fjac[4, 2] = q[4]*U_dot2 + press*U_dot2 + U*p_dot[2]

  Fjac[1, 3] = q[1]*U_dot3
  Fjac[2, 3] = q[2]*U_dot3 + dir[1]*p_dot[3]
  Fjac[3, 3] = U + q[3]*U_dot3 + dir[2]*p_dot[3]
  Fjac[4, 3] = q[4]*U_dot3 + press*U_dot3 + U*p_dot[3]

  Fjac[1, 4] = 0
  Fjac[2, 4] = dir[1]*p_dot[4]
  Fjac[3, 4] = dir[2]*p_dot[4]
  Fjac[4, 4] = U + U*p_dot[4]

  return nothing

end

function calcPressure_diff{Tsol}(params::ParamType{2},
                            q::AbstractArray{Tsol,1}, p_dot::AbstractVector{Tsol} )
  # calculate pressure for a node
  # q is a vector of length 4 of the conservative variables

  t1 = 1/(q[1]*q[1])
  t2 = q[2]*q[2]
  t3 = q[3]*q[3]

  p_dot[1] = (params.gamma_1)*( 0.5*(t2*t1 + t3*t1))
  p_dot[2] = -(params.gamma_1)*(q[2]/q[1])
  p_dot[3] = -(params.gamma_1)*(q[3]/q[1])
  p_dot[4] = params.gamma_1
  return  (params.gamma_1)*(q[4] - 0.5*(t2 + t3)/q[1])
end



