mutable struct ParamType{Tdim, T}
  res_vals::Array{T, 1}
  res_vals2::Array{T, 1}
  q_vals::Array{T, 1}
  v_vals2::Array{T, 1}
  flux_vals1::Array{T, 1}
  sat_vals::Array{T, 1}
  euler_fluxjac::Array{T, 2}
  nrm::Array{T, 1}
  gamma::Float64
  gamma_1::Float64

  function ParamType{Tdim, T}(numDofPerNode) where {Tdim, T}
    res_vals = zeros(T, numDofPerNode)
    res_vals2 = zeros(T, numDofPerNode)
    q_vals = zeros(T, numDofPerNode)
    v_vals2 = zeros(T, numDofPerNode)
    flux_vals1 = zeros(T, numDofPerNode)
    sat_vals = zeros(T, numDofPerNode)
    euler_fluxjac = zeros(T, numDofPerNode, numDofPerNode)
    nrm = zeros(T, Tdim)
    gamma = 1.4
    gamma_1 = 0.4

    return new(res_vals, res_vals2, q_vals, v_vals2, flux_vals1, sat_vals, euler_fluxjac, nrm, gamma, gamma_1)
  end
end



@doc """
### EulerEquationMod.RoeSolver
  This calculates the Roe flux for boundary conditions at a node. The inputs
  must be in *conservative* variables.

  **Inputs**
   * params : ParamType
   * q  : conservative variables of the fluid
   * qg : conservative variables of the boundary
   * aux_vars : vector of all auxiliary variables at this node
   * nrm : scaled face normal vector (x-y space, outward normal of the face q
        lives on)
   * use_efix: 1 = use entropy fix, 0 = do not use entropy fix (integer)

  **Outputs**
   * flux : vector to populate with solution

  Aliasing restrictions:  none of the inputs can alias params.res_vals1,
                          params.res_vals2, params.q_vals, params.flux_vals1, or
                          params.sat


"""->
function RoeSolver{Tmsh, Tsol, Tres}(params::ParamType{2},
                                     q::AbstractArray{Tsol,1},
                                     qg::AbstractArray{Tsol, 1},
                                     aux_vars::AbstractArray{Tres, 1},
                                     nrm::AbstractArray{Tmsh,1},
                                     flux::AbstractArray{Tres, 1},
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
  fac = d1_0/q[1]
  uL = q[2]*fac; vL = q[3]*fac;
  phi = d0_5*(uL*uL + vL*vL)
  HL = gamma*q[4]*fac - gami*phi # Total enthalpy, H = e + 0.5*(u^2 + v^2) + p/rho,
                                 # where e is the internal energy per unit mass

  # The right side of the Roe solver comprises the boundary conditions
  fac = d1_0/qg[1]
  uR = qg[2]*fac; vR = qg[3]*fac;
  phi = d0_5*(uR*uR + vR*vR)
  HR = gamma*qg[4]*fac - gami*phi # Total Enthalpy

  # Averaged states
  sqL = sqrt(q[1])
  sqR = sqrt(qg[1])
  fac = d1_0/(sqL + sqR)
  u = (sqL*uL + sqR*uR)*fac
  v = (sqL*vL + sqR*vR)*fac
  H = (sqL*HL + sqR*HR)*fac

  dq = params.v_vals2 # zeros(Tsol, 4)
  for i=1:length(dq)
    dq[i] = q[i] - qg[i]
  end
  sat = params.sat_vals
#  calcSAT(params, nrm, dq, sat, u, v, H, use_efix)
  
  sat[1] = q[1]*dq[1]
  sat[2] = u*dq[2]
  sat[3] = v*dq[3]
  sat[4] = H*dq[4]

  euler_flux = params.flux_vals1
  # calculate Euler flux in wall normal directiona
  # because edge numbering is rather arbitary, any memory access is likely to
  # be a cache miss, so we recalculate the Euler flux
  v_vals = params.q_vals
  nrm2 = params.nrm
  nrm2[1] = nx   # why are we assigning to nrm2?
  nrm2[2] = ny

#  convertFromNaturalToWorkingVars(params, q, v_vals)
  calcEulerFlux(params, q, aux_vars, nrm2, euler_flux)

  for i=1:4  # ArrayViews does not support flux[:] = .
    flux[i] = (sat_fac*sat[i] + euler_flux[i])
  end

  return nothing

end # ends the function eulerRoeSAT

function calcEulerFlux{Tmsh, Tsol, Tres}(params::ParamType{2},
                      q::AbstractArray{Tsol,1},
                      aux_vars::AbstractArray{Tres, 1},
                      dir::AbstractArray{Tmsh},  F::AbstractArray{Tsol,1})
# calculates the Euler flux in a particular direction at a point
# eqn is the equation type
# q is the vector (of length 4), of the conservative variables at the point
# aux_vars is the vector of auxiliary variables at the point
# dir is a vector of length 2 that specifies the direction
# F is populated with the flux (is a vector of length 4)
# 2D  only


  press = calcPressure(params, q)
#  press = getPressure(aux_vars)
#  press = @getPressure(aux_vars)
  U = (q[2]*dir[1] + q[3]*dir[2])/q[1]
  F[1] = q[1]*U
  F[2] = q[2]*U + dir[1]*press
  F[3] = q[3]*U + dir[2]*press
  F[4] = (q[4] + press)*U

  return nothing

end

function calcPressure{Tsol}(params::ParamType{2},
                            q::AbstractArray{Tsol,1} )
  # calculate pressure for a node
  # q is a vector of length 4 of the conservative variables

  return  (params.gamma_1)*(q[4] - 0.5*(q[2]*q[2] + q[3]*q[3])/q[1])
end



