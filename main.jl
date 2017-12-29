include("roe.jl")
#include("DualMod.jl")

function main()
  q = [2.0, 3.0, 4.0, 7.0]
  qg = [2.0, 3.0, 4.0, 7.0]
  aux_vars = Float64[0.0]
  nrm = Float64[2.0, 2.0]
  F = zeros(Float64, 4)
  params = ParamType{2, Float64}(4)
  RoeSolver(params, q, qg, aux_vars, nrm, F)
  println("F = \n", F)
end

main()
