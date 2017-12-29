include("roe.jl")
#include("DualMod.jl")
include("DualMod_static.jl")

function main()
  q = [2.0, 3.0, 4.0, 7.0]
  qg = [2.0, 3.0, 4.0, 7.0]
  aux_vars = Float64[0.0]
  nrm = Float64[2.0, 2.0]
  F = zeros(Float64, 4)
  params = ParamType{2, Float64}(4)

  numEl = 10000
  res = zeros(4, numEl)
  @time for i=1:numEl
    q[2] =+ 0.01
    RoeSolver(params, q, qg, aux_vars, nrm, F)

    for j=1:4
      res[j, i] = F[j]
    end
  end
#  println("F = \n", F)
end

function main2()
  q = Dual{4}[2.0, 3.0, 4.0, 7.0]
  qg = Dual{4}[2.0, 3.0, 4.0, 7.0]
  aux_vars = Dual{4}[0.0]
  nrm = Dual{4}[2.0, 2.0]
  F = zeros(Dual{4}, 4)
  params = ParamType{2, Dual{4}}(4)

  numEl = 10000
  res = zeros(Dual{4}, 4, numEl)
  @time for i=1:numEl
    q[2] =+ 0.01
    RoeSolver(params, q, qg, aux_vars, nrm, F)

    for j=1:4
      res[j, i] = F[j]
    end
  end
#  println("F = \n", F)
end

main()
main()
main2()
main2()
