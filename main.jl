include("roe.jl")
#include("DualMod.jl")
include("DualMod_static.jl")
include("roe_diff.jl")

global const NDUAL = 8
global const numEl = 100000
# flux evaluation
function main()
  q = [2.0, 3.0, 4.0, 7.0]
  qg = [2.0, 3.0, 4.0, 7.0]
  aux_vars = Float64[0.0]
  nrm = Float64[2.0, 2.0]
  F = zeros(Float64, 4)
  params = ParamType{2, Float64}(4)

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

# dual number evaluation
function main2()
  q = Dual{NDUAL}[2.0, 3.0, 4.0, 7.0]
  qg = Dual{NDUAL}[2.0, 3.0, 4.0, 7.0]
  aux_vars = Dual{NDUAL}[0.0]
  nrm = Dual{NDUAL}[2.0, 2.0]
  F = zeros(Dual{NDUAL}, 4)
  params = ParamType{2, Dual{NDUAL}}(4)

  res = zeros(Dual{NDUAL}, 4, numEl)
  @time for i=1:numEl
    q[2] =+ 0.01
    RoeSolver(params, q, qg, aux_vars, nrm, F)

    for j=1:4
      res[j, i] = F[j]
    end
  end

  return res
#  println("F = \n", F)
end

# complex numbers
function main3()
  q = Complex128[2.0, 3.0, 4.0, 7.0]
  qg = Complex128[2.0, 3.0, 4.0, 7.0]
  aux_vars = Complex128[0.0]
  nrm = Float64[2.0, 2.0]
  F = zeros(Complex128, 4)
  params = ParamType{2, Complex128}(4)

  resL = zeros(4, 4, numEl)
  resR = zeros(4, 4, numEl)
  h = 1e-20
  pert = Complex128(0, h)
  @time for i=1:numEl
    q[2] =+ 0.01

    for j=1:4
      q[j] += pert
      RoeSolver(params, q, qg, aux_vars, nrm, F)

      for k=1:4
        resL[k, j, i] = imag(F[k])/h
      end

      q[j] -= pert
    end

    for j=1:4
      qg[j] += pert
      RoeSolver(params, q, qg, aux_vars, nrm, F)

      for k=1:4
        resR[k, j, i] = imag(F[k])/h
      end

      qg[j] -= pert
    end
  end  # end loop i

  return resL, resR
end

# hand-written AD
function main4()
  q = [2.0, 3.0, 4.0, 7.0]
  qg = [2.0, 3.0, 4.0, 7.0]
  aux_vars = Float64[0.0]
  nrm = Float64[2.0, 2.0]
  F_dotL = zeros(Float64, 4, 4)
  F_dotR = zeros(Float64, 4, 4)
  params = ParamType{2, Float64}(4)

  resL = zeros(4, 4, numEl)
  resR = zeros(4, 4, numEl)
  @time for i=1:numEl
    q[2] =+ 0.01
    RoeSolver_diff(params, q, qg, aux_vars, nrm, F_dotL, F_dotR)

    for j=1:4
      for k=1:4
        resL[k, j, i] = F_dotL[k, j]
        resR[k, j, i] = F_dotR[k, j]
      end
    end
  end

  return resL, resR
end

println("Flux evaluation:")
main()
main()
println("Dual number {$NDUAL} evaluation:")
main2()
main2()
println("Complex number evaluation:")
main3()
main3()
println("Hand-written AD evaluation:")
main4()
main4()


println("\nVerifying AD")
res3L, res3R = main3()
res4L, res4R = main4()

#=
println("\n----- Comparing left results -----")
println("res3L = \n", res3L)
println("res4L = \n", res4L)
println("diff = \n", res3L - res4L)
println("diffnorm = \n", vecnorm(res3L - res4L))

println("\n----- Comparing right results -----")
println("res3R = \n", res3R)
println("res4R = \n", res4R)
println("diff = \n", res3R - res4R)
println("diffnorm = \n", vecnorm(res3R - res4R))
=#

@assert maximum(abs.(res3L - res4L)) < 1e-14
@assert maximum(abs.(res3R - res4R)) < 1e-14
