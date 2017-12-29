#module DualMod

import Base: +, -, *, /
import Base: >, >=, <, <=, ==, !=
import Base: sin, cos, exp, sqrt
import Base: abs, dot, zero
import Base: convert

export AbstDual, Dual
export >, >=, <, <=, ==, !=, +, -, *, /
export sin, cos, exp, sqrt, abs
export dot, initFrom

abstract type AbstDual <: Number end

struct Dual{DIM} <: AbstDual
    v::Float64
    d::Array{Float64,1}

    function Dual{DIM}() where {DIM}
        v = 0.0
        d = zeros(Float64, DIM)
        new(v, d)
    end
    function Dual{DIM}(v::Real) where {DIM}
        d = zeros(Float64, DIM)
        new(v, d)
    end
    function Dual{DIM}(v::Float64, 
                       d::AbstractArray{Float64, 1}) where {DIM}
        new(v, d)
    end
    function Dual{DIM}(v::Float64, 
                       i::Integer) where {DIM}
        d = zeros(Float64, DIM)
        d[i] = 1.0
        new(v, d)
    end
end


function convert{DIM}(::Type{Dual{DIM}}, v::Real)
    return Dual{DIM}(v)
end
function >{DIM}(d1::Dual{DIM}, d2::Dual{DIM})
    return d1.v > d2.v
end
function >{DIM}(d1::Dual{DIM}, s::Real)
    return d1.v > s
end
function >{DIM}(s::Real, d1::Dual{DIM})
    return s > d1.v
end

function >={DIM}(d1::Dual{DIM}, d2::Dual{DIM})
    return d1.v >= d2.v
end
function >={DIM}(d1::Dual{DIM}, s::Real)
    return d1.v >= s
end
function >={DIM}(s::Real, d1::Dual{DIM})
    return s >= d1.v
end

function <{DIM}(d1::Dual{DIM}, d2::Dual{DIM})
    return d1.v < d2.v
end
function <{DIM}(d1::Dual{DIM}, s::Real)
    return d1.v < s
end
function <{DIM}(s::Real, d1::Dual{DIM})
    return s < d1.v
end

function <={DIM}(d1::Dual{DIM}, d2::Dual{DIM})
    return d1.v <= d2.v
end
function <={DIM}(d1::Dual{DIM}, s::Real)
    return d1.v <= s
end
function <={DIM}(s::Real, d1::Dual{DIM})
    return s <= d1.v
end

function !={DIM}(d1::Dual{DIM}, d2::Dual{DIM})
    return d1.v != d2.v
end
function !={DIM}(d1::Dual{DIM}, s::Real)
    return d1.v != s
end
function !={DIM}(s::Real, d1::Dual{DIM})
    return s != d1.v
end

function =={DIM}(d1::Dual{DIM}, d2::Dual{DIM})
    return d1.v == d2.v
end
function =={DIM}(d1::Dual{DIM}, s::Real)
    return d1.v == s
end
function =={DIM}(s::Real, d1::Dual{DIM})
    return s == d1.v
end

function +{DIM}(d1::Dual{DIM}, d2::Dual{DIM})
    v = d1.v + d2.v
    d = Array{Float64}(DIM)
    for i = 1 : DIM
        d[i] = d1.d[i] + d2.d[i]
    end
    return Dual{DIM}(v, d)
end
function +{DIM}(d1::Dual{DIM}, s::Real)
    v = d1.v + s
    d = zeros(Float64, DIM)
    for i = 1 : DIM
        d[i] = d1.d[i]
    end
    return Dual{DIM}(v, d)
end
function +{DIM}(s::Real, d1::Dual{DIM})
    v = d1.v + s
    d = zeros(Float64, DIM)
    for i = 1 : DIM
        d[i] = d1.d[i]
    end
    return Dual{DIM}(v, d)
end

function -{DIM}(s::Real, d1::Dual{DIM})
    v = s - d1.v
    d = zeros(Float64, DIM)
    for i = 1 : DIM
        d[i] = -d1.d[i]
    end
    return Dual{DIM}(v, d)
end

function -{DIM}(d1::Dual{DIM}, s::Real)
    v = d1.v - s
    d = zeros(Float64, DIM)
    for i = 1 : DIM
        d[i] = d1.d[i]
    end
    return Dual{DIM}(v, d)
end

function -{DIM}(d1::Dual{DIM}, d2::Dual{DIM})
    v = d1.v - d2.v
    d = Array{Float64}(DIM)
    for i = 1 : DIM
        d[i] = d1.d[i] - d2.d[i]
    end
    return Dual{DIM}(v, d)
end

function *{DIM}(d1::Dual{DIM}, d2::Dual{DIM})
    v = d1.v * d2.v

    d = Array{Float64}(DIM)
    for i = 1 : DIM
        d[i] = d1.d[i]*d2.v + d1.v*d2.d[i]
    end
    return Dual{DIM}(v, d)
end

function *{DIM}(d1::Dual{DIM}, s::Real)
    v = d1.v * s

    d = Array{Float64}(DIM)
    for i = 1 : DIM
        d[i] = d1.d[i]*s
    end
    return Dual{DIM}(v, d)
end

function *{DIM}(s::Real, d1::Dual{DIM})
    v = d1.v * s

    d = Array{Float64}(DIM)
    for i = 1 : DIM
        d[i] = d1.d[i]*s
    end
    return Dual{DIM}(v, d)
end
function /{DIM}(d1::Dual{DIM}, d2::Dual{DIM})
    v = d1.v/d2.v
    d = Array{Float64}(DIM)
    v2 = 1./(d2.v * d2.v)
    for i = 1 : DIM
        d[i] = (d1.d[i] * d2.v - d1.v * d2.d[i]) * v2
    end
    return Dual{DIM}(v, d)
end

function /{DIM}(d1::Dual{DIM}, s::Real)
    v = d1.v/s
    d = Array{Float64}(DIM)
    for i = 1 : DIM
        d[i] = d1.d[i] / s
    end
    return Dual{DIM}(v, d)
end

function /{DIM}(s::Real, d1::Dual{DIM})
    v = s / d1.v
    d = Array{Float64}(DIM)
    v2 = -1./(d1.v * d1.v)
    for i = 1 : DIM
        d[i] = d1.d[i] * v2
    end
    return Dual{DIM}(v, d)
end

function sin{DIM}(d1::Dual{DIM})
    d = Array{Float64}(DIM)
    v = sin(d1.v)
    dsin = cos(d1.v)
    for i = 1 : DIM
        d[i] = d1.d[i] * dsin
    end
    return Dual{DIM}(v, d)
end

function cos{DIM}(d1::Dual{DIM})
    d = Array{Float64}(DIM)
    v = cos(d1.v)
    dcos = -sin(d1.v)
    for i = 1 : DIM
        d[i] = d1.d[i] * dcos
    end
    return Dual{DIM}(v, d)
end

function exp{DIM}(d1::Dual{DIM})
    d = Array{Float64}(DIM)
    v = exp(d1.v)
    dexp = v
    for i = 1 : DIM
        d[i] = d1.d[i] * dexp
    end
    return Dual{DIM}(v, d)
end

function sqrt{DIM}(d1::Dual{DIM})
    # @assert(d1.v >= 0.0)
    d = Array{Float64}(DIM)
    v = sqrt(d1.v)
    dsqrt = 0.5/v
    for i = 1 : DIM
        d[i] = d1.d[i] * dsqrt
    end
    return Dual{DIM}(v, d)
end

function zero{DIM}(::Type{Dual{DIM}})
    return Dual{DIM}()
end

function abs{DIM}(d1::Dual{DIM})
    if d1.v >= 0.0
        return d1
    else
        return 0.0 - d1
    end
end

function dot{DIM}(v1::AbstractArray{Dual{DIM}, 1},
                  v2::AbstractArray{Dual{DIM}, 1})
    val = Dual{DIM}()
    len = length(v1)
    for i = 1 : len
        val += v1[i] * v2[i]
    end
    return val
end

function dot{DIM, T}(v1::AbstractArray{Dual{DIM}, 1},
                     v2::AbstractArray{T, 1})
    val = Dual{DIM}()
    len = length(v1)
    for i = 1 : len
        val += v1[i] * v2[i]
    end
    return val
end

function dot{DIM, T}(v1::AbstractArray{T, 1},
                     v2::AbstractArray{Dual{DIM}, 1})
    val = Dual{DIM}()
    len = length(v1)
    for i = 1 : len
        val += v1[i] * v2[i]
    end
    return val
end


#
# This function helps unify the APIs in evaluation of residual.
# 
#
function initFrom{Tval<:Number}(::Type{Tval},
                                arr::AbstractArray{Float64, 1})
    return arr
end
function initFrom{Tval<:Number}(::Type{Tval},
                                arr::AbstractArray{Float64, 2})
    return arr
end
function initFrom{Tval<:AbstDual}(::Type{Tval},
                                  arr::AbstractArray{Float64, 1})
    DIM = UInt8(length(arr))
    darr = Array{Dual{DIM}}(DIM)
    for i = 1 : DIM
        darr[i] = Dual{DIM}(arr[i], i)
    end
    return darr
end
function initFrom{Tval<:AbstDual}(::Type{Tval},
                                  arr::AbstractArray{Float64, 2})
    dim1 = size(arr, 1)
    dim2 = size(arr, 2)
    DIM = UInt8(dim1 * dim2)
    darr = Array{Dual{DIM}}(dim1, dim2)
    idx = 1
    for j = 1 : dim2
        for i = 1 : dim1
            darr[i,j] = Dual{DIM}(arr[i, j], idx)
            idx += 1
        end
    end
    return darr
end

function initFrom{Tval<:Number}(::Type{Tval},
                  arr1::AbstractArray{Float64, 2},
                  arr2::AbstractArray{Float64, 2})
    return arr1, arr2
end

function initFrom{Tval<:AbstDual}(::Type{Tval}, 
                                  arr1::AbstractArray{Float64, 2},
                                  arr2::AbstractArray{Float64, 2})
    dim1 = size(arr1, 1)
    dim2 = size(arr1, 2)
    # DIM = 2 * dim1 * dim2
    DIM = UInt8(2*dim1 * dim2)
    darr1 = Array{Dual{DIM}}(dim1, dim2)
    darr2 = Array{Dual{DIM}}(dim1, dim2)
    idx = 1
    offset = dim1 * dim2
    for j = 1 : dim2
        for i = 1 : dim1
            darr1[i,j] = Dual{DIM}(arr1[i, j], idx)
            darr2[i,j] = Dual{DIM}(arr2[i, j], idx + offset)
            idx += 1
        end
    end
    return darr1, darr2
end
#end # end of mudule DualMod
