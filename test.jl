using LinearAlgebra
using JLD2
using InteractiveUtils

@load "Samples_Prog_1.jld2"

function lafact(A, b)
    A = hcat(A, b)
    l, u, p = LinearAlgebra.lu(A, check=false)
end

@time lafact(rand(3, 3), rand(3, 1))

for sample in Samples
    A = sample.A
    b = sample.b
    @time lafact(A, b)
end

include("EE5609TeamID20.jl")

for sample in Samples
    A = sample.A
    b = sample.b
    # @time rankconsistency(A, b)
    u, r, c = rankconsistencyTeamID20(A, b)
    u = convert(Matrix{Int64}, u)
    r = convert(Int64, r)
    println(u == sample.U, ' ', r == sample.r, ' ', c == sample.consistent)
end





