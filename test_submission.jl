include("EE5609TeamID20.jl") # including the file from Team number 07
m = 10 # must work for any m <= 500
n = 5 # must work for any n <= 500
A = rand(0:3,(m,n)) # this is a Matrix{Int64} array with entries 0,1,2,3
b = rand(0:3,m) # b::Vector{Int64}

while true
    A = rand(0:3,(m,n)) # this is a Matrix{Int64} array with entries 0,1,2,3
    b = rand(0:3,m) # b::Vector{Int64}
    U, r, consistent = rankconsistencyTeamID20(A,b)
    if consistent
        println(typeof(U))
        println(typeof(r))
        println(consistent)
        println(U)
        println(A)
        println(b)
        break
    end
end