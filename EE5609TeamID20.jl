const inv_map::Vector{Int8} = [1, 3, 2]

@inline function -(x::Int8)
    return x
end

@inline function +(x::Vector{Int8}, y::Vector{Int8})
    @inbounds @simd for i in 1:length(y)
        x[i] ⊻= y[i]
    end
    return x
end

@inline function -(x::Int8, y::Int8)
    return x ⊻ y
end

@inline function *(x::Int8, y::Int8)
    return ifelse(x==0 || y==0 , convert(Int8, 0), ifelse(x==1, y, ifelse(y==1, x, ifelse(x!=y, convert(Int8, 1), ifelse(x == 2, convert(Int8, 3), convert(Int8, 2))))))
end

@inline function inv(x::Int8)
    return inv_map[x]
end

@inline function +(x::Int64, y::Int64)
    return Base.:+(x, y)
end

@inline function -(x::Vector{Int8}, y::Vector{Int8})
    @inbounds @simd for i in 1:length(x)
        x[i] ⊻= y[i]
    end
    return x
end

@inline function *(x::Vector{Int8}, y::Vector{Int8})
    @inbounds @simd for i in 1:length(x)
        x[i] *= y[i]
    end
    return x
end

@inline function *(x::Vector{Int8}, y::Int8)
    @inbounds @simd for i in 1:length(x)
        x[i] *= y
    end
    return x
end


function rankconsistencyTeamID20(A, b)
    A = convert(Matrix{Int8}, A)
    b = convert(Vector{Int8}, b)
    A = hcat(A, b)
    m, n = size(A)
    rank = 0
    @inbounds begin
        for row = 1:Base.:-(m, 1)
            col = min(row, n)
            maxel, maxidx = findmax(A[row:end, col])
            if iszero(maxel)
                for i in Base.:+(col, 1):n
                    maxel, maxidx = findmax(A[row:end, i])
                    if !iszero(maxel)
                        col = i
                        break
                    end
                end
            end 

            if iszero(maxel)
                # All elements are zeros (no pivot); echelon form obtained
                # If prev row has an element only in b, the system is inconsistent
                if A[rank, Base.:-(n, 1)] == 0 && A[rank, n] != 0
                    return convert(Matrix{Int64}, A), convert(Int64, rank), false
                end
                return convert(Matrix{Int64}, A), convert(Int64, rank), true
            end
            
            # Pivot element found, increment rank
            rank = Base.:+(rank, 1)

            # index of row in matrix (findmax gives index of row in matrix slice, not in whole matrix)
            maxidx = Base.:-(Base.:+(maxidx, row), 1)

            if maxidx != row
                # Interchange
                @simd for i = 1:n
                    tmp = A[row,i]
                    A[row,i] = A[maxidx,i]
                    A[maxidx,i] = tmp
                end
            end

            # Update all rows below current row
            inverse = inv(A[row, col])
            @simd for i = Base.:+(row, 1):m
                A[i, :] = A[i, :] - A[row, :]*(A[i, col] * inverse)
            end
        end
    end
    # If last row is nnz
    rank_aug = rank
    if !iszero(findmax(A[m, 1:n])[1])
        rank_aug = Base.:+(rank, 1) 
    end
    if !iszero(findmax(A[m, 1:Base.:-(n, 1)])[1])
        rank = Base.:+(rank, 1) 
    end
    
    return convert(Matrix{Int64}, A[:, 1:Base.:-(n, 1)]), convert(Int64, rank), rank >= rank_aug
end

