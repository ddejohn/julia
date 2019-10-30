using LinearAlgebra


# decompose A into LU, with partial pivoting
function LU_decomp(A::Array{Float64,2})
    n = size(A,1)
    id = [1. for _ in 1:n]
    tx_mats = []

    pivot!(A)

    for i in 1:n
        tx = diagm(id)
        for j in i+1:n
            tx[j,i] = -A[j,i]/A[i,i]
        end
        A = tx*A
        push!(tx_mats, tx)
    end
    
    A, inv(reduce(*, reverse(tx_mats)))
end


# Gaussian elimination
function ref!(A::Array{Float64,2})
    rows = size(A,1)
    pivot!(A)

    for i in 1:rows
        for j in i+1:rows
            p = -A[j,i]/A[i,i]
            A[j,:] += p*A[i,:]
        end
        if 0 in diag(A)
            pivot!(A)
        end
    end
end


# Perform partial pivoting on matrix A
function pivot!(A::Array{Float64,2})
    for col in 1:size(A,1)-1
        swap!(A, col, (col-1) + argmax(abs.(A[col:end, col])))
    end

    for i in 2:size(A,1)
        if A[i,i] == A[i-1,i]
            swap!(A, i, i+1)
        end
    end
end


# swap rows i and j of matrix A
function swap!(A::Array{Float64,2}, i::Int, j::Int)
    for k in 1:size(A,2)
        A[i,k], A[j,k] = A[j,k], A[i,k]
    end
end


# solve a system using backward substitution
function backward_sub(A, b)
    num_rows = size(A,1)
    xx = zeros(num_rows)
    for i in num_rows:-1:1
        xi = mapreduce(j -> A[i,j]*xx[j], +, i+1:num_rows, init=0)
        xx[i] = (b[i]-xi)/A[i,i]
    end
    return xx
end


# solve a system using forward substitution
function forward_sub(A, b)
    num_rows = size(A,1)
    xx = zeros(num_rows)
    for i in 1:num_rows
        xi = mapreduce(j -> A[i,j]*xx[j], +, 1:i-1, init=0)
        xx[i] = (b[i]-xi)/A[i,i]
    end
    return xx
end


# print matrix
function print_mat(A)
    println()
    for r in 1:size(A,1)
        println(A[r,:])
    end
    println()
end


function gauss_solve(A, b)
    A = hcat(A,b)
    ref!(A)
    A, b = A[:,1:end-1], A[:, end]
    return backward_sub(A, b)
end


# decompose A into LU and solve for b
function LU_solve(A, b)
    U, L = LU_decomp(A)
    return backward_sub(U, forward_sub(L, b))
end;