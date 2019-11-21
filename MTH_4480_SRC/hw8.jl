using LinearAlgebra



# decompose A into LU, with partial pivoting
function LU_decomp(A::Array{Float64,2})
    n = size(A,1)

    # identity matrix
    id = [1. for _ in 1:n]

    # list for storing transformation matrices
    tx_mats = []

    # perform partial pivoting on A first
    pivot!(A)

    for i in 1:n
        # initialize the transformation matrix
        tx = diagm(id)
        
        # for each row below a pivot, record the "pivoting scalar" and place
        # it in the corresponding position in the transformation matrix
        for j in i+1:n
            tx[j,i] = -A[j,i]/A[i,i]
        end
        
        # perform the transformation
        A = tx*A
        
        # record the transformation matrix
        push!(tx_mats, tx)
    end
    
    # return U and L where L is the inverse of
    # the product of the transformation matrices
    A, inv(reduce(*, reverse(tx_mats)))
end


# Gaussian elimination
function ref!(A::Array{Float64,2})
    rows = size(A,1)
    pivot!(A)

    # zero out the entries in the rows below each pivot
    for i in 1:rows
        for j in i+1:rows
            # "pivoting scalar" which will cancel the entries below the pivot
            p = -A[j,i]/A[i,i]
            A[j,:] += p*A[i,:]
        end
    end
end


# Perform partial pivoting on matrix A
function pivot!(A::Array{Float64,2})
    # search for the largest absolute value in each pivot column
    # row swap so that this value is in the pivot position
    for col in 1:size(A,1)-1
        swap!(A, col, (col-1) + argmax(abs.(A[col:end, col])))
    end

    # check if row reduction will result in
    # a zero on the diagonal, repivot if needed
    for i in 2:size(A,1)
        if A[i,i] == A[i-1,i]
            swap!(A, i, i+1)
        end
    end
end


# swap rows i and j of matrix A
function swap!(A::Array{Float64,2}, i::Int, j::Int)
    # swap two rows element-wise (faster than slicing)
    for k in 1:size(A,2)
        A[i,k], A[j,k] = A[j,k], A[i,k]
    end
end


# solve a system using backward substitution
function backward_sub(A, b)
    num_rows = size(A,1)
    
    # initializes an empty solution vector
    xx = zeros(num_rows)
    for i in num_rows:-1:1
        # solves for each variable using the
        # previously solved variables on which it depends
        xi = mapreduce(j -> A[i,j]*xx[j], +, i+1:num_rows, init=0)
        xx[i] = (b[i]-xi)/A[i,i]
    end
    return xx
end


# solve a system using forward substitution
function forward_sub(A, b)
    num_rows = size(A,1)
    
    # initializes an empty solution vector
    xx = zeros(num_rows)
    for i in 1:num_rows
        # solves for each variable using the
        # previously solved variables on which it depends
        xi = mapreduce(j -> A[i,j]*xx[j], +, 1:i-1, init=0)
        xx[i] = (b[i]-xi)/A[i,i]
    end
    return xx
end


# Perform gaussian elimination and solve for b using backward substitution
function gauss_solve(A, b)
    # augment A with b
    A = hcat(A,b)
    
    # reduce to echelon form
    ref!(A)
    
    # split A and b again
    A, b = A[:,1:end-1], A[:, end]

    return backward_sub(A, b)
end


# compare gaussian elimination vs LU decomposition
# vs Julia's built-in matrix methods
function compare(A, b)
    @time LU = LU_solve(A, b)
    @time gauss = gauss_solve(A, b)
    @time built_in = A\b

    println("\nLU:       ", LU)
    println("gauss:    ", gauss)
    println("built-in: ", built_in)
    println()
end


# decompose A into LU and solve for b
function LU_solve(A, b)
    U, L = LU_decomp(A)
    return backward_sub(U, forward_sub(L, b))
end;