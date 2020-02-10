function ols(X::AbstractMatrix, y::AbstractVector)
    X = hcat(ones(size(X,1)), X)
    w = (X'X)\(X'y)
    return x -> hcat(ones(size(x,1)), x)*w
end