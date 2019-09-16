function ndd(; data)
    n = length(data)
    # initialize nxm matrix
    m = zeros(Float64, n, n)
    
    # unzips list of (x,y) coords into separate lists
    # set first column of DD matrix to f(x) for x in xx
    xx, m[:,1] = map(y -> map(x -> x[y], data), 1:2)

    # generate lower triangular NDD matrix
    for i in 2:n
        for j in 2:i
            m[i,j] = (m[i,j-1]-m[i-1,j-1])/(xx[i]-xx[i-j+1])
        end
    end

    # generate Newton polynomials
    function poly(i, x)
        p = 1
        for j in 1:i-1
            p *= (x-xx[j]) # (x - xj) for j from 1 to i-1
        end
        return p
    end

    # return a lambda in x: 
    # sum the diagonal entries as coefficients
    # on the generated Newton polynomial
    return x-> sum(m[i,i]*poly(i, x) for i in 1:n)
end


data = [(-0.1, 1.81818), (0., 2.), (0.2, 2.5), (0.3, 2.85714), (0.35, 3.07692)]

f = ndd(data=data)
println(f(0.4))