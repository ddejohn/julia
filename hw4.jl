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

    return m[1,1] + sum(m[i,i]*(x -> x-xx[j] for j in 1:n-1) for i in 2:n)
end



# data = [(1,3), (2, 3), (1,1), (0,0), (5,2)]

# xx,yy = points(data)

# print(xx)
# print(yy)

data = [(-0.1, 1.81818), (0., 2.), (0.2, 2.5), (0.3, 2.85714)]

f = ndd(data=data)
print(f(2))