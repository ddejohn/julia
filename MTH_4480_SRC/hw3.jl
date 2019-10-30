function newton(; f, df, xn, max_iter=100, tol=1e-4)
    n=1
    # Print initial guess xn and f(xn)
    # println("x0 = $(xn)\nf(x0) = $(f(xn))\n")

    while n <= max_iter
        # Newton's method: set current x to previous x minus the ratio 
        # of f and f' evaluated at the previous x
        x = xn - f(xn)/df(xn)

        # Print x and f(x) on each iteration
        # println("x$(n) = $(x)\nf(x$(n)) = $(f(x))\n")
        
        # If f(x) is within 'tol' of zero, we've found our root
        if abs(x-xn) < tol || abs(f(x)) < tol
            xn = x
            break
        end
        
        # Else set xn to current x and keep going
        n += 1
        xn = x
    end

    # If we've performed the maximum number of iterations
    # AND we're still not within the specified tolerance, we failed
    if (n == max_iter) & (abs(f(xn)) > tol)
        out = "failed to find root after $(max_iter) iterations.\n\n"
    else
        out = "root found near x = $(xn) in $(n) iterations.\n\n"
    end

    println(out)
end


function secant(; f, x0, x1, max_iter=100, tol=1e-4)
    n=1

    # Print initial guess xn and f(xn)
    # println("x0 = $(x0)\nx1 = $(x1)\n")

    while n <= max_iter
        # Secant method: set current x to previous x minus the ratio 
        # of f and f' approximated with the previous *two* x-values
        fx0, fx1 = f(x0), f(x1)
        x = x1 - fx1*(x1-x0)/(fx1-fx0)

        # If f(x) is within 'tol' of zero, we've found our root
        if abs(x-x1) < tol || abs(f(x)) < tol
            x1 = x
            break
        end
        
        # Else set xn to current x and keep going
        n += 1
        x0 = x1
        x1 = x
    end

    # If we've performed the maximum number of iterations
    # AND we're still not within the specified tolerance, we failed
    if (n == max_iter) & (abs(f(x1)) > tol)
        out = "failed to find root after $(max_iter) iterations.\n\n"
    else
        out = "root found near x = $(x1) in $(n) iterations.\n\n"
    end

    println(out)
end


# The functions and their derivatives for exercise 1a-d
# This is a list of lists, where each inner list
# contains f(x) and its derivative, both expressed as
# anonymous functions (like lambdas from Python)
funcs = [
    [x -> x^3 - 2x^2 - 5,        x -> 3x^2 - 4x],
    [x -> x^3 + 3x^2 - 1,        x -> 3x^2 + 6x],
    [x -> x - cos(x),            x -> 1 + sin(x)],
    [x -> x - 0.8 - 0.2sin(x),   x -> 1 - 0.2cos(x)]
]

# Initial guesses
xis = [
    1, -2.5, 1, 1
]

for (func, x0) in zip(funcs, xis)
    f, df = func                    # unzip funcs into f(x) and f'(x)
    newton(f=f, df=df, xn=x0)       # call Newton's method with initial guess x0
end


# Functions for secant method expressed as anonymous functions
funcs = [
    x -> x^3 - 2x^2 - 5,
    x -> x^3 + 3x^2 - 1,
    x -> x - cos(x),
    x -> x - 0.8 - 0.2sin(x)
]

# Initial guesses
xis = [
    [1, 1.1],
    [-2.5, -2.4],
    [1, 0.9],
    [1, 0.9]
]

for (func, xi) in zip(funcs, xis)
    x1, x0 = xi                     # unzip xns into x1 and x0
    secant(f=func, x1=x1, x0=x0)    # call secant method with initial guess xn
end;

f = x -> x^3 - 2x^2 - 5
df = x -> 3x^2 - 4x

x0 = 1.0
x1 = 1.01

newton(f=f, df=df, xn=x0)
secant(f=f, x0=x0, x1=x1)