setprecision(64)

include("algebra.jl")

data = Array{BigFloat, 2}( [
    234 231 321 213; 
    123 234 352 302; 
    342 413 643 521; 
    354 432 642 523])

# data[1, :] : wektor współrzędnych x satelit
# data[2, :] : wektor współrzędnych y satelit
# data[3, :] : wektor współrzędnych z satelit
# data[4, :] : wektor czasów zegara satelit

# data[:, i] : dane z i-tej satelity

C = 1 # dla testów

solution₁ = Array{BigFloat}([BigFloat("3.406027414052647875887849459335517363700416540639706155009830367636869268148760385343e+02"),
                            BigFloat("2.245621396967134977305526534214474645554267237451765234641071107647602855199079076453e+02"),
                            BigFloat("8.802209804658630277006504596604842085109012249373217265353711792757437615214584985008e+02"),
                            BigFloat("9.119971653230215237445022093317229764207208762060580107951832552971171535934449105973e+02")])

function f(x::Array{BigFloat})
    return ((-data .+ x) .^ 2)' * [1, 1, 1, -C]
end

JacobianArray = Array{Function, 2}(4, 4)

for i = 1:4
    for j = 1:3
        JacobianArray[i, j] = function(x) return 2*x[j] - 2*data[j, i] end
    end
    JacobianArray[i, 4] = function(x) return 2*C*data[4, i] - 2*C*x[4] end
end

function Jacobian(JacobianArray::Array{Function, 2}, X::Array{BigFloat})
    sol = Array{BigFloat, 2}(4, 4)
    for i = 1:4
        for j = 1:4
            sol[i, j] = JacobianArray[i, j](X)
        end
    end
    return sol
end

function NewtonMethod(X::Array{BigFloat}, ϵ::BigFloat)
    iteration_count = 1
    while true 
        @show X - solution₁
        δ = Jacobian(JacobianArray, X) \ f(X)
        X = X - δ
        if norm(δ) < ϵ
            return X
        end
    end
end

x = Array{BigFloat, 1}(3)
x = Array{BigFloat, 1}([1, 2, 3, 4])
# [223.905, 61.0486, 413.446, 258.905]
result = NewtonMethod(Array{BigFloat, 1}([1, 1, 1, 0]), BigFloat("0.0000000000001"))

#space = getSpaceFromPoints(data[1,:], data[2,:], data[3,:], data[4,:])
#result2 = reflect(result, space)
#@show f(result)
#@show f(result2)