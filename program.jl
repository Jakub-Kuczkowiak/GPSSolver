setprecision(32)

data = Array{BigFloat, 2}( [
    234 231 321 213; 
    123 234 352 302; 
    342 413 643 521; 
    354 432 642 523])

@show data[1, :]

# data[1, :] : wektor współrzędnych x satelit
# data[2, :] : wektor współrzędnych y satelit
# data[3, :] : wektor współrzędnych z satelit
# data[4, :] : wektor czasów zegara satelit

# data[:, i] : dane z i-tej satelity

C = 1 # dla testów

function f(i::Int64, X::Array{BigFloat})
    return (x[1] - data[1, i]) * (x[1] - data[1, i]) + 
            (x[2] - data[2, i]) * (x[2] - data[2, i]) +
            (x[3] - data[3, i]) * (x[3] - data[3, i]) - 
            C * (data[4, i] - x[4]) * (data[4, i] - x[4])
end

function F(X::Array{BigFloat})
    return [f(j, X) for j = 1:4]
end

JacobianArray = Array{Function, 2}(4, 4)

for j = 1:3
    for i = 1:4 JacobianArray[j, i] = function(x) return 2*x[j] - 2*data[j, i] end end
end
for i = 1:4 JacobianArray[4, i] = function(x) return 2*C*data[4, i] - 2*C*x[4] end end

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
    δ = - inv(Jacobian(JacobianArray, X)) * F(X)
    while norm(δ) > ϵ        
        X = X + δ
        δ = - inv(Jacobian(JacobianArray, X)) * F(X)
        @show X
    end
    return X + δ
end

x = Array{BigFloat, 1}(3)
x = Array{BigFloat, 1}([1, 2, 3, 4])

@show NewtonMethod(Array{BigFloat, 1}([1, 1, 1, 1]), BigFloat(1));