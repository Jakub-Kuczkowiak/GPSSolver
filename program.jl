setprecision(32)

data = Array{BigFloat, 2}( [234 543 213 4732; 132 756 326 5648; 321 564 435 3543; 46832 63434 42683 83429])

# data[1, :] : wektor współrzędnych x satelit
# data[2, :] : wektor współrzędnych y satelit
# data[3, :] : wektor współrzędnych z satelit
# data[4, :] : wektor czasów zegara satelit

C = 3000000

function f(i::Int64, X::Array{BigFloat})
    return (x[1] - data[1, i]) * (x[1] - data[1, i]) + 
            (x[2] - data[2, i]) * (x[2] - data[2, i]) +
            (x[3] - data[3, i]) * (x[3] - data[3, i]) - 
            C * (data[4, i] - x[4]) * (data[4, i] - x[4])
end

function F(X::Array{BigFloat})
    return [f(j, X) for j = 1:4]
end

Jacobian = Array{Function, 2}(4, 4)

for j = 1:3
    for i = 1:4 Jacobian[j, i] = function(x) return 2*x[j] - 2*data[j, i] end end
end
for i = 1:4 Jacobian[4, i] = function(x) return 2*C*data[4, i] - 2*C*x[4] end end

x = Array{BigFloat, 1}(3)
x = Array{BigFloat, 1}([1, 2, 3, 4])

@show Jacobian[1, 1](234)