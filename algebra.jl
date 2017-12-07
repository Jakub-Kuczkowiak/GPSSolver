function getSpaceFromPoints(a::Array{BigFloat}, b::Array{BigFloat}, c::Array{BigFloat}, d::Array{BigFloat})
    P = hcat(a, b, c, d)'
    v = copy(P[1, :])
    P = P .- v'
    x = P[:, 1]
    y = P[:, 2]
    z = P[:, 3]
    t = P[:, 4]
    A = y[2] * z[3] * t[4] + z[2] * t[3] * y[4] + t[2] * y[3] * z[4] - y[2] * z[4] * t[3] - z[2] * y[3] * t[4] - t[2] * z[3] * y[4]
    B = x[2] * t[3] * z[4] + z[2] * x[3] * t[4] + t[2] * z[3] * x[4] - x[2] * t[4] * z[3] - z[2] * t[3] * x[4] - t[2] * x[3] * z[4]
    C = x[2] * y[3] * t[4] + y[2] * t[3] * x[4] + t[2] * x[3] * y[4] - x[2] * t[3] * y[4] - y[2] * t[4] * x[3] - t[2] * x[4] * y[3]
    D = x[2] * z[3] * y[4] + y[2] * x[3] * z[4] + z[2] * y[3] * x[4] - x[2] * z[4] * y[3] - y[2] * x[4] * z[3] - z[2] * y[4] * x[3]
    E = - v[1] * A - v[2] * B - v[3] * C - v[4] * D
    return [A, B, C, D, E]
end

function reflect(x::Array{BigFloat}, s::Array{BigFloat})
    u = - sum(s[1:end-1] .* x + s[end]) / sum( s[1:end-1] .* s[1:end-1])
    return (x + 2 * u * s[1:end-1])
end

# @show getSpaceFromPoints(
#    Array{BigFloat}([1, 2, 3, 4]),
#    Array{BigFloat}([2, 4, 4, 5]), 
#    Array{BigFloat}([3, 4, 5, 6]), 
#    Array{BigFloat}([5, 6, 4, 3]))