setprecision(256)

# geneartor testów
function gen_test(satellites)
    # wektor [x; y; z; t] - (x, y, z): nasza pozycja, t: błąd zegara
    X = vcat(normalize(Array{BigFloat}(2 * rand(3) - 1) .^ 3), Array{BigFloat}(randn(1)) .^ 3)
    data = Array{BigFloat, 2}(4, satellites)
    for i = 1:satellites
        v = normalize(Array{BigFloat}(2 * rand(3) - 1) .^ 3)
        data[:, i] = X + vcat(v, norm(v) / C)
    end
    return data, X
end

n = 8 # liczba satelit

data, solution₁ = gen_test(n)

# data[1, :] : wektor współrzędnych x satelitów
# data[2, :] : wektor współrzędnych y satelitów
# data[3, :] : wektor współrzędnych z satelitów
# data[4, :] : wektor czasów zegara satelitów

# data[:, i] : dane z i-tej satelity
# solution₁ = poszukiwane rozwiązanie (niekoniecznie jedyne)

# wyświetlanie wyników
function prnt(x::Array{BigFloat}, fx::Array{BigFloat})
    @printf "x - sol: [%.4e, %.4e, %.4e, %.4e]\tf(x): [%.4e, %.4e, %.4e, %.4e]\n" x[1] x[2] x[3] x[4] fx[1] fx[2] fx[3] fx[4]
end

# Funkcja F jak w sprawozdaniu
# rozważamy dane tylko z 4 pierwszych satelitów
function f4(x::Array{BigFloat})
    return ((-data[:, 1:4] .+ x) .^ 2)' * [1, 1, 1, -C*C]
end


# Funkcja F jak w sprawozdaniu
# bierze pod uwagę dane ze wszystkich satelitów
function f(x::Array{BigFloat})
    return ((-data .+ x) .^ 2)' * [1, 1, 1, -C*C]
end

# Jakobian jak w sprawozdaniu
# dla 4 pierwszych satelitów
∂f4arr = Array{Function, 2}(4, 4)
for i = 1:4
    for j = 1:3
        ∂f4arr[i, j] = function(x) return 2*x[j] - 2*data[j, i] end
    end
    ∂f4arr[i, 4] = function(x) return 2*C*C*data[4, i] - 2*C*C*x[4] end
end

function ∂f4(x::Array{BigFloat})
    return map((f) -> f(x), ∂f4arr)
end


# Jakobian jak w sprawozdaniu
# dla wszystkich satelitów
∂farr = Array{Function, 2}(n, 4)
for i = 1:n
    for j = 1:3
        ∂farr[i, j] = function(x) return 2*x[j] - 2*data[j, i] end
    end
    ∂farr[i, 4] = function(x) return 2*C*C*data[4, i] - 2*C*C*x[4] end
end

function ∂f(x::Array{BigFloat})
    return map((f) -> f(x), ∂farr)
end

# Metoda Newtona - Rhapsona, klasyczna: dla 4 satelitów
function NewtonMethod(f, ∂f, ϵ; x=Array{BigFloat}([0, 0, 0, 0]), max_iter=1000, log=true, solution=solution₁)
    i = 0
    while norm(f(x)) > ϵ && i < max_iter
        δ = ∂f(x) \ f(x)
        x = x - δ
        i += 1
        if log prnt(abs.(x - solution), abs.(f(x))) end
    end
    return x
end

# ważona metoda Newtona - Rhapsona: dla dowolnej liczby satelitów
function NewtonMethod(f, ∂f, ϵ; x=Array{BigFloat}([0, 0, 0, 0]), max_iter=1000, log=true)
    i = 0
    while norm(f(x)) > ϵ && i < max_iter
        δ = (∂f(x)' * ∂f(x)) \ ∂f(x)' * f(x)
        x = x - δ
        i += 1
        if log prnt(abs.(x - solution₁), abs.(f(x))) end
    end
    return x
end

# Metoda najmniejszych kwadratów
function S(x::Array{BigFloat})
    return [sum(4 * f(x) .* (x[i] - data[i, :])) for i = 1:4] .* [1, 1, 1, -C*C]
end

∂Sarr = Array{Function, 2}(4, 4)
for i = 1:4
    for j = 1:4
        ∂Sarr[i, j] = (i == 4 || j == 4) ? 
            function(x) return - C * C * 8 * sum((x[j] - data[j, :]) .* (x[i] - data[i, :])) end :
            function(x) return 8 * sum((x[j] - data[j, :]) .* (x[i] - data[i, :])) end

    end
end

for i = 1:3
    ∂Sarr[i, i] = function(x) return 4 * sum(f(x) + (2 * (data[i, :] - x[i]) .* (data[i, :] - x[i]))) end
end

∂Sarr[4, 4] = function(x) return - 4 * C * C * sum(f(x) - (2 * C * C * (data[4, :] - x[4]) .* (data[4, :] - x[4]))) end

# Jakobian funkcji błędu S
function ∂S(x::Array{BigFloat})
    return map((f) -> f(x), ∂Sarr)
end