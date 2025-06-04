using QuadGK, Roots, Plots

# Same ItemParams and functions as before...
struct ItemParams
    η::Float64
    α::Float64
    n::Float64
    h::Float64
    ω::Float64
    p::Float64
    c::Float64
    v::Float64
end

T0 = 1.0
A = 0.0

items = [
    ItemParams(10.0, 3.0, 1.5, 2.0, 5.0, 20.0, 10.0, 1.0),
    ItemParams(15.0, 2.5, 0.8, 1.5, 4.0, 25.0, 12.0, 1.2)
]

μ(i::ItemParams) = i.α * i.η / (i.α - 1)

function S0(i::ItemParams)
    β = i.h / (i.h + i.ω)
    thresh = i.n / (i.α + i.n)
    if thresh <= β
        return i.η * ((i.α + i.n) * i.ω / (i.α * (i.h + i.ω)))^(1 / i.n)
    else
        return i.η * (i.n * (i.h + i.ω) / (i.h * (i.α + i.n)))^(1 / i.α)
    end
end

function Si_lambda(i::ItemParams, λ)
    if λ > i.ω / i.v
        return 0.0
    end
    rhs = (i.h + λ * i.v) / (i.h + i.ω)
    thresh = i.n / (i.α + i.n)
    if thresh <= rhs
        return i.η * ((i.α + i.n) * (i.ω - λ * i.v) / (i.α * (i.h + i.ω)))^(1 / i.n)
    else
        return i.η * (i.n * (i.h + i.ω) / ((i.h + λ * i.v) * (i.α + i.n)))^(1 / i.α)
    end
end

function g(λ, W)
    sum(Si_lambda(i, λ) * i.v for i in items) - W
end

function inventory_cost(i::ItemParams, S::Float64, λ::Float64)
    μi = μ(i)
    if S <= i.η
        J = 0.0
    else
        J = i.η^i.α / (i.α - 1) * S^(-i.α + 1) + S - i.α * i.η / (i.α - 1)
    end
    C = (i.h + i.ω) * i.n / (i.n + 1) * J +
        i.ω * i.n / (i.n + 1) * (μi - S) -
        λ * i.v * S / (i.n + 1)
    return C
end

function solve_for_W(W)
    S0s = [S0(i) for i in items]
    total_unconstrained = sum(i.v * S for (i, S) in zip(items, S0s))

    if total_unconstrained <= W
        λ = 0.0
        Sstars = S0s
    else
        λ_upper = maximum(i.ω / i.v for i in items)
        λ = find_zero(λ -> g(λ, W), (0.0, λ_upper), Bisection())
        Sstars = [Si_lambda(i, λ) for i in items]
    end

    total_cost = sum(inventory_cost(i, S, λ) for (i, S) in zip(items, Sstars))
    total_profit = sum((i.p - i.c) * μ(i) / T0 for i in items) - total_cost

    return λ, Sstars, total_cost, total_profit
end

# === Sweep W and Plot === #
W_vals = range(10.0, 100.0, length=100)
profits = Float64[]
costs = Float64[]
S_plot = [[] for _ in 1:length(items)]

for W in W_vals
    λ, Sstars, cost, profit = solve_for_W(W)
    push!(profits, profit)
    push!(costs, cost)
    for (k, S) in enumerate(Sstars)
        push!(S_plot[k], S)
    end
end

# === Plot === #
plot(W_vals, profits, lw=2, xlabel="Warehouse Capacity W", ylabel="Expected Profit", title="Profit vs Storage Capacity", legend=false)

plot_S = plot(title="Optimal Inventory Level Sᵢ vs W", xlabel="W", ylabel="Sᵢ")
for (i, Si_vals) in enumerate(S_plot)
    plot!(plot_S, W_vals, Si_vals, lw=2, label="Item $i")
end

plot(plot_S, layout=(1, 1))

