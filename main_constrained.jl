using QuadGK
using Roots
import Pkg
Pkg.add("Plots")

# ========== Inputs ========== #
struct ItemParams
    η::Float64
    α::Float64
    n::Float64
    h::Float64
    ω::Float64
    p::Float64
    c::Float64
    v::Float64  # storage per unit
end

T0 = 1.0
A = 0.0
W = 50.0  # warehouse capacity

items = [
    ItemParams(10.0, 3.0, 1.5, 2.0, 5.0, 20.0, 10.0, 1.0),
    ItemParams(15.0, 2.5, 0.8, 1.5, 4.0, 25.0, 12.0, 1.2)
]

# ========== Demand Mean ========== #
μ(i::ItemParams) = i.α * i.η / (i.α - 1)

# ========== Unconstrained S₀ (Eq. 38) ========== #
function S0(i::ItemParams)
    β = i.h / (i.h + i.ω)
    thresh = i.n / (i.α + i.n)
    if thresh <= β
        return i.η * ((i.α + i.n) * i.ω / (i.α * (i.h + i.ω)))^(1 / i.n)
    else
        return i.η * (i.n * (i.h + i.ω) / (i.h * (i.α + i.n)))^(1 / i.α)
    end
end

# ========== Constrained S(λ) (Eq. 42) ========== #
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

# ========== g(λ) Function (Eq. 30) ========== #
function g(λ)
    sum(Si_lambda(i, λ) * i.v for i in items) - W
end

# ========== Expected Inventory Cost (Eq. 33) ========== #
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

# ========== Final Solver ========== #
function solve_with_constraints()
    # Step 1: try unconstrained
    S0s = [S0(i) for i in items]
    total_space = sum(i.v * S for (i, S) in zip(items, S0s))
    if total_space <= W
        println("Unconstrained optimal solution is feasible.")
        λstar = 0.0
        Sstars = S0s
    else
        # Step 2: root-finding
        λ_upper = maximum(i.ω / i.v for i in items)
        λstar = find_zero(g, (0.0, λ_upper), Bisection())
        Sstars = [Si_lambda(i, λstar) for i in items]
    end

    println("\n=== Optimal Results ===")
    total_cost = 0.0
    total_profit = 0.0
    for (j, i) in enumerate(items)
        S = Sstars[j]
        cost = inventory_cost(i, S, λstar)
        profit = (i.p - i.c) * μ(i) / T0 - cost
        println("Item $(j):")
        println("  S*: ", round(S, digits=4))
        println("  Cost: ", round(cost, digits=4))
        println("  Profit: ", round(profit, digits=4))
        total_cost += cost
        total_profit += profit
    end
    println("Total Cost: ", round(total_cost, digits=4))
    println("Total Profit: ", round(total_profit, digits=4))
end

# ========== Run Solver ========== #
solve_with_constraints()
