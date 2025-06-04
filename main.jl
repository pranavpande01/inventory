using QuadGK  # For numerical integration
using Roots   # For root-finding

# ========== Input Parameters ========== #
struct ItemParams
    η::Float64      # Scale parameter of Pareto
    α::Float64      # Shape parameter of Pareto
    n::Float64      # Power demand index
    h::Float64      # Holding cost
    ω::Float64      # Backlogging cost
    p::Float64      # Selling price
    c::Float64      # Purchasing price
end

T0 = 1.0           # Scheduling period
A = 0.0            # Replenishing cost (set to nonzero if needed)

# Define items
items = [
    ItemParams(10.0, 3.0, 1.5, 2.0, 5.0, 20.0, 10.0),
    ItemParams(15.0, 2.5, 0.8, 1.5, 4.0, 25.0, 12.0)
]

# ========== Helper Functions ========== #

function μ(i::ItemParams)
    return i.α * i.η / (i.α - 1)
end

function S0(i::ItemParams)
    β = i.h / (i.h + i.ω)
    threshold = i.n / (i.α + i.n)

    if threshold <= β
        return i.η * ((i.α + i.n) * i.ω / (i.α * (i.h + i.ω)))^(1/i.n)
    else
        return i.η * (i.n * (i.h + i.ω) / (i.h * (i.α + i.n)))^(1/i.α)
    end
end

function expected_inventory_cost(i::ItemParams, S::Float64)
    if S <= i.η
        J = 0.0
    else
        J = i.η^i.α / (i.α - 1) * S^(-i.α + 1) + S - i.α * i.η / (i.α - 1)
    end

    μi = μ(i)
    C = (i.h + i.ω) * i.n / (i.n + 1) * J +
        i.ω * i.n / (i.n + 1) * (μi - S)
    return C
end

function expected_profit(i::ItemParams, S::Float64)
    μi = μ(i)
    revenue = (i.p - i.c) * μi / T0
    cost = expected_inventory_cost(i, S)
    return revenue - cost
end

# ========== Main Routine ========== #

for (j, item) in enumerate(items)
    println("Item $(j):")
    S_star = S0(item)
    cost = expected_inventory_cost(item, S_star)
    profit = expected_profit(item, S_star)

    println("  Optimal Inventory Level S*: ", round(S_star, digits=4))
    println("  Expected Cost: ", round(cost, digits=4))
    println("  Expected Profit: ", round(profit, digits=4))
end

