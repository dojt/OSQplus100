using JSON
using DataFrames
using Plots

# Load JSON file
function load_json(filepath)
    open(filepath, "r") do file
        return JSON.parse(file)
    end
end

# Process JSON data into a DataFrame
function process_data(json_data)
    data = DataFrame(n=Int[], delta=Float64[], ub=Float64[], lb=Float64[], diff=Float64[], mu=Float64[], alpha=Float64[])
    for entry in json_data
        n = entry["n"]
        delta = entry["delta"]
        ub = entry["data"]["bounds"]["ub"]
        lb = entry["data"]["bounds"]["lb"]
        diff = ub - lb
        mu = (ub + lb) / 2
        alpha = ub/delta
        push!(data, (n, delta, ub, lb, diff, mu, alpha))
    end

    return data
end

# Fit c values for diff = c * delta
function fit_coefficients(data, y_col)
    results = DataFrame(n=Int[], c=Float64[])
    for n in unique(data.n)
        subset = filter(row -> row.n == n, data)
        x = subset.delta
        y = subset[!, y_col]
        c = sum(x .* y) / sum(x .* x)  # Simple linear regression
        push!(results, (n, c))
    end
    return results
end

# Load data
file_path = "fixed-all-compare.json"  # Replace with your JSON file path
json_data = load_json(file_path)

# Process data
df = process_data(json_data)

# Fit coefficients for diff and mu
diff_coefficients = fit_coefficients(df, :diff)
mu_coefficients = fit_coefficients(df, :mu)

# Print tables
println("Difference Coefficients:")
println(diff_coefficients)

println("\nMidpoint Coefficients:")
println(mu_coefficients)

# Plot diff vs. delta for each n
plt = plot(
    xscale=:log10, yscale=:log10,
    xlabel="Delta (log scale)", ylabel="Difference (log scale)",
    title="Difference Between Upper and Lower Bounds vs. Delta"
)

for n in unique(df.n)
    subset = filter(row -> row.n == n, df)
    plot!(subset.delta, subset.diff, marker=:circle, label="n=$n")
end
plt = plot!(legend=:topright)

plt2 = plot(
    xscale=:log10, # yscale=:log10,
    xlabel="Delta (log scale)", ylabel="Alpha",
    title="Alpha vs Delta"
)

for n in unique(df.n)
    subset = filter(row -> row.n == n, df)
    sorted_indices = sortperm(subset.delta)
    delta_sorted = subset.delta[sorted_indices]
    alpha_sorted = subset.alpha[sorted_indices]
    plot!(delta_sorted, alpha_sorted, marker=:circle, label="n=$n")
end
plt2 = plot!(legend=:topright)



#savefig("difference_vs_delta.png")  # Save plot as PNG
