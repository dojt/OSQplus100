include("ifpp.jl")

using CSV
using DataFrames

using Statistics: mean

using Printf
using Test

using DOT_NiceMath            # ⋅=* etc
using DOT_NiceMath.Numbers64  # ℝ=Float64, ℤ=Int64, etc

using .Intrinsic_False_Positive_Probability_Bounds
ifpp_lowerbound(;n::ℤ , m::ℤ) = Intrinsic_False_Positive_Probability_Bounds.lowerbound(;n=big(n),α=2/big(m))
ifpp_upperbound(;n::ℤ , m::ℤ) = Intrinsic_False_Positive_Probability_Bounds.upperbound(;n=big(n),α=2/big(m))


# Tabulate some values:

let
    println("n\t1/α\tlb\tub\t")
    for n = 3 : 10
        for m = 2 : 8
            α = 2/m
            lb = ifpp_lowerbound(;n,m)
            ub = ifpp_upperbound(;n,m)

            @printf("%2d\t%5.3f\t%10.3e\t%10.3e\n", n, α, lb, ub)
        end
    end
end

#    n       1/α    lb              ub
#     3      1.000    2.857e-01       3.810e-01
#     3      0.667    3.571e-02       1.607e-01
#     3      0.500    0.000e+00       1.129e-01
#     3      0.400    0.000e+00       9.301e-02
#     3      0.333    0.000e+00       8.229e-02
#     3      0.286    0.000e+00       7.562e-02
#     3      0.250    0.000e+00       7.108e-02
#     4      1.000    2.381e-01       4.571e-01
#     4      0.667    2.381e-02       1.446e-01
#     4      0.500    1.190e-03       9.030e-02
#     4      0.400    0.000e+00       6.975e-02
#     4      0.333    0.000e+00       5.925e-02
#     4      0.286    0.000e+00       5.293e-02
#     4      0.250    0.000e+00       4.874e-02
#     5      1.000    1.935e-01       6.194e-01
#     5      0.667    1.498e-02       1.470e-01
#     5      0.500    5.760e-04       8.156e-02
#     5      0.400    1.200e-05       5.907e-02
#     5      0.333    0.000e+00       4.816e-02
#     5      0.286    0.000e+00       4.183e-02
#     5      0.250    0.000e+00       3.774e-02
#     6      1.000    1.541e-01       6.882e-01
#     6      0.667    8.961e-03       1.225e-01
#     6      0.500    2.560e-04       6.042e-02
#     6      0.400    4.000e-06       4.102e-02
#     6      0.333    3.572e-08       3.211e-02
#     6      0.286    0.000e+00       2.711e-02
#     6      0.250    0.000e+00       2.396e-02
#     7      1.000    1.207e-01       7.839e-01
#     7      0.667    5.144e-03       1.046e-01
#     7      0.500    1.058e-04       4.588e-02
#     7      0.400    1.197e-06       2.920e-02
#     7      0.333    7.875e-09       2.194e-02
#     7      0.286    3.076e-11       1.802e-02
#     7      0.250    0.000e+00       1.560e-02
#     8      1.000    9.341e-02       8.300e-01
#     8      0.667    2.856e-03       8.310e-02
#     8      0.500    4.121e-05       3.239e-02
#     8      0.400    3.268e-07       1.933e-02
#     8      0.333    1.529e-09       1.394e-02
#     8      0.286    4.343e-12       1.113e-02
#     8      0.250    7.539e-15       9.436e-03
#     9      1.000    7.160e-02       8.841e-01
#     9      0.667    1.546e-03       6.638e-02
#     9      0.500    1.530e-05       2.300e-02
#     9      0.400    8.277e-08       1.287e-02
#     9      0.333    2.665e-10       8.909e-03
#     9      0.286    5.311e-13       6.914e-03
#     9      0.250    6.639e-16       5.743e-03
#    10      1.000    5.452e-02       9.109e-01
#    10      0.667    8.200e-04       5.129e-02
#    10      0.500    5.465e-06       1.580e-02
#    10      0.400    1.972e-08       8.284e-03
#    10      0.333    4.252e-11       5.508e-03
#    10      0.286    5.757e-14       4.156e-03
#    10      0.250    4.997e-17       3.381e-03

let α=0.4
    for n = 4 : 36
        lb = ifpp_lowerbound(;n,m=round(ℤ,2/α))
        println("n = $n, \tlb = $lb")
    end
end


using Plots
plotly()

function my_plot_lb_ub()
    # Define ranges for n and m
    n_max = 65

    m_values = [ 3,4,5,6,7,8 ]
    n_values = 3:n_max

    # Compute corresponding α values
    α_values = 2 ./ m_values

    # Initialize matrices to store lower and upper bounds
    lb_matrix = zeros(Float64, length(n_values), length(m_values))
    ub_matrix = zeros(Float64, length(n_values), length(m_values))

    # Compute lb and ub for each combination of n and m
    println("Computing IFPP bounds for n = 3 to $n_max and gven m's...")
    for (i, n) in enumerate(n_values)
        for (j, m) in enumerate(m_values)
            α = α_values[j]

            lb = ifpp_lowerbound(;n, m) |> Float64
            ub = ifpp_upperbound(;n, m) |> Float64

            @assert isfinite(lb)
            @assert isfinite(ub)

            lb_matrix[i, j] = lb
            ub_matrix[i, j] = ub
        end
    end
    println("Computation completed.")

    # Create a list to store individual plots
    plots_list = []

    # Iterate over each α to create separate plots
    for j in 1:length(m_values)
        m = m_values[j]
        α = α_values[j]

        # Filter n_values and corresponding lb/ub values for the current α
        lbvalid_n_indices = findall(n -> m ≤ n ≤ n_max, n_values) # Condition for positive lower bound
        ubvalid_n_indices = findall(n -> 3 ≤ n ≤ n_max, n_values) # Superset
        lbfiltered_n_values = n_values[lbvalid_n_indices]
        ubfiltered_n_values = n_values[ubvalid_n_indices]
        lb_plot = Float64.(lb_matrix[lbvalid_n_indices, j]) # Filter lb_matrix
        ub_plot = Float64.(ub_matrix[ubvalid_n_indices, j]) # Filter ub_matrix

        min_y = min( minimum(lb_plot), minimum(ub_plot), 1e-12) # Find minimum across both datasets and 1e-12 floor
        min_exponent = round(Int, log10(min_y))


        # Create a plot for the current α
        p = scatter(
            lbfiltered_n_values, lb_plot, # Use filtered_n_values here
            label = "",       # Label box takes too much space
            xlabel = "𝑛",
            xticks = 8: 8 :n_max,
            yticks = 10.0 .^( 0: -10 :min_exponent),
            ylabel = "Probability",
            title = "𝛼 = $(round(α, digits=3))",
            yscale = :log10,
            linewidth = 2,
            marker = (:circle, 2),
            markerstrokewidth = 0,
            legend = :topright
        )

        # Add the upper bound to the same plot (using filtered_n_values)
        scatter!(
            p,
            ubfiltered_n_values, ub_plot,  # Use filtered_n_values here
            label = "",       # Label box takes too much space
            linewidth = 2,
            marker = (:diamond, 2),
            markerstrokewidth = 0
        )

        hline!([1e-12], color = :lightgray, linestyle = :dash, label="")

        push!(plots_list, p)
    end

    # Determine the layout based on the number of α values
    num_plots = length(plots_list)
    cols = 2  # Number of columns in the layout
    rows = ceil(Int, num_plots / cols)

    # Combine all individual plots into a single figure with subplots
    combined_plot = plot(plots_list..., layout = (rows, cols), size = (1200, 400 * rows), bottom_margin = 15Plots.mm)

    # Display the combined plot
    display(combined_plot)

    # Optionally, save the plot to a file
    # savefig(combined_plot, "ifpp_bounds_plots.png")

    return plots_list
end

#------------------------------------------------------------------------------
# Compare with data from numerical-stochastic simulations
#------------------------------------------------------------------------------

function plot_Probability_of_Basis_Accept_ESTIMATE(csv_filename::String)
    # Read the CSV file into a DataFrame
    data = CSV.read(csv_filename, DataFrame)

    # Print all column names (for debugging)
    println("Column names in the data:")
    println(names(data))

    # Find unique values of n and α
    n_values = data.n |> unique |> sort
    α_values = unique(data[!, Symbol("α")]) |> sort

    # Filter α_values to only include those where 2/α is an integer
    valid_α_values = filter(α -> isinteger(2/α), α_values)

    println("values of n:        ", n_values...)
    println("Values of α:        ", α_values...)
    println("Valid α values:     ", valid_α_values...)
    println("Discarded α values: ", setdiff(α_values, valid_α_values)...)

    plots_list = []  # Initialize a list to store plots

    for α in valid_α_values
        println("Creating plot for α = $α") # Debugging output

        filtered_data = filter(row -> row[Symbol("α")] == α, data)

        # Compute statistics for plotting
        stats = combine(groupby(filtered_data, :n),
                        :n => first,
                        Symbol("Probability of Basis-Accept, ESTIMATE") => mean => :mean_prob,
                        Symbol("Probability of Basis-Accept, ESTIMATE") => minimum => :min_prob,
                        Symbol("Probability of Basis-Accept, ESTIMATE") => maximum => :max_prob
                        )

        # Prepare data for plotting
        n_values_plot = stats.n
        mean_prob = stats.mean_prob
        min_prob = stats.min_prob
        max_prob = stats.max_prob

        # Get corresponding m value for α
        m = round(Int, 2 / α)

        # --- Consistent n_values and filtering for bounds ---
        n_max = maximum(data.n)
        n_values_for_bounds = 3:n_max

        lb_values = [ ifpp_lowerbound.(; n,m)   for n in n_values_for_bounds ]
        ub_values = [ ifpp_upperbound.(; n,m)   for n in n_values_for_bounds ]

        lbvalid_n_indices = findall(n -> m ≤ n ≤ n_max, n_values_for_bounds)
        ubvalid_n_indices = findall(n -> 3 ≤ n ≤ n_max, n_values_for_bounds)
        lbfiltered_n_values = n_values_for_bounds[lbvalid_n_indices]
        ubfiltered_n_values = n_values_for_bounds[ubvalid_n_indices]
        lb_plot = Float64.(lb_values[lbvalid_n_indices])
        ub_plot = Float64.(ub_values[ubvalid_n_indices])

        min_y = min(minimum(lb_plot), minimum(ub_plot), 1e-12)
        min_exponent = round(Int, log10(min_y))

        # --- Plotting ---
        plt = scatter(
            lbfiltered_n_values, lb_plot,
            label = "", xlabel = "𝑛", xticks = 8:8:n_max,
            yticks = 10.0 .^ (0:-10:min_exponent), ylabel = "Probability",
            title = "α = $(round(α, digits=3))", yscale = :log10,
            linewidth = 2, marker = (:circle, 2), markerstrokewidth = 0,
            legend = :topright
        )

        scatter!(plt, ubfiltered_n_values, ub_plot, label = "",
                 linewidth = 2, marker = (:diamond, 2), markerstrokewidth = 0)

        scatter!(plt, n_values_plot, mean_prob, yerror = (mean_prob .- min_prob, max_prob .- mean_prob),
                 label = "", marker = (:hline, 4), markerstrokewidth = 1)

        hline!([1e-12], color = :lightgray, linestyle = :dash, label = "")

        push!(plots_list, plt)
    end

    # Determine the layout based on the number of α values
    num_plots = length(plots_list)
    cols = 2  # Number of columns in the layout
    rows = ceil(Int, num_plots / cols)

    combined_plot = plot(plots_list..., layout = (rows, cols), size = (1200, 400 * rows), bottom_margin = 15Plots.mm)

    display(combined_plot)  # Display the combined plot


    return plots_list
end

# Example usage: create and display plots for all valid α values in "output.csv"
plots = plot_Probability_of_Basis_Accept_ESTIMATE("output.csv")
plot(plots..., layout=(length(plots), 1)) # Example layout, adjust as needed

function plot_estimates(csv_filename::String)
    # Read and preprocess data (same as in plot_Probability_of_Basis_Accept_ESTIMATE)
    data = CSV.read(csv_filename, DataFrame)

    # Print all column names (for debugging)
    println("Column names in the data:")
    println(names(data))

    # Find unique values of n and α
    n_values = data.n |> unique |> sort
    α_values = unique(data[!, Symbol("α")]) |> sort

    # Filter α_values to only include those where 2/α is an integer
    valid_α_values = filter(α -> isinteger(2/α), α_values)

    n_values_to_plot = [7, 8, 9, 10]  # n values for individual plots
    α_range = 0.51: -0.001 :0.24  # α range for horizontal axis
    plots_list = []

    for n in n_values_to_plot
        plt = plot(xlabel = "α", ylabel = "Probability", yscale = :log10, title = "n = $n", legend=:topright) # Initialize plot

        # Plot upper bounds for all α in the range
        ub_values =
            [
                Intrinsic_False_Positive_Probability_Bounds.upperbound(;n=big(n),α=big(α)) |> Float64
                for α in α_range
            ]
#        plot!(plt, collect(α_range), ub_values, label = "", color=:red)

        # α values for error bars and lower bounds
        α_error_bars = [#=1/2,=# 3/8 #=, 1/3, 1/4=#]
        α_lower_bounds = [1/2, 1/3, 1/4]

        # Plot estimates with error bars
        for α in α_error_bars
            filtered_data = filter(row -> row.n == n && row[Symbol("α")] == α, data)
            if !isempty(filtered_data) # Check if data exists for this (n, α) combination
                mean_prob = mean(filtered_data[!, Symbol("Probability of Basis-Accept, ESTIMATE")])
                min_prob = minimum(filtered_data[!, Symbol("Probability of Basis-Accept, ESTIMATE")])
                max_prob = maximum(filtered_data[!, Symbol("Probability of Basis-Accept, ESTIMATE")])
                scatter!(plt, [α], [mean_prob], yerror = [(mean_prob - min_prob, max_prob - mean_prob)],
                         label = "",
                         marker = (:hline, 4), markerstrokewidth = 1, color=:green)
            else
               println("Warning: No data found for n = $n and α = $α")
            end
        end

        # Plot lower bounds
        for α in α_lower_bounds
            lb = ifpp_lowerbound(; n = n, m = round(Int, 2 / α))
            scatter!(plt, [α], [lb], label = "", marker = (:circle, 2), markerstrokewidth = 0, linewidth=2, color=:blue)
        end

        push!(plots_list, plt)

    end


    # Combined plot (adjust layout as needed)
    combined_plot = plot(plots_list..., layout = (length(n_values_to_plot), 1), size = (800, 300 * length(n_values_to_plot)))
    display(combined_plot)

    return plots_list, combined_plot

end

plots = plot_estimates("output.csv")
