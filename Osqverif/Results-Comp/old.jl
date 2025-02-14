🔝 old ↓ new


#------------------------------------------------------------------------------
# Compare with data from numerical-stochastic simulations
#------------------------------------------------------------------------------
# File paths
input_csv_file = "output.csv"

# Read the CSV file into a DataFrame
data = CSV.read(input_csv_file, DataFrame)

# Display the first few rows to confirm correct reading
# println("Loaded data:")
# println(first(data, 5))

# Print all column names
println("Column names in the data:")
println(names(data))

# Sorting data into structures for plotting
# Convert independent variables to appropriate data structures
n_values = unique(data.n)
α_values = unique(data[!, Symbol("α")])
iterations = unique(data._iter)

println("Unique n values: ", n_values)
println("Unique α values: ", α_values)
println("Unique iteration numbers: ", iterations)

# Example: Sorting data by n and α for analysis
sorted_data = sort(data, [:n, Symbol("α"), :_iter])
# println("Data sorted by n, α, and iteration number:")
# println(first(sorted_data, 5))

# Data structures for plotting (example)
# Group data by (n, α) for potential plotting
grouped_data = groupby(data, [:n, Symbol("α")])
# for group in grouped_data
#     println("Group (n = $(group[!, :n][1]), α = $(group[!, Symbol("α")][1])):")
#     println(first(group, 5))
# end


function plot_Probability_of_Basis_Accept_ESTIMATE_powerof2(;α::ℝ)
        @assert α ∈ [1/2, 1/4, 1/8]  "Only negative powers of 2 allowed for α"

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
    basis_accept_fraction_values = basis_accept_fraction(;α).(n_values_plot)

    # Create the candlestick plot
    plt = plot(n_values_plot, mean_prob, ribbon=(mean_prob .- min_prob, max_prob .- mean_prob),
               label="Probability of Basis-Accept, ESTIMATE", xlabel="n", ylabel="Probability", legend=:topright)
    plt = plot!(n_values_plot, basis_accept_fraction_values, label="Basis Accept Fraction", linestyle=:dash)

    println("Candlestick plot for α = 1/2 has been created.")
    return plt
end

pt = plot_Probability_of_Basis_Accept_ESTIMATE_powerof2(;α=1/2)



-----------------------------------------------------------------------------------------------------------------

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

    println("Values of n:        ", n_values)
    println("Values of α:        ", α_values)
    println("Valid α values:     ", valid_α_values)
    println("Discarded α values: ", setdiff(α_values, valid_α_values))

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

    # Find unique values of n and α
    n_values = data.n |> unique |> sort
    α_values = unique(data[!, Symbol("α")]) |> sort

    # Filter α_values to only include those where 2/α is an integer
    valid_α_values = filter(α -> isinteger(2/α), α_values)

    # Print all column names (for debugging)
    println("Column names in the data:")
    println(names(data))
    println("Values of n:        ", n_values)
    println("Values of α:        ", α_values)

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
        α_error_bars = [1/2, 3/8 , 1/3, 1/4]
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
