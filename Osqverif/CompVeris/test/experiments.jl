# CompVeris/test/experiments.jl
# (c) Dirk Oliver Theis
# Open source license: CC0

using DOT_NiceMath
using DOT_NiceMath.Numbers64

using CompVeris
using CompVeris.Verif: Decision_t, Accept, Reject
using CompVeris.Verif: DFE, MYNS
using CompVeris.Stabilizer_Foolers: Scenario_t, FALSE_POSITIVE, FALSE_NEGATIVE

using JSON
using Printf

using Plots

using Test

# ========================================================================================
# Plot results
# ========================================================================================

function make_result_plots(results)
    (;false_accepts,false_rejects) = results

    # Create separate plots for false accepts and false rejects
    p_accepts = plot(title="False Accepts", xlabel="m", ylabel="False Accepts %", legend=:topright)
    p_rejects = plot(title="False Rejects", xlabel="m", ylabel="False Rejects %", legend=:topright)


    methods = collect(keys(false_accepts)) ∪ collect(keys(false_rejects)) |> sort |> unique!


    for method in methods
        data = false_accepts[method]
        x = collect(keys(data))
        y = collect(values(data)) * 100 # percent %
        scatter!(p_accepts, x, y, label=method, markersize=1, markerstrokewidth=0)
    end

    for method in methods
        data = false_rejects[method]
        x = collect(keys(data))
        y = collect(values(data)) * 100 # percent %
        scatter!(p_rejects, x, y, label=method, markersize=1, markerstrokewidth=0)
    end

    return (false_accepts = p_accepts,
            false_rejects = p_rejects)
end

function make_margin_plot(results)
    margins = results.margin  # Extract margin data

    if isempty(margins)
        println("No margin data found.")
        return nothing # Or handle the case as needed
    end

    num_instances = length(margins)

    breakpoints = [margins[i].om_breakpoint for i in 1:num_instances]
    expvals = [margins[i].om_expval for i in 1:num_instances]

    p = plot(title="Margin for Min-of-Means", xlabel="Data point", ylabel="Value", legend=:topright)
    scatter!(p, 1:num_instances, breakpoints,
             label="",#"Breakpoint",
             markersize=1, markerstrokewidth=1, markershape=:hline)
    scatter!(p, 1:num_instances, expvals,
             label="",#"Expval",
             markersize=1, markerstrokewidth=0)

    return p #return plot object
end

# ========================================================================================
# Test
# ========================================================================================

compare()
results = compare(; δ=1e-1, n=2, num_instances=2, numerators_of_m = 1.0: -0.1 :0.1, reps=5)
plots = make_result_plots(results)

full_set__compare(DENSITY_DIAGONAL_INSTANCES, n_range=2:3, d_range=2:3, num_instances=2, reps=2)

full_set__compare_bad_states(DENSITY_DIAGONAL_INSTANCES_SMARTBAD ;
                             n_range=3:7, d_range=5:6, num_instances=2, reps=2)

# ========================================================================================
# Profiling
# ========================================================================================

using Profile
using ProfileView

@profile compare(; δ=1e-2, n=5, num_instances=1, numerators_of_m=9:-1:7, reps=10)

#
# Either:
Profile.print()
#
# Or:
ProfileView.view()
#

# ========================================================================================
# Running it
# ========================================================================================

results = compare(; δ=1e-2, n=5, num_instances=10,
                  numerators_of_m = 0.01: 0.01 :1.0, reps=100)
plots = make_result_plots(results)
display(plot(plots.false_accepts, plots.false_rejects, layout = (2, 1)))
# savefig("false_accepts_rejects.png")  # Optionally save to a file

# ========================================================================================

results = compare(; δ=1e-2, n=5, num_instances=100, numerators_of_m = 4.0: 1.0 :4.0, reps=1000,   keep_RNGs=true)
bounds = CompVeris.gap(results.margin, results.expval_broken)
println("midpoint: ",(bounds.ub+bounds.lb)/2, " margin: ",(bounds.ub-bounds.lb)/2)

display(make_margin_plot(results))

plots = make_result_plots(results)
display(plot(plots.false_accepts, plots.false_rejects, layout = (2, 1)))


# ========================================================================================
# Indep UAR Stabilizers
# ========================================================================================


using Plots
plotly()
#using StatsBase
#using StatsPlots

using DOT_NiceMath
using DOT_NiceMath.Numbers64

using CompVeris
using CompVeris.Stabilizer_Foolers

Supp(f) = begin L=length(f) ; [ ξ for ξ=1:L if f[ξ] != 0 ] end

function slack(p ::Vector{ℝ} , y ::ℤ ; N, α)
    α/2  −  ∑((count_ones(y & x) % 2) * p[x] for x in 1:N-1)
end

function accept_rows(p ::Vector{ℝ};n,α)
    N = 2^n ; @assert N-1 == length(p)
    [ y for y = 1:N-1 if slack(p,y;N,α) ≥ -1e-10 ]
end

function unif_p( Ξ ;n)
    N = 2^n
    p = zeros(ℝ, N-1)
    𝑝 = 1/length(Ξ)
    for ξ ∈ Ξ
        p[ξ] = 𝑝
    end
    return p
end

function basdens(acceptrows ::Vector{ℤ} ; n, repet) ::ℝ
    N = 2^n
    L = length(acceptrows)
    n ≤ L || return 0.0

    count = 0
    for _r = 1:repet
        B = sample(1:L, n; replace=false)
        for i=1:n
            B[i] = acceptrows[B[i]]
        end
        if Stabilizer_Foolers._rank(B;n) == n
            count +=1
        end
    end
    return count / repet
end

function basis_accept_fraction(;n)
    # 𝐧 ::BigInt = n
    𝟓 ::BigInt = 5
    𝟒 ::BigInt = 4
    𝟑 ::BigInt = 3
    𝟐 ::BigInt = 2
    (  𝟓^n − 4⋅𝟒^n + 6⋅𝟑^n − 4⋅𝟐^n + 1  )/( (𝟐^n-1)⋅(𝟐^n-2)⋅(𝟐^n-4)⋅(𝟐^n-8) )
end

plt =
    let # reduced
        plt = nothing # plot(;xlabel="Slack")
        minute = minutes = Stabilizer_Foolers.min

        for n = 4:10
            N = 2^n
            α = 0.5 # 2/n ⋅ 1.0
            println("""_____________________________________________________________________________________
                       Start: n=$(n), 2/n=$(2/n), α=$(α)""")

            (;opt,slacks) =
                reduced_false_positive_analysis(
                    ; n, α,
                    time_limit       = max(1,n-7) * 10minutes,
                    perturb          = 1,
                    # fix_Id           = false, # default: yes, iff α ≥ 2/n
                    basis_repeats = 1_000_000 )

            println("""_____________________________________________________________________________________
                       Finished n=$(n), 2/n=$(2/n), α=$(α)""")

            p = abs.(round.(opt.p .* (1<<16))) ./ (1<<16)
            @show Supp(p)
            @show sort(unique(p))
            @show p

            if plt !== nothing
                slacks *= -1 ; sort!(slacks)
                vals = unique( slacks )
                pdf  = [ count(s->(s==x), slacks)/length(slacks)  for x ∈ vals ]
                cdf  = cumsum(pdf)
                #@show vals
                #@show cdf
                plot!(plt,vals,cdf
                      ; seriestype=:steppost, label="n=$(n),α-2/n=$(α-2/n)")
            end
        end
        plt
    end



let
    for n = 4:4
        α = 2/n ⋅ 1.1
        @show α
        minute = minutes = Stabilizer_Foolers.min

        for d = 12 : 12
            δ = min( 1/5 , (1/2)^d )
            false_positive_analysis( ; n, δ,
                                     ε        = α⋅δ,
                                     fid      = 1-δ,
                                     time_limit = 2minutes,
                                     skip_basis_stuff = false)
        end
    end
end

let n=8
    minute = minutes = Stabilizer_Foolers.min
    #for δ = min(1/3,2/n) : -0.1: 1e-3
    δ = 1.52587890625
    while (δ > 1e-6)  δ /= 10
        false_negative_analysis( ; n, δ,
                                 time_limit = 1minute,
                                 skip_basis_stuff = false)
    end
end


# ========================================================================================
# Working With Bases
# ========================================================================================


num_bases(;dim::Int64) = begin @assert UInt64(dim)^2<64 ; prod( (UInt64(2)^dim-2^j) for j=0:dim-1 ) |> Int64 end
num_bases(dim::Int) = begin 𝟐=BigInt(2); prod( (𝟐^dim-𝟐^j) for j=0:dim-1 ) end
