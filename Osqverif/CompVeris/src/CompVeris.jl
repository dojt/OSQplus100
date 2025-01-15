# CompVeris/src/CompVeris.jl
# (c) Dirk Oliver Theis
# Open source license: CC0

module CompVeris

# ====================================================================================================
# Global helpers
# ====================================================================================================

fourier_coeff(x::UInt64,y::UInt64)  = ( count_ones(x & y) % 2 == 1 ? -1.0 : +1.0 )

# ====================================================================================================
# Include sub-modules
# ====================================================================================================

include("./density_operators.jl")
include("./instances.jl")
include("./measure.jl")
include("./verif.jl")

include("./stabilizer_foolers.jl")

# ====================================================================================================
# Exports
# ====================================================================================================

export Density_Operators, Instances, Measure, Verif
export Stabilizer_Foolers
export compare, full_set__compare, full_set__compare_bad_states
export Instance_type_t,
    DENSITY_OPERATOR_INSTANCES, DENSITY_DIAGONAL_INSTANCES, DENSITY_DIAGONAL_INSTANCES_SMARTBAD

# ====================================================================================================
# Full processes
# ====================================================================================================

using .Verif: Decision_t, Accept, Reject
using .Verif: DFE, MYNS
using .Stabilizer_Foolers: Scenario_t, FALSE_POSITIVE, FALSE_NEGATIVE

using .Density_Operators: ℂHermitian_t
using .Measure:           Operator_or_Diagonal_t

using LinearAlgebra: Hermitian

using JSON
using Printf

using DOT_NiceMath
using DOT_NiceMath.Numbers64


function intrinsic_falsepositive_12(;n::ℤ)
    @assert n ≥ 1

    function basis_accept_fraction()
        𝟓 ::BigInt = 5
        𝟒 ::BigInt = 4
        𝟑 ::BigInt = 3
        𝟐 ::BigInt = 2
        (  𝟓^n − 4⋅𝟒^n + 6⋅𝟑^n − 4⋅𝟐^n + 1  )/( (𝟐^n-1)⋅(𝟐^n-2)⋅(𝟐^n-4)⋅(𝟐^n-8) )
    end

    if n ≤ 3
        return (lower_bound=nothing, upper_bound=0.0)
    end

    ret = (lower_bound = basis_accept_fraction(),
           upper_bound = -1.0                       )

    if n==4
        ret.upper_bound = (2⋅3⋅4)/Stabilizer_Foolers.numbas(4) # Basis Worst-Case Bound
    else
        @assert :to==:do # TODO: Implement exact bound from proof of
                         # Theorem thm:false-positive-asymptotic
                         # "Asymptotic Analysis of False-Positives"
    end
end



const Margin_t = @NamedTuple{ om_breakpoint::ℝ, om_expval::ℝ }

@enum Instance_type_t DENSITY_OPERATOR_INSTANCES DENSITY_DIAGONAL_INSTANCES DENSITY_DIAGONAL_INSTANCES_SMARTBAD


function handle_broken_state(scen ::Scenario_t, ϱ::Operator_or_Diagonal_t; n::ℤ,δ::ℝ)

    filename =
        if scen==FALSE_NEGATIVE
            Printf.@sprintf("broken-goodstate-%02u-%#a.bin",n,δ)
        else # scen==FALSE_POSITIVE
            Printf.@sprintf("broken-badstate-%02u-%#a.bin",n,δ)
        end

    open(filename, "w") do io
        if ϱ isa Vector{ℝ}
            write(io, ϱ) # Write the diagonal directly
        elseif ϱ isa ℂHermitian_t
            write(io, ϱ.data)  # Write the data field of the Hermitian matrix
        else
            @assert false "You shouldn't be here!!"
        end
    end
end

function compare_good_states(inst_type :: Instance_type_t
                             ;
                             δ                = 1e-2,
                             n                = 4,
                             num_instances    = 1,
                             numerators_of_m  = 1:2,
                             reps             = 10,
                             keep_RNGs        = false)

    @assert δ > 1e-9 "Are you nuts?!!"

    all_false_rejects = Dict{String,Dict{Int,Float64}}()
    margin            = Vector{ Margin_t }()
    expval_broken     = false

    keep_RNGs ||  Instances.reset()

    m_set = collect(Int(round(numerator_of_m / δ)) for numerator_of_m in numerators_of_m) |> sort |> unique!

    println()

    for _ii = 1:num_instances

        ϱ_good ::Operator_or_Diagonal_t =
            if     inst_type == DENSITY_OPERATOR_INSTANCES
                ψ₀ = zeros(ℂ, 2^n)
                ψ₀[1] = 1   # Set the first element to 1 to represent state |0...0⟩
                Instances.make_rnd_dens_op(; n, ψ₀, ε=δ)    # basis independent
            elseif inst_type == DENSITY_DIAGONAL_INSTANCES
                Instances.make_rnd_dens_diag(Instances.Dumb(); n, ε=δ) # basis dependent: only works for ψ₀=|0...0⟩
            else
                @assert false "You shouldn't be here!"
            end

        for m in m_set

            wrong_counts = Dict{String,Int}()

            println("Random GOOD states:  \tδ = $(δ), m = $(m)")

            keep_RNGs || Measure.reset()
            keep_RNGs || MYNS.reset()
            keep_RNGs || DFE.reset()

            for _iter = 1:reps

                let
                    list = MYNS.verify( ; n, δ, ϱ=ϱ_good,
                                        num_shots_per_stab=m)
                    for (name,stuff) in list
                        get!(wrong_counts, name, 0)
                        (;decision,witness,expval) = stuff
                        wrong = decision != Accept
                        if name=="Min-of-Means"
                            breakpoint = 1-3δ
                            if wrong
                                wrong_counts[name] += 1
                                if (expval < breakpoint)
                                    println("👉❗️❗️❗️ expval=$(expval), breakpoint=$(breakpoint), witness=$(witness)")
                                    expval_broken = true
                                    handle_broken_state(FALSE_NEGATIVE,ϱ_good;n,δ)
                                end
                            end
                            push!(margin, (om_breakpoint=1-breakpoint,om_expval=1-expval) )
                        end
                    end
                end

                let
                    dfe_ret = DFE.verify(; n, δ, ϱ=ϱ_good,
                                         num_stabs = m, num_shots_per_stab = 1)
                    get!(wrong_counts, "DFE(m,1)", 0)
                    dfe_ret.decision==Accept || ( wrong_counts["DFE(m,1)"] += 1 )
                end
                let
                    dfe_ret = DFE.verify(; n, δ, ϱ=ϱ_good,
                                         num_stabs = m, num_shots_per_stab = n)
                    get!(wrong_counts, "DFE(m,n)", 0)
                    dfe_ret.decision==Accept || ( wrong_counts["DFE(m,n)"] += 1 )
                end
                let
                    dfe_ret = DFE.verify(; n, δ, ϱ=ϱ_good,
                                         num_stabs = m, num_shots_per_stab = n)
                    get!(wrong_counts, "DFE(n,m)", 0)
                    dfe_ret.decision==Accept || ( wrong_counts["DFE(n,m)"] += 1 )
                end

            end #^ for (reps)

            for (name, count) in wrong_counts
                get!(
                    get!(all_false_rejects, name, Dict{Int,Float64}()),
                    m,
                    0.0)
                all_false_rejects[name][m] += count / (num_instances ⋅ reps)
            end

        end #^ for m (numerators)

    end #^ for (num_instances)

    println("\nCompleted.")

    return (false_rejects=all_false_rejects, margin, expval_broken)
end #^ compare_good_states()


function compare_bad_states(inst_type :: Instance_type_t
                            ;
                            δ                = 1e-2,
                            n                = 4,
                            num_instances    = 1,
                            numerators_of_m  = 1:2,
                            reps             = 10,
                            keep_RNGs        = false)

    @assert δ > 1e-9 "Are you nuts?!!"

    all_false_accepts = Dict{String,Dict{Int,Float64}}()
    margin            = Vector{ Margin_t }()
    expval_broken     = false

    keep_RNGs || Instances.reset()

    m_set = collect(Int(round(numerator_of_m / δ)) for numerator_of_m in numerators_of_m) |> sort |> unique!

    println()

    for _ii = 1:num_instances

        ϱ_bad ::Operator_or_Diagonal_t =
            if     inst_type == DENSITY_OPERATOR_INSTANCES
                ψ₀ = zeros(ℂ, 2^n)
                ψ₀[1] = 1    # Set the first element to 1 to represent state |0...0⟩
                Instances.make_rnd_dens_op(; n, ψ₀, ε=5δ)     # basis independent
            elseif inst_type == DENSITY_DIAGONAL_INSTANCES    # basis dependent: only works for ψ₀=|0...0⟩
                Instances.make_rnd_dens_diag(Instances.Dumb(); n, ε=5δ)
            elseif inst_type == DENSITY_DIAGONAL_INSTANCES_SMARTBAD
                β = 3δ # ❗ breakpoint  ❗ Change this when the breakpoint changes  ❗ ❗ ❗ ❗ ❗ ❗ ❗ ❗ ❗
                Instances.make_rnd_dens_diag(Instances.Smart(); n, ε=5δ, α=β/5δ)
            else
                @assert false "You shouldn't be here!"
            end

        for m in m_set

            wrong_counts = Dict{String,Int}()

            println("Random BAD states: \tδ = $(δ), m = $(m)")

            keep_RNGs || Measure.reset()
            keep_RNGs || MYNS.reset()
            keep_RNGs || DFE.reset()

            for _iter = 1:reps

                let
                    list = MYNS.verify( ; n, δ, ϱ=ϱ_bad,
                                        num_shots_per_stab=m)
                    for (name,stuff) in list
                        get!(wrong_counts, name, 0)
                        (;decision,witness,expval) = stuff
                        wrong = decision != Reject
                        if name=="Min-of-Means"
                            breakpoint = 1-3δ
                            if wrong
                                wrong_counts[name] += 1
                                if (expval > breakpoint)
                                    println("👉❗️❗️❗️ expval=$(expval), breakpoint=$(breakpoint), witness=$(witness)")
                                    expval_broken = true
                                    handle_broken_state(FALSE_POSITIVE,ϱ_bad;n,δ)
                                end
                            end
                            push!(margin, (om_breakpoint=1-breakpoint,om_expval=1-expval) )
                        end
                    end
                end

                let
                    dfe_ret = DFE.verify(; n, δ, ϱ=ϱ_bad,
                                         num_stabs = m, num_shots_per_stab = 1)
                    get!(wrong_counts, "DFE(m,1)", 0)
                    dfe_ret.decision==Reject || ( wrong_counts["DFE(m,1)"] += 1 )
                end
                let
                    dfe_ret = DFE.verify(; n, δ, ϱ=ϱ_bad,
                                         num_stabs = m, num_shots_per_stab = n)
                    get!(wrong_counts, "DFE(m,n)", 0)
                    dfe_ret.decision==Reject || ( wrong_counts["DFE(m,n)"] += 1 )
                end
                let
                    dfe_ret = DFE.verify(; n, δ, ϱ=ϱ_bad,
                                         num_stabs = n, num_shots_per_stab = m)
                    get!(wrong_counts, "DFE(n,m)", 0)
                    dfe_ret.decision==Reject || ( wrong_counts["DFE(n,m)"] += 1 )
                end

            end #^ for (reps)


            for (name, count) in wrong_counts
                get!(
                    get!(all_false_accepts, name, Dict{Int,Float64}()),
                    m,
                    0.0)
                all_false_accepts[name][m] += count / (num_instances ⋅ reps)
            end

        end #^ for m (numerators)

    end #^ for (num_instances)

    println("\nCompleted.")

    return (false_accepts=all_false_accepts, margin, expval_broken)
end #^ compare_good_states()

function compare(inst_type :: Instance_type_t
                ; kwargs...)

    good_results = compare_good_states(inst_type;kwargs...)
    bad_results = compare_bad_states(inst_type;kwargs...)

    # Combine margins.  Since we are just concatenating them, we can do it directly.
    margin = vcat(good_results.margin, bad_results.margin)

    # Check if expval_broken is true in either result
    expval_broken = good_results.expval_broken || bad_results.expval_broken

    return (false_accepts = bad_results.false_accepts, 
            false_rejects = good_results.false_rejects,
            margin, expval_broken)

end #^ compare()

function gap(margins ::Vector{ Margin_t }, expval_broken)

    if expval_broken
        println("❗️❗️❗️ Expectation values are broken: Some are on the wrong side of the breakpoint❗️❗️❗️ ")
        return (ub=NaN,lb=NaN)
    end

    ub = +Inf
    lb = -Inf
    for m in margins
        if m.om_expval > m.om_breakpoint
            ub = min(ub, m.om_expval)
        else
            lb = max(lb, m.om_expval)
        end
    end
    return (ub=ub,lb)
end #^ gap

function full_set__compare(inst_type :: Instance_type_t
                           ;
                           n_range       = 2:6,
                           d_range       = 2:19,
                           num_instances = 32,
                           reps          = 1024 )
    @assert inst_type ∈ [DENSITY_OPERATOR_INSTANCES, DENSITY_DIAGONAL_INSTANCES] """
    `SMARTBAD`instances only work with `compare_bad_states` functions.
    """

    for n = n_range
        for d = d_range
            δ = (1/2)^d
            filename = Printf.@sprintf("compare-%02u-%#a.json", n,δ)

            results =
                compare(inst_type; δ, n, num_instances, reps,
                        keep_RNGs=true,
                        numerators_of_m = δ:δ) # numbers of individual shots don't play a role

            (;margin,expval_broken) = results
            bds = gap(margin, results.expval_broken)
            open(filename, "w") do io
                JSON.print(io, Dict("expval_broken" => expval_broken, "bounds" => bds))
            end
        end
    end
end


function full_set__compare_bad_states(inst_type :: Instance_type_t
                                      ;
                                      n_range       = 2:6,
                                      d_range       = 2:19,
                                      num_instances = 32,
                                      reps          = 1024)
    for n = n_range
        for d = d_range
            δ = (1/2)^d
            filename = Printf.@sprintf("compareBad-%02u-%#a.json", n,δ)

            results =
                compare_bad_states(inst_type; δ, n, num_instances, reps,
                                   keep_RNGs=true,
                                   numerators_of_m = δ:δ) # numbers of individual shots don't play a role

            (;margin,expval_broken) = results
            bds = gap(margin, results.expval_broken)
            open(filename, "w") do io
                JSON.print(io, Dict("expval_broken" => expval_broken, "bounds" => bds))
            end
        end
    end
end

end # module CompVeris
#EOF
