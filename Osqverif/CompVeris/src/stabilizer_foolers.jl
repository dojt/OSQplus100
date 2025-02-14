# CompVeris/src/density_operators.jl
# (c) Dirk Oliver Theis
# Open source license: CC0
module Stabilizer_Foolers

export false_positive_analysis, false_negative_analysis
export reduced_false_positive_analysis

export lpfeas_stabilizer_fooler, lpfeas_stabilizer_fooler_g

using HiGHS
const THE_OPTIMIZER = HiGHS.Optimizer

using DOT_NiceMath
using DOT_NiceMath.Numbers64

# import LinearAlgebra
# using  Statistics: mean
# using Distributions

using JuMP
using JuMP

using Mods
using StatsBase: sample
using LinearAlgebraX: rankx, nullspacex, detx
import Combinatorics

import ..fourier_coeff
fourier_coeff(x ::I, y ::I) where{I<:Integer} = fourier_coeff(UInt64(x),UInt64(y))

const Opt_ℝ = Union{ℝ,Nothing}
const min   = 60.0
@enum Scenario_t FALSE_POSITIVE FALSE_NEGATIVE

function mip_stabilizer_fooler( ;
                                n        ::ℤ,
                                ε        ::ℝ,
                                fid      ::Opt_ℝ      = nothing,
                                p        ::Opt_ℝ      = nothing,
                                scenario ::Scenario_t = FALSE_POSITIVE,
                                time_limit =10min)
    @assert p===nothing || fid===nothing
    @assert p!==nothing || fid!==nothing
    @assert p  ===nothing || 0 ≤ p   ≤ 1
    @assert fid===nothing || 0 ≤ fid ≤ 1

    @assert 0 ≤ ε  ≤ 1
    @assert n ≥ 1

    N = 2^n

    model = Model(HiGHS.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "time_limit", time_limit)

    @variable(model, f[0:N-1] ≥ 0)
    @variable(model, big[1:N-1], Bin)

    @constraint(model, ∑(f[x] for x in 0:N-1) == 1)

    if fid !== nothing
        if scenario == FALSE_POSITIVE
            println("mip_stabilizer_fooler: Maximizing probability of ACCEPT for upper-bounded fidelity (false positive scenario)")
            @objective(model, Max, ∑(big[y] for y in 1:N-1))
            @constraint(model, f[0] ≤ fid)
            for y in 1:N-1
                @constraint(model,
                            ∑(   fourier_coeff(x,y) ⋅ f[x] for x in 0:N-1 )
                            ≥
                            (1-ε) ⋅ big[y]
                           )
            end
        elseif scenario == FALSE_NEGATIVE
            println("mip_stabilizer_fooler: Maximizing probability of REJECT for lower-bounded fidelity (false negative scenario)")
            @objective(model, Min, ∑(big[y] for y in 1:N-1))
            @constraint(model, f[0] ≥ fid)
            for y in 1:N-1
                @constraint(model,
                            ∑(   fourier_coeff(x,y) ⋅ f[x] for x in 0:N-1 )
                            ≤
                            1 − ε⋅(1-big[y])
                           )
            end
        else
            @assert false "This shoudn't be reachable!"
        end
    elseif p !== nothing
        if scenario == FALSE_POSITIVE
            println("mip_stabilizer_fooler: Minimizing fidelity for lower-bounded probability of ACCEPT (false positive scenario)")
            @objective(model,  Min, f[0])
            @constraint(model, ∑( big[y] for y = 1:N-1 ) ≥ p⋅(N-1))
            for y in 1:N-1
                @constraint(model,
                            ∑(   fourier_coeff(x,y) ⋅ f[x] for x in 0:N-1 )
                            ≥
                            (1-ε) ⋅ big[y]
                           )
            end
        elseif scenario == FALSE_NEGATIVE
            println("mip_stabilizer_fooler: Maximizing fidelity for upper-bounded probability of ACCEPT (false negative scenario)")
            @objective(model,  Max, f[0])
            @constraint(model, ∑( big[y] for y = 1:N-1 ) ≤ p⋅(N-1))
            for y in 1:N-1
                @constraint(model,
                            ∑(   fourier_coeff(x,y) ⋅ f[x] for x in 0:N-1 )
                            ≤
                            1 − ε⋅(1-big[y])
                           )
            end
        else
            @assert false "This shoudn't be reachable!"
        end
    else
        @assert false "This shoudn't be reachable!"
    end


    optimize!(model)

    if termination_status(model) != MOI.OPTIMAL
        @show termination_status(model)
    end

    return (status       = termination_status(model),
            primal_bound = objective_value(model),
            dual_bound   = objective_bound(model),
            f            = JuMP.value.(f),
            big          = JuMP.value.(big)                 )
end

function lpfeas_stabilizer_fooler( ;
                                   n        ::ℤ,
                                   ε        ::ℝ,
                                   fid      ::Opt_ℝ,
                                   scenario ::Scenario_t = FALSE_NEGATIVE,
                                   time_limit =10min)
    @assert scenario == FALSE_NEGATIVE   "This makes sense only for the false-negative scenario"
    @assert 0 ≤ fid ≤ 1

    @assert 0 ≤ ε  ≤ 1
    @assert n ≥ 1

    N = 2^n

    println("lpfeas_stabilizer_fooler: Testing for single big y for lower-bounded fidelity (false negative scenario)")

    returnlist = []

    for y=1:N-1
        model = Model(HiGHS.Optimizer)
        set_silent(model)
        set_optimizer_attribute(model, "time_limit", time_limit)

        @variable(model, f[0:N-1] ≥ 0)

        @objective(model, Min,   ∑( fourier_coeff(x,y) ⋅ f[x]   for x in 0:N-1 )    )

        @constraint(model, ∑(f[x] for x in 0:N-1) == 1)

        @constraint(model, f[0] ≥ fid)

        optimize!(model)

        if termination_status(model) != MOI.OPTIMAL
            @show termination_status(model)
        end

        print("\nf̂($(y)) = \t$(objective_value(model))")

        if objective_value(model) ≤ 1 - 0.9ε  ||  termination_status(model) != MOI.OPTIMAL
            print("\t❗️")
            push!(returnlist,
                    (y            = y,
                     status       = termination_status(model),
                     primal_bound = objective_value(model),
                     dual_bound   = objective_bound(model),
                     f            = JuMP.value.(f)              )
                   )
        end
    end #^ for

    return returnlist
end

function lpfeas_stabilizer_fooler_g( ;
                                     n        ::ℤ,
                                     δ        ::ℝ,
                                     time_limit = 10min,
                                     dumpsols = false)
    @assert 0 < δ < 1
    @assert n ≥ 1

    N = 2^n

    returnlist = []

    for y = 1:N-1
        model = Model(HiGHS.Optimizer)
        set_silent(model)
        set_optimizer_attribute(model, "time_limit", time_limit)

        @variable(model, g[0:N-1])  # g can be negative at 0

        @objective(model, Min, ∑(fourier_coeff(x, y) * g[x] for x in 0:N-1))

        if dumpsols
            # name constrants
            @constraint(model, g[0] ≥ -δ, base_name = "0l")
            @constraint(model, g[0] ≤  0, base_name = "0u")
            @constraint(model, ∑( g[x] for x in 0:N-1 ) == 0, base_name = "sum")
            @constraint(model, [x in 1:N-1], g[x] ≥ 0, base_name = "nng")
        else
            @constraint(model, -δ ≤ g[0] ≤ 0)
            @constraint(model, [x in 1:N-1], g[x] ≥ 0 )
            @constraint(model, ∑( g[x] for x in 0:N-1 )   == 0)
        end

        optimize!(model)

        status = termination_status(model)
        primal_bound = objective_value(model)

        print("\nĝ($y) = \t$primal_bound")

        if primal_bound < -2δ || status != MOI.OPTIMAL
            print("\t❗️")
            push!(returnlist, (y            = y,
                               status       = status,
                               primal_bound = primal_bound,
                               dual_bound   = objective_bound(model),
                               g            = JuMP.value.(g)))
        end

        if dumpsols
            dual_values = Dict()
            for (con_name, con_obj) in list_of_constraint_types(model)
                dual_values[con_name] = dual.(all_constraints(model, con_name, con_obj))
            end

            dual_values["0l"] = dual(constraint_by_name(model, "0l"))
            dual_values["0u"] = dual(constraint_by_name(model, "0u"))
            dual_values["sum"] = dual(constraint_by_name(model, "sum"))

            #
            println("\nPrimal:")
            println(JuMP.value.(g))
            println("Dual:")
            println(dual_values)
        end

    end

    println()
    return returnlist
end


# --------------------------------------------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------------------------------------------

function estimate_bigprb(f :: Vector{ℝ},
                         ;
                         n     ::ℤ,
                         ε     ::ℝ,
                         repet ::Int,
                         eps   ::ℝ = 1e-9 ) ::ℝ

    N = 2^n
    @assert size(f) == (N,)
    @assert 0 ≤ ε ≤ 1

    count = 0
    for _r = 1:repet
        y = rand(1:N-1)
        if ∑( fourier_coeff(y,x) ⋅ f[1+ x] for x in 0:N-1 ) ≥ 1-ε − eps
            count +=1
        end
    end
    return count / repet
end

function estimate_acceptprb(f :: Vector{ℝ},
                            ;
                            n     ::ℤ,
                            ε     ::ℝ,
                            repet ::Int,
                            eps   ::ℝ = 1e-9) ::ℝ

    N = 2^n
    @assert size(f) == (N,)
    @assert 0 ≤ ε ≤ 1

    B ::Vector{Int64} = Vector{Int64}(undef,n)

    count = 0
    for _r = 1:repet
        for j=1:n
            B[j] = rand(1:N-1)
        end
        if minimum(
                     ∑( fourier_coeff(b,x) ⋅ f[1+ x]   for x in 0:N-1 )
                     for b ∈ B
                  ) ≥ 1-ε − eps
            count +=1
        end
    end
    return count / repet
end

function print_stats(;n,δ,ε,opt,fid,p)
    N = 2^n
    println("n                                                 ", n)
    println("δ                                                 ", δ)
    println("ε                                                 ", ε)
    println("Deterministic worst case fidelity for that ε:     ", 1-n⋅ε/2)
    println("Optimization termination status:                  ", opt.status)
    println("Optimization primal bound:                        ", opt.primal_bound)
    println("Optimization dual bound:                          ", opt.dual_bound)
    println("Fidelity:                                         ", opt.f[0]  , " (compare to $fid)")
    println("Probability that a random query is big (rounded): ", ∑(round.(opt.big))/(N-1), " (compare to $(p!=nothing ? ceil(p⋅(N-1))/(N-1) : p))")
    println("Probability that a random query is big:           ", ∑(opt.big)/(N-1))

    f = ℝ[ opt.f[x] for x=0:N-1 ]
    𝐸 = estimate_bigprb( f
                         ; n,ε,
                         repet = 1_000_000 )
    println("Probability that a random query is big, ESTIMATE: ", 𝐸)

    println("Probability of Accept:                            ", ( ∑(opt.big)/(N-1) )^n )
    println("Probability of Accept (primal bound):             ", ( ∑(opt.primal_bound)/(N-1) )^n )
    println("Probability of Accept (dual bound):               ", ( ∑(opt.dual_bound)/(N-1) )^n )

    𝐸 = estimate_acceptprb( f
                            ; n,ε,
                            repet = 1_000_000 )
    println("Probability of Accept, ESTIMATE:                  ", 𝐸 )

    println("Probability of Reject:                            ", 1-( ∑(opt.big)/(N-1) )^n )
end

# ==============================================================================================================
# § Test Against Random *Bases*
# ==============================================================================================================

function _rand_fullrk_mod2mtx(n :: ℤ) ::Matrix{Mod{2}}

    M ::Matrix{Mod{2}} = zeros(Mod{2},n,n)

    while rankx(M) < n
        M = rand(Mod{2},n,n)
    end
    return M
end

function _sample_basis!(list  ::Vector{Int64}
                        ;
                        n     ::ℤ             ) ::Vector{Int64}
    ;                                                            ; @assert 1 ≤ n < 63
    N    = Int64(2)^n
    oh   = Int64(0)
    one  = Int64(1)

    resize!(list, n)

    M = _rand_fullrk_mod2mtx(n)

    for ℓ = 1:n
        x = oh;
        for k = 1:n
            if M[k,ℓ] == 1
                x |=  ( 1<<(k-1) )
            end
        end
        list[ℓ] = x
    end

    return list
end

function _rank(L::Vector{ℤ}; n) ::ℤ
    m = length(L)
    if m==0   return 0   end

    M ::Matrix{Mod{2}} = zeros(Mod{2},n,m)
    for ℓ = 1:m
        y = L[ℓ]
        for k = 1:n
            M[k,ℓ] = ( y >> (k-1) )&1
        end
    end
    return rankx(M)
end


function estimate_basis_acceptprb(f :: Vector{ℝ},
                                  ;
                                  n     ::ℤ,
                                  ε     ::ℝ,
                                  repet ::Int,
                                  eps   ::ℝ = 1e-9) ::ℝ

    N = 2^n
    @assert size(f) == (N,)
    @assert 0 ≤ ε ≤ 1

    B ::Vector{Int64} = Int64[]

    count = 0
    for _r = 1:repet
        _sample_basis!(B;n)
        if minimum(
                     ∑( fourier_coeff(b,x) ⋅ f[1+ x]   for x in 0:N-1 )
                     for b ∈ B
                  ) ≥ 1-ε − eps
            count +=1
        end
    end
    return count / repet
end

function false_positive_analysis( ;
                                  n        = 2,
                                  δ        = 0.1, # ≤ 2/n
                                  ε        = 3δ,
                                  fid      = 1-5δ, # nothing
                                  p        = nothing,
                                  time_limit = 1min,
                                  skip_basis_stuff = true)
    scenario = FALSE_POSITIVE
    N   = 2^n

    println("\n👉𝐋𝐄𝐓'𝐒 𝐆𝐎!!!")
    println("====================================================================================")

    opt = mip_stabilizer_fooler(; n,
                                scenario,
                                ε, p, fid,
                                time_limit)

    print("\n")
    println("------------------------------------------------------------------------------------")
    print_stats(;n,δ,ε,opt,fid,p)

    if ! skip_basis_stuff
        f = ℝ[ opt.f[x] for x=0:N-1 ]
        𝐸 = estimate_basis_acceptprb( f
                                      ; n,ε,
                                      repet = 1_000_000 )
        println("Probability of Basis-Accept, ESTIMATE:            ", 𝐸 )
    else
        println("Probability of Basis-Accept, ESTIMATE:            ", "(skipped)")
    end

    return opt
end

function false_negative_analysis( ;
                                  n        = 2,
                                  δ        = 0.1, # ≤ 2/n
                                  ε        = 3δ,
                                  fid      = 1-δ,
                                  p        = nothing,
                                  time_limit = 1min,
                                  skip_basis_stuff = true)
    scenario = FALSE_NEGATIVE
    N   = 2^n

    println("\n👉𝐋𝐄𝐓'𝐒 𝐆𝐎!!!")
    println("====================================================================================")

    opt = mip_stabilizer_fooler(; n,
                                scenario,
                                ε, p, fid,
                                time_limit =1min)

    print("\n")
    println("------------------------------------------------------------------------------------")
    print_stats(;n,δ,ε,opt,fid,p)

    if ! skip_basis_stuff
        f = ℝ[ opt.f[x] for x=0:N-1 ]
        𝐸 = estimate_basis_acceptprb( f
                                      ; n,ε,
                                      repet = 1_000_000 )
        println("Probability of Basis-Accept, ESTIMATE:            ", 𝐸 )
    else
        println("Probability of Basis-Accept, ESTIMATE:            ", "(skipped)")
    end

    return opt
end

# ===============================================================================================
#-------------------------------------------------------------------
# reduced stuff
#-------------------------------------------------------------------


const RndHypergr_t  = @NamedTuple{w::ℤ, p::ℝ}
const Fix_Stuff_t   = Union{ Bool ,  ℤ ,  RndHypergr_t }


function reduced_fooler(; n::Int, α::Float64,
                        fix_Id     ::Fix_Stuff_t     = false,
                        perturb    ::ℤ               = -1,
                        time_limit                   = 10min )
    @assert n ≥ 1
    @assert α ≥ 0

    N = 2^n

    model = Model(HiGHS.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "time_limit", time_limit)

    @variable(model, p[1:N-1] ≥ 0)
    @variable(model, acceptrow[1:N-1], Bin)


    @constraint(model, ∑(p[x] for x in 1:N-1) == 1)                   # Probability distribution

    c =
        if perturb ≥ 0
            [ rand( perturb : perturb+1 )    for y = 1 : N-1 ]
        else
            [ 1                              for y = 1 : N-1 ]
        end
    @objective(model, Max, ∑( c[y]⋅acceptrow[y]  for y in 1:N-1))

    for y in 1:N-1
        @constraint(model,
                    ∑((count_ones(y & x) % 2) * p[x] for x in 1:N-1)
                    ≤ 1 - acceptrow[y] * (1 - α / 2)
                   )
    end

    function fix_this(S)
        s = ∑(  1<<(j-1) for j ∈ S  )
        @constraint(model, acceptrow[ s ] == 1)
    end


    if fix_Id isa Bool || fix_Id isa ℤ
        for w = 1:fix_Id                                              # yes, this works even if fix_Id is a bool
            for S in Combinatorics.combinations(1:n,w)
                fix_this(S)
            end
        end
    else
        @assert fix_Id isa RndHypergr_t   "Why are you here?!!"
        @assert fix_Id.w ≥ 2
        @assert 0 < fix_Id.p ≤ 1
        let w=1
            for S in Combinatorics.combinations(1:n,w)
                fix_this(S)
            end
        end
        let w = fix_Id.w,
            p = fix_Id.p
            for S in Combinatorics.combinations(1:n,w)
                if rand() ≤ p
                    fix_this(S)
                end
            end
        end
    end

    optimize!(model)

    if termination_status(model) != MOI.OPTIMAL
        @show termination_status(model)
    end

    return (status       = termination_status(model),
            primal_bound = objective_value(model),
            dual_bound   = objective_bound(model),
            p            = JuMP.value.(p),
            acceptrow    = JuMP.value.(acceptrow))
end

function slack_distrib(p; n, α)
    @assert n ≥ 1
    N = 2^n
    @assert (N - 1,) == size(p)
    @assert 0 ≤ α

    slacks = ℝ[]  ; sizehint!(slacks,N-1)

    for y = 1 : N-1
        push!(slacks,

              α/2 − ∑((count_ones(y & x) % 2) * p[x] for x in 1:N-1)

             )
    end
    sort!(slacks)
    return slacks
end

function estimate_basis_reduced(p; n, α, repet, eps=1e-10, verbose::Bool=false)
    @assert n ≥ 1
    N = 2^n
    @assert (N - 1,) == size(p)
    @assert 0 ≤ α

    B = Int64[]

    count = 0
    for _r = 1:repet
        _sample_basis!(B; n)
        m = maximum(
                   ∑((count_ones(b & x) % 2) * p[x] for x in 1:N-1)
                   for b in B
               )
        if m ≤ α / 2 +eps
            count += 1
        end
        if verbose && m > α/2 && m ≈ α/2
            @show α/2-m
        end
    end
    return count / repet
end

numbas(dim::Int) = begin 𝟐=BigInt(2); prod( (𝟐^dim-𝟐^j) for j=0:dim-1 ) end

function basdens(𝐸::ℝ,a::ℤ;n)
    res=𝐸⋅numbas(n)
    for j=0:n-1
        res /= (a-j)
    end
    return res
end

function basdens(acceptrows ::Vector{ℤ} ; n, repet) ::ℝ
    N = 2^n
    L = length(acceptrows)
    @assert repet ≥ 1
    n ≤ L || return 0.0

    count = 0
    for _r = 1:repet
        B = sample(1:L, n; replace=false)
        for i=1:n
            B[i] = acceptrows[B[i]]
        end
        if _rank(B;n) == n
            count +=1
        end
    end
    return count / repet
end

function reduced_false_positive_analysis(;
                                         n  :: ℤ = 2,
                                         α  :: ℝ = 1.0,
                                         time_limit        ::Float64 = 1min,
                                         basis_repeats     ::Int     = 0,
                                         perturb           ::ℤ       = -1,
                                         fix_Id            ::Fix_Stuff_t  = ( α ≥ 2/n) )
    @assert n ≥ 1
    N = 2^n

    println("\n👉𝐋𝐄𝐓'𝐒 𝐆𝐎 (reduced)!!!")
    println("====================================================================================")

    opt = reduced_fooler(;n, α, time_limit, perturb, fix_Id)

    num_ambigrows = count(0.01 < opt.acceptrow[y] < 0.99 for y=1:N-1)
    acceptrowlist = [ y for y=1:N-1 if opt.acceptrow[y] > 0.5 ]
    num_accept    = length(acceptrowlist)

    print("\n")
    println("------------------------------------------------------------------------------------")
    println("n                                       ", n)
    println("α                                       ", α)
    println("Optimization termination status:        ", opt.status)
    println("Optimization primal bound:              ", opt.primal_bound)
    println("Optimization dual bound:                ", opt.dual_bound)
    println("Number of ambiguous rows:               ", num_ambigrows)
    println("Number of selected rows:                ", ∑(round.(opt.acceptrow)))
    println("Fraction of selected rows:              ", ∑(opt.acceptrow) / (N - 1))
    println("Fraction of selected rows (rounded):    ", ∑(round.(opt.acceptrow)) / (N - 1))
    println("Number of accepting rows:               ", num_accept)
    println("Rank of set of accepting rows:          ", _rank(acceptrowlist;n))

    if basis_repeats > 0
        𝐸 = estimate_basis_reduced(opt.p; n, α, repet=basis_repeats)
        println("Probability of Basis-Accept, ESTIMATE:            ", 𝐸)
        println("Basis density in accepting rows from accept prob: ", Float64( basdens(𝐸,num_accept;n) ))
        println("Basis density in accepting rows:                  ", Float64( basdens(acceptrowlist;n,repet=1_000_000) ))
    else
        println("Probability of Basis-Accept, ESTIMATE:            ", "(skipped)")
    end

    if 1+α ≈ 1+2/n
        println("WARNING: α ≈ 2/n.  The Basis-Accept prob should be $(factorial(n)/numbas(n))")
    elseif α < 2/n
        println("WARNING: α < 2/n.  The Basis-Accept prob should be 0.")
    elseif α > 2/n
        println("α > 2/n: The Basis-Accept prob should be at least $(Float64(factorial(n)/numbas(n)))")
    end

    return ( opt    = opt,
             slacks = slack_distrib(opt.p;n,α) )
end


end #^ module Stabilizer_Foolers
