# CompVeris/src/instances.jl
# (c) Dirk Oliver Theis
# Open source license: CC0

"""
Module `Instances`

## Exports:
 * Function `make_rnd_dens_op` — Creates a random density operator close to a given pure state.

## Not exported, but useful:
 * Function `reset()` — resets the RNG
"""
module Instances
    export make_rnd_dens_op

    using DOT_NiceMath
    using DOT_NiceMath.Numbers64 # ℤ=Int64, ℝ=Float64, ℂ=ComplexF64, ⋅=*
    using Random

    using JuMP
    using CSDP

    using LinearAlgebra: Hermitian

    using ..Density_Operators: ℂHermitian_t, make_sdp
    import ..fourier_coeff

    const THE_INST_SEED = 1234
    const MY_RNG_r      = Ref(  Xoshiro(THE_INST_SEED)  )

    function reset()
        MY_RNG_r[] = Xoshiro(THE_INST_SEED)
    end

    # =============================================================================================================
    # Random density operator
    # =============================================================================================================

    function make_rnd_dens_op(;n::ℤ, ψ₀ ::Vector{ℂ}, ε ::ℝ) ::ℂHermitian_t

        N = 2^n ;                    @assert 1 ≤ n ≤ 24
        ;                            @assert size(ψ₀) == (N,)
        @assert ε ≥ 0

        # Constraint: Ensure closeness to the pure state
        LHS_constraint = Hermitian(ψ₀ ⋅ ψ₀')
        constraint = (LHS = LHS_constraint, rhs = 1 - ε)

        # Random Hermitian cost matrix
        C = Hermitian(randn(MY_RNG_r[], ℂ, N, N))

        # Use make_sdp to find the optimal density operator
        My_Optimizer = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)
        model, cplx = make_sdp(My_Optimizer; n, C, equations=[constraint])

        optimize!(model)

        # Extract the complex density matrix
        ϱ_real = value.(model[:ϱ])
        ϱ_complex = [cplx(ϱ_real, k, ℓ) for k in 1:N, ℓ in 1:N]
        return Hermitian(ϱ_complex)
    end #^ make_rnd_dens_op()

    # =============================================================================================================
    # Random diagonal
    # =============================================================================================================

    abstract type   Kind_Of_Rnd_Dens_Diag_t end
    struct Dumb  <: Kind_Of_Rnd_Dens_Diag_t end
    struct Smart <: Kind_Of_Rnd_Dens_Diag_t end


    function make_rnd_dens_diag(::Dumb ; n::ℤ, ε ::ℝ) ::Vector{ℝ}

        N = 2^n ;                    @assert 2 ≤ n ≤ 24
        @assert ε ≥ 0

        ϱ_diag    =   rand(MY_RNG_r[], ℝ,N)
        denom     =   ∑( ϱ_diag[2:N] )
        ϱ_diag   .*= ε/denom
        ϱ_diag[1] =  1-ε

        return ϱ_diag

    end #^ make_rnd_dens_diag(::Dumb)


    using HiGHS
    const THE_OPTIMIZER = HiGHS.Optimizer

    function _lp_worst_case!(ϱ_diag; n::ℤ, ε ::ℝ, α ::ℝ, σ::ℝ = 0.1)
        N = 2^n

        model = Model(THE_OPTIMIZER)
        set_silent(model)

        @variable(  model, ϱ[1:N] ≥ 0)

        @objective( model, Min, σ ⋅ ∑(randn()⋅ϱ[1+ x ]  for x=1:N-1 ) )

        @constraint(model, ∑(ϱ)      == 1)
        @constraint(model, ϱ[1+ 0 ]  ≤  1-ε)

        for j = 0:n-1
            @constraint(model,
                        ∑( fourier_coeff(x,2^j) ⋅ ϱ[1 + x]  for x=0:N-1 )        ≥       1 - α⋅ε
                        )
        end

        optimize!(model)

        if termination_status(model) != MOI.OPTIMAL
            @show termination_status(model)
        end

        resize!(ϱ_diag, N)
        ϱ_diag .= value.(ϱ)
        ;
    end

    function make_rnd_dens_diag(::Smart ; n::ℤ, ε ::ℝ, α ::ℝ = NaN) ::Vector{ℝ}

        N = 2^n ;                    @assert 2 ≤ n ≤ 24
        @assert ε ≥ 0
        @assert α !== NaN "\n❗ You need to give the `α` parameter for the 'Smart' version of"
                          "\n❗ `CompVeris.Instances.make_rnd_dens_diag()`.\n"
        @assert 0 < α < 1

        ϱ_diag = ℝ[]
        _lp_worst_case!(ϱ_diag ; n,ε,α)

        return ϱ_diag

    end #^ make_rnd_dens_diag(::Smart)

end #^ module Instances
#EOF
