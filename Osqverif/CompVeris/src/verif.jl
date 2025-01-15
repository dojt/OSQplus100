# CompVeris/src/verif.jl
# (c) Dirk Oliver Theis
# Open source license: CC0

"""
Module `Verif`

## Exports:
 * Enum-type `Decision_t` — can hold the value `Accept` or `Reject`.
 * Sub-module `DFE` — Direct Fidelity Estimation
 * Sub-module `MYNS` — My New Shit


### The `verify` function
The `verify` function in each sub-module returns a named tuple consisting of the
`decision` (`Accept` or `Reject`) as the first entry, and a `witness` as the
second.  The type of `witness` depends on the method.


All `verify` functions have some keyword arguments in common:
```
    function verify( ;
                     n                         :: ℤ,
                     δ                         :: ℝ,
                     ϱ                         :: ℂHermitian_t,
                     #= more keyword args =#                       )
```
* `n` — number of qubits
* `ϱ` — output of the circuit
* `δ` — expected level of noise in ϱ

## Not exported, but useful in the sub-modules:

The submodules don't export anything, but each contains this function:

 * Function `reset()` — resets the RNG

E.g., `DFE.reset()` resets the RNG local to the sub-module `DFE`.

### ❗ 👉 The `verify` function in the MYNS sub-module...

... currently returns a `Dict` of named tuples (assigned to strings) — not just
a single named tuple.  The reason is that I'm trying several methods and want to
compare them on the same randomness.
"""
module Verif
export Decision_t, DFE, MYNS

@enum Decision_t :: Bool begin
    Reject = false
    Accept = true
end


using DOT_NiceMath
using DOT_NiceMath.Numbers64 # ℤ=Int64, ℝ=Float64, ℂ=ComplexF64, ⋅=*

import ..Measure
using  ..Measure: Operator_or_Diagonal_t

const THE_VERIF_SEED = 5678

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Direct Fidelity Estimation
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


"""
Module `DFE`

Implements Direct Fidelity Estimation, stabilizer state variant.
"""
module DFE
    using Random

    using DOT_NiceMath
    using DOT_NiceMath.Numbers64 # ℤ=Int64, ℝ=Float64, ℂ=ComplexF64, ⋅=*

    import ..THE_VERIF_SEED
    import ..Operator_or_Diagonal_t
    import ..Decision_t, ..Accept, ..Reject

    using ..Measure


    const MY_RNG_r  = Ref(  Xoshiro(THE_VERIF_SEED)  )

    function reset()
        MY_RNG_r[] = Xoshiro(THE_VERIF_SEED)
    end

    const Ret_t = @NamedTuple{ decision::Decision_t, witness::ℝ }

    function verify( ;
                     n                  :: ℤ,
                     δ                  :: ℝ,
                     ϱ                  :: Operator_or_Diagonal_t,
                     num_stabs          :: ℤ = n,
                     num_shots_per_stab :: ℤ = 1024              )  :: Ret_t
        rng     = MY_RNG_r[]
        N       = 2^n

        expvals = ℝ[]
        Z_list  = rand(rng, UInt64, num_stabs)
        measure!(expvals,ϱ ; n, Z_list)

        # the estimate:
        Σ = ∑(
                sample_many(num_shots_per_stab, expvals[stab] )
                for stab = 1:num_stabs
             )

        𝐸 = Σ / (num_stabs⋅num_shots_per_stab)

        return (decision = ( 𝐸 < 1-3δ ? Reject : Accept ) ,
                witness  = 𝐸                                )
    end
end #^ module DFE

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# My New Shit
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module MYNS
    using Random

    using DOT_NiceMath
    using DOT_NiceMath.Numbers64 # ℤ=Int64, ℝ=Float64, ℂ=ComplexF64, ⋅=*

    import ..THE_VERIF_SEED
    import ..Operator_or_Diagonal_t
    import ..Decision_t, ..Accept, ..Reject

    using ..Measure

    using LinearAlgebraX: rankx, nullspacex, detx
    using Statistics: mean
    using Mods

    const MY_RNG_r  = Ref(  Xoshiro(THE_VERIF_SEED)  )

    function reset()
        MY_RNG_r[] = Xoshiro(THE_VERIF_SEED)
    end

    #---------------------------------------------------------------------------------------------------------
    #
    # Helpers
    #
    #---------------------------------------------------------------------------------------------------------

    function _rand_fullrk_mod2mtx(;n :: ℤ) ::Matrix{Mod{2}}

        M ::Matrix{Mod{2}} = zeros(Mod{2},n,n)

        while rankx(M) < n
            M = rand(Mod{2},n,n)
        end
        return M
    end

    function _sample_basis!(list  ::Vector{UInt64}
                            ;
                            n     ::ℤ             ) ::Vector{UInt64}
        ;                        ; @assert 1 ≤ n < 63
        N    = UInt64(2)^n
        oh   = UInt64(0)
        one  = UInt64(1)

        resize!(list, n)

        M = _rand_fullrk_mod2mtx(;n)

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

    #---------------------------------------------------------------------------------------------------------
    #
    # Main work: `verify()`
    #
    #---------------------------------------------------------------------------------------------------------

    const Ret_t        = @NamedTuple{ decision::Decision_t, witness::ℝ, expval ::ℝ }
    const Actual_Ret_t = Dict{String,Ret_t}

    # const STORAGE = Vector{Int8}(undef, 1)

    function verify( ;
                     n                  :: ℤ,
                     δ                  :: ℝ,
                     ϱ                  :: Operator_or_Diagonal_t,
                     num_shots_per_stab :: ℤ  =1024    ) ::Actual_Ret_t

        rng     = MY_RNG_r[]
        N       = 2^n

        expvals = ℝ[]
        Z_list  = zeros(UInt64, n)
        _sample_basis!(Z_list;n)

        measure!(expvals,ϱ ; n, Z_list)

        # for b = 1:n
        #     for _shot = 1:num_shots_per_stab
        #         M = sample( expvals[stab] )
        #         # What should I do with it?!??
        #     end
        # end

        𝒮 = Vector{Int64}(undef, n)
        for stab = 1:n
            𝑓 = expvals[stab]
            𝒮[stab] = sample_many(num_shots_per_stab, 𝑓 )
        end

        return_this = Actual_Ret_t()

        let 𝐸 = ∑(
                     𝒮[stab]
                     for stab = 1:n
                 ) / (n⋅num_shots_per_stab)

            return_this["Average-Everything"] =
                (  decision = ( 𝐸 < 1-3δ ? Reject : Accept ) ,
                   witness  = 𝐸,
                   expval   = mean(expvals)                   )
        end
        let 𝑚𝑖𝑛 = minimum(  𝒮[stab]/num_shots_per_stab
                            for stab=1:n
                         )

            return_this["Min-of-Means"] =
                (  decision = ( 𝑚𝑖𝑛 < 1-3δ ? Reject : Accept ) ,
                   witness  = 𝑚𝑖𝑛,
                   expval   = minimum(expvals)                 )
        end

        𝒮 = nothing
        return return_this
    end

end #^ module MYNS


end #^ module Verif
        # 𝒮 = Vector{ℝ}(undef, n)
        # for stab = 1:n                                  # column major
        #     𝑓 = expvals[stab]
        #     𝒮[stab] = ∑( Int64( sample(𝑓) )  for _shot = 1:num_shots_per_stab  ) / num_shots_per_stab
        # end
        #
        # return_this = Actual_Ret_t()
        #
        # let 𝐸 = mean( 𝒮 )
        #
        #     return_this["Average-Everything"] =
        #         (  decision = ( 𝐸 < 1-3δ ? Reject : Accept ) ,
        #            witness  = 𝐸                              )
        # end
        # let 𝑚𝑖𝑛 = minimum( 𝒮 )
        #
        #     return_this["Min-of-Means"] =
        #         (  decision = ( 𝑚𝑖𝑛 < 1-3δ ? Reject : Accept ) ,
        #            witness  = 𝑚𝑖𝑛                              )
        # end
#EOF
