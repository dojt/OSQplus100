# CompVeris/src/measure.jl
# (c) Dirk Oliver Theis
# Open source license: CC0

"""
Module `Measure`

## Exports:
 * Function `measure!` — Simulates measurements of Pauli-Z operators on a quantum state.
 * Function `sample` —  Samples a measurement outcome (±1) based on a given expectation value.

## Not exported, but useful:
 * Function `reset()` — resets the RNG
"""
module Measure
    export measure!, sample, sample_many

    using DOT_NiceMath
    using DOT_NiceMath.Numbers64 # ℤ=Int64, ℝ=Float64, ℂ=ComplexF64, ⋅=*
    using Random
    using Distributions

    import ..fourier_coeff
    using ..Density_Operators: ℂHermitian_t
    const Operator_or_Diagonal_t = Union{ ℂHermitian_t , Vector{ℝ} }


    const THE_MEASURE_SEED = 91827
    const MY_RNG_r         = Ref(  Xoshiro(THE_MEASURE_SEED)  )

    function reset()
        MY_RNG_r[] = Xoshiro(THE_MEASURE_SEED)
    end

    function _binomial(rng, m::Int64, p::Float64) ::Int64
        dist = Binomial(m, clamp(p,0.0,1.0))
        return rand(rng, dist)
    end

    sample(                𝜇 ::ℝ) ::Int8  = let 𝑝₊ = (1+𝜇)/2 ; Int8( 2⋅(rand(MY_RNG_r[]) ≤ 𝑝₊ ? 1 : 0) −1 ) end
    sample_many(m ::Int64, 𝜇 ::ℝ) ::Int64 = let 𝑝₊ = (1+𝜇)/2 ;       2⋅_binomial(MY_RNG_r[], m,𝑝₊)     −m   end

    #
    # Main Work
    #

    make_z(a::Integer ;  n) ::Vector{ℝ}      = ℝ[ fourier_coeff(UInt64(x),UInt64(a))   for x = UInt64(0) : UInt64(2^n-1) ]
    make_z!(vec::Vector{ℝ}, a::Integer;  n)  = begin # fast version
        for x = UInt(0) : UInt(2^n-1)
            vec[1+ x] = fourier_coeff(x,UInt64(a))
        end
        nothing;
    end

   expval(a::UInt64, ϱ ::Vector{ℝ}) ::ℝ =
       begin
           N = UInt64( length(ϱ) )
           ∑( fourier_coeff(x,a)⋅ϱ[1+ x]    for x = UInt64(0) : N-1 )
       end


    """
    Function
        measure!(res_list::Vector{T}; n::ℤ, ϱ::ℂHermitian_t, Z_list::Vector{UInt64}) where {T <: Union{Int8, ℝ}}

    Simulates measurements of Pauli-Z operators on a quantum state represented by a density matrix.

    This function takes a density matrix `ϱ` and a list of Pauli-Z operators `Z_list` to measure.
    It simulates the measurements and stores the results in `res_list`.

    # Arguments
    - `res_list`: A pre-allocated vector to store the measurement results. 
                   The type of the elements determines the type of measurement:
                   - `Vector{Int8}`:  Simulates projective measurements and stores ±1 for each shot.
                   - `Vector{ℝ}`:     Calculates and stores the expectation values.
    - `n`: An integer representing the number of qubits.
    - `ϱ`: The density matrix representing the quantum state.
    - `Z_list`: A vector of integers, each representing a Pauli-Z operator to measure.
                Each integer is a bitstring, where each bit corresponds to a qubit,
                with 1 indicating a Z operator on that qubit and 0 indicating the identity.

    # Details
    The function uses the following approach:
    1. For each Pauli-Z operator in `Z_list`:
        - Calculate the expectation value of the operator in the state `ϱ`.
        - If `res_list` is a `Vector{Int8}`:
            - Simulate a projective measurement by drawing a random number and comparing it to the expectation value.
            - Store +1 or -1 in `res_list` based on the measurement outcome.
        - If `res_list` is a `Vector{ℝ}`:
            - Store the calculated expectation value directly in `res_list`.

    # Example Usage
    ```julia
    using Density_Operators, Measure

    n = 2  # Number of qubits
    ϱ = Density_Operators.make_ϱ!(Hermitian(Zeros(2^n,2^n)),0b00;n=n) # |00⟩
    res_list = Int8[]
    Measure.measure!(res_list; n=n, ϱ=ϱ, Z_list=[0b01, 0b10])  # Measures Z on qubit 1 and 2
    println(res_list)  # Example output: [1, 1]

    res_list = Float64[]
    Measure.measure!(res_list; n=n, ϱ=ϱ, Z_list=[0b01, 0b10])  # Measures Z on qubit 1 and 2
    println(res_list)  # Example output: [1.0, 1.0]
    ```
    """
    function measure!(res_list :: Vector{Int8},
                      ϱ        :: ℂHermitian_t
                      ;
                      n        :: ℤ,
                      Z_list   :: Vector{UInt64})

        N   = 2^n ;      @assert 1 ≤ n ≤ 24

        empty!(   res_list)
        sizehint!(res_list, length(Z_list))

        Z = make_z(UInt64(0);n)
        for bit_Z in Z_list
            make_z!(Z,bit_Z;n)
            expval = ∑( Z[x]⋅ϱ[x,x]  for x=1:N ) |> ℜ
            res    = sample(expval)
            push!(res_list, res )
        end
    end #^ measure!()

    function measure!(res_list :: Vector{Int8},
                      ϱ        :: Vector{ℝ}
                      ;
                      n        :: ℤ,
                      Z_list   :: Vector{UInt64})

        N   = 2^n ;      @assert 1 ≤ n ≤ 24

        empty!(   res_list)
        sizehint!(res_list, length(Z_list))

        for bit_Z in Z_list
            res    = sample( expval(bit_Z,ϱ) )
            push!(res_list, res )
        end
    end #^ measure!()

    function measure!(res_list :: Vector{ℝ},
                      ϱ        :: ℂHermitian_t
                      ;
                      n        :: ℤ,
                      Z_list   :: Vector{UInt64})

        N   = 2^n ;      @assert 1 ≤ n ≤ 24

        empty!(   res_list)
        sizehint!(res_list, length(Z_list))

        Z = make_z(UInt64(0);n)
        for bit_Z in Z_list
            make_z!(Z,bit_Z;n)
            expval = ∑( Z[x]⋅ℜ( ϱ[x,x] )  for x=1:N )
            push!(res_list, expval )
        end
    end #^ measure!()

    function measure!(res_list :: Vector{ℝ},
                      ϱ        :: Vector{ℝ}
                      ;
                      n        :: ℤ,
                      Z_list   :: Vector{UInt64})

        N   = 2^n ;      @assert 1 ≤ n ≤ 24

        empty!(   res_list)
        sizehint!(res_list, length(Z_list))

        for bit_Z in Z_list
            push!(res_list, expval(bit_Z,ϱ) )
        end
    end #^ measure!()

end #^ module Measure
#EOF
