# CompVeris/src/density_operators.jl
# (c) Dirk Oliver Theis
# Open source license: CC0
module Density_Operators

export make_sdp, make_ϱ

using DOT_NiceMath
using DOT_NiceMath.Numbers64

using LinearAlgebra: Hermitian, tr, issymmetric

using JuMP


const ℂHermitian_t = Hermitian{ℂ,Matrix{ℂ}}
const ℝSymmetric_t = Hermitian{ℝ,Matrix{ℝ}}

"""
Type const `Constraint_t = @NamedTuple{ LHS :: ℂHermitian_t, rhs ::ℝ }`

Type for a constraint, consisting of a Hermitian LHS matrix and a real scalar RHS.
"""
const Constraint_t = @NamedTuple{ LHS :: ℂHermitian_t, rhs ::ℝ }

"""
Function
    make_sdp(optimizer; n, C, equations=Constraint_t[], greaterthans=Constraint_t[])

Creates a JuMP optimization model for optimizing over density operators.

This function takes a complex Hermitian matrix `C` representing the objective function
and lists of complex Hermitian constraints to construct a JuMP optimization model.
It converts the complex-valued problem into a real-valued SDP that can be solved
using JuMP's SDP solvers.

# Arguments
- `optimizer`: The JuMP optimizer to use (e.g., `Mosek.Optimizer`).
- `n`: An integer representing the number of qubits (the dimension of the density matrix is 2^n).
- `C`: A complex Hermitian matrix representing the objective function.
- `equations`: A vector of `Constraint_t` tuples representing equality constraints LHS⋅ϱ == rhs.
- `greaterthans`: A vector of `Constraint_t` tuples representing inequality constraints LHS⋅ϱ ≥ rhs.

# Returns
- A named tuple `(model, cplx)` where:
  - `model`: The JuMP optimization model.
  - `cplx`: A function to retrieve complex entries from the optimized real-valued matrix.

# Example Usage
```julia
using Density_Operators
using JuMP, MosekTools

# Define problem parameters
n = 2  # Number of qubits
C = Hermitian([1 0; 0 -1])  # Example objective function

# Define constraints (optional)
constraint1 = (LHS = Hermitian([0 1; 1 0]), rhs = 0.5)
equations = [constraint1]

# Create the JuMP model
model, cplx = make_sdp(Mosek.Optimizer; n=n, C=C, equations=equations)

# Solve the optimization problem
optimize!(model)

# Retrieve the optimized density matrix
ρ_real = value.(model[:ϱ])
ρ_complex = [cplx(ρ_real, k, ℓ) for k in 1:2^n, ℓ in 1:2^n]
```
"""
function make_sdp( optimizer
                   ;
                   n            :: ℤ,
                   C            :: ℂHermitian_t,
                   equations    :: Vector{ Constraint_t } = Constraint_t[],
                   greaterthans :: Vector{ Constraint_t } = Constraint_t[]  ) ::@NamedTuple{model::Model,cplx::Function}

    N = 2^n ;                    @assert 1 ≤ n ≤ 24
    ;                            @assert size(C) == (N,N)
    M = 2N

    function cplx(ϱ, k::I,ℓ::I) ::ℂ where{I <: Integer}
        @assert (M,M)==size(ϱ) ; @assert 1≤k≤N ;@assert 1≤ℓ≤N
        return ϱ[2k-1,2ℓ-1] + 𝒊⋅ϱ[2k,2ℓ-1]
    end

    # Convert complex matrices to real matrices
    C_real = zeros(ℝ, M, M)
    for i in 1:N, j in 1:N
        C_real[2i-1:2i, 2j-1:2j] = [real(C[i,j]) -imag(C[i,j]); imag(C[i,j]) real(C[i,j])]
    end

    @assert size(C) == (N, N)
    @assert all(size(eq.LHS) == (N, N) for eq in equations)
    @assert all(size(gt.LHS) == (N, N) for gt in greaterthans)

    # Store real matrices directly
    equations_real = Vector{Tuple{ℝSymmetric_t, ℝ}}() 
    sizehint!(equations_real, length(equations))
    for eq in equations
        LHS_real = zeros(ℝ, M, M)
        for i in 1:N, j in 1:N
            LHS_real[2i-1:2i, 2j-1:2j] = [real(eq.LHS[i,j]) -imag(eq.LHS[i,j]); imag(eq.LHS[i,j]) real(eq.LHS[i,j])]
        end
        @assert issymmetric(LHS_real) "BUG!! LHS_real is not symmetric!"
        push!(equations_real, (Hermitian(LHS_real), 2⋅eq.rhs)) 
    end

    greaterthans_real = Vector{Tuple{ℝSymmetric_t, ℝ}}()
    sizehint!(greaterthans_real, length(greaterthans))
    for gt in greaterthans
        LHS_real = zeros(ℝ, M, M)
        for i in 1:N, j in 1:N
            LHS_real[2i-1:2i, 2j-1:2j] = [real(gt.LHS[i,j]) -imag(gt.LHS[i,j]); imag(gt.LHS[i,j]) real(gt.LHS[i,j])]
        end
        @assert issymmetric(LHS_real) "BUG!! LHS_real is not symmetric!"
        push!(greaterthans_real, (Hermitian(LHS_real), 2⋅gt.rhs))
    end

    model = Model(optimizer)

    @variable(  model, ϱ[1:M,1:M], PSD)
    @constraint(model, tr(ϱ) == 2) # Trace should be 2 for the real-valued matrix

    # Enforce shapes of 2x2 blocks:         x -y         for z = x+iy
    #                                       y  x
    for k = 1:N
        for ℓ = 1:N
            @constraint(model, ϱ[ 2k-1, 2ℓ-1 ] ==  ϱ[ 2k  , 2ℓ ])
            @constraint(model, ϱ[ 2k  , 2ℓ-1 ] == −ϱ[ 2k-1, 2ℓ ])
        end
    end

    # Use the real-valued matrices in the objective and constraints
    @objective( model, Min, tr(C_real ⋅ ϱ)) # Objective function should be a scalar, use trace

    for (LHS_real, rhs) in equations_real
        @constraint(model, tr(LHS_real ⋅ ϱ) == rhs)
    end

    for (LHS_real, rhs) in greaterthans_real
        @constraint(model, tr(LHS_real ⋅ ϱ) ≥ rhs) 
    end

    return (model=model, cplx=cplx)
end


"""
Function
    make_ϱ!(ϱ ::Matrix{ℂ},  x₀::T  ;  n::ℤ)    where T<:Unsigned

Create a density operator for the pure state |x₀⟩:
the a computational basis state with a single non-zero entry at position `x₀`.
"""
function make_ϱ!(ϱ ::ℂHermitian_t,  x₀::T  ;  n::ℤ) :: ℂHermitian_t    where T<:Unsigned
    N = 2^n
    ϱ.data .= zeros(ℂ, N,N)
    ϱ.data[1+ x₀,1+ x₀] = 1
    return ϱ
end

end #^ module Density_Operators
# EOF
