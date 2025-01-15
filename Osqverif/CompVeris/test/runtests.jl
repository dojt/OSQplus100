# CompVeris/test/runtests.jl
# (c) Dirk Oliver Theis
# Open source license: CC0

using Test

using JuMP
using CSDP
using LinearAlgebra: tr, Hermitian

using Printf

using DOT_NiceMath
using DOT_NiceMath.Numbers64

a::ℝ ⪅ b::ℝ = ( a - 1e-4 ≤ b )

using CompVeris
using CompVeris.Verif: Decision_t, Accept, Reject
using CompVeris.Stabilizer_Foolers: FALSE_POSITIVE, FALSE_NEGATIVE
using CompVeris.Stabilizer_Foolers: fourier_coeff

@testset "Package CompVeris" verbose=true begin


@testset "Sub-Module `Density_Operators`" verbose=true begin
    @testset "make_sdp" verbose=true begin
        @testset "n=1" verbose=true begin
            n = 1
            C = Hermitian(ℂ[1 0; 0 -1])

            My_Optimizer = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

            let # Test case 1: No constraints
                model, cplx = Density_Operators.make_sdp(My_Optimizer; n=n, C)
                optimize!(model)
                ϱ_real = value.(model[:ϱ])
                ϱ_complex = [cplx(ϱ_real, k, ℓ) for k in 1:2^n, ℓ in 1:2^n]
                @test isapprox(tr(ϱ_complex), 1.0, atol=1e-5) # Check trace
                @test isapprox(ϱ_complex, ℂ[0 0; 0 1], atol=1e-5) # Check solution
            end

            let # Test case 2: With equality constraint
                constraint1 = (LHS = Hermitian(ℂ[0 1; 1 0]), rhs = 0.5)
                equations = [constraint1]
                model, cplx = Density_Operators.make_sdp(My_Optimizer; n=n, C, equations)
                optimize!(model)
                ϱ_real = value.(model[:ϱ])
                ϱ_complex = [cplx(ϱ_real, k, ℓ) for k in 1:2^n, ℓ in 1:2^n]
                @test isapprox(tr(ϱ_complex), 1.0, atol=1e-5) # Check trace
                @test isapprox(tr(constraint1.LHS ⋅ ϱ_complex), constraint1.rhs, atol=1e-5) # Check constraint
            end

            let # Test case 3: With inequality constraint
                constraint2 = (LHS = Hermitian(ℂ[0 1; 1 0]), rhs = 0.4)
                greaterthans = [constraint2]
                model, cplx = Density_Operators.make_sdp(My_Optimizer; n=n, C, greaterthans)
                optimize!(model)
                ϱ_real = value.(model[:ϱ])
                ϱ_complex = [cplx(ϱ_real, k, ℓ) for k in 1:2^n, ℓ in 1:2^n]
                @test isapprox(tr(ϱ_complex), 1.0, atol=1e-5) # Check trace
                # @test ℑ( tr(constraint2.LHS ⋅ ϱ_complex) ) +1 ≈ 1
                @test ℜ( tr(constraint2.LHS ⋅ ϱ_complex) ) ≥ constraint2.rhs - 1e-5 # Check constraint
            end
        end
        @testset "n=1:3" verbose=true begin
            My_Optimizer = optimizer_with_attributes(CSDP.Optimizer, MOI.Silent() => true)

            for n in 1:3  # Test for 1, 2, and 3 qubits
                N = 2^n
                # Example objective function: try to make a diagonal matrix
                C = Hermitian(ℂ[ k==ℓ ? 1 : 0  for k=1:N, ℓ=1:N ])

                let # Test case 1: No constraints
                    model, cplx = Density_Operators.make_sdp(My_Optimizer; n=n, C)
                    optimize!(model)
                    ϱ_real = value.(model[:ϱ])
                    ϱ_complex = [cplx(ϱ_real, k, ℓ) for k in 1:N, ℓ in 1:N]
                    @test isapprox(tr(ϱ_complex), 1.0, atol=1e-5) # Check trace

                    # For this C, the solution should be diagonal
                    @test all( isapprox(ϱ_complex[k,ℓ], 0, atol=1e-4) 
                                for k=1:N for ℓ=1:N if k≠ℓ ) 
                end

                let # Test case 2: With equality constraint
                    # Constraint:  <1|ϱ|1> = 0.5
                    constraint1 = (LHS = Hermitian(ℂ[ k==ℓ==1 ? 1 : 0  for k=1:N, ℓ=1:N ]), rhs = 0.5)
                    equations = [constraint1]
                    model, cplx = Density_Operators.make_sdp(My_Optimizer; n=n, C, equations)
                    optimize!(model)
                    ϱ_real = value.(model[:ϱ])
                    ϱ_complex = [cplx(ϱ_real, k, ℓ) for k in 1:N, ℓ in 1:N]
                    @test isapprox(tr(ϱ_complex), 1.0, atol=1e-5) # Check trace
                    @test isapprox(tr(constraint1.LHS ⋅ ϱ_complex), constraint1.rhs, atol=1e-5) # Check constraint
                end

                let # Test case 3: With inequality constraint
                    # Constraint:  <1|ϱ|1> ≥ 0.4
                    constraint2 = (LHS = Hermitian(ℂ[ k==ℓ==1 ? 1 : 0  for k=1:N, ℓ=1:N ]), rhs = 0.4)
                    greaterthans = [constraint2]
                    model, cplx = Density_Operators.make_sdp(My_Optimizer; n=n, C, greaterthans)
                    optimize!(model)
                    ϱ_real = value.(model[:ϱ])
                    ϱ_complex = [cplx(ϱ_real, k, ℓ) for k in 1:N, ℓ in 1:N]
                    @test isapprox(tr(ϱ_complex), 1.0, atol=1e-5) # Check trace
                    @test ℜ( tr(constraint2.LHS ⋅ ϱ_complex) ) ≥ constraint2.rhs - 1e-5 # Check constraint
                end
            end 
        end
    end
end


@testset "Sub-Module `Measure`" verbose=true begin
    @testset "measure!" verbose=true begin
        n = 2
        N = 2^n
        ϱ = Density_Operators.make_ϱ!(Hermitian(zeros(ℂ,N,N)), 0b00;n) # |00⟩
        res_list = Int8[]

        # Test measuring Z0 (should always return +1)
        Measure.measure!(res_list, ϱ; n, Z_list=UInt64[0])
        @test all(res_list .== 1)

        # Test measuring Z1 (should always return +1)
        Measure.measure!(res_list, ϱ; n, Z_list=UInt64[1])
        @test all(res_list .== 1)

        # Test measuring Z2 (should always return +1)
        Measure.measure!(res_list, ϱ; n, Z_list=UInt64[2])
        @test all(res_list .== 1)

        # Test measuring Z3 (should always return +1)
        Measure.measure!(res_list, ϱ; n, Z_list=UInt64[3])
        @test all(res_list .== 1)

        # Test with a different pure state |11⟩
        ϱ = Density_Operators.make_ϱ!(ϱ,0b11 ; n=n)
        Measure.measure!(res_list, ϱ; n, Z_list=UInt64[0b00])
        @test all(res_list .== 1)
        Measure.measure!(res_list, ϱ; n, Z_list=UInt64[0b01])
        @test all(res_list .== -1)
        Measure.measure!(res_list, ϱ; n, Z_list=UInt64[0b10])
        @test all(res_list .== -1)
        Measure.measure!(res_list, ϱ; n, Z_list=UInt64[0b11])
        @test all(res_list .== 1)

        # Test with a mixed state
        ϱ = 0.5 * Density_Operators.make_ϱ!(ϱ,                       0b00 ; n) +
            0.5 * Density_Operators.make_ϱ!(Hermitian(zeros(ℂ,N,N)), 0b11 ; n)
        Measure.reset() # Reset RNG for consistent results
        Measure.measure!(res_list, ϱ; n, Z_list=UInt64[0b00, 0b11])
        @test res_list == [1, 1]

        # Test expectation values with the mixed state
        res_list = ℝ[]
        Measure.measure!(res_list, ϱ; n, Z_list=UInt64[0b00, 0b01, 0b10, 0b11])
        @test res_list ≈ [1.0, 0.0, 0.0, 1.0] atol=1e-6
    end
end

@testset "Sub-Module Verif" verbose=true begin
    @testset "Sub-sub module DFE" verbose=true begin
        @testset "Simple runs" verbose=true begin
            for n = 1:5
                N = 2^n
                ϱ = Density_Operators.make_ϱ!(Hermitian(zeros(ℂ,N,N)), UInt(0);n) # |0...0⟩
                Measure.reset() # Reset RNG for consistent results
                Verif.DFE.reset()
                (;decision,witness) = Verif.DFE.verify( ; n,  δ=1e-8, ϱ,
                                                        num_stabs=1, num_shots_per_stab=1024)
                @test witness ≈ 1.0 atol=1e-6
                @test decision == Accept
            end
        end
        @testset "Random runs" verbose=true begin
            for e = 1:3
                δ = 10.0^-e
                for n = 1:5
                    Instances.reset()
                    Measure.reset() # Reset RNG for consistent results
                    Verif.DFE.reset()
                    ψ₀ = zeros(ℂ, 2^n)
                    ψ₀[1] = 1  # Set the first element to 1 for |0...0⟩
                    ϱ = Instances.make_rnd_dens_op(;n, ψ₀, ε=3δ)
                    @assert δ > 1e-7 "Are you nuts?!!"
                    m = Int(ceil( 100/δ ))
                    (;decision,witness) = Verif.DFE.verify( ; n, δ, ϱ,
                                                            num_stabs=m, num_shots_per_stab=1)
                    @test 1-5δ ≤ witness ≤ 1-δ
                    @test !( witness < 1-3δ )  || decision==Reject
                    @test !( witness ≥ 1-3δ )  || decision==Accept
                end
            end
        end
        @testset "Random accepts" verbose=true begin
            for e = 1:3
                δ = 10.0^-e
                for n = 1:5
                    N = 2^n
                    ψ₀ = zeros(ℂ, 2^n)
                    ψ₀[1] = 1  # Set the first element to 1 for |0...0⟩
                    Instances.reset()
                    Measure.reset() # Reset RNG for consistent results
                    Verif.DFE.reset()
                    ϱ = Instances.make_rnd_dens_op(;n, ψ₀, ε=δ)
                    @assert δ > 1e-7 "Are you nuts?!!"
                    m = Int(ceil( 100/δ ))
                    (;decision,witness) = Verif.DFE.verify( ; n, δ, ϱ,
                                                            num_stabs=m, num_shots_per_stab=1)
                    @test decision==Accept
                end
            end
        end
        @testset "Random rejects" verbose=true begin
            for e = 1:3
                δ = 10.0^-e
                for n = 1:5
                    N = 2^n
                    ψ₀ = zeros(ℂ, 2^n)
                    ψ₀[1] = 1  # Set the first element to 1 for |0...0⟩
                    Instances.reset()
                    Measure.reset() # Reset RNG for consistent results
                    Verif.DFE.reset()
                    ϱ = Instances.make_rnd_dens_op(;n, ψ₀, ε=5δ)
                    @assert δ > 1e-7 "Are you nuts?!!"
                    m = Int(ceil( 100/δ ))
                    (;decision,witness) = Verif.DFE.verify( ; n, δ, ϱ,
                                                            num_stabs=m, num_shots_per_stab=1)
                    @test decision==Reject
                end
            end
        end

        @testset "Sub-sub module DFE with Diagonal Instances" verbose=true begin  # New testset
            for e = 1:3
                δ = 10.0^-e
                for n = 2:5 # n ≥ 2
                    Instances.reset()
                    Measure.reset()
                    Verif.DFE.reset()


                    ϱ_diag ::Vector{ℝ} = Instances.make_rnd_dens_diag(Instances.Dumb();n=n, ε=3δ)
                    @test all( ϱᵢ ≥ 0  for ϱᵢ in ϱ_diag  )
                    @test ∑( ϱ_diag ) ≈ 1
                    @assert δ > 1e-7 "Are you nuts?!!"
                    m = Int(ceil(100 / δ))
                    (; decision, witness) = Verif.DFE.verify(; n, δ, ϱ=ϱ_diag,
                                                             num_stabs=m, num_shots_per_stab=1)
                    @test 1 - 5δ ≤ witness ≤ 1 - δ  # The witness should be close to 1-ε
                    @test !(witness < 1 - 3δ) || decision == Reject # Below threshold should be rejected
                    @test !(witness ≥ 1 - 3δ) || decision == Accept  # Above or equal to threshold should be accepted
                end
            end
        end

        @testset "Random accepts with diagonal instances" verbose=true begin
             for e = 1:3
                δ = 10.0^-e
                for n = 2:5 # n ≥ 2
                    Instances.reset()
                    Measure.reset() # Reset RNG for consistent results
                    Verif.DFE.reset()
                    ϱ_diag = Instances.make_rnd_dens_diag(Instances.Dumb();n, ε=δ)
                    @test all( ϱᵢ ≥ 0  for ϱᵢ in ϱ_diag  )
                    @test ∑( ϱ_diag ) ≈ 1

                    @assert δ > 1e-7 "Are you nuts?!!"
                    m = Int(ceil( 100/δ ))
                    (;decision,witness) = Verif.DFE.verify( ; n, δ, ϱ=ϱ_diag, num_stabs=m, num_shots_per_stab=1)
                    @test decision==Accept
                end
            end
        end

        @testset "Random rejects with diagonal instances" verbose=true begin
             for e = 1:3
                δ = 10.0^-e
                for n = 2:5 # n ≥ 2
                    Instances.reset()
                    Measure.reset() # Reset RNG for consistent results
                    Verif.DFE.reset()
                    ϱ_diag = Instances.make_rnd_dens_diag(Instances.Dumb();n, ε=5δ)
                    @test all( ϱᵢ ≥ 0  for ϱᵢ in ϱ_diag  )
                    @test ∑( ϱ_diag ) ≈ 1

                    @assert δ > 1e-7 "Are you nuts?!!"
                    m = Int(ceil( 100/δ ))
                    (;decision,witness) = Verif.DFE.verify( ; n, δ, ϱ=ϱ_diag, num_stabs=m, num_shots_per_stab=1)
                    @test decision==Reject
                end
            end
        end
    end #^ testset Sub-sub module DFE
end

@testset "Sub-module `Stabilizer_Foolers`" verbose=true begin
    n = 2
    N = 2^n
    δ = 0.01

    min = Stabilizer_Foolers.min

    @testset "FALSE_POSITIVE with fid" begin
        ε = 3δ
        fid = 1 - δ
        p = nothing
        scenario = FALSE_POSITIVE
        opt =
            redirect_stdout(devnull) do
                Stabilizer_Foolers.mip_stabilizer_fooler(;n, ε, fid, p, scenario, time_limit=1min)
            end
        @test opt.f[0] ≤ fid
        @test all(sum(fourier_coeff(x, y) * opt.f[x] for x in 0:N-1) ≥ (1 - ε) * opt.big[y] -1e-9   for y in 1:N-1)
        # Add a test that checks the objective is maximized, perhaps comparing to a naive solution.  This requires more problem-specific knowledge.
    end

    @testset "FALSE_NEGATIVE with fid" begin
        ε = 5δ
        fid = 1 - δ
        p = nothing
        scenario = FALSE_NEGATIVE
        opt =
            redirect_stdout(devnull) do
                Stabilizer_Foolers.mip_stabilizer_fooler(;n, ε, fid, p, scenario, time_limit=1min)
            end
        @test opt.f[0] ≥ fid
        @test all(sum(fourier_coeff(x, y) * opt.f[x] for x in 0:N-1) ≤ 1 − ε*(1-opt.big[y]) for y in 1:N-1)
        # Add a test that checks the objective is minimized.
    end

    @testset "FALSE_POSITIVE with p" begin
        ε = 5δ
        fid = nothing
        p = 0.9
        scenario = FALSE_POSITIVE
        opt =
            redirect_stdout(devnull) do
                Stabilizer_Foolers.mip_stabilizer_fooler(;n, ε, fid, p, scenario, time_limit=1min)
            end
        @test sum(opt.big) ≥ p * (N-1)
        @test all(sum(fourier_coeff(x, y) * opt.f[x] for x in 0:N-1) ≥ (1 - ε) * opt.big[y] for y in 1:N-1)
        # Add a test that checks the objective (f[0]) is minimized.
    end


    @testset "FALSE_NEGATIVE with p" begin
        ε = 3δ
        fid = nothing
        p = 0.9
        scenario = FALSE_NEGATIVE
        opt =
            redirect_stdout(devnull) do
                Stabilizer_Foolers.mip_stabilizer_fooler(;n, ε, fid, p, scenario, time_limit=1min)
            end
        @test sum(opt.big) ≤ p * (N-1)
        @test all(sum(fourier_coeff(x, y) * opt.f[x] for x in 0:N-1) ≤ 1 − ε*(1-opt.big[y]) for y in 1:N-1)
        # Add a test that checks the objective (f[0]) is maximized.
    end


    # Edge case tests
    redirect_stdout(devnull) do
        @test_throws AssertionError Stabilizer_Foolers.mip_stabilizer_fooler(;n, ε=1.0)  # ε out of range
        @test_throws AssertionError Stabilizer_Foolers.mip_stabilizer_fooler(;n, ε=0.0, fid=1.1) # fid out of range
        @test_throws AssertionError Stabilizer_Foolers.mip_stabilizer_fooler(;n, ε=0.0, p=1.1)   # p out of range
        @test_throws AssertionError Stabilizer_Foolers.mip_stabilizer_fooler(;n, ε=0.0, fid=0.5, p=0.5) # both fid and p provided
        @test_throws AssertionError Stabilizer_Foolers.mip_stabilizer_fooler(;n, ε=0.0) # neither fid nor p provided.
    end
    redirect_stdout(devnull) do
        Stabilizer_Foolers.false_positive_analysis();
        Stabilizer_Foolers.false_negative_analysis();
    end

    # reduced stuff

    @testset "reduced_false_positive_analysis" begin
        n = 2
        α = 0.5  # Example α value
        (;opt,slacks) = redirect_stdout(devnull) do
            Stabilizer_Foolers.reduced_false_positive_analysis(; n, α, time_limit=61.0) # Short time limit for testing
        end

        # Check if the optimization ran successfully (adjust tolerance as needed)
        @test opt.status == MOI.OPTIMAL || opt.status == MOI.TIME_LIMIT

        # Check that the sum of probabilities is 1 (within a tolerance)
        @test isapprox(sum(opt.p), 1.0, atol=1e-6)

        # Basic check on the constraint (choose a specific y for this example)
        y = 1
        left_side = sum((count_ones(y & x) % 2) * opt.p[x] for x in 1:2^n - 1)
        right_side = 1 - opt.acceptrow[y] * (1 - α / 2) # As implemented in reduced_fooler.jl
        @test left_side ≤ right_side + 1e-6 # Added tolerance
    end

end

@testset "Function handle_broken_state" verbose=true begin

    @testset "With density operator" verbose=true begin
    # Test setup
    n = 1
    δ = 0.1
    ϱ_dummy    = Hermitian( [1.0+2.0im 3.0-4.0im; 5.0+6.0im 7.0-8.0im] ) # Complex test matrix
    ϱ_with_nan = Hermitian( [NaN 1.0+2.0im; 7.0-8.0im Inf-Inf*im]      )


    mktempdir() do tmpdir
        cd(tmpdir) do

            @testset "FALSE_NEGATIVE scenario" begin
                CompVeris.handle_broken_state(FALSE_NEGATIVE, ϱ_dummy; n=n, δ=δ)
                expected_filename = Printf.@sprintf("broken-goodstate-%02u-%#a.bin", n, δ)
                @test isfile(expected_filename)

                ϱ_read =
                    open(expected_filename) do io
                        vec = Vector{ℂ}(undef,2^n⋅2^n)
                        read!(io, vec)
                        Hermitian( reshape(vec, (2^n,2^n)) )
                    end

                @test isequal(ϱ_read,ϱ_dummy) # takes care of NaNs

                CompVeris.handle_broken_state(FALSE_NEGATIVE, ϱ_with_nan; n=n, δ=δ)

                ϱ_read =
                    open(expected_filename) do io
                        vec = Vector{ℂ}(undef,2^n⋅2^n)
                        read!(io, vec)
                        Hermitian( reshape(vec, (2^n,2^n)) )
                    end

                @test isequal(ϱ_read, ϱ_with_nan)  # Test with NaN, Inf
            end

            @testset "FALSE_POSITIVE scenario" begin
                CompVeris.handle_broken_state(FALSE_POSITIVE, ϱ_dummy; n=n, δ=δ)
                expected_filename = Printf.@sprintf("broken-badstate-%02u-%#a.bin", n, δ)
                @test isfile(expected_filename)

                ϱ_read =
                    open(expected_filename) do io
                        vec = Vector{ℂ}(undef,2^n⋅2^n)
                        read!(io, vec)
                        Hermitian( reshape(vec, (2^n,2^n)) )
                    end

                @test isequal(ϱ_read,ϱ_dummy)

                CompVeris.handle_broken_state(FALSE_POSITIVE, ϱ_with_nan; n=n, δ=δ)
                ϱ_read =
                    open(expected_filename) do io
                        vec = Vector{ℂ}(undef,2^n⋅2^n)
                        read!(io, vec)
                        Hermitian( reshape(vec, (2^n,2^n)) )
                    end
                @test isequal(ϱ_read, ϱ_with_nan) # Test with NaN, Inf
             end

        end
    end
    end #^ testset for density operators

    @testset "With Diagonal Vector" verbose=true begin
        n = 1
        δ = 0.1
        diag_dummy    = [1.0, 2.0]       # Test diagonal
        diag_with_nan = [NaN, Inf] # Diagonal with NaN and Inf

        mktempdir() do tmpdir
            cd(tmpdir) do
                @testset "FALSE_NEGATIVE scenario" begin
                    CompVeris.handle_broken_state(FALSE_NEGATIVE, diag_dummy; n=n, δ=δ)
                    expected_filename = Printf.@sprintf("broken-goodstate-%02u-%#a.bin", n, δ)
                    @test isfile(expected_filename)

                    diag_read = open(expected_filename) do io
                        vec = Vector{ℝ}(undef,2^n)
                        read!(io, vec)
                        vec
                    end

                    @test isequal(diag_read, diag_dummy)


                    CompVeris.handle_broken_state(FALSE_NEGATIVE, diag_with_nan; n=n, δ=δ)
                    diag_read = open(expected_filename) do io
                        vec = Vector{ℝ}(undef,2^n)
                        read!(io, vec)
                        vec
                    end

                    @test isequal(diag_read, diag_with_nan) # Test with NaN, Inf
                end

                @testset "FALSE_POSITIVE scenario" begin
                    CompVeris.handle_broken_state(FALSE_POSITIVE, diag_dummy; n=n, δ=δ)
                    expected_filename = Printf.@sprintf("broken-badstate-%02u-%#a.bin", n, δ)
                    @test isfile(expected_filename)

                    diag_read = open(expected_filename) do io
                        vec = Vector{ℝ}(undef,2^n)
                        read!(io, vec)
                        vec
                    end

                    @test isequal(diag_read, diag_dummy)

                    CompVeris.handle_broken_state(FALSE_POSITIVE, diag_with_nan; n=n, δ=δ)
                    diag_read = open(expected_filename) do io
                        vec = Vector{ℝ}(undef,2^n)
                        read!(io, vec)
                        vec
                    end
                    @test isequal(diag_read, diag_with_nan)  # Test with NaN, Inf
                end
            end
        end
    end #^ testset With Diagonal Vector

end #^ testset for `handle_broken_state`

@testset "Functions full_set__compare*" verbose=true begin
    opt =
        redirect_stdout(devnull) do
            full_set__compare(DENSITY_OPERATOR_INSTANCES;
                              n_range=2:3, d_range=2:3, num_instances=2, reps=2)
        end
    @test true
    opt =
        redirect_stdout(devnull) do
            full_set__compare(DENSITY_DIAGONAL_INSTANCES;
                              n_range=2:3, d_range=2:3, num_instances=2, reps=2)
        end
    @test true
    opt =
        redirect_stdout(devnull) do
            full_set__compare_bad_states(DENSITY_OPERATOR_INSTANCES;
                                         n_range=2:3, d_range=2:3, num_instances=2, reps=2)
        end
    @test true
    opt =
        redirect_stdout(devnull) do
            full_set__compare_bad_states(DENSITY_DIAGONAL_INSTANCES ;
                                         n_range=2:3, d_range=2:3, num_instances=2, reps=2)
        end
    @test true
    opt =
        redirect_stdout(devnull) do
            full_set__compare_bad_states(DENSITY_DIAGONAL_INSTANCES_SMARTBAD ;
                                         n_range=2:3, d_range=2:3, num_instances=2, reps=2)
        end
    @test true
end

end #^ testset CompVeris
#EOF
