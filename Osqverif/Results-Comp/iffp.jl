module Intrinsic_False_Positive_Probability_Bounds

    using Test

    using DOT_NiceMath             # ⋅=* etc
    using DOT_NiceMath.NumbersBig  # ℝ=BigFloat, ℤ=BigInt, etc

    𝟏 ::ℤ = 1
    𝟐 ::ℤ = 2
    𝟑 ::ℤ = 3

    function _matrix_counting(; n::ℤ, m::ℤ) ::ℤ
        @assert 1 ≤ m < n

        bino(n::ℤ, k::ℤ) ::BigInt =
            let 𝐧 = big(n), 𝐤 = big(k)
                factorial(𝐧) ÷ factorial(𝐤) ÷ factorial(𝐧-𝐤)
            end

        return (-1)^m ⋅ ∑(   bino(m,j) ⋅ (-1)^j ⋅ (𝟏+j)^n   for j = 0 : m   )
    end

    function lowerbound(; n::ℤ, α::ℝ) ::ℝ

        ;                                   @assert 1 ≤ n < 256
        ;                                   @assert 0 < α ≤ 2
        m =  round(ℤ, 𝟐/α )            ;    @assert m ≈ 2/α     # 2/α must be integer

        if m > n
            return 0.0 |> ℝ
        end

        denominator ::ℤ = prod(  𝟐^n - 𝟐^j    for j = 0 : m-1  )
        numerator   ::ℤ =
            if m < n
                _matrix_counting(;n,m)
            elseif m == n
                factorial(big(n))
            end

        return numerator / denominator

    end #^ lowerbound()

    function upperbound(;n::ℤ, α::ℝ) ::BigFloat
        ;                                   @assert 1 ≤ n < 256
        ;                                   @assert 0 < α ≤ 2

        @assert n ≥ 3

        α < 2 || return big(1.0)

        r(m::ℤ) ::BigFloat = begin @assert m>0 ; 𝟏/(  𝟐^(m+1) - 2  ) end

        η(k::ℤ ; α::ℝ,n::ℤ) =
            let ret = ( 1 − 2r(max(k,n-k)) )/( 𝟐 − α )   ;            @assert ret!=0 ; @assert ret > 0
                ret
            end

        return prod(  η(k ; α,n )   for k = 0 : n-1  )
    end #^ upperbound()


    @testset "Testing ifpp_upper/lower-bound" verbose=true begin

        for m = 𝟏:32
            α = 2/m

            for n = 𝟑 : 32
                lb = lowerbound(;n,α)
                ub = upperbound(;n,α)

                @test isfinite(lb)
                @test isfinite(ub)

                @test lb == 1 == ub   ||   lb < ub
            end
        end

    end #^ testset

end #^ module Intrinsic_False_Positive_Probability_Bounds
