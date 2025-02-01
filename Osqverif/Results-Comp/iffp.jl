# Results-Comp/iffp.jl
module Intrinsic_False_Positive_Probability_Bounds

    using Test

    using DOT_NiceMath             # `⋅`, `∑`=sum, `∏`=prod, etc
    using DOT_NiceMath.NumbersBig  # `ℝ`=BigFloat, `ℤ`=BigInt, etc


    function _matrix_counting(; n::ℤ, m::ℤ) ::ℤ
        @assert 1 ≤ m < n

        binomial(n::ℤ, k::ℤ) ::ℤ = factorial(n) ÷ factorial(k) ÷ factorial(n-k)

        return (-1)^m ⋅ ∑(   binomial(m,j) ⋅ (-1)^j ⋅ (1+j)^n   for j = 0 : m   )
    end

    function lowerbound(; n::ℤ, α::ℝ) ::ℝ
        ;                                   @assert 0 < α ≤ 2
        m =  round(ℤ, 2/α )            ;    @assert m ≈ 2/α     # 2/α must be integer

        if m > n
            return 0.0 |> ℝ
        end

        denominator ::ℤ = ∏(  2^n - 2^j    for j = 0 : m-1  )
        numerator   ::ℤ =
            if m < n
                _matrix_counting(;n,m)
            else # m == n
                factorial(n)
            end

        return numerator / denominator

    end #^ lowerbound()

    function upperbound(;n::ℤ, α::ℝ) ::ℝ
        @assert 0 < α ≤ 2
        @assert n ≥ 3

        α < 2 || return 1.0 |> ℝ

        r(m::ℤ) ::ℝ = 1/(  2^(m+1) - 2  )

        η(k::ℤ ; α::ℝ,n::ℤ) =
            let ret = ( 1 − 2r(max(k,n-k)) )/( 2 − α )
                ret
            end

        return ∏(  η(k ; α,n )   for k = 0 : n-1  )
    end #^ upperbound()


    @testset "Testing ifpp_upper/lower-bound" verbose=true begin

        𝟏 ::ℤ = 1
        𝟑 ::ℤ = 3

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

#EOF
