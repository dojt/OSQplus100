# (c) Dirk Oliver Theis

using DOT_NiceMath
using DOT_NiceMath.Numbers64

using CompVeris
using CompVeris.Verif: Decision_t, Accept, Reject
using CompVeris.Verif
using CompVeris.Measure


# Thm.  If everybody in a fixed basis is at least 1-δ, then the fidelity is at least 1-n⋅δ/2.
# Q.    If everybody in a random basis is at least 1-δ, then the fidelity is at least >> 1-n⋅δ/2??


function tryit_good( ;
                     δ                = 1e-2,
                     n                = 4,
                     num_instances    = 1,
                     num_bases        = 1)

    𝑚𝑎𝑥 = -Inf

    ψ₀ = zeros(ℂ, 2^n)
    ψ₀[1] = 1  # Set the first element to 1 for |0...0⟩

    expvals = ℝ[]
    Z_list  = zeros(UInt64, n)

    for _ii = 1:num_instances

        ϱ = Instances.make_rnd_dens_op(;n, ψ₀, ε=δ)  # good state


        for _bb = 1:num_bases
            MYNS._sample_basis!(Z_list;n)
            measure!(expvals ; n,ϱ, Z_list)

            𝑚𝑖𝑛 = minimum( expvals )
            𝑡ℎ𝑚 = 1-n⋅(1-𝑚𝑖𝑛)/2

            𝑚𝑎𝑥 = max(𝑡ℎ𝑚,𝑚𝑎𝑥)

        end

    end

@show 1-δ, 𝑚𝑎𝑥

    return  𝑚𝑎𝑥
end


# Test
tryit_good(;δ=1e-2, n=2, num_instances=2, num_bases=2)




# Runs

@time tryit_good(;δ=1e-2, n=6, num_instances=50, num_bases=200)
