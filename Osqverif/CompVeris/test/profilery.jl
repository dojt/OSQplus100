# CompVeris/test/runtests.jl
# (c) Dirk Oliver Theis
# Open source license: CC0

using Profile
using ProfileView  # Optional, for graphical representation

using DOT_NiceMath
using DOT_NiceMath.Numbers64

using CompVeris


function run_some_DFEs()

    for e = 1:6
        ε = 10.0^-e
        for n = 1:5
            N = 2^n
            ψ₀ = zeros(ℂ, 2^n)
            ψ₀[1] = 1  # Set the first element to 1 for |0...0⟩
            ϱ = Instances.make_rnd_dens_op(; n, ψ₀, ε=3ε)
            Measure.measure_reset() # Reset RNG for consistent results
            m = Int(ceil(100 / ε))
            E = Verif.Verif_DFE.verify(; n, ϱ, num_stabs=m, num_shots_per_stab=1)
            if 1-5ε ≤ E ≤ 1-ε
                println("n=$(n), ε=$(ε): Test passed with $(1-5ε) ≤ E=$(E) ≤ $(1-ε)")
            else
                println("Test failed:  Oops $(1-5ε) ≤ E=$(E) ≤ $(1-ε)")
            end
        end #^ for n
    end #^ for e(i.e., ε)

end #^ fn

@profile run_some_DFEs()

@show "Are we done?"

# Print the profiling results
Profile.print()

# Optional: Use ProfileView for a graphical representation
ProfileView.view()

#EOF












