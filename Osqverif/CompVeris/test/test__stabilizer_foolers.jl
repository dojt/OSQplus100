using DOT_NiceMath
using DOT_NiceMath.Numbers64

using CompVeris
using CompVeris.Stabilizer_Foolers

const minute  = Stabilizer_Foolers.min
const minutes = Stabilizer_Foolers.min

Supp(f) =  begin
               L=length(f)
               [ ξ for ξ=1:L if f[ξ] != 0 ]
           end

for n = 4:10
    #for α = 1/4: 1/8 : 1/2
    let α = 1/3
        for _iter = 1:10

            println("#######################################################################################")
            println("#######################################################################################")
            println("vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv")
            @show n,α,_iter

            println("""_____________________________________________________________________________________
                       Start: n=$(n), 2/n=$(2/n), α=$(α)""")

            (;opt,slacks) =
                reduced_false_positive_analysis(
                    ; n, α,
                    time_limit       = max(1,n-6) * 30minutes,
                    perturb          = 1,
                    # fix_Id           = false, # default: yes, iff α ≥ 2/n
                    basis_repeats = 1_000_000 )

            println("""_____________________________________________________________________________________
                       Finished n=$(n), 2/n=$(2/n), α=$(α)""")

            p = abs.(round.(opt.p .* (1<<16))) ./ (1<<16)
            @show Supp(p)
            @show sort(unique(p))
            @show p

            println("^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")
            flush(stdout)
            println()
        end
    end
end
