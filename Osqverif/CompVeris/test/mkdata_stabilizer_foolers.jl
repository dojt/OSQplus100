using DOT_NiceMath
using DOT_NiceMath.Numbers64

using CompVeris
using CompVeris.Stabilizer_Foolers

const minute  = Stabilizer_Foolers.min
const minutes = Stabilizer_Foolers.min

# MIP

for n = 2:10
    for d = 3:19
        δ = min( (1/2)^d  , 2/n)

        false_negative_analysis( ; n, δ,
                                 time_limit = 180minutes,
                                 skip_basis_stuff = false)
        false_positive_analysis( ; n, δ,
                                 time_limit = 180minutes,
                                 skip_basis_stuff = false)
    end
end


# LP-end


let totalresultlist = []
    println()
    for n = 9:9
        for d = 3:19
            δ = min( (1/2)^d  , 2/n)

            let ε          = 3δ,
                fid        = 1-δ,
                time_limit = 10minutes

                println("1-ε = $(1-ε)")
                resultlist = lpfeas_stabilizer_fooler(;n,ε,fid,time_limit)
                println()
                append!(totalresultlist,resultlist)
            end
        end
    end
    @show totalresultlist
    totalresultlist
end


# LP-after shift

let totalresultlist = []
    println()
    for n = 16:16
        for d = 20:20
            δ = min( (1/2)^d  , 2/n)

            let time_limit = 10minutes

                println("\n-3δ = $(-3δ)")
                resultlist = lpfeas_stabilizer_fooler_g(;n,δ,time_limit)
                append!(totalresultlist,resultlist)
            end
        end
    end
    @show totalresultlist
    totalresultlist
end



let totalresultlist = []
    println()
    for n = 2:2
        for d = 4:4
            δ = min( (1/2)^d  , 2/n)

            let time_limit = 10minutes

                println("\n-3δ = $(-3δ)")
                resultlist = lpfeas_stabilizer_fooler_g(;n,δ,time_limit, dumpsols=true)
                append!(totalresultlist,resultlist)
            end
        end
    end
    @show totalresultlist
    totalresultlist
end
