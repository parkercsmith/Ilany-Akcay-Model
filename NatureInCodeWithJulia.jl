#!/usr/bin/env julia

using Plots
using StatsBase

mutable struct frequencies
    p::Float64
    q::Float64
    pp::Float64
    pq::Float64
    qq::Float64

    function frequencies(P::Float64)
        p = P
        q = 1 - p
        pp = p * p
        qq = q * q
        pq = p * q * 2
        new(p,q,pp,pq,qq)
    end
end

mutable struct plotter
    pArray::Array{Float64, 1}
    qArray::Array{Float64, 1}
    ppArray::Array{Float64, 1}
    pqArray::Array{Float64, 1}
    qqArray::Array{Float64, 1}

    function plotter()
        pArray = []
        qArray = []
        ppArray = []
        pqArray = []
        qqArray = []
        new(pArray, qArray, ppArray, pqArray, qqArray)
    end
end

frequencyStuff = frequencies(.6)
plottingStuff = plotter()

function runGens(frequs::frequencies, plots::plotter)
    freqs = frequs
    firstAlleles = []
    for(i) in 1:1000
        if(rand()<freqs.p)
            push!(firstAlleles, 0)
        else
            push!(firstAlleles, 1)
        end
    end
    secondAlleles = []
    for(i) in 1:1000
        if(rand()<freqs.p)
            push!(secondAlleles, 0)
        else
            push!(secondAlleles, 1)
        end
    end
    for(i) in 1:10000
        newFirst = []
        newSecond = []
        for(ii) in 1:1000

            mateID = Int64(round(rand()*1000+0.5))

            mateID = 0
            mateWeights = weights([1,2,3,4,4,3,2,1])
            if(ii==1)
                mateID = sample([997,998,999,1000,2,3,4,5], mateWeights)
            elseif(ii==1000)
                mateID = sample([996,997,998,999,1,2,3,4], mateWeights)
            elseif(ii==999)
                mateID = sample([995,996,997,998,999,1,2,3], mateWeights)
            elseif(ii==998)
                mateID = sample([994,995,996,997,998,999,1,2], mateWeights)
            elseif(ii==997)
                mateID = sample([993,994,995,996,997,998,999,1], mateWeights)
            elseif(ii==2)
                mateID = sample([998,999,1000,1,3,4,5,6], mateWeights)
            elseif(ii==3)
                mateID = sample([999,1000,1,2,4,5,6,7], mateWeights)
            elseif(ii==4)
                mateID = sample([1000,1,2,3,5,6,7,8], mateWeights)
            else
                mateID = sample([ii-4,ii-3,ii-2,ii-1,ii+1,ii+2,ii+3,ii+4], mateWeights)
            end


            if(rand()<.5)
                push!(newFirst, firstAlleles[ii])
            else
                push!(newFirst, firstAlleles[mateID])
            end
            if(rand()<.5)
                push!(newSecond, secondAlleles[ii])
            else
                push!(newSecond, secondAlleles[mateID])
            end
        end

        firstAlleles = copy(newFirst)
        secondAlleles = copy(newSecond)

        pSum = 0.0
        for(ii) in 1:1000
            if(firstAlleles[ii]==0)
                pSum+=1.0
            end
            if(secondAlleles[ii]==0)
                pSum+=1.0
            end
        end
        pSum/=2000
        freqs = frequencies(pSum)

        push!(plots.pArray, freqs.p)
        push!(plots.qArray, freqs.q)
        push!(plots.ppArray, freqs.pp)
        push!(plots.pqArray, freqs.pq)
        push!(plots.qqArray, freqs.qq)
    end
    #=plot(plots.pArray, label = "p")
    plot!(plots.qArray, label = "q")
    plot!(plots.ppArray, label = "pp")
    plot!(plots.pqArray, label = "pq")
    plot!(plots.qqArray, label = "qq")=#
end

runGens(frequencyStuff, plottingStuff)
