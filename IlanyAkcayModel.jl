population = []

popSize = 100
let population = population
    for(i) in 1:popSize
        push!(population, [])
    end
end

function indexOf(pop::Array{Any, 1}, parent::Int64, child::Int64)
    for(n) in 1:length(pop)
        if(pop[n] === pop[parent][child])
            return n
        end
    end
end

numInitEdges = round(1.5*length(population))
let population = population
    repairedLinks = 0
    for(i) in 1:numInitEdges
            firstID = Int(round(rand()*length(population)+.5))
            secondID = firstID
            while(secondID == firstID)
                secondID = Int(round(rand()*length(population)+.5))
                for(j) in 0:length(population[firstID])
                    if(j>0)
                        if(population[firstID][j] === population[secondID])
                            secondID = firstID
                        end
                    end
                end
            end
            if(firstID != secondID)
                push!(population[firstID], population[secondID])
                push!(population[secondID], population[firstID])
            else
                println("Link rebuilt")
            end
    end
end

using Plots
function linkBarChart(pop::Array{Any, 1})
    freqData = []
    for(i) in 1:(length(pop)-1)
        counter = 0
        for(ii) in 1:length(pop)
            if(length(pop[ii]) >= i)
                counter = counter + 1
            end
        end
        push!(freqData, counter)
    end
    zeroCounter = 0
    for(i) in 1:length(pop)
        if(length(pop[i])==0)
            zeroCounter = zeroCounter + 1
        end
    end
    pushfirst!(freqData, (zeroCounter+freqData[1]))
    #println(freqData)
    plotAdjust = length(pop)/15
    linkBars = bar(0:1:100, freqData)
    linkBars
end


#println("New")

#running generations to stabilize (goal is 20n)
numGens = 200

for(g) in 1:numGens
    spliceID = Int(round(rand()*length(population)+.5))
    spliceContents = []
    #print("Cell to be removed has length of $(length(population[spliceID])): ")
    for(c) in 1:length(population[spliceID])
        push!(spliceContents, indexOf(population, spliceID, c))
    end
    #println(spliceContents)
    for(s) in 1:length(population)
        for(ss) in 1:length(population[s])
            if(population[s][ss] === population[spliceID])
                #println("$(s) loses its $(ss) connection")
                splice!(population[s], ss)
                break
            end
        end
    end
    #individuals who know the individual that will soon die lose their connections to the latter
    splice!(population, spliceID)
    #one individual dies at random


    momIndex = Int(round(rand()*length(population)+.5))
    pN = .7
    pR = .3
    #variables for inheritance

    push!(population, [population[momIndex]])
    push!(population[momIndex], population[length(population)])
    #forms baby with connection to mother

    for(n) in 1:length(population[momIndex])
        if(rand() < pN)
            push!(population[length(population)], population[momIndex][n])
            push!(population[n], population[length(population)])
        end
    end
    #adds links with inheritance from mom's links

    for(r) in 1:length(population)
        momLink = false
        for(n) in 1:length(population[momIndex])
            if(population[r] === population[momIndex][n])
                momLink = true
            end
        end
        if(!momLink)
            if(rand() < pR)
                push!(population[length(population)], population[r])
                push!(population[r], population[length(population)])
            end
        end
    end
    #forms random connections with nodes not connected to mother

end

#println("poplength: $(length(population))")
linkBarChart(population)
