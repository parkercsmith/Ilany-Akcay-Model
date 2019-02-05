using StatsBase

population = []
#represents all nodes in population
edgeMatrix = []
#matrix; holds all edges in population where row and
#col are indices of nodes

global meanCoopRatio = 0.0

mutable struct Node
    strategy::Int8 #1 if cooperator, 0 if defector
    payoff::Float64
    fitness::Float64
    ID::Int64
end

popSize = 100
numInitialCoops = popSize/2
#population = collect(1:popSize)
let numInitialCoops = numInitialCoops
    for(i) in 1:popSize
        if(numInitialCoops > 0)
            push!(population, Node(1, 0, 1, i))
            numInitialCoops -= 1
        else
            push!(population, Node(0, 0, 1, i))
        end
    end
end
#initialize a population of popSize nodes; each node
#contains a strategy 1/0 representing coop/defect,
#as well as a default payoff of 0, fitness of 1, and ID


edgeMatrix = zeros(Int8, popSize, popSize)
#initialize an empty edge matrix

#=function indexOf(pop::Array{Any, 1}, parent::Int64, child::Int64)
    for(n) in 1:length(pop)
        if(pop[n] === pop[parent][child])
            return n
        end
    end
end=#
#OBSOLETE; allowed recursive arrays to identify
#their own index

numInitEdges = Int(round(1.5*popSize))
let population = population
    for(i) in 1:numInitEdges
            firstID = Int(round(rand()*length(population)+.5))
            secondID = firstID
            while(secondID == firstID)
                secondID = Int(round(rand()*length(population)+.5))
                if(edgeMatrix[firstID, secondID] == 1 || edgeMatrix[secondID, firstID] ==1)
                    secondID = firstID
                    continue
                end
            end
            edgeMatrix[firstID, secondID] = 1
            edgeMatrix[secondID, firstID] = 1
    end
end
#identifies a number of initial edges based on
#popSize, creates numInitEdges distinct pairs of
#distinct indices, and changes value of corresponding
#edgeMatrix placeholders to 1

using Plots
function linkBarChart()
    #=
    freqData = []
    for(i) in 1:popSize
        counter = 0
        for(ii) in 1:popSize
            degreeMarker = 0
            for(iii) in 1:popSize
                if(edgeMatrix[ii, iii] == 1)
                    degreeMarker += 1
                end
            end
            if(degreeMarker >= i)
                counter += 1
            end
        end
        push!(freqData, counter)
    end
    #Examines each row of edgeMatrix, counts present
    #edges, and forms cumulative frequency data
    #based on the count for every possible degree.

    #=zeroCounter = 0
    for(i) in 1:popSize
        noEdges = true
        for(ii) in 1:popSize
            if(edgeMatrix[i, ii] == 1)
                noEdges = false
            end
        end
        if(noEdges)
            zeroCounter += 1
        end
    end=#
    #Above segment counts number of nodes with no
    #edges. This is unnecessary for a cumulative
    #frequency graph like the current one, but
    #is useful to have if a normal frequency graph
    #is ever used. Below push statement should also
    #be altered in this case.

    pushfirst!(freqData, popSize#=(zeroCounter+freqData[1])=#)
    #println(freqData)

    degData = []
    for(d) in 1:length(freqData)-1
        push!(degData, (freqData[d]-freqData[d+1]))
    end
    degSum = 0
    for(d) in 1:length(degData)
        degSum += (d * degData[d])
    end
    meanDegree = Int64(round(degSum/popSize))
    println("Mean Degree: $(meanDegree)")
    =#

    #=
    adjacencySum = 0
    #SCI FAIR Counts individuals two away from root node
    #that also share an edge to the root node; looks for
    #"triangles" of 1-edge length
    for(i) in 1:popSize
        firstNeighbors = []
        secondNeighbors = 0
        secondAdjacencyCounter = 0
        for(ii) in 1:popSize
            if(edgeMatrix[i, ii] == 1)
                push!(firstNeighbors, ii)
            end
        end
        for(j) in 1:length(firstNeighbors)
            for(jj) in 1:popSize
                if(edgeMatrix[firstNeighbors[j], jj] == 1)
                    secondNeighbors += 1
                    if(edgeMatrix[i, jj] == 1)
                        secondAdjacencyCounter += 1
                    end
                end
            end
        end
        adRat = 0
        if(secondNeighbors>0)
            adRat = Float64(secondAdjacencyCounter/secondNeighbors)
        end
        adjacencySum += adRat
        GC.gc()
    end

    meanAdjacencyRatio = Float64(adjacencySum/popSize)
    println("Mean Adjacency Ratio at a = 2: $(meanAdjacencyRatio)")
    =#

    #calculating the frequency of cooperation
    coopCount = 0.0
    for(i) in 1:popSize
        if(population[i].strategy == 1)
            coopCount += 1.0
        end
    end
    coopRatio = coopCount/popSize
    global meanCoopRatio
    meanCoopRatio += coopRatio
    #println("Frequency of Cooperation: $(coopRatio)")
    #=
    linkBars = bar(0:1:popSize, freqData)
    linkBars
    =#
end


#println("New")

#numGens is the number of generations to run
#should be 500
numGens = 500
pN = .4
pR = .03
benefit = 2.0
synergism = 0.0
cost = 0.5
delta = 0.1
mu = .001
#important variables for inheritance
#pN is the probability that infant nodes form
#edges with their parent's neighbors
#pR is the probability that infant nodes form
#edges with random strangers in the population
#benefit is the total fitness added to neighbors of
#cooperators
#synergism is the additional benefit gained by two
#cooperator who simultaneously cooperate
#cost is the fitness decrease for a cooperator
#delta is the strength of selection in the population
#mu is the probability an infant mutates to the strategy
#its parent does not possess

function runGens(gens::Int64)
    for(g) in 1:(gens * popSize)
        spliceID = Int(round(rand()*popSize+.5))
        #identifies the ID of the node that dies
        for(s) in 1:popSize
            for(ss) in 1:popSize
                if( s == spliceID || ss == spliceID )
                    edgeMatrix[s, ss] = 0
                end
            end
        end
        population[spliceID].payoff = 0
        population[spliceID].fitness = 1
        #individuals who know the individual that will soon die lose their connections to the latter

        fitSum = 0.0
        fitnesses = []
        for(i) in 1:popSize
            fitSum += population[i].fitness
        end
        for(i) in 1:popSize
            push!(fitnesses, population[i].fitness/fitSum)
        end
        fitWeights = weights(Array{Float64, 1}(fitnesses))
        momIndex = sample(1:popSize, fitWeights)

        #=
        tempPop = copy(population)
        fitted = false
        momIndex = tempPop[1].ID
        #fitSum -= population[momIndex].fitness
        while(!fitted)
            if(rand() < (population[momIndex].fitness/fitSum))
                fitted = true
                if(rand() > mu)
                    population[spliceID].strategy = population[momIndex].strategy
                else
                    population[spliceID].strategy = ((population[momIndex].strategy) * (-1)) + 1
                end
            else
                popfirst!(tempPop)
                momIndex = tempPop[1].ID
                fitSum = 0.0
                for(p) in 1:length(tempPop)
                    fitSum += tempPop[p].fitness
                end
                #fitSum -= population[momIndex].fitness
            end
        end
        =#

        #selects node to birth a new node
        #FITNESS NOW IMPLEMENTED

        if(rand() > mu)
            population[spliceID].strategy = population[momIndex].strategy
        else
            population[spliceID].strategy = ((population[momIndex].strategy) * (-1)) + 1
        end
        edgeMatrix[spliceID, momIndex] = 1
        edgeMatrix[momIndex, spliceID] = 1
        #infant node forms an edge to its parent

        for(i) in 1:popSize
            momNeighbor = false
            if(edgeMatrix[momIndex, i] == 1)
                momNeighbor = true
            end
            #separates neighbor nodes from stranger nodes

            if( i != spliceID && edgeMatrix[spliceID, i] == 0)
                if(momNeighbor && rand() < pN)
                    edgeMatrix[i, spliceID] = 1
                    edgeMatrix[spliceID, i] = 1
                elseif(!momNeighbor && rand() < pR)
                    edgeMatrix[i, spliceID] = 1
                    edgeMatrix[spliceID, i] = 1
                end
            end
            #rules out self and mother from population,
            #then randomly forms edges with neighbors
            #and strangers separately at frequencies of
            #pN and pR respectively.
        end

        #cooperation occurs
        for(i) in 1:popSize
            if(population[i].strategy == 1)
                beneficiaryNodes = []
                coopBeneficiaries = []
                for(ii) in 1:popSize
                    if(edgeMatrix[i, ii] == 1)
                        push!(beneficiaryNodes, ii)
                        if(population[ii].strategy == 1)
                            push!(coopBeneficiaries, ii)
                        end
                    end
                end
                for(b) in 1:length(beneficiaryNodes)
                    population[beneficiaryNodes[b]].payoff += (benefit/length(beneficiaryNodes))
                end
                for(c) in 1:length(coopBeneficiaries)
                    coopDegree = 0
                    for(cc) in 1:popSize
                        if(edgeMatrix[coopBeneficiaries[c], cc] == 1)
                            coopDegree += 1
                        end
                    end
                    population[coopBeneficiaries[c]].payoff += synergism/(coopDegree*length(beneficiaryNodes))
                end
                population[i].payoff -= cost
            end
        end

        #calculate fitness based on payoff from above
        for(i) in 1:popSize
                population[i].fitness = (1.0 + delta) ^ population[i].payoff
                population[i].payoff = 0
        end

        GC.gc()
        if((g > 100 * popSize) && (g % popSize == 0))
            linkBarChart()
        end
    end
end

for(x) in 1:100
    global meanCoopRatio
    #println(popSize)
    runGens(numGens)
    meanCoopRatio = meanCoopRatio/400.0
    println("Mean Cooperation Ratio: $(meanCoopRatio)")
    meanCoopRatio = 0.0
end
#println("poplength: $(length(population))")
