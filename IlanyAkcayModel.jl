population = []
#represents all nodes in population
edgeMatrix = []
#matrix; holds all edges in population where row and
#col are indices of nodes

popSize = 100
let population = population
    for(i) in 1:popSize
        push!(population, i)
    end
end
#initialize a population of popSize nodes; each node
#contains only its index

for(p) in 1:popSize
    push!(edgeMatrix, [])
    for(pp) in 1:popSize
        push!(edgeMatrix[p], 0)
    end
end
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

numInitEdges = round(1.5*length(population))
let population = population
    repairedLinks = 0
    for(i) in 1:numInitEdges
            firstID = Int(round(rand()*length(population)+.5))
            secondID = firstID
            while(secondID == firstID)
                secondID = Int(round(rand()*length(population)+.5))
                if(edgeMatrix[firstID][secondID] == 1 || edgeMatrix[secondID][firstID] ==1)
                    secondID = firstID
                end
            end
            edgeMatrix[firstID][secondID] = 1
            edgeMatrix[secondID][firstID] = 1
    end
end
#identifies a number of initial edges based on
#popSize, creates numInitEdges distinct pairs of
#distinct indices, and changes value of corresponding
#edgeMatrix placeholders to 1

using Plots
function linkBarChart()
    freqData = []
    for(i) in 1:(length(edgeMatrix))
        counter = 0
        for(ii) in 1:length(edgeMatrix)
            degreeMarker = 0
            for(iii) in 1:length(edgeMatrix[ii])
                if(edgeMatrix[ii][iii] == 1)
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
    for(i) in 1:length(edgeMatrix)
        noEdges = true
        for(ii) in 1:length(edgeMatrix[i])
            if(edgeMatrix[i][ii] == 1)
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
    linkBars = bar(0:1:100, freqData)
    linkBars
end


#println("New")

#Running generations to stabilize
#Goal is 20*popSize
numGens = 2000

for(g) in 1:numGens
    spliceID = Int(round(rand()*length(population)+.5))
    #identifies the ID of the node that dies
    for(s) in 1:length(edgeMatrix)
        for(ss) in 1:length(edgeMatrix[s])
            if( s == spliceID || ss == spliceID )
                edgeMatrix[s][ss] = 0
            end
        end
    end
    #individuals who know the individual that will soon die lose their connections to the latter

    momIndex = Int(round(rand()*length(population)+.5))
    #selects node to birth a new node
    pN = .7
    pR = .3
    #variables for inheritance
    #pN is the probability that infant nodes form
    #edges with their parent's neighbors
    #pR is the probability that infant nodes form
    #edges with random strangers in the population

    edgeMatrix[spliceID][momIndex] = 1
    edgeMatrix[momIndex][spliceID] = 1
    #infant node forms an edge to its parent

    for(i) in 1:length(edgeMatrix)
        momNeighbor = false
        if(edgeMatrix[momIndex][i] == 1)
            momNeighbor = true
        end
        #separates neighbor nodes from stranger nodes

        if( i != spliceID && edgeMatrix[spliceID][i] == 0)
            if(momNeighbor && rand() < pN)
                edgeMatrix[i][spliceID] = 1
                edgeMatrix[spliceID][i] = 1
            elseif(!momNeighbor && rand() < pR)
                edgeMatrix[i][spliceID] = 1
                edgeMatrix[spliceID][i] = 1
            end
        end
        #rules out self and mother from population,
        #then randomly forms edges with neighbors
        #and strangers separately at frequencies of
        #pN and pR respectively.
    end
    GC.gc()
end

#println("poplength: $(length(population))")
linkBarChart()
