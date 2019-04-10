#!/usr/bin/env julia

using StatsBase
using Plots
using ArgParse
using FileIO
using JLD2

#Node structure represents individuals in population
mutable struct Node
    strategy::Int8 #1 if cooperator, 0 if defector
    payoff::Float64
    fitness::Float64
    ID::Int64
    pN::Float64
    pR::Float64
end

#globalstuff structure exists in one instance; contains all universal variables
mutable struct globalstuff
    popSize::Int64 #number of individuals in population
    edgeMatrix::Array{Int64, 2} #matrix of binary Ints; 1 where connection between row and col index is present
    meanCoopRatio::Float64 #summed over 500 generations, then calculated after core function
    meanProbNeighbor::Float64
    meanProbRandom::Float64
    meanDegree::Float64
    meanCoopDegree::Float64
    meanDefDegree::Float64
    meanCoopDefDistance::Float64
    population::Array{Node, 1} #contains 100 Node structs
    numInitialCoops::Int64 #number of cooperators at initialization of globalstuff; used once per run
    numInitEdges::Int64 #number of edges in network at initialization of globalstuff; used once per run
    numGens::Int64 #generations per run; currently last 400 add into meanCoopRatio
    pN::Float64 #probability that infant connects to parent's neighbor
    pR::Float64 #probability that infant connects to random Node
    #pB is assumed to be 1 in this simulation
    benefit::Float64 #payoff to be distributed amongst cooperator's neighbors
    synergism::Float64 #payoff for mutual cooperators; inversely related to degrees
    cost::Float64 #payoff lost by cooperators
    delta::Float64 #strength of selection
    mu::Float64 #chance of strategy mutation in infant
    cLink::Float64 #payoff lost upon connection of two nodes

    #generic constructor takes no parameters; uses same globals for each run
    #CHANGE VARIABLES HERE
    function globalstuff(ben::Float64, cL::Float64)
        #basic variables for measurement
        popSize = 100
        meanProbNeighbor = 0.0
        meanProbRandom = 0.0
        meanDegree = 0.0
        meanCoopRatio = 0.0
        meanCoopDegree = 0.0
        meanDefDegree = 0.0
        meanCoopDefDistance = 0.0

        #initialization factors (may need to be adjusted)
        numInitialCoops = popSize/2
        numInitEdges = Int(round(1.5*popSize))

        #specific details of simulation
        numGens = 100000
        pN = .75 #NOT IN USE; see Node pN
        pR = .0001 #NOT IN USE; see Node pR
        benefit = ben
        synergism = 0.0
        cost = -0.5
        delta = 0.5
        mu = .01
        cLink = cL

        #defines population as a Node Array with 100 spots, then fills it with cooperators at beginning and defectors afterwards
        population = Array{Node, 1}(undef, popSize)
        for(i) in 1:popSize
            population[i] = Node(0, 0, 1, i, 0.75, 0.0001)#adjust pn and pr here
            if(numInitialCoops > 0)
                population[i].strategy = 1
                numInitialCoops -= 1
            end
        end

        #empties edgeMatrix, then finds pairs of distinct indices to connect
        edgeMatrix = zeros(Int64, popSize, popSize)
        for(i) in 1:numInitEdges
                firstID = Int(round(rand()*popSize+.5))
                secondID = firstID
                while(secondID == firstID)
                    secondID = Int(round(rand()*popSize+.5))
                    if(edgeMatrix[firstID, secondID] == 1 || edgeMatrix[secondID, firstID] ==1)
                        secondID = firstID
                        continue
                    end
                end
                edgeMatrix[firstID, secondID] = 1;
                edgeMatrix[secondID, firstID] = 1;
        end

        #constructs globalstuff structure
        new(popSize, edgeMatrix, meanCoopRatio, meanProbNeighbor, meanProbRandom, meanDegree, meanCoopDegree, meanDefDegree, meanCoopDefDistance, population, numInitialCoops, numInitEdges, numGens, pN, pR, benefit, synergism, cost, delta, mu, cLink)
    end
end

#measurement function; name and contents regularly change
#CURRENTLY: counts cooperators and adds the frequency of cooperation to meanCoopRatio
function countCoops(over::globalstuff)
    pNTotal = 0.0
    for(i) in 1:over.popSize
        pNTotal += over.population[i].pN
    end
    pNTotal /= over.popSize
    over.meanProbNeighbor += pNTotal
    #= THIS SEGMENT CHANGES FOR EACH TYPE OF DATA
    pRTotal = 0.0
    for(i) in 1:over.popSize
        pRTotal += over.population[i].pR
    end
    pRTotal /= over.popSize
    over.meanProbRandom += pRTotal

    coopTotal = 0.0
    for(i) in 1:over.popSize
        if(over.population[i].strategy == 1)
            coopTotal += 1.0
        end
    end
    coopTotal /= over.popSize
    over.meanCoopRatio += coopTotal

    degTotal = 0.0
    coopDegTotal = 0.0
    defDegTotal = 0.0
    for(i) in 1:over.popSize
        degCounter = 0
        for(ii) in 1:over.popSize
            if(over.edgeMatrix[i, ii] == 1)
                degCounter += 1
            end
        end
        degTotal += degCounter
        if(over.population[i].strategy == 1)
            coopDegTotal += degCounter
        else
            defDegTotal += degCounter
        end
    end
    degTotal /= over.popSize
    defDegTotal /= over.popSize
    coopDegTotal /= over.popSize
    over.meanDegree += degTotal
    over.meanCoopDegree += coopDegTotal
    over.meanDefDegree += defDegTotal

    distanceTotal = 0.0
    for(i) in 1:over.popSize
        found = false
        usualSuspects = zeros(Int64, over.popSize * over.popSize)
        susCounter = 1
        distCount = 1
        for(ii) in 1:over.popSize
            if(over.edgeMatrix[i, ii] == 1)
                usualSuspects[susCounter] = ii
                susCounter += 1
            end
        end
        while(!found)
            dC = distCount
            if(dC == 1)
                for(s) in 1:length(usualSuspects)
                    if(usualSuspects[s] != 0)
                        if(over.population[usualSuspects[s]].strategy != over.population[i].strategy)
                            found = true
                        end
                    end
                end
            else
                for(s) in 1:length(usualSuspects)
                    if(usualSuspects[s] != 0)
                        for(ii) in 1:over.popSize
                            if(over.edgeMatrix[usualSuspects[s], ii] == 1)
                                usualSuspects[susCounter] = ii
                                susCounter += 1
                            end
                        end
                    end
                end
                dC -= 1
            end
        end
        distanceTotal += distCount
    end
    distanceTotal /= over.popSize
    over.meanCoopDefDistance += distanceTotal
    =#

    #vestigial bar chart formatting (may be used later)
    #=
    linkBars = bar(0:1:popSize, freqData)
    linkBars
    =#
end

#core function; conducts 100 birth-death events per generation with selection for high payoff
#accepts globalstuff struct as a parameter
function runGens(over::globalstuff)
    for(g) in 1:(over.numGens * over.popSize)

        #randomly selects a node to die and resets its connections, payoff and fitness
        spliceID = Int(round(rand()*over.popSize+.5))
        over.edgeMatrix[spliceID, :] .= 0
        over.edgeMatrix[:, spliceID] .= 0
        over.population[spliceID].payoff = 0
        over.population[spliceID].fitness = 1

        #calculates relative fitnesses of nodes, creates a weight vector, and samples a parent with selection weight
        fitSum = 0.0
        fitnesses = zeros(over.popSize)
        for(i) in 1:over.popSize
            fitSum += over.population[i].fitness
        end
        for(i) in 1:over.popSize
            fitnesses[i] = over.population[i].fitness/fitSum
        end
        fitWeights = weights(Array{Float64, 1}(fitnesses))
        momIndex = sample(1:over.popSize, fitWeights)

        #mutates 1 in 1000 infants, then connects infant to parent
        if(rand() > over.mu)
            over.population[spliceID].pN = over.population[momIndex].pN
        else
            mutAddend = randn()/100.0
            over.population[spliceID].pN = over.population[momIndex].pN + mutAddend
            if(over.population[spliceID].pN > 1)
                over.population[spliceID].pN = 1
            elseif(over.population[spliceID].pN < 0)
                over.population[spliceID].pN = 0
            end
            #over.population[spliceID].pN = rand()
        end
        if(rand() > over.mu)
            over.population[spliceID].pR = over.population[momIndex].pR
        else
            mutAddend = randn()/100.0
            over.population[spliceID].pR = over.population[momIndex].pR + mutAddend
            if(over.population[spliceID].pR > 1)
                over.population[spliceID].pR = 1
            elseif(over.population[spliceID].pR < 0)
                over.population[spliceID].pR = 0
            end
        end
        if(rand() > over.mu)
            over.population[spliceID].strategy = over.population[momIndex].strategy
        else
            over.population[spliceID].strategy = ((over.population[momIndex].strategy) * (-1)) + 1
        end
        over.edgeMatrix[spliceID, momIndex] = 1
        over.edgeMatrix[momIndex, spliceID] = 1

        #formation of edges from infant to various nodes following social inheritance
        for(i) in 1:over.popSize

            #identifies parent neighbors in population, then connects to neighbors with prob pN and strangers with prob pR
            momNeighbor = false
            if(over.edgeMatrix[momIndex, i] == 1)
                momNeighbor = true
            end
            if( i != spliceID && over.edgeMatrix[spliceID, i] == 0)
                if(momNeighbor && rand() < over.population[spliceID].pN)
                    over.edgeMatrix[i, spliceID] = 1
                    over.edgeMatrix[spliceID, i] = 1
                elseif(!momNeighbor && rand() < over.population[spliceID].pR)
                    over.edgeMatrix[i, spliceID] = 1
                    over.edgeMatrix[spliceID, i] = 1
                end
            end
        end

        #cooperation occurs at every cooperator node
        for(i) in 1:over.popSize
            if(over.population[i].strategy == 1)

                #counts degree of cooperator and mutual cooperators
                benCount = 0
                coopCount = 0
                for(ii) in 1:over.popSize
                    if(over.edgeMatrix[i, ii] == 1)
                        benCount+=1
                        if(over.population[ii].strategy == 1)
                            coopCount+=1
                        end
                    end
                end

                #initializes then fills arrays of indices of beneficiary neighbors and mutually cooperating neighbors
                beneficiaryNodes = collect(1:benCount)
                benCount = 1
                coopBeneficiaries = collect(1:coopCount)
                coopCount = 1
                for(ii) in 1:over.popSize
                    if(over.edgeMatrix[i, ii] == 1)
                        beneficiaryNodes[benCount] = ii
                        benCount += 1
                        if(over.population[ii].strategy == 1)
                            coopBeneficiaries[coopCount] = ii
                            coopCount+=1
                        end
                    end
                end

                #gives payoff to beneficiaries relative to degree of cooperator as well as synergistic benefits for mutual cooperators
                benCount -= 1
                coopCount -= 1
                for(b) in 1:(benCount)
                    over.population[beneficiaryNodes[b]].payoff += (over.benefit/(benCount))
                end
                for(c) in 1:(coopCount)
                    coopDegree = 0
                    for(cc) in 1:over.popSize
                        if(over.edgeMatrix[coopBeneficiaries[c], cc] == 1)
                            coopDegree += 1
                        end
                    end
                    over.population[coopBeneficiaries[c]].payoff += over.synergism/(coopDegree*(benCount))
                end

                #subtracts cost from identified cooperator node's payoff
                over.population[i].payoff -= over.cost
            end

            #nodes pay for all of their edges at cLink
            for(ii) in 1:over.popSize
                if(over.edgeMatrix[i, ii] == 1 || over.edgeMatrix[ii, i] == 1)
                    over.population[i].payoff -= over.cLink
                end
            end
        end

        #calculates fitnesses based on payoffs per Akcay then resets payoffs
        for(i) in 1:over.popSize
                over.population[i].fitness = (1.0 + over.delta) ^ over.population[i].payoff
                over.population[i].payoff = 0
        end

        #counts cooperators in last 400 generations after 100 birth-death events
        if((g > 200 * over.popSize) && (g % over.popSize == 0))
            countCoops(over)
        end
    end
end

argTab = ArgParseSettings(description = "arguments and stuff, don't worry about it")
@add_arg_table argTab begin
    "--cLink"
        arg_type = Float64
        default = 0.0
end
parsedArgs = parse_args(ARGS, argTab)
costLink = parsedArgs["cLink"]

for(b) in 0:1:10 #edit here
    #=benWeights = zeros(1)
    benWeights[b] = 1
    benWeights = weights(Array{Float64, 1}(benWeights))
    currBenefit = sample([3.0] , benWeights)=#
    currBenefit = b
    currBenefit = Float64(currBenefit)
    #replicates data for 100 simulations
    finalMeanPN = 0.0
    benVal = 0.0
    repSims = 10
    for(x) in 1:repSims

        #initializes globalstuff structure with generic constructor
        overlord = globalstuff(currBenefit, costLink)
        benVal = overlord.benefit

        #checks efficiency of simulation while running it
        runGens(overlord)

        #divides meanCooperationRatio by last 400 generations to get a true mean, then outputs
        overlord.meanProbNeighbor = overlord.meanProbNeighbor/80000.0
        #= THIS SEGMENT CHANGES FOR EACH TYPE OF DATA
        overlord.meanProbRandom = overlord.meanProbRandom/80000.0
        overlord.meanDegree = overlord.meanDegree/80000.0
        overlord.meanCoopDegree = overlord.meanCoopDegree/80000.0
        overlord.meanDefDegree = overlord.meanDefDegree/80000.0
        overlord.meanCoopRatio = overlord.meanCoopRatio/80000.0
        overlord.meanCoopDefDistance = overlord.meanCoopDefDistance/80000.0
        =#

        if(x==1)
        #    println("Simulation at popSize = $(overlord.popSize) and cLink = $(overlord.cLink)")
        end
        finalMeanPN += overlord.meanProbNeighbor
        #= THIS SEGMENT CHANGES FOR EACH TYPE OF DATA
        finalMeanPR += overlord.meanProbRandom
        finalMeanDegree += overlord.meanDegree
        finalMeanCoopDegree += overlord.meanCoopDegree
        finalMeanDefDegree += overlord.meanDefDegree
        finalMeanDistance += overlord.meanCoopDefDistance
        finalMeanCoopRatio += overlord.meanCoopRatio
        =#
    end
    finalMeanPN /= 10.0
    #= THIS SEGMENT CHANGES FOR EACH TYPE OF DATA
    finalMeanPR /= 10.0
    finalMeanDegree /= 10.0
    finalMeanCoopDegree /= 10.0
    finalMeanDefDegree /= 10.0
    finalMeanDistance /= 10.0
    finalMeanCoopRatio /= 10.0
    =#

    #=
    popSizeStr = "$(currPopSize)"
    while(length(popSizeStr)<4)
        popSizeStr = "0" * popSizeStr
    end
    =#
    save("expDataPN_CL$(costLink)_B$(benVal).jld2", "parameters", [costLink, benVal], "meanPN", finalMeanPN#=, "meanPR", finalMeanPR, "meanDegree", finalMeanDegree, "meanDefectorDegree", finalMeanDefDegree, "meanCooperatorDegree", finalMeanCoopDegree, "meanDistanceFromDefToCoop", finalMeanDistance, "meanCooperationRatio", finalMeanCoopRatio=#)
end
