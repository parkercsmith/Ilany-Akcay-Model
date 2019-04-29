#!/usr/bin/env julia

using JLD2
using ArgParse
using StatsBase
#using FileIO
#using Plots

mutable struct NetworkParameters

    #measurement data
    meanCoopRatio::Float64
    meanProbNeighbor::Float64
    meanProbRandom::Float64
    meanDegree::Float64
    meanCoopDegree::Float64
    meanDefDegree::Float64
    meanCoopDefDistance::Float64

    #Node characteristics
    popPN::Array{Float64, 1}
    popPR::Array{Float64, 1}
    popStrategies::Array{Int64, 1}
    popPayoff::Array{Float64, 1}
    popFitness::Array{Float64, 1}

    #structural variables
    numGens::Int64
    popSize::Int64
    edgeMatrix::Array{Int64, 2}

    #social variables
    cost::Float64
    benefit::Float64
    synergism::Float64
    linkCost::Float64
    mu::Float64
    delta::Float64

    function NetworkParameters(b::Float64, cL::Float64)

        numGens = 100000
        popSize = 100
        popPN = zeros(Float64, popSize)
        popPN[:] .= 0.5
        popPR = zeros(Float64, popSize)
        popPR[:] .= 0.0001
        popStrategies = zeros(Int64, popSize)
        popStrategies[2:2:popSize] .= 1
        popFitness = zeros(Float64, popSize)
        popFitness[:] .= 1.0

        edgeMatrix = zeros(Int64, 100, 100)
        cost = 1.0
        synergism = 0.0
        benefit = b
        linkCost = cL
        mu = .01
        delta = 0.5

        new(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, popPN, popPR, popStrategies, zeros(Float64, popSize), popFitness, numGens, popSize, edgeMatrix, cost, benefit, synergism, linkCost, mu, delta)
    end
end

#measurement functions
function coopRatio(network::NetworkParameters)
    coopCount = 0
    for(i) in 1:network.popSize
        if(network.popStrategies[i] == 1)
            coopCount += 1
        end
    end
    coopCount /= network.popSize
    network.meanCoopRatio += coopCount
end

function probNeighbor(network::NetworkParameters)
    pNTotal = 0.0
    for(i) in 1:network.popSize
        pNTotal += network.popPN[i]
    end
    pNTotal /= network.popSize
    network.meanProbNeighbor += pNTotal
end

function probRandom(network::NetworkParameters)
    pRTotal = 0.0
    for(i) in 1:network.popSize
        pRTotal += network.popPR[i]
    end
    pRTotal /= network.popSize
    network.meanProbRandom += pRTotal
end

function degrees(network::NetworkParameters)
    cooperatorsPresent = 0
    degTotal = 0
    coopDegTotal = 0
    defDegTotal = 0
    for(i) in 1:network.popSize
        degCounter = 0
        for(ii) in 1:network.popSize
            if(network.edgeMatrix[i, ii] == 1)
                degCounter += 1
            end
        end
        degTotal += degCounter
        if(network.population[i].strategy == 1)
            coopDegTotal += degCounter
            cooperatorsPresent += 1
        else
            defDegTotal += degCounter
        end
    end
    degTotal /= network.popSize
    defDegTotal /= (network.popSize - cooperatorsPresent)
    coopDegTotal /= cooperatorsPresent
    network.meanDegree += degTotal
    network.meanCoopDegree += coopDegTotal
    network.meanDefDegree += defDegTotal
end

function distance(network::NetworkParameters)
    distanceTotal = 0.0
    for(i) in 1:network.popSize
        found = false
        usualSuspects = zeros(Int64, network.popSize)
        susCounter = 1
        distCount = 0
        for(ii) in 1:network.popSize
            if(network.edgeMatrix[i, ii] == 1)
                usualSuspects[susCounter] = ii
                susCounter += 1
            end
        end
        while(!found)
            distCount += 1
            dC = distCount
            if(dC == 1)
                for(s) in 1:length(usualSuspects)
                    if(usualSuspects[s] != 0)
                        if(network.popStrategies[usualSuspects[s]] != over.popStrategies[i])
                            found = true
                        end
                    end
                end
            else
                for(s) in 1:length(usualSuspects)
                    if(usualSuspects[s] != 0)
                        for(ii) in 1:network.popSize
                            if(over.edgeMatrix[usualSuspects[s], ii] == 1 && !(ii in usualSuspects))
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
    distanceTotal /= network.popSize
    network.meanCoopDefDistance += distanceTotal
end

#evolution functions
function death(network::NetworkParameters)
    deadID = sample(1:network.popSize)
    network.edgeMatrix[deadID, :] .= 0
    network.edgeMatrix[:, deadID] .= 0
    network.popPayoff[deadID] = 0
    network.popFitness[deadID] = 1
    deadID
end

function findMom(network::NetworkParameters)
    fitSum = sum(network.popFitness)
    fitnesses = zeros(over.popSize)
    for(i) in 1:network.popSize
        fitnesses[i] = network.popFitness[i]/fitSum
    end
    fitWeights = weights(Array{Float64, 1}(fitnesses))
    momIndex = sample(1:network.popSize, fitWeights)
    momIndex
end

function birth(network::NetworkParameters, child::Int64, parent::Int64)
    network.popStrategies[child] = network.popStrategies[parent]
    if(rand()<network.mu)
        network.popStrategies[child] -= 1
        network.popStrategies[child] *= -1
    end
    network.popPN[child] = network.popPN[parent]
    if(rand()<network.mu)
        network.popPN[child] += randn()/100
        if(network.popPN[child] > 1.0)
            network.popPN[child] = 1.0
        elseif(network.popPN[child] < 0.0)
            network.popPN[child] = 0.0
        end
    end
    network.popPR[child] = network.popPR[parent]
    if(rand()<network.mu)
        network.popPR[child] += randn()/100
        if(network.popPR[child] > 1.0)
            network.popPR[child] = 1.0
        elseif(network.popPR[child] < 0.0)
            network.popPR[child] = 0.0
        end
    end
    network.edgeMatrix[parent, child] = 1
    network.edgeMatrix[child, parent] = 1

    for(i) in 1:network.popSize
        if(i != child && network.edgeMatrix[i, child] == 0)
            if(network.edgeMatrix[i, parent] == 1)
                if(rand() < network.popPN[child])
                    network.edgeMatrix[i, child] = 1
                    network.edgeMatrix[child, i] = 1
                end
            else()
                if(rand() < network.popPR[child])
                    network.edgeMatrix[i, child] = 1
                    network.edgeMatrix[child, i] = 1
                end
            end
        end
    end
end

function cooperate(network::NetworkParameters)
    for(i) in 1:network.popSize
        if(network.popStrategies[i] == 1)
            B = network.benefit/getDegree(network, i)
            D = network.synergism/getDegree(network, i)

            for(ii) in 1:network.popSize
                if(network.edgeMatrix[i, ii] == 1)
                    network.popPayoff[ii] += B
                    if(network.popStrategies[ii] == 1)
                        network.popPayoff[i] += D/getDegree(network, ii)
                        network.popPayoff[ii] += D/getDegree(network, ii)
                    end
                end
            end

            network.popPayoff[i] -= network.cost
        end
        network.popPayoff[i] -= network.linkCost * getDegree(network, i)
    end
end

function resolveFitnesses(network::NetworkParameters)
    for(i) in 1:network.popSize
        network.popFitness[i] = (1.0 + network.delta) ^ network.popPayoff[i]
    end
    network.popPayoff[:] .= 0.0
end

function getDegree(network::NetworkParameters, node::Int64)
    sum(network.edgeMatrix[node, :])
end

function run(CL::Float64, BEN::Float64)
    #=
    [1]finalMeanPN
    [2]finalMeanPR
    [3]finalMeanDegree
    [4]finalMeanCoopDegree
    [5]finalMeanDefDegree
    [6]finalMeanDistance
    [7]finalMeanCoopRatio
    =#
    dataArray = zeros(7)
    repSims = 10
    for(x) in 1:repSims

        #initializes globalstuff structure with generic constructor
        network = NetworkParameters(BEN, CL)

        #checks efficiency of simulation while running it
        for(g) in 1:(network.numGens * network.popSize)

            childID = death(network)
            parentID = findMom(network)
            birth(network, childID, parentID)
            cooperate(network)
            resolveFitnesses(network)

            if(g > (network.numGens * network.popSize / 5) && (g % network.popSize) == 0)
                coopRatio(network)
                probNeighbor(network)
                probRandom(network)
                degrees(network)
                distance(network)
            end
        end

        #divides meanCooperationRatio by last 400 generations to get a true mean, then outputs
        network.meanProbNeighbor /= 80000.0
        network.meanProbRandom /= 80000.0
        network.meanDegree /= 80000.0
        network.meanCoopDegree /= 80000.0
        network.meanDefDegree /= 80000.0
        network.meanCoopRatio /= 80000.0
        network.meanCoopDefDistance /= 80000.0

        dataArray[1] += overlord.meanProbNeighbor
        dataArray[2] += overlord.meanProbRandom
        dataArray[3] += overlord.meanDegree
        dataArray[4] += overlord.meanCoopDegree
        dataArray[5] += overlord.meanDefDegree
        dataArray[6] += overlord.meanCoopDefDistance
        dataArray[7] += overlord.meanCoopRatio
    end
    dataArray[:] ./= Float64(repSims)
    save("expData_CL$(CL)_B$(BEN).jld2", "parameters", [CL, BEN], "meanPN", dataArray[1], "meanPR", dataArray[2], "meanDegree", dataArray[3], "meanCooperatorDegree", dataArray[4], "meanDefectorDegree", dataArray[5], "meanDistanceFromDefToCoop", dataArray[6], "meanCooperationRatio", dataArray[7])
end

argTab = ArgParseSettings(description = "arguments and stuff, don't worry about it")
@add_arg_table argTab begin
    "--cLink"
        arg_type = Float64
        default = 0.0
end
parsedArgs = parse_args(ARGS, argTab)
currCostLink = parsedArgs["cLink"]
for(b) in 0:10
    currBenefit = Float64(b)
    run(currCostLink, currBenefit)
end
