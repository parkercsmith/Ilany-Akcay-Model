#!/usr/bin/env julia

using FileIO
using JLD2
using Plots
pyplot()

dataTypeList = ["Coop","Degs","PN","PR","Dist"]
dataMatrix = zeros(Float64, 11, 11, 9)
#[1] = linking cost
#[2] = benefit value
#[3] = cooperation ratio
#[4] = mean degree
#[5] = mean cooperator degree
#[6] = mean defector degree
#[7] = pN
#[8] = pR
#[9] = mean distance from cooperator to defector


for (c) in 0:10
    for (b) in 0:10
        for(i) in 1:5
            currentDict = load("expData"*dataTypeList[i]*"_CL$(Float64(c)/10.0)_B$(b).0.jld2")

            if(i==1)
                dataMatrix[c+1, b+1, 1] = currentDict["parameters"][1]
                dataMatrix[c+1, b+1, 2] = currentDict["parameters"][2]
                dataMatrix[c+1, b+1, 3] = currentDict["meanCooperationRatio"]
            end
            if(i==2)
                dataMatrix[c+1, b+1, 4] = currentDict["meanDegree"]
                dataMatrix[c+1, b+1, 5] = currentDict["meanCooperatorDegree"]
                dataMatrix[c+1, b+1, 6] = currentDict["meanDefectorDegree"]
            end
            if(i==3)
                dataMatrix[c+1, b+1, 7] = currentDict["meanPN"]
            end
            if(i==4)
                dataMatrix[c+1, b+1, 8] = currentDict["meanPR"]
            end
            if(i==5)
                dataMatrix[c+1, b+1, 9] = currentDict["meanDistanceFromDefToCoop"]
            end
        end
    end
end

xs = [string("cLink ", i) for i = 0.0:0.1:1.0]
ys = [string("B ", i) for i = 0:10]
#PN
z = dataMatrix[1:11,1:11,3]
heatmap(xs, ys, z, aspect_ratio = 1)
