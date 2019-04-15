#!/usr/bin/env julia

using FileIO
using JLD2
using Plots
pyplot()

dataTypeList = ["Coop","Degs","PN","PR","Dist"]
dataMatrix = zeros(Float64, 9, 121)
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
                dataMatrix[1][(c*11)+b] = currentDict["parameters"][1]
                dataMatrix[2][(c*11)+b] = currentDict["parameters"][2]
                dataMatrix[3][(c*11)+b] = currentDict["meanCooperationRatio"]
            end
            if(i==2)
                dataMatrix[4][(c*11)+b] = currentDict["meanDegree"]
                dataMatrix[5][(c*11)+b] = currentDict["meanCooperatorDegree"]
                dataMatrix[6][(c*11)+b] = currentDict["meanDefectorDegree"]
            end
            if(i==3)
                dataMatrix[7][(c*11)+b] = currentDict["meanPN"]
            end
            if(i==4)
                dataMatrix[8][(c*11)+b] = currentDict["meanPR"]
            end
            if(i==5)
                dataMatrix[9][(c*11)+b] = currentDict["meanDistanceFromDefToCoop"]
            end
        end
    end
end

xs = [string("cLink ", i) for i = 0.0:0.1:1.0]
ys = [string("B ", i) for i = 0:10]
#PN
z = float(0:1)
heatmap(xs, ys, z, aspect_ratio = 1)
