using DataFrames
using Printf
using Statistics
using StatsBase
using Random
using Flux
using LsqFit
using GLM
using MLJLinearModels
using LinearAlgebra

hartreeToWavenumberConversion::Float64 = 219474.63

inputFile::Vector{String} = readlines()
include("PotentialEnergyModel.jl")
include("ReadInput.jl")
include("FittingModels.jl")
using .PotentialEnergyModel
using .ReadInput
using .FittingModels


println("Parsing input...")
inputBlocks::Vector{Vector{String}} = readBlocks(inputFile)
model::potentialEnergyModel = readModel(inputBlocks[1])
println()
grid::DataFrame = readGrid(inputBlocks[2], model.numberOfModes)
println()
println("End of Input")
println()
molecule::String = model.molecule
include("$(molecule).jl")


println("Minimum in the potential energy:")
minimumEnergy::Float64 = minimum(grid[:, :energy])
@printf("%12.8f %s\n", minimumEnergy, "hartree")
minimumEnergy *= hartreeToWavenumberConversion
@printf("%12.8f %s\n", minimumEnergy, "cm-1")
grid[:, :energy] .= grid[:, :energy].*hartreeToWavenumberConversion
numberOfPoints::Int64 = size(grid)[1]
println()
# println(grid[:, :energy])

# Weight factors by Partridge and Schwenke
weightOfPoints::Vector{Float64} = computeWeightOfPoint.(grid[:, :energy] .- minimumEnergy)
weightOfPoints = weightOfPoints ./ sum(weightOfPoints)

setSymmetryOperations(defineSymmetryOperations())


gridPoints::Matrix{Float64} = zeros(numberOfPoints, model.numberOfModes)
for i in 1:numberOfPoints
    gridPoints[i, :] = grid[i, :coordinates]
end
# oldParameters::Vector{Float64} = [model.equilibriumParameters; model.linearParameters]
println()
println("Beginning fit...")
newParameters::Vector{Float64} = LsqCurveFit(computePotentialOfPoints, computeDerivativesOfPoints, gridPoints, [model.equilibriumParameters; model.linearParameters], grid[:, :energy], weightOfPoints)
model = updatePotentialEnergyModel(model, newParameters)
println("Done!")
println()

grid = computePotentialOfGrid(grid, [model.equilibriumParameters; model.linearParameters])
println(minimum(grid[:, :potential]))
# println(model.linearParameters[1])
println(minimum(grid[:, :potential]) - model.linearParameters[1])
grid[:, :potential] .= grid[:, :potential] .- model.linearParameters[1]
grid[:, :energy] .= grid[:, :energy] .- model.linearParameters[1]
grid[:, :error] .= grid[:, :energy] .- grid[:, :potential]
# grid[:, :potential] .= grid[:, :potential] .- minimum(grid[:, :potential])
# grid[:, :error] = grid[:, :energy] - grid[:, :potential]
println("Grid of final energies:")
println()
printComputedEnergies(grid, minimumEnergy)
println()
print("Printing new potential parameters...")
println()
printModel(model)
print("End of Fit")
# computeDerivativesOfPoints(grid[:, :coordinates], [model.equilibriumParameters; model.linearParameters])

grid[:, :xi] .= generateXiCoordinates.(grid[:, :localModes])


# expansionCoefficients = expansionCoefficients[expansionCoefficients.expansionOrder .< 4, :]
numberOfParameters::Int64 = size(expansionCoefficients)[1]
println("Invariance of current potential...")
checkPotentialForInvariance(grid, expansionCoefficients, symmetryOperations)
println("Initializing fit coordinates...")
@time xiPowers::Matrix{Float64} = setupFitVariables(grid, symmetryOperations, expansionCoefficients)
# @time grid::DataFrame = setupFitVariables(grid, symmetryOperations, expansionCoefficients)
println("Done!")


# Ridge Regression 
# potentialEnergyModel = @formula(E ~ .)
# reducedGrid::DataFrame = hcat(grid[:, :E], grid[:, 5:end])
# normalisedGrid::DataFrame
# for i in 1:numberOfParameters + 1
#     normalisedGrid[:, i] .= (reducedGrid[:, i] .- mean(reducedGrid[:, i]))./std(reducedGrid[:, i])
# end 
# lambda::Float64 = 1.0
# println("Begin fitting...")
# @time potentialFit = lm(potentialEnergyModel, reducedGrid, RidgeReg(lambda = lambda))
# println("Done!")

# fittedParameters::Vector{Float64} = coef(potentialFit)
# error = coeftable(potentialFit).se
# fittedPotential::Vector{Float64} = predict(potentialEnergyModel, reducedGrid)
# grid[:, :fittedEnergies] .= fittedPotential
# grid[:, :obsMinusCalc] .= grid[:, :E] .- grid[:, :fittedEnergies]

# Least Squares Fit
energies::Vector{Float64} = convert(Vector, grid[:, :E])
expansionParameters::Vector{Float64} = convert(Vector, expansionCoefficients[:, :coefficientValue])
weights::Vector{Float64} = convert(Vector, grid[:, :weight])
function potentialEnergyModel(xiPowers::Matrix{Float64}, expansionParameters::Vector{Float64})::Vector{Float64}
    numberOfPoints::Int64 = size(xiPowers)[1]
    potential::Vector{Float64} = zeros(numberOfPoints)
    for i in 1:numberOfPoints
        potential[i] = dot(xiPowers[i, :], expansionParameters)
    end
    return potential
end

function derivatives(xiPowers::Matrix{Float64}, expansionParameters::Vector{Float64})::Matrix{Float64}
    return xiPowers
end

normalisedEnergies::Vector{Float64} = (energies .- mean(energies))./std(energies)
normalisedXiPowers::Matrix{Float64} = zeros(numberOfPoints, numberOfParameters)
for i in 1:numberOfParameters
    normalisedXiPowers[:, i] = (xiPowers[:, i] .- mean(xiPowers[:, i]))./std(xiPowers[:, i]) 
end
println("Begin fitting...")

Random.seed!(123)

# @time potentialFit = curve_fit((xiPowers, expansionParameters) -> potentialEnergyModel(xiPowers, expansionParameters),
#     (xiPowers, expansionParameters) -> derivatives(xiPowers, expansionParameters),
#     xiPowers, energies, weights, expansionParameters)
fittedParameters::Vector{Float64} = LsqCurveFit(potentialEnergyModel, derivatives, xiPowers, expansionParameters, energies, weights)
# fittedParameters::Vector{Float64} = expansionCoefficients[:, 3] 
# @time fittedParameters::Vector{Float64} = TikhonovRegularisation(xiPowers, expansionParameters, energies)
println("Done!")

# println(potentialFit)
# println(size(expansionParameters))
# println(size(fittedParameters))

# fittedParameters::Vector{Float64} = potentialFit.param
fittedPotential::Vector{Float64} = potentialEnergyModel(xiPowers, fittedParameters)
# fittedPotential::Vector{Float64} = potentialEnergyModel(xiPowers[:, 2:end], fittedParameters[2:end]) .+ fittedParameters[1]
grid[:, :fittedEnergies] .= fittedPotential # .- fittedParameters[1]
grid[:, :obsMinusCalc] .= grid[:, :E] .- grid[:, :fittedEnergies]
println("Parameters of fit:")
displayResult::String = ""
for i in 1:size(fittedParameters)[1]
    global displayResult = ""
    displayResult = displayResult*"$(expansionCoefficients[i, 1])"
    displayResult = displayResult*"  "
    powers = expansionCoefficients[i, 2]
    for j in 1:numberOfModes
        displayResult = displayResult*"$(powers[j])"*" "
    end
    displayResult = displayResult*"  "
    # displayResult = displayResult*"$(expansionCoefficients[i, 3])"*"  "
    displayResult = displayResult*"$(fittedParameters[i])"
    # displayResult = displayResult*"$(sigma[i])"
    println(displayResult)
end
println("Displaying energies at grid")
# for i in 1:size(grid)[1]
#     global displayResult = ""
#     gridCoordinates = grid[i, 1]
#     for j in 1:numberOfModes
#         displayResult = displayResult*"$(gridCoordinates[j])"*" "
#     end
#     displayResult = displayResult*" "
#     displayResult = displayResult*"$(grid[i, :E])"
#     displayResult = displayResult*"  "
#     displayResult = displayResult*"$(grid[i, :fittedEnergies])"
#     displayResult = displayResult*"  "
#     displayResult = displayResult*"$(grid[i, :obsMinusCalc])"
#     println(displayResult)
# end
println(grid)
println("Total rms:")
println(sqrt(mean(grid[:, :obsMinusCalc].^2)))
println("rms for energies below 10000 cm-1:")
println(sqrt(mean(filter(row -> row.E <= 10000, grid)[:, :obsMinusCalc].^2)))
println("rms for energies below 5000 cm-1:")
println(sqrt(mean(filter(row -> row.E <= 5000, grid)[:, :obsMinusCalc].^2)))
println("rms for energies below 1000 cm-1:")
println(sqrt(mean(filter(row -> row.E <= 1000, grid)[:, :obsMinusCalc].^2)))
println("Error in parameters: ")
println(error)
