using ..PotentialEnergyModel

global equilibriumParametersInFitGlobal::Vector{Int64} = model.equilibriumParametersInFit
global linearParametersInFitGlobal::Vector{Int64} = model.linearParametersInFit
global powersGlobal::Matrix{Int64} = model.powers
global symmetryOperationsGlobal::Array{Float64}

function setParametersGlobal(model::potentialEnergyModel)
    global equilibriumParametersInFitGlobal = model.equilibriumParametersInFit
    global linearParametersInFitGlobal = model.linearParametersInFit
    global powersGlobal = model.powers
end

function defineSymmetryOperations(case::String="C3v(M)")::Array{Float64}
    symmetryOperations::Array{Float64} = zeros(Float64, 6, 6, 6)
    if case =="C3v(M)"
        symmetryOperations[1, :, :] = Matrix(1I, 6, 6)
        # (123)
        symmetryOperations[2, :, :] = [
            0 1 0 0 0 0; 
            0 0 1 0 0 0; 
            1 0 0 0 0 0;
            0 0 0 -1/2 sqrt(3)/2 0;
            0 0 0 -sqrt(3)/2 -1/2 0;
            0 0 0 0 0 1;            
            ]
        # (132)
        symmetryOperations[3, :, :] = [
            0 0 1 0 0 0;
            1 0 0 0 0 0;
            0 1 0 0 0 0;
            0 0 0 -1/2 -sqrt(3)/2 0;
            0 0 0 sqrt(3)/2 -1/2 0;
            0 0 0 0 0 1;
        ]
        # (12)*
        symmetryOperations[4, :, :] = [
            0 1 0 0 0 0;
            1 0 0 0 0 0;
            0 0 1 0 0 0;
            0 0 0 -1/2 sqrt(3)/2 0;
            0 0 0 sqrt(3)/2 1/2 0;
            0 0 0 0 0 1;
        ]
        # (23)*
        symmetryOperations[5, :, :] = [    
            1 0 0 0 0 0;
            0 0 1 0 0 0;
            0 1 0 0 0 0;
            0 0 0 1 0 0;
            0 0 0 0 -1 0;
            0 0 0 0 0 1;
        ]
        # (13)*
        symmetryOperations[6, :, :] = [  
            0 0 1 0 0 0;
            0 1 0 0 0 0;
            1 0 0 0 0 0;
            0 0 0 -1/2 -sqrt(3)/2 0;
            0 0 0 -sqrt(3)/2 1/2 0;
            0 0 0 0 0 1;
        ]
    elseif case == "D3h(M)"
        symmetryOperations = zeros(Float64, 3, 6, 6)
        symmetryOperations[1, :, :] = Matrix(1I, 6, 6)
        # (23)
        symmetryOperations[2, :, :] = [
            1 0 0 0 0 0;
            0 0 1 0 0 0;
            0 1 0 0 0 0;
            0 0 0 1 0 0;
            0 0 0 0 1 0;
            0 0 0 0 0 -1;        
            ]
        # (123)
        symmetryOperations[3, :, :] = [
            0 1 0 0 0 0;
            0 0 1 0 0 0;
            1 0 0 0 0 0;
            0 0 0 -1/2 sqrt(3)/2 0;
            0 0 0 -sqrt(3)/2 -1/2 0;
            0 0 0 0 0 1;
            ]
    end
    return symmetryOperations
end

function setSymmetryOperations(symmetryOperations::Array{Float64})
    global symmetryOperationsGlobal = symmetryOperations
end

function generateXiMatrix(grid::DataFrame, model::potentialEnergyModel, symmetryOperations::Array{Float64})::Matrix{Float64}
    convertToRadians::Float64 = pi/180.00
    numberOfGridPoints::Int64 = size(grid)[1]
    degreesOfFreedom::Int64 = model.numberOfModes
    numberOfLinearParameters::Int64 = size(model.linearParameters)[1]
    numberOfSymmetryOperations::Int64 = size(symmetryOperations)[1]
    equilibriumBondLength::Float64 = model.equilibriumParameters[1]
    morseParameter::Float64 = model.equilibriumParameters[2]
    
    xiMatrix::Matrix{Float64} = zeros(numberOfGridPoints, numberOfLinearParameters)

    for i in 1:numberOfGridPoints
        coordinates::Vector{Float64} = grid[1, i]
        xiCoordinates::Vector{Float64} = zeros(degreesOfFreedom)
        xiCoordinates[1:3] = 1.0 .- exp.(-morseParameter.*(coordinates[1:3] .- equilibriumBondLength))
        xiCoordinates[4:5] = convertToRadians.*coordinates[4:5]
        xiCoordinates[6] = 1.0 - sin(coordinates[6]*convertToRadians)
        for j in 1:numberOfLinearParameters
            powers::Vector{Int64} = model.powers[j, :]
            for k in 1:numberOfSymmetryOperations
                xiMatrix[i, j] += prod((symmetryOperations[k, :, :]*xiPowers).^powers)
            end
            xiMatrix /= numberOfSymmetryOperations
        end
    end
    return xiMatrix
end

function potential(coordinates::Vector{Float64}, parameters::Vector{Float64})::Float64
    convertToRadians::Float64 = pi/180.00
    degreesOfFreedom::Int64 = size(coordinates)[1]
    numberOfEquilibriumParameters::Int64 = size(equilibriumParametersInFitGlobal)[1]
    numberOfLinearParameters::Int64 = size(linearParametersInFitGlobal)[1]
    symmetryOperations::Array{Float64} = symmetryOperationsGlobal
    numberOfSymmetryOperations::Int64 = size(symmetryOperations)[1]
    equilibriumBondLength::Float64 = parameters[1]
    morseParameter::Float64 = parameters[2]
    powersMatrix::Matrix{Int64} = powersGlobal
    
    xi::Vector{Float64} = zeros(degreesOfFreedom)
    xi[1:3] = 1.0 .- exp.(-morseParameter.*(coordinates[1:3] .- equilibriumBondLength))
    xi[4:5] = convertToRadians.*coordinates[4:5]
    xi[6] = 1.0 - sin(coordinates[6]*convertToRadians)
    potential::Float64 = 0
    for i in 1:numberOfLinearParameters
        powers::Vector{Int64} = powersMatrix[i, :]
        xiNew = 0
        for j in 1:numberOfSymmetryOperations
            xiTransform::Float64 = prod((symmetryOperations[j, :, :]*xi).^powers)
            xiNew += xiTransform
        end
        xiNew /= numberOfSymmetryOperations
        potential += parameters[i+numberOfEquilibriumParameters]*xiNew
    end
    return potential
end

function computePotentialOfPoints(gridPoints::Matrix{Float64}, parameters::Vector{Float64})::Vector{Float64}
    numberOfGridPoints::Int64 = size(gridPoints)[1]
    potentialEnergy::Vector{Float64} = zeros(numberOfGridPoints)
    for i in 1:numberOfGridPoints
        potentialEnergy[i] = potential(gridPoints[i, :], parameters)
    end
    return potentialEnergy
end

function computePotentialOfGrid(grid::DataFrame, parameters::Vector{Float64})::DataFrame
    numberOfGridPoints::Int64 = size(grid)[1]
    potentialEnergy::Vector{Float64} = zeros(numberOfGridPoints)
    for i in 1:numberOfGridPoints
        potentialEnergy[i] = potential(grid[i, 1], parameters)
    end
    grid[:, :potential] = potentialEnergy
    return grid
end

function computeDerivativesAtPoint(coordinates::Vector{Float64}, parameters::Vector{Float64})::Vector{Float64}
    convertToRadians::Float64 = pi/180.00
    degreesOfFreedom::Int64 = size(coordinates)[1]
    equilibriumParametersInFit::Vector{Int64} = equilibriumParametersInFitGlobal
    numberOfEquilibriumParameters::Int64 = size(equilibriumParametersInFit)[1]
    linearParametersInFit::Vector{Int64} = linearParametersInFitGlobal
    numberOfLinearParameters::Int64 = size(linearParametersInFitGlobal)[1]
    symmetryOperations::Array{Float64} = symmetryOperationsGlobal
    numberOfSymmetryOperations::Int64 = size(symmetryOperations)[1]
    equilibriumParametersInFit::Vector{Int64}
    equilibriumBondLength::Float64 = parameters[1]
    morseParameter::Float64 = parameters[2]
    powersMatrix::Matrix{Int64} = powersGlobal
    
    xi::Vector{Float64} = zeros(degreesOfFreedom)
    xi[1:3] = 1.0 .- exp.(-morseParameter.*(coordinates[1:3] .- equilibriumBondLength))
    xi[4:5] = convertToRadians.*coordinates[4:5]
    xi[6] = 1.0 - sin(coordinates[6]*convertToRadians)
    derivatives::Vector{Float64} = zeros(size(parameters)[1])
    for i in 1:numberOfLinearParameters
        powers::Vector{Int64} = powersMatrix[i, :]
        for j in 1:numberOfSymmetryOperations
            transformedXi::Vector{Float64} = symmetryOperations[j, :, :]*xi
            transformedCoordinates::Vector{Float64} = symmetryOperations[j, :, :]*coordinates
            termToAdd::Float64 = 0
            if equilibriumParametersInFit[1] == 1
                termToAdd = (-powers[1]*morseParameter*exp(morseParameter*(transformedCoordinates[1] - equilibriumBondLength))*transformedXi[1]^(powers[1] - 1)*transformedXi[2]^powers[2]*transformedXi[3]^powers[3]*transformedXi[4]^powers[4]*transformedXi[5]^powers[5]*transformedXi[6]^powers[6] +
                -powers[2]*morseParameter*exp(morseParameter*(transformedCoordinates[2] - equilibriumBondLength))*transformedXi[1]^powers[1]*transformedXi[2]^(powers[2] - 1)*transformedXi[3]^powers[3]*transformedXi[4]^powers[4]*transformedXi[5]^powers[5]*transformedXi[6]^powers[6] +
                -powers[3]*morseParameter*exp(morseParameter*(transformedCoordinates[3] - equilibriumBondLength))*transformedXi[1]^powers[1]*transformedXi[2]^powers[2]*transformedXi[3]^(powers[3]-1)*transformedXi[4]^powers[4]*transformedXi[5]^powers[5]*transformedXi[6]^powers[6])*parameters[i+numberOfEquilibriumParameters]
                if isnan(termToAdd) == false
                    derivatives[1] += termToAdd
                end
            end
            if equilibriumParametersInFit[2] == 1
                termToAdd = (powers[1]*(transformedCoordinates[1] - equilibriumBondLength)*exp(morseParameter*(transformedCoordinates[1] - equilibriumBondLength))*transformedXi[1]^(powers[1] - 1)*transformedXi[2]^powers[2]*transformedXi[3]^powers[3]*transformedXi[4]^powers[4]*transformedXi[5]^powers[5]*transformedXi[6]^powers[6] +
                powers[2]*(transformedCoordinates[2] - equilibriumBondLength)*exp(morseParameter*(transformedCoordinates[2] - equilibriumBondLength))*transformedXi[1]^powers[1]*transformedXi[2]^(powers[2] - 1)*transformedXi[3]^powers[3]*transformedXi[4]^powers[4]*transformedXi[5]^powers[5]*transformedXi[6]^powers[6] +
                powers[3]*(transformedCoordinates[3] - equilibriumBondLength)*exp(morseParameter*(transformedCoordinates[3] - equilibriumBondLength))*transformedXi[1]^powers[1]*transformedXi[2]^powers[2]*transformedXi[3]^(powers[3]-1)*transformedXi[4]^powers[4]*transformedXi[5]^powers[5]*transformedXi[6]^powers[6])*parameters[i+numberOfEquilibriumParameters]
                if isnan(termToAdd) == false
                    derivatives[2] += termToAdd
                end
            end
            if linearParametersInFit[i] == 1
                derivatives[i+numberOfEquilibriumParameters] += prod(transformedXi.^powers)
            end
        end
    end
    derivatives = derivatives ./ numberOfSymmetryOperations
    # println(derivatives)
    return derivatives
end

function computeDerivativesOfPoints(gridPoints::Matrix{Float64}, parameters::Vector{Float64})::Matrix{Float64}
    numberOfGridPoints::Int64 = size(gridPoints)[1]
    numberOfParameters::Int64 = size(parameters)[1]
    jacobian::Matrix{Float64} = zeros(numberOfGridPoints, numberOfParameters)
    for i in 1:numberOfGridPoints
        jacobian[i, :] = computeDerivativesAtPoint(gridPoints[i, :], parameters)
    end
    return jacobian
end