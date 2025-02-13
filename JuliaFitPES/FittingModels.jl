module FittingModels
    using Printf, DataFrames, LsqFit
    export computeWeightOfPoint, printComputedEnergies, LsqCurveFit

    function computeWeightOfPoint(energy::Float64, energyThreshold=10000.0::Float64)::Float64
        weight::Float64 = (tanh(−0.0006*(energy - energyThreshold)) + 1.002002002)/2.002002002
        if energy > 10000.0
            weight = weight/(0.0001*energy)
        else
            weight = weight/(0.0001*10000.0)
        end
        return weight
    end

    function printComputedEnergies(grid::DataFrame, minimumEnergy::Float64)
        numberOfGridPoints::Int64 = size(grid)[1]
        degreesOfFreedom::Int64 = size(grid[1, 1])[1]
        for i in 1:numberOfGridPoints
            coordinates::Vector{Float64} = grid[i, 1]
            for j in 1:degreesOfFreedom
                @printf("%12.8f   ", coordinates[j])
            end
            @printf("%12.8f   ", grid[i, :energy])
            @printf("%12.8f   ", grid[i, :potential])
            @printf("%12.8f\n", grid[i, :error])
        end
    end
    function LsqCurveFit(potentialEnergyModel::Function, derivatives::Function, gridOfPoints::Matrix{Float64}, parameters::Vector{Float64}, energies::Vector{Float64}, weights::Vector{Float64})::Vector{Float64}
        @time potentialFit = curve_fit((gridOfPoints, parameters) -> potentialEnergyModel(gridOfPoints, parameters),
        (gridOfPoints, parameters) -> derivatives(gridOfPoints, parameters),
        gridOfPoints, energies, weights, parameters)
        return potentialFit.param
    end
end