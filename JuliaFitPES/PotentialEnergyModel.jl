module PotentialEnergyModel
    using Printf
    export potentialEnergyModel, updatePotentialEnergyModel, printModel

    mutable struct potentialEnergyModel
        molecule::String
        numberOfModes::Int64
        equilibriumParametersInFit::Vector{Int64}
        equilibriumParameters::Vector{Float64}
        fourierType::Vector{String}
        linearParametersInFit::Vector{Int64}
        linearParameters::Vector{Float64}
        powers::Matrix{Int64}

        function potentialEnergyModel(molecule::String, numberOfModes::Int64,
            equilibriumParametersInFit::Vector{Int64}, equilibriumParameters::Vector{Float64},
            fourierType::Vector{String}, linearParametersInFit::Vector{Int64}, linearParameters::Vector{Float64}, 
            powers::Matrix{Int64})
            new(molecule::String, numberOfModes::Int64,
            equilibriumParametersInFit::Vector{Int64}, equilibriumParameters::Vector{Float64},
            fourierType::Vector{String}, linearParametersInFit::Vector{Int64}, linearParameters::Vector{Float64}, 
            powers::Matrix{Int64})
        end
    end
    function updatePotentialEnergyModel(model::potentialEnergyModel, newParameters::Vector{Float64})::potentialEnergyModel
        numberOfEquilibriumParameters::Int64 = size(model.equilibriumParameters)[1]
        newEquilibriumParameters::Vector{Float64} = newParameters[1:numberOfEquilibriumParameters]
        newLinearParameters::Vector{Float64} = newParameters[numberOfEquilibriumParameters+1:end]
        newModel::potentialEnergyModel = potentialEnergyModel(model.molecule, model.numberOfModes,
        model.equilibriumParametersInFit, newEquilibriumParameters,
        model.fourierType, model.linearParametersInFit, newLinearParameters, 
        model.powers)
        return newModel
    end
    function printModel(model::potentialEnergyModel)
        numberOfEquilibriumParameters::Int64 = size(model.equilibriumParameters)[1]
        numberOfLinearParameters::Int64 = size(model.linearParameters)[1]
        for i in 1:numberOfEquilibriumParameters
            @printf("              %1.0f   %12.8f  \n", model.equilibriumParametersInFit[i], model.equilibriumParameters[i])
        end
        for i in 1:numberOfLinearParameters
            @printf("%s ", model.fourierType[i])
            for j in 1:model.numberOfModes
                @printf(" %1.0f ", model.powers[i, j])
            end
            @printf("  %1.0f ", model.linearParametersInFit[i])
            @printf("         %12.8f \n", model.linearParameters[i])
        end
    end
end