module ReadInput
    using ..PotentialEnergyModel, DataFrames
    export readBlocks, readModel, readGrid 

    function readBlocks(inputFile::Vector{String})::Vector{Vector{String}}
        println("model")
        keywords::Vector{String} = ["model", "grid"]
        addToBlock::Bool = false
        inputBlocks::Vector{Vector{String}} = []
        for keyword in keywords
            newBlock::Vector{String} = []
            for line in inputFile
                if lowercase(line) == lowercase(keyword)
                    addToBlock = true
                elseif lowercase(line) == "end"
                    addToBlock = false
                end
                if addToBlock
                    push!(newBlock, line)
                end
            end
            push!(inputBlocks, newBlock)
        end
        return inputBlocks
    end

    function readModel(modelBlock::Vector{String})::potentialEnergyModel
        moleculeLine::String = modelBlock[2]
        println(moleculeLine)
        molecule::String = split(moleculeLine, r"\s+")[2]
        modesLine::String = modelBlock[3]
        println(modesLine)
        numberOfModes::Int64 = parse(Int64, split(modesLine, r"\s+")[2])

        equilibriumParametersLine::String = modelBlock[4]
        println(equilibriumParametersLine)
        numberOfEquilibriumParameters::Int64 = parse(Int64, split(equilibriumParametersLine, r"\s+")[3])
        equilibriumParametersInFit::Vector{Int64} = zeros(numberOfEquilibriumParameters)
        equilibriumParameters::Vector{Float64} = zeros(numberOfEquilibriumParameters)        
        for i in 5:4+numberOfEquilibriumParameters
            println(modelBlock[i])
            equilibriumParameterSplit::Array{SubString{String}} = split(modelBlock[i], r"\s+")
            equilibriumParametersInFit[i-4] = parse(Int64, equilibriumParameterSplit[2])
            equilibriumParameters[i-4] = parse(Float64, equilibriumParameterSplit[3])
        end

        numberOfLinearParameters::Int64 = size(modelBlock)[1] - (4 + numberOfEquilibriumParameters)
        fourierType::Vector{String} = []
        linearParametersInFit::Vector{Int64} = zeros(numberOfLinearParameters)
        linearParameters::Vector{Float64} = zeros(numberOfLinearParameters)
        powers::Matrix{Int64} = zeros(numberOfLinearParameters, numberOfModes)
        for i in 5+numberOfEquilibriumParameters:4+numberOfEquilibriumParameters+numberOfLinearParameters
            println(modelBlock[i])
            parameterSplit::Array{SubString{String}} = split(modelBlock[i], r"\s+")
            push!(fourierType, parameterSplit[1])
            powers[i-(4+numberOfEquilibriumParameters), :] = parse.(Int64, parameterSplit[2:1+numberOfModes])
            linearParametersInFit[i-(4+numberOfEquilibriumParameters)] = parse(Int64, parameterSplit[end - 1])
            linearParameters[i-(4+numberOfEquilibriumParameters)] = parse(Float64, parameterSplit[end])
        end

        println("end")
        model::potentialEnergyModel = potentialEnergyModel(molecule, numberOfModes, equilibriumParametersInFit, 
            equilibriumParameters, fourierType, linearParametersInFit, linearParameters, powers)
        return model
    end

    function readGrid(gridBlock::Vector{String}, numberOfModes::Int64)::DataFrame
        println(gridBlock[1])
        numberOfGridPoints::Int64 = size(gridBlock)[1] - 1
        gridPoint::Vector{Vector{Float64}} = []
        energies::Vector{Float64} = zeros(numberOfGridPoints)
        pointNumber::Vector{Int64} = zeros(numberOfGridPoints)
        for i in 2:numberOfGridPoints+1
            println(gridBlock[i])
            gridBlockSplit = split(gridBlock[i], r"\s+")
            push!(gridPoint, parse.(Float64, gridBlockSplit[1:numberOfModes]))
            energies[i - 1] = parse(Float64, gridBlockSplit[end-1])
            pointNumber[i - 1] = parse(Int64, gridBlockSplit[end])
        end
        println("end")
        grid::DataFrame = DataFrame(coordinates=gridPoint, energy=energies, point=pointNumber)
        return grid
    end
end