
function _CheckContiguousTypes(typeDict::Dict{Int64, Element})
    n = length(typeDict)
    for t in 1:n
        haskey(typeDict, t) || error("typeDict keys must be 1:$(n), missing $(t).")
    end
    return n
end

function InitConstantsByType(typeDict::Dict{Int64, Element}, parameters::Parameters)
    n = _CheckContiguousTypes(typeDict)
    V_upterm = zeros(Float64, n, n)
    a_U = zeros(Float64, n, n)
    E_m = zeros(Float64, n)
    S_e_upTerm = zeros(Float64, n, n)
    S_e_downTerm = zeros(Float64, n, n)
    x_nl = zeros(Float64, n, n)
    a = zeros(Float64, n, n)
    Q_nl = zeros(Float64, n, n)
    Q_loc = zeros(Float64, n, n)
    qMax = zeros(Float64, n, n)
    sigma = zeros(Float64, n)
    log_info("")
    log_info("Vibration σ for each type:")
    for p in 1:n
        radius_p, mass_p, Z_p, _, _, α_p, β_p = TypeToProperties(p, typeDict)
        for t in 1:n
            radius_t, _, Z_t, _, _, _, _ = TypeToProperties(t, typeDict)
            V_upterm[p, t] = BCA.ConstantFunctions.V_upterm(Z_p, Z_t)
            a_U[p, t] = BCA.ConstantFunctions.a_U(Z_p, Z_t)
            S_e_upTerm[p, t] = BCA.ConstantFunctions.S_e_upTerm(p, Z_p, Z_t, mass_p, α_p)
            x_nl[p, t] = BCA.ConstantFunctions.x_nl(p, Z_p, Z_t, β_p)
            a[p, t] = BCA.ConstantFunctions.a(Z_p, Z_t)
            Q_nl[p, t] = BCA.ConstantFunctions.Q_nl(Z_p, Z_t, parameters.pMax)
            Q_loc[p, t] = BCA.ConstantFunctions.Q_loc(Z_p, Z_t)
            qMax[p, t] = radius_p + radius_t
        end
        E_m[p] = BCA.ConstantFunctions.E_m(Z_p, mass_p)
        sigma[p] = TemperatureToSigma(parameters.temperature, parameters.DebyeTemperature, mass_p)
        log_info("  Type $(p): σ = $(round(sigma[p]; digits=3)) Å")
    end
    return ConstantsByType(PairTable(V_upterm), PairTable(a_U), TypeTable(E_m),
                           PairTable(S_e_upTerm), PairTable(S_e_downTerm), PairTable(x_nl),
                           PairTable(a), PairTable(Q_nl), PairTable(Q_loc), PairTable(qMax),
                           TypeTable(sigma))
end


function InitθτFunctions(parameters::Parameters, constantsByType::ConstantsByType)
    typeDict = parameters.typeDict
    n = _CheckContiguousTypes(typeDict)
    θTable = Matrix{ΘτInterpolation}(undef, n, n)
    τTable = Matrix{ΘτInterpolation}(undef, n, n)
    log_separator()
    log_info("Loading θ and τ functions...")
    for type_p in 1:n
        for type_t in 1:n
            mass_p = typeDict[type_p].mass
            mass_t = typeDict[type_t].mass
            θInterpolation, τInterpolation = θτFunctions(mass_p, mass_t, type_p, type_t, constantsByType, parameters)
            θTable[type_p, type_t] = θInterpolation
            τTable[type_p, type_t] = τInterpolation
            log_debug("  $(parameters.typeDict[type_p].name) → $(parameters.typeDict[type_t].name) loaded")
        end
    end
    log_success("All θ and τ functions initialized")
    log_separator()
    return InterpTable(θTable), InterpTable(τTable)
end


function θτFunctions(mass_p::Float64, mass_t::Float64, type_p::Int64, type_t::Int64, constantsByType::ConstantsByType, parameters::Parameters)
    E_p_axis = Float64[]
    p_axis = Float64[]
    θMatrix = Matrix{Float64}(undef, 0, 0)
    τMatrix = Matrix{Float64}(undef, 0, 0)
    try
        E_p_axis, p_axis, θMatrix, τMatrix = LoadθτData(type_p, type_t, parameters)
    catch
        EPowerRange = parameters.EPowerRange    
        pPowerRange = parameters.pPowerRange
        nE = length(EPowerRange)
        np = length(pPowerRange)
        θMatrix = Array{Float64, 2}(undef, nE, np)
        τMatrix = Array{Float64, 2}(undef, nE, np)
        N = length(EPowerRange)
        @showprogress @threads for i in 1:N
            E_p_power = EPowerRange[i]
            E_p = 10.0^E_p_power
            for (j, p_power) in enumerate(pPowerRange)
                p = 10.0^p_power
                θ, τ = BCA.θτ(E_p, mass_p, mass_t, type_p, type_t, p, constantsByType)
                θMatrix[i, j] = θ
                τMatrix[i, j] = τ
            end
        end
        E_p_axis = collect(EPowerRange)
        p_axis = collect(pPowerRange)    
        SaveθτData(type_p, type_t, θMatrix, τMatrix, E_p_axis, p_axis, parameters)
    end
    # interpolate
    θFunction = interpolate((E_p_axis, p_axis), θMatrix, Gridded(Linear()))
    τFunction = interpolate((E_p_axis, p_axis), τMatrix, Gridded(Linear()))
    return θFunction, τFunction
end
