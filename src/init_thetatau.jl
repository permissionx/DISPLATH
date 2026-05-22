
function InitConstantsByType(typeDict::Dict{Int64, Element}, parameters::Parameters)
    V_upterm = Dict{Vector{Int64}, Float64}()
    a_U = Dict{Vector{Int64}, Float64}()
    E_m = Dict{Int64, Float64}()
    S_e_upTerm = Dict{Vector{Int64}, Float64}()
    S_e_downTerm = Dict{Vector{Int64}, Float64}()
    x_nl = Dict{Vector{Int64}, Float64}()
    a = Dict{Vector{Int64}, Float64}()
    Q_nl = Dict{Vector{Int64}, Float64}()  
    Q_loc = Dict{Vector{Int64}, Float64}()
    types = keys(typeDict)
    qMax = Dict{Vector{Int64}, Float64}()
    sigma = Dict{Int64, Float64}()
    log_info("")
    log_info("Vibration σ for each type:")
    for p in types
        radius_p, mass_p, Z_p, _, _, α_p, β_p = TypeToProperties(p, typeDict)
        for t in types
            radius_t, _, Z_t, _, _, _, _ = TypeToProperties(t, typeDict)
            V_upterm[[p,t]] = BCA.ConstantFunctions.V_upterm(Z_p, Z_t)
            a_U[[p,t]] = BCA.ConstantFunctions.a_U(Z_p, Z_t)
            S_e_upTerm[[p,t]] = BCA.ConstantFunctions.S_e_upTerm(p, Z_p, Z_t, mass_p, α_p)
            x_nl[[p,t]] = BCA.ConstantFunctions.x_nl(p, Z_p, Z_t, β_p)
            a[[p,t]] = BCA.ConstantFunctions.a(Z_p, Z_t)
            Q_nl[[p,t]] = BCA.ConstantFunctions.Q_nl(Z_p, Z_t, parameters.pMax)
            Q_loc[[p,t]] = BCA.ConstantFunctions.Q_loc(Z_p, Z_t)
            qMax[[p,t]] = radius_p + radius_t
        end
        E_m[p] = BCA.ConstantFunctions.E_m(Z_p, mass_p)
        sigma[p] = TemperatureToSigma(parameters.temperature, parameters.DebyeTemperature, mass_p)
        log_info("  Type $(p): σ = $(round(sigma[p]; digits=3)) Å")
    end
    return ConstantsByType(V_upterm, a_U, E_m, S_e_upTerm, S_e_downTerm, x_nl, a, Q_nl, Q_loc, qMax, sigma)
end


function InitθτFunctions(parameters::Parameters, constantsByType::ConstantsByType)
    typeDict = parameters.typeDict
    θFunctions = Dict{Vector{Int64}, Function}()
    τFunctions = Dict{Vector{Int64}, Function}()
    log_separator()
    log_info("Loading θ and τ functions...")
    for type_p in keys(typeDict)
        for type_t in keys(typeDict)
            mass_p = typeDict[type_p].mass
            mass_t = typeDict[type_t].mass
            θInterpolation, τInterpolation = θτFunctions(mass_p, mass_t, type_p, type_t, constantsByType, parameters)
            θFunctions[[type_p, type_t]] = (E_p, p) -> θInterpolation(E_p, p)
            τFunctions[[type_p, type_t]] = (E_p, p) -> τInterpolation(E_p, p)
            log_debug("  $(parameters.typeDict[type_p].name) → $(parameters.typeDict[type_t].name) loaded")
        end
    end
    log_success("All θ and τ functions initialized")
    log_separator()
    return θFunctions, τFunctions
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
