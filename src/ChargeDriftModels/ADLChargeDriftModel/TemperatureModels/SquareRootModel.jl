@doc raw"""

mutable struct SquareRootModel{T <: SSDFloat}

Values needed to parametrize the longitudinal drift velocity temperature dependence as a function of the electric field strength.

## Background information 

The parameterization for the longitudinal drift velocity temperature dependence, ``v_l``, as a function of the electric 
field strength, ``E``, was proposed by [M.Ali Omar, Lino Reggiani](https://www.sciencedirect.com/science/article/pii/0038110187900633):
```math
v_l(T) = v_l(T_0)f(T)/f(T_0) 
```
where v_l(T_0) is a known velocity at reference temperature ``T_0``. The temperature dependence, ``f(T)``,  is parametrized as:
```math
f(T) =  V_s(T)\frac{E/E_c(T)}{(1 + (E/E_c(T))^2)^{1/2}}
```
where ``E_c(T) = V_s(T)/\mu_0(T)`` and  ``V_s(T) = B\tanh^{1/2}(\theta/2T)`` and  ``\mu_0(T) = A/T^P``.
The four parameters, ``A`` (``\mu_0`` at ``T=1``K), ``P``, ``B``, ``\theta``, are different
for electrons and holes. The work was confined to the fields along the <100> axis. Here the same values are assumed for the <111> axis. 
    
## Fields
* `temperature::T`: Actual temperature of the crystal, ``T`` in the parameterization shown above.
* `reftemperature::T`: Temperature at which the reference velocities to be scaled were measured. ``T_0`` in the parameterization shown above.
* `A_e::T`: Parameter ``A`` for electrons in the parameterization shown above.
* `P_e::T`: Parameter ``P`` for electrons in the parameterization shown above.
* `B_e::T`: Parameter ``B`` for electrons in the parameterization shown above.
* `theta_e::T`: Parameter ``\theta`` for electrons in the parameterization shown above.
* `A_h::T`: Parameter ``A`` for holes in the parameterization shown above.
* `P_h::T`: Parameter ``P`` for holes in the parameterization shown above.
* `B_h::T`: Parameter ``B`` for holes in the parameterization shown above.
* `theta_h::T`: Parameter ``\theta`` for holes in the parameterization shown above.

"""

mutable struct SquareRootModel{T <: SSDFloat} <: AbstractTemperatureModel{T}

    #temperatures needed for the drift velocity rescaling
    temperature::T      # actual temperature of the crystal
    reftemperature::T   # temperature at which the reference values were measured

    #fit parameters for electrons (e) and holes (h) for <100> and <111> axes
    A_e::T
    P_e::T
    B_e::T
    theta_e::T

    A_h::T
    P_h::T
    B_h::T
    theta_h::T
end

function SquareRootModel(config_file::Dict; T::Type{<:AbstractFloat} = Float32)::SquareRootModel
    return SquareRootModel{T}(config_file::Dict)
end

function SquareRootModel{T}(config_file::Dict; temperature::Union{Missing, T} = missing)::SquareRootModel where T <: SSDFloat
    config_file["temperature_dependence"]["model"] != "SquareRoot" ? error() : nothing
    ismissing(temperature) ? temperature = config_file["temperature_dependence"]["temperature"] : nothing
    m = SquareRootModel{T}(
        temperature, #temperature = config_file["temperature_dependence"]["temperature"]
        config_file["drift"]["velocity"]["temperature"],

        config_file["temperature_dependence"]["parameters"]["e100"]["A"],
        config_file["temperature_dependence"]["parameters"]["e100"]["P"],
        config_file["temperature_dependence"]["parameters"]["e100"]["B"],
        config_file["temperature_dependence"]["parameters"]["e100"]["theta"],

        config_file["temperature_dependence"]["parameters"]["h100"]["A"],
        config_file["temperature_dependence"]["parameters"]["h100"]["P"],
        config_file["temperature_dependence"]["parameters"]["h100"]["B"],
        config_file["temperature_dependence"]["parameters"]["h100"]["theta"]
    )
end

function scale_to_given_temperature(Emag::T, m::SquareRootModel{T})::NTuple{4,T} where T <: SSDFloat

    function f(A::T, P::T, B::T, theta::T)::T where T <: SSDFloat
        mu_0_1 = A * m.reftemperature^P
        v_s_1 = B * sqrt(tanh(0.5 * theta / m.reftemperature))
        E_c_1 = v_s_1 / mu_0_1
        mu_0_2 = A * m.temperature^P
        v_s_2 = B * sqrt(tanh(0.5 * theta / m.temperature))
        E_c_2 = v_s_2 / mu_0_2;
        return (v_s_2 * (Emag/E_c_2) / sqrt(1.0 + (Emag/E_c_2)^2)) / (v_s_1 * (Emag/E_c_1) / sqrt(1.0 + (Emag/E_c_1)^2))
    end

    scale_e100 = f(m.A_e, m.P_e, m.B_e, m.theta_e)
    scale_e111 = f(m.A_e, m.P_e, m.B_e, m.theta_e)
    scale_h100 = f(m.A_h, m.P_h, m.B_h, m.theta_h)
    scale_h111 = f(m.A_h, m.P_h, m.B_h, m.theta_h)

    return scale_e100, scale_e111, scale_h100, scale_h111
end

print(io::IO, tm::SquareRootModel{T}) where {T <: SSDFloat} = print(io, "SquareRootModel{$T}")
function println(io::IO, tm::SquareRootModel{T}) where {T <: SSDFloat}
    println("\n________SquareRootModel________")
    println("V_s(T) = Btanh^{1/2}(theta/2T) and mu_0(T) = A/T^P\n")
    println("---Temperature settings---")
    println("Crystal temperature:  \t$(tm.temperature)")
    println("Reference temperature:\t$(tm.reftemperature)\n")

    println("---Fitting parameters---")
    println("   \te100      \th100")
    println("A \t$(tm.A_e)   \t$(tm.A_h)")
    println("P \t$(tm.P_e)   \t$(tm.P_h)")
    println("B \t$(tm.B_e)   \t$(tm.B_h)")
    println("theta \t$(tm.theta_e)   \t$(tm.theta_h)")
end
