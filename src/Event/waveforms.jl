const WaveformSamples{T<:RealQuantity} = AbstractVector{T}
const TimeAxis{T<:RealQuantity} = AbstractVector{T}

export RDWaveform

"""
    RDWaveform
Represents a radiation detector signal waveform.
Fields:
* `time`: time axis, typically a range
* `value`: sample values
"""
struct RDWaveform{
    T<:RealQuantity,U<:RealQuantity,
    TV<:TimeAxis{T},UV<:WaveformSamples{U},
}
    time::TV
    value::UV
end


RDWaveform{T,U,TV,UV}(wf::RDWaveform) where {T,U,TV,UV} = RDWaveform{T,U,TV,UV}(wf.time, wf.value)
Base.convert(::Type{RDWaveform{T,U,TV,UV}}, wf::RDWaveform) where {T,U,TV,UV} = RDWaveform{T,U,TV,UV}(wf)



function _unitful_axis_label(axis_label::AbstractString, units::Unitful.Unitlike)
    u_str = units == NoUnits ? "" : "$units"
    u_str2 = replace(u_str, "μ" => "u")
    isempty(u_str2) ? axis_label : "$axis_label [$u_str2]"
end

function prep_for_plotting(X::AbstractVector{<:RealQuantity}, axis_label::AbstractString)
    X_units = unit(eltype(X))
    X_unitless = X_units == NoUnits ? X : ustrip.(X)
    Array(X_unitless), _unitful_axis_label(axis_label, X_units), X_units
end

@recipe function f(wf::RDWaveform)
    X, X_label = prep_for_plotting(wf.time, "t")
    Y, Y_label = prep_for_plotting(wf.value, "sample value")

    seriestype = get(plotattributes, :seriestype, :line)

    if seriestype in (:line, :scatter)
        @series begin
            seriestype := seriestype
            xlabel --> X_label
            ylabel --> Y_label
            (X, Y)
        end
    else
        error("seriestype $seriestype not supported")
    end
end