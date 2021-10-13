# This should be in RadiationDetectorSignals.jl
@recipe function f(wvs::Vector{Union{Missing, RadiationDetectorSignals.RDWaveform}})
    @series begin
        RadiationDetectorSignals.RDWaveform[wv for wv in skipmissing(wvs)]
    end
end