@recipe function f(wvs::CustomIDVector{<:RadiationDetectorSignals.RDWaveform})
    @series begin
        label --> wvs.idx'
        linecolor --> wvs.idx'
        unitformat --> :slash
        [wv for wv in wvs.data]
    end
end