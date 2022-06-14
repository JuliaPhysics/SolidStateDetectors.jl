@recipe function f(wvs::CustomIDArray{RadiationDetectorSignals.RDWaveform})
    @series begin
        wvs.data
    end
end