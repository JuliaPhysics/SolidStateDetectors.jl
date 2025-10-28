"""
    plot_electron_and_hole_contribution(evt::Event{T}, sim::Simulation{T}, contact_id::Int)
    
Plots the waveform as well as the electron and hole contribution to the waveform 
of a [`Contact`](@ref) with a given `contact_id` of an [`Event`](@ref).

## Arguments
* `evt::Event{T}`: [`Event`](@ref) in which the charges have already been drifted.
* `sim::Simulation{T}`: [`Simulation`](@ref) which defines the setup in which the charges in `evt` were drifted.
* `contact_id::Int`: The `id` of the [`Contact`](@ref) for which the waveform should be split into electron and hole contribution.

## Keywords 
* `n_samples::Int`: Number of samples with which the waveforms will be plotted. The default is the number of samples of the original waveform.

## Example
```julia 
using Plots
using SolidStateDetector
T = Float32

simulation = Simulation{T}(SSD_examples[:InvertedCoax])
simulate!(simulation)
event = Event([CartesianPoint{T}(0.02,0.01,0.05)])
simulate!(event, simulation)

contact_id = 1
plot_electron_and_hole_contribution(evt, sim, contact_id, n_samples = 300)
```

!!! note 
    The drift paths in `evt` need to be calculated using [`drift_charges!`](@ref) before calling this function.
    
!!! note 
    This method requires to load the package [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
    
See also [`get_electron_and_hole_contribution`](@ref).
"""
function plot_electron_and_hole_contribution end

@userplot Plot_electron_and_hole_contribution
@recipe function f(gdd::Plot_electron_and_hole_contribution; linewidth = 2, n_samples = missing)
    
    sim::Simulation = gdd.args[2]
    T = get_precision_type(sim)
    evt::Event{T} = gdd.args[1]
    contact_id::Int = gdd.args[3]
    
    ismissing(n_samples) ? n_samples = length(evt.waveforms[contact_id].signal) : nothing

    # check in which units the waveforms are calibrated
    signal_unit::Unitful.Units = unit(evt.waveforms[contact_id].signal[end])
    wf::NamedTuple{(:electron_contribution, :hole_contribution), <:Tuple{RDWaveform, RDWaveform}} = get_electron_and_hole_contribution(evt, sim, contact_id; signal_unit)
    
    unitformat --> :slash
    yguide --> signal_unit isa Unitful.Units{<:Any, Unitful.Charge} ? "Charge" : "Energy"

    @series begin
        linecolor := :red 
        linewidth --> linewidth
        label --> "Electron contribution"
        add_baseline_and_extend_tail(wf.electron_contribution, 0, n_samples)
    end
    
    @series begin
        linecolor := :green 
        linewidth --> linewidth
        label --> "Hole contribution"
        add_baseline_and_extend_tail(wf.hole_contribution, 0, n_samples)
    end
    
    @series begin
        linecolor --> contact_id
        linewidth --> linewidth
        seriesalpha --> 0.5
        label --> "Sum"
        add_baseline_and_extend_tail(evt.waveforms[contact_id], 0, n_samples)
    end
end


# This should be in RadiationDetectorSignals.jl
@recipe function f(wvs::Vector{Union{Missing, RadiationDetectorSignals.RDWaveform}})
    @series begin
        RadiationDetectorSignals.RDWaveform[wv for wv in skipmissing(wvs)]
    end
end