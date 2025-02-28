---
title: "End-to-end simulations of semiconductor detectors using the `Geant4.jl` extension of `SolidStateDetectors.jl`"
tags:
  - Julia
  - physics
  - solid state detectors
  - electric field calculation
  - pulse shape simulation
  - charge trapping
  - Geant4
authors:
  - name: Felix Hagemann
    orcid: 0000-0001-5021-3328
    corresponding: true
    # equal-contrib: true
    affiliation: "1, 2"
  - name: Julian Henzler
    orcid: 0009-0000-6977-8687
    # equal-contrib: true
    affiliation: 1
  - name: Benedikt Nagler
    orcid: 0009-0004-0906-4286
    # equal-contrib: true
    affiliation: 1
  - name: Oliver Schulz
    orcid: 0000-0002-4200-5905
    # equal-contrib: true
    affiliation: 1
affiliations:
 - name: Max Planck Institut für Physik, Boltzmannstr. 8, 85748 Garching, Germany
   index: 1
   ror: "0079jjr10"
 - name: Space Sciences Laboratory, University of California, Berkeley, 7 Gauss Way, Berkeley CA 94720, USA
   index: 2
   ror: "048400679"
date: 30 May 2026
bibliography: paper.bib
---

# Summary

`SolidStateDetectors.jl` is an open-source Julia package for calculating electric fields, simulating charge transport and generating waveforms in semiconductor detectors.
Since its initial publication in 2021, it has been significantly extended to enhance the accuracy of the charge transport simulations, including the simulation of charge cloud effects, charge trapping effects and reduced charge collection near contacts.

This paper presents the new `Geant4.jl` extension that incorporates Monte Carlo event generation into the simulation pipeline.
Based on the Julia wrapper of Geant4, the extension generates realistic input event distributions using the same `SolidStateDetectors.jl` detector object used for both field calculation and charge transport simulations.
This approach eliminates the need to manually maintain separate descriptions of the detector setup for event and signal generation.
By reducing the complexity of setting up Geant4 simulations and eliminating potential inconsistencies between geometry implementations in different tools, this extension makes end-to-end simulations from an initial particle source to detector waveforms more accessible.


# Statement of Need

Solid state detectors are widely used in fundamental physics research, primarily due to their excellent energy resolution and high detection efficiency.
These detectors operate through the direct conversion of radiation energy into a measurable electronic signal: 
a particle depositing energy in the semiconductor creates electron-hole pairs, which drift in an applied electric field and induce charges on the contacts that are read out as the detector signal.
A precise understanding of this detector response is critical for modern physics experiments. 
Simulated waveforms are used to develop and test pulse shape analysis routines, to assist with detector calibrations, to explain features observed in experimental data and, recently, to train machine learning models for event classification.

`SolidStateDetectors.jl` is an open-source Julia software package to calculate complex 3D electric fields and simulate charge transport in semiconductor detectors.
Building on its initial publication [@Abt:2021], it has been expanded to include advanced charge transport effects such as charge cloud diffusion and self-repulsion, charge trapping [@Boggs:2023] and reduced charge collection in lithium-drifted contact layers [@Zhang:2026]. 

\newpage

Previously, simulating the detector response to a particle source using `SolidStateDetectors.jl` still required an additional component: 
the ability to generate realistic event distributions in the detector from that given source.
This task was traditionally handled by the Geant4 toolkit [@Geant4:2003], a widely used Monte Carlo framework written in C++ to simulate particle-matter interactions.
However, Geant4 requires specifying the detector simulation setup either directly in code or using a GDML configuration file,
while `SolidStateDetectors.jl` uses its own YAML configuration file format which carries additional information required for charge transport and waveform simulations.
Users therefore had to switch between programming languages, maintain two independent representations of the same detector setup and keep them consistent by hand.
This approach was not only inefficient but also prone to errors.

The `Geant4.jl` extension presented in this paper solves this problem: 
it automatically generates a GDML description of the detector geometry from the `Simulation` object at runtime and passes it to `Geant4.jl`, the Julia wrapper of Geant4.
In this way, Monte Carlo event generation and signal simulation run entirely from Julia, driven by the same YAML configuration file.
The software presented here provides physicists working on solid state detectors with a powerful tool to simulate and characterize the detector response to a particle source.



# State of the Field

In the current landscape for semiconductor detector simulations, users typically have to compromise between geometric flexibility and Geant4 integration.
For macroscopic high-purity germanium detectors, common choices include highly optimized C and C++ tools like `fieldgen`/`siggen` [@Siggen:2021] and the AGATA detector library [@ADL:2016].
While their charge transport algorithms are exceptionally fast, their native internal field solvers are built heavily around specific detector geometries.
To model arbitrary 3D geometries using the AGATA detector library, users must calculate the electric fields and weighting potentials in external software like SIMION and import them as cubic grids.
Even when paired with Geant4 frameworks, the particle tracking and charge transport simulations are separate steps using independent geometry representations.
For silicon and diamond detectors, the ROOT-based Weightfield2 [@Weightfield2:2015] handles rapid signal simulation exceptionally well, but is restricted entirely to 2D cross-sections.

Full 3D end-to-end simulations require coupling charge transport simulations with Monte Carlo particle trackers like Geant4.
Tools like the ROOT-based KDetSim [@KDetSim:2016] and the Python framework RASER [@RASER:2026] take this approach.
Within this category, Allpix^2^ [@Allpix2:2018] is a prominent and widely adopted framework that features tight C++ Geant4 integration and is specifically designed for microscopic silicon pixel trackers.
While Allpix^2^ can generate simple analytic field profiles, it lacks a native solver for complex 3D electric fields and relies on imported field maps from external, often commercial, TCAD software.

`SolidStateDetectors.jl` addresses these limitations by providing a native environment for both field calculation and Monte Carlo integration.
Written entirely in Julia, this open-source package calculates complex 3D electric fields using constructive solid geometry, allowing for arbitrary 3D geometries and removing the dependency on external field solvers. 
The `Geant4.jl` extension unifies the simulation pipeline: a single YAML configuration file drives the field calculation, charge transport, and the automatic generation of the GDML configuration file for Geant4.
This integrated approach ensures that all simulation stages rely on the same detector representation, preventing geometric mismatches between event and signal generation.

\newpage 

# Software Design

The `Geant4.jl` functionality in `SolidStateDetectors.jl` is implemented using Julia's package extension mechanism.
A package extension is a code module inside a host package that is only loaded when a specific trigger package is also present in the Julia session.
The trigger package is declared as a *weak* dependency, meaning that it does not need to be installed by users who do not need the extended functionality.

Geant4 is a large installation with its own dependencies.
Adding the Julia wrapper `Geant4.jl` as a weak dependency of `SolidStateDetectors.jl` allows users to calculate electric fields, simulate charge transport and generate detector signals without forcing them to install Geant4 if realistic input event distributions are not required.
Julia's multiple dispatch allows the extension to attach Geant4-specific methods to existing `SolidStateDetectors.jl` types, keeping a clean separation between the code in the extension and the base package.

Geometries in `SolidStateDetectors.jl` are defined via a YAML configuration file, see \autoref{appendix:config_file} for an example.
All volumes are built using *Constructive Solid Geometry* (CSG) from primitive volumes by combining them through boolean operations into a binary geometry tree [@CSG:1989].

Geometries in Geant4 are also built using CSG and can be defined via a GDML configuration file.
The `Geant4.jl` extension of `SolidStateDetectors.jl` generates the GDML file at runtime from the same `Simulation` object that is used for electric field, charge transport and waveform simulations.
The extension traverses the `SolidStateDetectors.jl` geometry tree recursively from the root downward, building the GDML representation bottom-up as required by Geant4.
The resulting GDML file contains the full geometrical description of the detector and its surroundings, along with the material properties of each component, and is passed directly to `Geant4.jl` to initialize a Geant4 application for Monte Carlo simulations.
This means all parts of the simulation pipeline are always based on the same geometry.


# Research Impact Statement

`SolidStateDetectors.jl` has been used across a variety of experimental setups involving semiconductor detectors, often in combination with Geant4 to generate input event distributions.

Originally developed to perform pulse shape simulations for germanium detectors searching for neutrinoless double-beta decay, its user-friendly interface, flexible geometry definition and modular architecture have allowed it to expand into other areas of detector research, including

- the development of position reconstruction algorithms in cross-strip germanium detectors [@COSI:2026; @Zhang:2025], 
- the characterization of radiation damage in germanium detectors [@Boggs:2025],
- the simulation of precise pulse shapes for pixelized silicon detectors [@Hayen:2023],
- the analysis of charge cloud effects in segmented silicon sensors [@Itoh:2025],
- the creation of pulse shape libraries for segmented planar germanium detectors [@Khandelwal:2026],
- the training of machine learning models for event discrimination [@CycleGAN:2026],
- the determination of detection limits in synchrotron measurements [@Goyal:2025],
- the design of novel detector geometries [@Alharbi:2026; @Alexander:2024].

Beyond its impact in research, the software has also been used for educational purposes, e.g.\ for teaching the working principle of semiconductor detectors during simulation tutorials.

The `Geant4.jl` extension presented in this paper builds on this foundation by integrating particle transport directly into the simulation pipeline.
With event generation, field calculation, charge transport simulation and signal generation unified within a single framework, we expect `SolidStateDetectors.jl` to grow further in the future.


# Example
\label{section:example}

`SolidStateDetectors.jl` is designed to provide an easy-to-use, modular end-to-end simulation workflow with only a few lines of code. 
As a demonstration, this section outlines the steps needed to simulate the detector response for the inverted-coaxial point contact germanium detector 50A to an uncollimated $^{228}$Th source, as reported by the @GERDA:2021.

The first step is to define the detector geometry and other simulation parameters in a YAML configuration file.
The configuration file `"GERDA50A.yaml"`, based on information from the paper, is listed in \autoref{appendix:config_file}.
The `Simulation` object is then initialized from this file:
```Julia
using SolidStateDetectors, Geant4, Unitful; T = Float32
sim = Simulation{T}("GERDA50A.yaml")
```

The next step is to define a particle source. 
`SolidStateDetectors.jl` currently supports defining a `MonoenergeticSource` emitting particles with a fixed energy, or an `IsotopeSource` simulating the decay chain of a radioactive isotope.
In this example, a $^{228}$Th isotope source is placed at $x = 15\,\text{cm}$, $y = 0\,\text{cm}$ and $z = 5\,\text{cm}$:
```Julia
source = IsotopeSource(90, 228, 0.0, 0.0, CartesianPoint(15u"cm",0u"cm",5u"cm")) 
```

In combination, the `Simulation` object and the particle source provide the information needed to generate the GDML file and to instantiate a `G4JLApplication` for running Geant4 simulations:
```Julia
app = G4JLApplication(sim, source)
events = run_geant4_simulation(app, 250_000)
```

The resulting variable `events` contains $250,000$ events with at least one hit in the detector and the corresponding hit positions and energies.
The spatial distribution of hits and the energy spectrum, convolved with Fano noise, are shown in \autoref{figure:event_distribution}.
This event distribution serves as input to the charge transport simulation.

![Event distribution induced by a $^{228}$Th isotope source placed at $x = 15\,\text{cm}$, $y = 0\,\text{cm}$ and $z = 5\,\text{cm}$ as generated using the `Geant4.jl` extension of `SolidStateDetectors.jl`: (a) distribution of the hit positions and (b) spectrum of the raw Geant hit energies convolved with Fano noise. In the energy spectrum, the double-escape peak (DEP), single-escape peak (SEP) and full-energy peak (FEP) associated with the $2614.5\,\text{keV}$ gamma emitted during $^{208}$Tl decay are labeled.\label{figure:event_distribution}](./figures/event_distribution_spectrum.pdf)


Further inputs to the charge transport and waveform simulation are the electric field and weighting potentials, which are calculated using `SolidStateDetectors.jl` via `simulate!`.
In this example, the net impurity concentration is set to a constant value of $-5\cdot10^9\,\text{cm}^{-3}$, where the negative sign denotes p-type impurities:
```Julia
sim.detector = SolidStateDetector(sim.detector,
    ConstantImpurityDensity{T}(-5e9u"cm^-3"))
simulate!(sim, refinement_limits = [0.2,0.1,0.05,0.02,0.01])
```
The refinement limits define the precision of the resulting potentials and fields.


The charge transport simulation is written in a modular way. 
`SolidStateDetectors.jl` offers a variety of charge drift and charge trapping models, and can also easily be extended to allow users to implement custom models.
In this example, the charge drift model from the AGATA detector library [@ADL:2016], scaled to a detector temperature of $95\,\text{K}$, and a charge trapping model based on drift length and thermal motion [@Boggs:2023] are used:
```Julia
sim.detector = SolidStateDetector(sim.detector,
    ADL2016ChargeDriftModel{T}(temperature = 95u"K"))
sim.detector = SolidStateDetector(sim.detector,
    BoggsChargeTrappingModel{T}(Dict("nσe-1"=>0.5e4u"cm", "nσh-1"=>1e4u"cm")))
```

Waveforms are simulated with charge cloud diffusion and self-repulsion enabled:
```Julia
wf = SolidStateDetectors.simulate_waveforms(events, sim,
    number_of_carriers = 40, number_of_shells = 1, Δt = 1u"ns", 
    diffusion = true, self_repulsion = true, max_nsteps = 3_000)
```

\autoref{figure:energy-shift} shows the region around the $^{208}$Tl double-escape peak at $1593\,\text{keV}$. The figure compares energy depositions simulated using `Geant4.jl` and the energies reconstructed from the simulated waveforms using `SolidStateDetectors.jl`, both convolved with Fano noise.
In this example, the presence of charge trapping in the simulation broadens the peak and shifts it to lower energies.
This demonstrates how charge transport effects can create spectral features which are not present in the raw Geant4 output.

![Zoom into the $^{208}$Tl double-escape peak at $1593\,\text{keV}$, simulated with the `Geant4.jl` extension of `SolidStateDetectors.jl` (green), as well as reconstructed from the waveforms after simulating with charge trapping (blue), both convolved with Fano noise.\label{figure:energy-shift}](./figures/energy_shift.pdf)


Many applications using solid state detectors also rely on pulse shape analysis, e.g.\ for event classification.
A standard pulse shape discrimination parameter used by the @GERDA:2022 to distinguish single-site from multi-site interactions is the ratio $A/E$, denoting the maximum slope of the waveform normalized to its amplitude.
\autoref{figure:aoe_correlations} depicts the correlation between this $A/E$ parameter and the $0.5$–$90\%$ rise time as determined from the simulated waveforms using the [`RadiationDetectorDSP.jl`](https://github.com/JuliaPhysics/RadiationDetectorDSP.jl) package.

![Correlation of $A/E$ with the $0.5$–$90\%$ rise time for simulated waveforms in the $^{208}$Tl double-escape peak at $1593\,\text{keV}$, for an uncollimated $^{228}$Th source placed at $x = 15\,\text{cm}$, $y = 0\,\text{cm}$ and $z = 5\,\text{cm}$.\label{figure:aoe_correlations}](./figures/AoE_risetime_correlations.pdf)


The simulation qualitatively predicts the existence of two main populations in the $A/E$ vs.\ rise time plane, similar to the ones reported by the @GERDA:2021.
These two populations are also visible in the marginalized distributions shown in \autoref{figure:aoe_histograms}.
Events closer to the point contact feature short drift paths with little diffusion, yielding higher $A/E$ values and shorter rise times.
Events in the outer volume drift longer and, therefore, expand more due to diffusion and self-repulsion, resulting in lower $A/E$ values and longer rise times.

![Marginalized (a) rise time and (b) $A/E$ distributions for the same events as shown in \autoref{figure:aoe_correlations}.\label{figure:aoe_histograms}](./figures/AoE_risetime_histograms.pdf)

This example demonstrates how straightforward end-to-end simulations have now become due to the integrated Geant4 support. 
It is not meant to perfectly reproduce published results by the @GERDA:2021. 
A good match between the simulated and the measured spectra would require precise knowledge of detector parameters like impurity gradients and geometry details like possible volume tapers. 
It would also require an exact electronics response model. 
Fine-tuning such simulation inputs is possible but well beyond the scope of this paper.

# Future Plans

In the future, several additional features are planned to further enhance the functionality of the `Geant4.jl` extension for `SolidStateDetectors.jl`:

**Support for Additional Event Generators**: 
Integrating additional event generators will allow simulations to cover more specialized scenarios, such as those required for studying rare processes like neutrinoless double-beta decay. 
A natural candidate is the [`DoubleBetaDecayGenerators.jl`](https://github.com/JuliaHEP/DoubleBetaDecayGenerators.jl) package from the JuliaHEP ecosystem [@JuliaHEP:2025]. 
This would allow the full simulation pipeline to be applied to rare-decay search scenarios without requiring external event generation tools.

**Integration of Electric Field Calculations**: 
Incorporating the electric field calculated by `SolidStateDetectors.jl` into the Geant4 simulation will enable the simulation to account for the influence of the internal field on charged particles moving inside the detectors. 
This will allow for more accurate modeling of the interaction between the detectors and radiation, as well as the resulting signal generation.

Together, these additions will broaden the range of supported physics applications while preserving the single configuration file setup.


# AI Usage Disclosure

No generative AI tools were used in the development of this software, the writing of this manuscript or the preparation of supporting materials.


# Acknowledgements

We thank Steven E. Boggs for fruitful discussions about the charge trapping model,
David Hervas Aguilar and Yuan-Ru Lin for their contributions to the `Geant4.jl` extension
and Ariana Pearson for testing the software.
We are especially grateful to Pere Mato, whose `Geant4.jl` Julia wrapper laid the essential foundation for the work presented in this paper.

\newpage

\appendix
\renewcommand{\thesection}{\Alph{section}}

# Appendix A: Example detector configuration file
\label{appendix:config_file}

This is the YAML configuration file for the GERDA inverted-coaxial detector 50A used for the simulations in the \nameref{section:example} section.
A schematic of the construction of the semiconductor and contact geometries is shown in \autoref{figure:csg}.
```
name: GERDA 50A             contacts:                  surroundings:
units:                        - material: HPGe           - name: Detector Holder
  length: mm                    id: 1                      material: Al
  angle: deg                    potential: 0V              potential: 0V
  potential: V                  geometry:                  id: 1
  temperature: K                  tube:                    geometry:
grid:                               r: 10                    tube:
  coordinates: cylindrical          h: 0                       r:
  axes:                       - material: HPGe                   from: 39
    r:                          id: 2                            to: 40
      to: 50                    potential: 3700V               h: 80.4
      boundaries: inf           geometry:                      origin:
    phi:                          union:                         z: 40.2
      from: 0                       - tube:              - name: Cryostat
      to: 0                             r:                 material: Al
      boundaries: periodic                from: 13         potential: 0V
    z:                                    to: 37.1         id: 2
      from: -10                         h: 0               geometry:
      to: 90                        - tube:                  union:
      boundaries: inf                   r:                     - tube:
medium: vacuum                            from: 37.1               r:
detectors:                                to: 37.1                   from: 47
  - semiconductor:                      h: 80.4                      to: 48 
      material: HPGe                    origin:                    h: 100
      temperature: 95K                    z: 40.2                  origin:
      geometry:                     - tube:                          z: 40
        difference:                     r:                     - tube:
          - tube:                         from: 5.25               r: 48
              r: 37.1                     to: 37.1                 h: 1
              h: 80.4                   h: 0                       origin:
              origin:                   origin:                      z: 89.5
                z: 40.2                   z: 80.4
          - tube:                   - tube:
              r: 5.25                   r:
              h: 42                       from: 5.25
              origin:                     to: 5.25
                z: 61.4                 h: 40
                                        origin: 
                                          z: 60.4
                                    - tube:
                                        r: 5.25
                                        h: 0
                                        origin:
                                          z: 40.4
```

![Schematic of the construction of the semiconductor of the GERDA 50A detector using a boolean difference of cylindrical primitive volumes, and of its contacts using boolean unions of cylindrical primitive volumes.\label{figure:csg}](./figures/csg_50a.pdf)


# References
