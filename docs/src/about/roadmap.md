# Roadmap

AlphaPEM is an active research project with several planned improvements and extensions.

## Current Development

**Version 2.0** (in development)
    - the transition from the original programming language to Julia, leveraging its high execution speed while
      maintaining a high-level language framework.
    - the abandonment of dictionary usage in favor of increased reliance on object-oriented programming.
    - the redesign of the AlphaPEM architecture so that the code is closer to industry standards.
    - the previous GUI has been replaced by a web-based interface.

## Planned Enhancements

### Near Term

**Spatial Modeling**
- 1D+1D thermal evolution (currently 1D in V2.0)
- 1D+1D+1D full-stack modeling with manifold channels

**Physics Refinements**
- More accurate auxiliary models
- ECSA (electrochemical surface area) degradation integration

### Medium Term

**Application Extensions**
- PEM electrolyzer support
- Degradation and lifetime prediction

**Calibration & Analysis**
- Sensitivity analysis enhancements
- Extended parameter identification workflows

### Long Term

- Advanced control strategies
- Real-time predictive models
- Hardware-in-the-loop integration

## Version History

- [V1.3](https://github.com/gassraphael/AlphaPEM/tree/65dd73ed306a054c80018447f7943b9d9f973ffb) - 2026.02.16 - This version of AlphaPEM includes: 
	- the addition of O2 flow to Pt particules which improves the modeling of overvoltage due to flooding at high curent densities.
		- the limiting liquid water saturation coefficient ($s_{lim}$) has been definitively removed, as this model replaces it.
	- the addition of liquid water flow inside the GC (with the sorption flow at the GDL/GC interface).
	- the spatial extension to 1D+1D (except thermal evolution which remains 1D for now).
- [V1.2](https://github.com/gassraphael/AlphaPEM/tree/b71f42878a186e17efeb7e97b5d7fb50d6e76827) - 2025.12.11 - This version of AlphaPEM includes: 
	- the addition of convective flow between the inlet, gas channels, and outlet of the cell, thereby removing the Pukrushpan equations (from Michigan University).
		- auxiliaries are temporarily removed, as they require reconstruction. 
	- the addition of the MPL to the simulated cell, in both the anode and cathode. 
	- effective diffusive flows for the dissolved water insided the CLs are introduced.
	- the addition of the open-source [ZSW GenStack](https://zenodo.org/records/14223364) as a calibrated fuel cell case study. 
- [V1.1](https://github.com/gassraphael/AlphaPEM/tree/11f07bd084a09cc6432f441b010d89d2a4229e4e) - 2025.08.18 - This version of AlphaPEM includes: 
	- the addition of heat transfer to the program, in cooperation with Pedro Affonso Nobrega (PERSEE, MINES PSL).
	- an improvement of the initial variable values: the algorithm waits for a given time to reach equilibrium, and then the experiment starts (step/pola/EIS).
	- the limiting liquid water saturation coefficient ($s_{lim}$) is temporarily removed for future refinement.
- [V1.0](https://github.com/gassraphael/AlphaPEM/tree/2b042c3d16d53fcd16779a5ffdc81eea75a9d012) - 2024.09.05 - This version of AlphaPEM corresponds to the one developed during Raphaël Gass's PhD from 2021 to 2024. 
	- It is based on a physics-based, one-dimensional (1D), dynamic, two-phase, and isothermal model.

## Contact & Collaboration

For inquiries about the roadmap or collaboration opportunities, contact us at: 
📧 [raphael.gass@univ-reunion.fr](mailto:raphael.gass@univ-reunion.fr)

AlphaPEM is an ongoing research effort. Newer versions may be available pending scientific publication.
Contact us to discuss potential collaboration.

---

**Roadmap last updated:** June 2026
