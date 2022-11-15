# Polarisation Correlation in Positron Emission Tomography

A project under Professor Pragya Das at the Indian Institute of Technology, Bombay.

The repository is divided into 
- **Simulation** which is composed of GATE macros. [GATE](https://gate.uca.fr/#/admin) is a simulation software for medical imaging and radiotherapy. This code generates ROOT files as output.
- **Simulation Results** consists of visuals of the simulation, and graphs obtained directly from the ROOT data.
- **Reconstruction** aims to include MATLAB code which will reconstruct the original image from simulated data.

Furthermore, a 4-month progress report is uploaded.

### Abstract

Medical imaging forms a crucial component of the non-invasive diagnostic procedures available in modern medicine. Thus, improving the quality of images while simultaneously trying to improve the safety of these procedures is the focus of a large body of research. In this project, I focus on one of the two broad categories of Nuclear Medical Imaging, PET scans, and attempt to understand how to improve their image quality using previously unexploited properties of the process causing the radiation.

Nuclear Medical Imaging takes advantage of the radioactive properties of certain materials and our ability to detect this radiation with high sensitivity. A radioactive material is injected into the body to be imaged, and detectors measure the radiation at different positions and angles. Then the source distribution is mathematically reconstructed. PET scans use a specific kind of radiation, annhilation photons. These are emitted in pairs, and hence they have various correlations in their properties, notably in their polarisations. I attempt to use this information to better identify true pairs from the detected photons and hence improve image quality.
