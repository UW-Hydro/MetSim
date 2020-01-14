---
title: 'MetSim: A Python package for estimation and disaggregation of meteoroligcal data'
tags:
    - Python
    - meteorology
    - hydrology
authors:
    - name: Andrew R. Bennett
      orcid: 0000-0002-7742-3138
      affiliation: 1
    - name: Joseph J. Hamman
      orcid: 0000-0001-7479-8439
      affiliation: 2
    - name: Bart Nijssen
      orcid: 0000-0002-4062-0322
      affiliation: 1
affiliations:
    - name: Department of Civil and Environmental Engineering, University of Washington
      index: 1
    - name: Climate and Global Dynamics Laboratory, National Center for Atmospheric Research
      index: 2
date: 08 January 2020
bibliography: paper.bib
---

# Summary

Hydrometeorological modeling is concerned with regions that have uncertain or unknown boundary conditions.
For example, incoming shortwave radiation, longwave radiation, and humidity have sparse and irregularly placed sensor locations with varying record lengths and observation intervals.
Further, even when such quantities are measured it is often at a daily resolution, while many environmental models require finer temporal resolution for simulation.
To provide the necessary data to solve the model equations in such circumstances we must be able to provide estimates for these quantities at the appropriate temporal resolution.
``MetSim`` is a Python package and standalone tool for the estimation of meteorological quantities at variable temporal resolutions that can address the issues described above.
``MetSim`` can be used to generate spatially distributed sub-daily timeseries of incoming shortwave radiation, outgoing longwave radiation, air pressure, specific humidity, relative humidity, vapor pressure, precipitation, and air temperature given daily timeseries of minimum temperature, maximum temperature, and precipitation.
We have based ``MetSim`` on methods from the Mountain Microclimate Simulator (``MTCLIM``) and the forcing preprocessor that was built into the Variable Infiltration Capacity (``VIC``) hydrological model version 4 [@Bohn:2013; @Thornton:1999; @Liang:1994].
``MetSim`` provides a modern workflow, building upon previous tools by improving performance, adding new IO routines, allowing for exact restarts, and providing an extensible architecture which can incorporate new features.

# Architecture and performance

``MetSim``'s architecture follows the common design of a model driver which coordinates high level operations and which delegates computation to several modules to do the bulk of the work.
The top level model driver handles all IO routines as well as job scheduling and parallelism, a schematic representation of ``MetSim``'s architecture is shown in figure 1.

![Figure 1: A schematic representation of the ``MetSim`` software flow](figure1.pdf)

``MetSim`` has three main computation modules for solar geometry, meteorological simulation, and temporal disaggregation.
The solar geometry module computes the daily potential radiation, daylength, transmittance of the atmosphere, and the fraction of daily radiation received at the top of atmosphere during each 30 second interval.
Computations are based on the algorithms described in @Whiteman:1986 as implemented in MTCLIM [@Thornton:1999].
The data from the solar geometry module is fed to the meteorology simulation module along with the input forcings.
``MetSim`` implements the estimation methods discussed in @Bohn:2013 and @Thornton:1999 to estimate the daily mean temperature, shortwave radiation, vapor pressure, and potential evapotranspiration.
If disaggregation to shorter time steps is configured, the data is passed from the meteorology simulation module to the disaggregation module.
@Bohn:2013 provides a further description and evaluation of these algorithms.

``MetSim`` implements several options for parallelism, which are primarily managed by the Dask library [@dask].
We explore ``MetSim``'s computational performance by conducting two scaling experiments.
Strong scaling experiments test how the total runtime is affected by adding processors for a fixed overall problem size.
Weak scaling experiments test how the total runtime is affected by adding processors proportional to the overall problem size.
All of the times reported for the scaling experiments were for a single year run at an hourly time step with default parameter and output settings.
For the strong scaling experiment we ran ``MetSim`` for one year at an hourly timestep over a domain of 6333 cells and ran using 2, 4, 8, 16, 32, 64, and 128 processors.
The time to complete each run is shown in figure 2 panel a.
The results show that scaling is nearly log-linear with the number of processors.

In the weak scaling experiment we ran ``MetSim`` using 2, 4, 8, 16, 32, 64, and 128 processors while varying the number of cells in the run to maintain a constant workload per processor.
We ran 125 cells per 2 processors, resulting in runs of 125, 250, 500, 1000, 2000, 4000, and 8000 cells, respectively.
The results of the weak scaling experiment are shown in panel b of figure 2.
Similarly to the strong scaling experiment, we see increasing penalties for adding additional processors.

![Figure 2: ``MetSim`` scaling performance](figure2.pdf)

# Applications & Related work

``MetSim`` has been used in several research applications predominantly for generating input to hydrologic models, though other applications are possible.
@Bohn:2019 extended the precipitation disaggregation component to include a new option which was shown to result in better streamflow predictions than the default method.
@Cheng:2020 used ``MetSim`` as a component of their modeling framework to explore how reservoirs affect stream temperatures, and how reservoir operations may be able to help mitigate the effects of climate change on warming stream temperatures.
The Climate Toolbox [@ClimateToolbox] uses ``MetSim`` to generate meteorological data as an intermediate step for developing hydrologic predictions.
``MetSim`` has many other possible uses and is developer-friendly enough for it to be extended to provide additional functionality.

# References
