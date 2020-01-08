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
---

# Summary

Hydrometeorological modeling is concerned with domains that have uncertain or unknown boundary conditions [@Beven:2012].
For example, incoming shortwave radiation, longwave radiation, and humidity have sparse and irregular sensor locations with varying record lengths and observation intervals.
Further, even when such quantities are measured it is often at a daily resolution, while many environmental models require finer temporal resolution for simulation.
To provide closure to the model equations in such circumstances we must be able to provide estimates for these quantities at the appropriate temporal resolution.
``MetSim`` is a software package and standalone tool for the estimation of sparsely observed meteorological quantities at variable temporal resolutions.
``MetSim`` can be used to generate spatially distributed sub-daily timeseries of incoming shortwave radiation, outgoing longwave radiation, air pressure, specific humidity, relative humidity, vapor pressure, precipitation, and air temperature given daily timeseries of minimum temperature, maximum temperature, and precipitation.
We have based ``MetSim`` on methods from the Mountain Microclimate Simulator (``MTCLIM``) and the forcing preprocessor that was built into the Variable Infiltration Capacity (``VIC``) hydrological model version 4 [@Bohn:2013; @Thornton:1999; @Liang:1994].
MetSim provides a modern workflow, building upon previous tools by improving performance, adding new IO routines, allowing for exact restarts, and providing an extensible architecture which can incorporate new features.


