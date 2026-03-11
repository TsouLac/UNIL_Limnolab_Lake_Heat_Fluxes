# Lake Morat / Murten Project


## General Project Information 

Lake Morat still suffers from past excessive eutrophication with elevated phosphate concentration resulting in anoxic condition in the hypolimnion in summer. A study from 2009 (Müller et al. 2009) suggested to decrease the phosphorus loading to 11 tons per year in order to prevent deep low oxygen concentration in Lake Morat. In this project we propose to evaluate such threshold based on: 
- (i) updated biogeochemical knowledge (see for instance Müller et al 2021),
- (ii) impact of ongoing climate change on the lake thermal structure (Råman Vinnå et al. 2020).

Finally we will assess how the exploitation of hydropower operation from Lake Schiffnen to Lake Morat will modify the physics and biogeochemistry of Lake Murten. We will specifically focus on the change in: 
- (i) thermal structure, 
- (ii) the change in Phosphate, and
- (iii) the change in deep oxygen concentration.

## LimnoLab Heat Fluxes

**Goal**

The goal of this LimnoLab is to quantify (i) the importance of the riverine throughflow heat transport compared to the atmospheric heat fluxes and (ii) how this contribution might change with the connection between Lake Schiffenen and Lake Morat.

Specifically we propose to:

* Estimate the heat fluxes (atmospheric and throughflow) over the period 1981 - 2022
* Estimate the relative importance of the throughflow for the period 1981 - 2022 and possible future connection between Lake Schiffenen and Lake Morat (using period 1981 - 2022 with an increased discharged into Lake Morat)

**(Proposed) Method**

* Provide a schematic of the problem including heat sources and sinks that illustrate the research question
* Prepare the dataset (this part is partly provided in this JupyterNotebook)
* Estimate the atmospheric heat flux (codes are provided)
* Estimate the river throughflow heat flux (codes are provided)
* Run the same calculation in case of a future connection between Schiffenen and Morat
* Evaluate the current contribution of the throughflow on the heat budget
* Evaluate the future contribution of Lake Schiffenen on the heat budget
* Discuss the role of the throughflow. You may consider the change in water residence time between the two scenario
 



### References

Fink, G., Schmid, M., Wahl, B., Wolf, T., & Wüest, A. (2014). Heat flux modifications related to climate‐induced warming of large European lakes. Water Resources Research, 50(3), 2072-2085.

Kiefer, I., Steinsberger, T., Wüest, A., & Müller, B. (2021). Netto-Ökosystemproduktion in Seen. Bestimmung aus Monitoring-Daten. Aqua & Gas, 101(4), 22-29.

Müller, B., Steinsberger, T., Schwefel, R., Gächter, R., Sturm, M., & Wüest, A. (2019). Oxygen consumption in seasonally stratified lakes decreases only below a marginal phosphorus threshold. Scientific Reports, 9(1) https://doi.org/10.1038/s41598-019-54486-3

Müller, B., & Schmid, M. (2009). Bilans du phosphore et de l'oxygène dans le lac de Morat. Eawag.

Piccolroaz, S., Zhu, S., Ladwig, R., Carrea, L., Oliver, S., Piotrowski, A. P., ... & Zhu, D. Z. (2024). Lake water temperature modeling in an Era of climate change: Data sources, models, and future prospects. Reviews of Geophysics, 62(1), e2023RG000816.

Raman Vinna, L., Wüest, A., Zappa, M., Fink, G., and Bouffard, D.: Tributaries affect the thermal response of lakes to climate change, Hydrol. Earth Syst. Sci., 22, 31–51, https://doi.org/10.5194/hess-22-31-2018, 2018

Råman Vinnå, L., Wüest, A., & Bouffard, D. (2017). Physical effects of thermal pollution in lakes. Water Resources Research, 53(5), 3968-3987.

Råman Vinnå, L., Medhaug, I., Schmid, M., & Bouffard, D. (2021). The vulnerability of lakes to climate change along an altitudinal gradient. Communications Earth & Environment, 2, 35 (10 pp.). https://doi.org/10.1038/s43247-021-00106-w

Schmid, M., & Read, J. (2022). Heat budget of lakes. Encyclopedia of Inland Waters, second edition, edited by: Mehner, T. and Tockner, K., Elsevier, 467-473.



### Time Frame
- Begin date: 1981
- End data: 2022

### Geographic coverage
Lake Morat [46.93 7.08]


## Data

### Meteorological data
Extracted from the MeteoSwiss station Neuchâtel (47.00 / 6.95)

https://www.meteosuisse.admin.ch/

Parameters:
- air temperature (hourly)
- solar radiation (hourly)
- wind speed (hourly)
- wind direction (hourly)
- relative humidity (hourly)

Data from the meteorological station from Neuchâtel (MeteoSwiss NEU) and can be found here:

 `./data/Forcing_Morat.dat `


### Inflow, Outflow
Extracted from FOEN https://www.hydrodaten.admin.ch/

Inflow: station Broye - Payerne, Caserne d 'aviation 2034

Outflow: Canal de la Broye - Sugiez 2447

Parameters:
- dischage (hourly)
- temperature (hourly)

The inflow temperature and discharge data are already formated as input file for the 1D hydrodyanmic model Simstrat
The discharge data of the inflow data can be found here

 `../data/Qin.dat ` 

The Temperature data of the inflow data can be found here

 `../data/Tin.dat `

We assume that the discharge of the outflow is equal to the discharge of the inflow


### Lake Monitoring 
Cantonal monitoring https://www.die3seen.ch/?lang=fr

Parameters:
- temperature

The lake temperature data provided here are already the output of a 1D hydrodynamic model (Simstrat):

 `../data/Results/T_out.dat `

Cantonal monitoring data can be found here:

 `../data/Temperature_Morat.csv `

They can be used to compare with the model but this is not a requirement here.

### Lake hypsometry
https://www.datalakes-eawag.ch/lakemorphology 

### other data/products
- 1D hydrothermal model https://simstrat.eawag.ch/LakeMurten
- hypsometry https://www.datalakes-eawag.ch/lakemorphology
- Future temperature trend Lake Morat: https://www.datalakes-eawag.ch/datadetail/16
- high frequency temperature Lake Morat (start Aug 2022): https://www.datalakes-eawag.ch/datadetail/956



## Code

[![License: MIT][mit-by-shield]][mit-by]

The easiest way to work with this Renku environment is to clone the Renku project

**Method**


**warning** You need to have [git](https://git-scm.com/downloads) and [git-lfs](https://git-lfs.github.com/) installed in order to successfully clone the repository.

- Clone the repository to your local machine using the command: 

 `git clone https://gitlab.renkulab.io/damien.bouffard/2024-limnolab-unil-lake-heat-fluxes.git`
 
 Note that the repository will be copied to your current working directory.

- Use Python 3 and install the requirements with:

 `pip install -r requirements.txt`

 The python version can be checked by running the command `python --version`. In case python is not installed or only an older version of it, it is recommend to install python through the anaconda distribution which can be downloaded [here](https://www.anaconda.com/products/individual). 

## Usage

### Adapt/Extend the Notebooks scripts

The python script can be updated in `notebook/`

## Folder structure


--data (here you find all the data needed)

--notebooks (here you find all the scripts)

--workflow (no need to open)

### Litterature

Relevant litterature for the project

### Notebook

Jupyter Notebook python scripts



## Authors
- Damien Bouffard, Eawag & UNIL
- Marie-Elodie Perga, UNIL


[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-g.svg?label=Data%20License
[mit-by]: https://opensource.org/licenses/MIT
[mit-by-shield]: https://img.shields.io/badge/License-MIT-g.svg?label=Code%20License
