# Meteorological factors and non-pharmaceutical interventions explain local differences in the spread of SARS-CoV-2 in Austria

Author: Katharina Ledebur, Complexity Science Hub Vienna [ledebur@csh.ac.at](ledebur@csh.ac.at)


A data-driven modelling approach based on an age-structured model that compares 116 Austrian regions to a suitably chosen control set of regions to explain variations in local transmission rates through a combination of meteorological factors, non-pharamaceutical interventions and mobility. 

Reference: 

## Contents
The repository contains two folders: ```Code``` contains the necessary script to run the simulation and arrive at the results and ```Figures``` containing the Figures used in the paper.
In ```Code``` you will find the Matlab scripts necessary to run to obtain the results and figures in thr paper. 
​
The ```Figures``` folder does not contain Figure 1 since it is a visual abstract.
​
## Requirements
To run the script you will need Matlab (R2017b).

### Data
- Cases

  Our study period ranges from July 1 2020 to May 15 2021. The daily confirmed Covid-19 cases of 116 districts in Austria were collected through the official Austrian 71 COVID-19 disease reporting system (EMS) [1]. A unique identification number was assigned to each positive Covid-19 test. People can be tested positive multiple times. We calculate daily numbers of positive tests per district and age group. In the model we employ the cumulative cases everyday for each one of the four age groups and the 116 districts. 
- Contact matrices

  To calculate the social mixing matrix entry for one of the four new age groups, we build the population-weighted sum over the corresponding social mixing matrix entries provided by Prem et al. (2017) [2]. The social mixing matrices for each federal state (used for the control sets) are again given by thepopulation-weighted sum over the individual districts.  
- Population

  The population size per district was available from the national statistics office and linked by the unique district ID.
- Input variable timeseries

  The eleven input variable timeseries for metereological factors, regional measures and mobility are available in the Zenodo repository.
- beta weights

  To make the analysis more robust with respect to potential variations in the duration of infectiousness, we performed simulations for a range of beta and calculated a weighted average of the effect sizes employing weights from Paul et al. [3]. The weights follow a gamma distribution with a maximum at 19d.

## References:
[1] https://www.ages.at/en/wissen-aktuell/publikationen/ epidemiologische-parameter-des-covid19-ausbruchs-oesterreich-2020/, accessed: 2021-08-04.

[2] Prem, K., Cook, A., and Jit, M.. Projecting social contact matrices in 152 countries using contact surveys and demographic data. PLOS Computational Biology, 13(9):e1005697, 10.1371/journal.pcbi.1005697 (2017)

[3] Paul, S., Lorin, E., Estimation of COVID-19 recovery and decease periods in Canada using delay model. Scientific Reports. 2021, 11, 23763, https://doi.org/10.1038/s41598-021-02982-w
