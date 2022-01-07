# Compartmental epidemic model to assess undocumented infections: applications to SARS-CoV-2 epidemics in Brazil - Datasets and Codes

This repository is part of the article "[]".

## Dictionaries

**Municipalities** :The files (a) *dictES.csv* and (b) *dictPR.csv*  yield some information about municipalities of (a) ES (B) PR states. These files have six columns:

1. *ID*: numeric key regarding calibration of confirmed cases time series
2. *ibgeID*: official code to identify the city
3. *name*: name of the city
4. *intermID*: official code of intermediate region to which the city belongs 
5. *imedID*: official code of immediate region to which the city belongs
6. *totPop2019*: population of the city estimated in 2019 

**Immediate and intermediate regions** The files (a) *dictImed.csv* and (b) *dictInterm.csv* yield some information about (a) Immediate and (b) Intermediate regions of PR and ES. These files have five columns:


1. *ID*: numeric key regarding calibration of confirmed cases time series
2. *imedID* or \verb|intermID|: official code to identify the region
3. *name*: name of the region
4. *state*: state to which the region belongs
5. *totPop2019*: population of the region estimated in 2019 


**States** The file *dictUF.csv* yield some information about PR and ES states. These files have five columns:

1. *ID*: numeric key regarding calibration of confirmed cases time series
2. *ibgeID*: official code to identify the state
3. *name*: name of the state
4. *uf*: abbreviation of the state's name
5. *totPop2019*: population of the state estimated in 2019 

## Time series

**Cases and deaths:** The files (a) *PR.csv*, (b) *ES.csv*, (c) *saopaulo.csv*  and (d) *manaus.csv* yield the time series of confirmed cases and deaths since April 1, 2020 for (a) All cities of PR state, (b) All cities of ES state, (c) São Paulo city and (d) Manaus city. These files have seven columns: 

1. *date*: date
2. *ibgeID*: official code to identify the city
3. *newCases*: new confirmed cases on that day 
4. *newDeaths*: new confirmed deaths on that day
5. *city*: name of the city
6. *totalCases*: accumulated cases 
7. *totalDeaths*: accumulated deaths 

**Calibration:** Within files (a) *imed.zip* and (b) *state.zip* we have the time series of accumulated cases and fatality ratio, aggregated for different geographical levels. In this, we have two types of files: *casesXX.dat* (XX refers to the calibrating IDs mentioned before) are accumulated cases while *lethXX.dat* are the daily fatalities).

## Calibration Code

The file *calibra.f90* is a program written in Fortran that executes the calibration algorithm described on Methods section of the main paper $1000$ times with different epidemiological parameters. This program has four inputs: the time series of accumulated cases and fatality, the initial date for calibration and the population of the region (state, city, etc). Besides that, this program has two output files: *epiQuantities.dat* and *hiddenCompart.dat*. The first has seven columns:
1. Days from the initial time 
2. Calibrated confirmed cases
3. Reference cases
4. Effective reproductive number
5. Fraction of susceptible population
6. Underreporting coefficient
7. Sample

On *hiddenCompart.dat*, we have time series for some compartments in the model: from left to right *S*, *E*, *A*, *I*, *C<sub>A</sub>* + *C<sub>I</sub>*, *R* + *R<sub>I</sub>* + *R<sub>A</sub>* + *D* and sample number.

## Python scripts and figures

**Calculation of underreporting coefficient**: the file *underreporting.ipynb* is a I-python script that calculates the underreporting coefficient starting from a time series of confirmed cases and deaths. At the end, it exhibits a graphic showing the evolution of this coefficient.

**Template for figures** The majority of figures in this work were generated with *matplotlib* and *seaborn* packages of *Python 3.7*. File *format_covid19br.mplstyle* contains the template (font family and sizes) for generating those figures and graphics.
