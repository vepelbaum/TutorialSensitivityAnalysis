# Tutorial: Sensitivity Analysis of Computational Models in Psychology

This repository contains the code and sampling matrices for “Sensitivity Analysis of Computational Models in Psychology: A Tutorial”. 

We present a practical workflow for conducting global sensitivity analysis (SA) of computational models in psychology. Using the General Escape Theory (GET) model of suicidal thoughts as an example, we demonstrate two approaches:

* Variance-based Sobol sensitivity analysis
* Distribution-based PAWN sensitivity analysis

The full script for the sensitivity analysis can be found in [Full_Analysis_GET.R](Full_Analysis_GET.R). To run this analysis, source [GET_model.R](GET_model.R) to load the GET model first. Computing the output for the sample matrices might take a while. The precompted matrices and their output can be found in the [samples](samples) folder.

### Reference:

Preprint: https://doi.org/10.17605/OSF.IO/HZ8AT

Note: All results in this repository are based on simulations of the GET model. The repository contains no empirical or participant data.
