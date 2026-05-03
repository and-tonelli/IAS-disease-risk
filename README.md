# Alien mammal introductions can reshape global viral sharing networks

[![Preprint](https://img.shields.io/badge/Publication-Preprint-D95C5C)](https://www.researchsquare.com/article/rs-7622438/v1)
[![App](https://img.shields.io/badge/VirAliNet-App-4A90D9)](https://and-tonelli.github.io/AlienViralSharing_app/)

## Overview

This repository contains the R scripts used to predict and map viral sharing between wild mammals, with a focus on how the establishment of alien species in novel ranges may reshape global viral sharing networks. The pipeline covers the construction of range overlap and viral sharing networks, dataset assembly, machine learning modelling, and the generation of figures and maps.

---

## Repository Structure

### Setup

* **`0_native_overlap.R`**
  Constructs the range overlap network between mammal species considering only their native distributions. Species pairs are connected if their native ranges overlap spatially.

* **`0_alien_native_overlap.R`**
  Constructs the overlap network between the established alien ranges of introduced mammals and the native ranges of all other mammal species. This captures novel contact opportunities created by alien introductions.

* **`0_Mammals_distances.R`**
  Computes pairwise Gower distances between mammal species using biological and foraging traits. These dissimilarity scores serve as predictors in the modelling pipeline.

* **`0_ViralSharing_Network.R`**
  Builds the observed viral sharing network by linking species pairs that share at least one known virus. This network forms the response variable for model training.

### Analysis

* **`1_Assemble_datasets.R`**
  Integrates viral sharing information, range overlap data, and pairwise distances into training and prediction datasets. Produces one dataset for observed (native range) species pairs and one for novel (alien–native) species pairs.

* **`2_Modelling.R`**
  Trains and tunes an XGBoost classifier to predict viral sharing probability between species pairs. Includes nested cross-validation for hyperparameter tuning and model evaluation, and outputs predicted sharing probabilities for both in-sample and out-of-sample pairs.

* **`3a_PartialDependence3D.R`**
  Generates three-dimensional partial dependence plots to visualise the joint effect of pairs of predictors on predicted viral sharing probability.

* **`3b_Variable_importance.R`**
  Assesses and plots the relative importance of predictors in the trained XGBoost model.

* **`4_FigureNetwork.R`**
  Builds the predicted viral sharing network and produces spatial maps of viral sharing.

---


## Software Requirements

R version and core packages

* **R:** version `[4.4.2]`
* **tidymodels:** version `[1.3.0]`
* **tidyverse:** version `[2.0.0]`
* **terra:** version `[1.8-86]`
