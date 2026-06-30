# Monte Carlo Simulation for the Maximum Estimator of a Uniform Distribution

This project presents a Monte Carlo simulation study of the maximum estimator for a Uniform(0, θ) distribution using **R**. The objective is to compare theoretical properties of the estimator with Monte Carlo approximations and evaluate the accuracy of the simulation under different experimental settings.

## Project Objectives

- Derive the theoretical bias, variance, and mean squared error (MSE) of the maximum estimator.
- Approximate these quantities using Monte Carlo simulation.
- Compare analytical and simulated results.
- Investigate the effect of sample size on estimation accuracy.
- Evaluate simulation error using 95% confidence intervals.
- Assess the robustness of the simulation across different values of θ.

## Methodology

The project follows these steps:

1. Implement the theoretical formulas for the bias, variance, and MSE of the maximum estimator.
2. Generate random samples from a Uniform(0, θ) distribution.
3. Perform Monte Carlo simulations with repeated sampling.
4. Compare theoretical and simulated quantities.
5. Analyze the effect of sample size on estimator performance.
6. Construct confidence intervals for simulation error.
7. Investigate the influence of different values of θ on the simulation results.

## Results

The simulation results closely match the theoretical values of the bias, variance, and mean squared error.

The analysis shows that:

- Monte Carlo estimates converge to the theoretical values as the sample size increases.
- Simulation errors remain within the corresponding 95% confidence intervals.
- The simulation framework provides stable and reliable estimates across different values of θ.

These findings are consistent with the theoretical behavior of the maximum estimator and illustrate the effectiveness of Monte Carlo simulation for validating statistical properties.


## Technologies

- R
- R Markdown
- Monte Carlo Simulation
- Statistical Computing
- Probability Theory
- Data Visualization

## Authors

**Irmak Özveren**  
