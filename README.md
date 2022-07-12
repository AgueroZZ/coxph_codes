# coxph_codes
Codes for Bayesian Inference for Cox Proportional Hazard Models with Partial Likelihoods, Nonlinear Covariate Effects and Correlated Observations

## Data

- **Section 4.1.1-4.1.2**: simulated datasets are generated in the R scripts.

- **Section 4.2**: dataset `kidney` in package `survival`.

- **Section 4.3**: dataset `Leuk` in package `INLA`.

## Scripts

The simulation/examples in the main paper are replicated in the following scripts. Specific points to pay attention to:
1. In each simulation/example, please compile the cpp template before run the R scripts.
2. Once the cpp file is compiled, please make sure the compiled file can be found in the same folder as the working directory.
3. All the required R packages will be loaded in the R script, please make sure they are installed in the system before proceed.
