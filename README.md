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

The results in the paper are reproduced as follows (please make sure each cpp file is compiled before running the R script in the corrsponding folder):

- **Section 4.1.1**:  
     - *script.R* runs the simulation and then *plot.R* does the plotting works in the main paper and in appendix B.
     - the simulation/plotting scripts of Appendix A and C can be found in the separate folders.

- **Section 4.1.2**:
	 - *script.R* runs the simulation and then *plot.R* does the plotting works in the main paper.


- **Section 4.2-4.3**:
	 - *script.R* \bold{both} runs the analysis and does the plotting works.
