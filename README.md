# BaTODO: A Novel Multi-Dose Randomized Basket Trial Design for Dose Testing and Optimization in Oncology Drug Development

R codes to implement BaTODO design for dose testing and optimization.

# Description

In the development of oncology drugs, dose testing and optimization are pivotal objectives in early-phase clinical trials. The judicious selection of an appropriate dose for further clinical studies can increase the probability of success in drug development and benefit patients in real-world drug usage. Basket trials enable concurrently evaluating multiple tumor types, thereby expediting the drug development process.
This paper proposes a novel multi-dose randomized basket trial design that considers the tradeoff of toxicity and efficacy endpoints. To facilitate efficient information sharing among tumor types and doses, we develop a novel statistical model that integrates the Bayesian hierarchical model and normal dynamic linear model. Besides, we utilize the Bayesian model averaging approach to address potential heterogeneity among tumor types.  Extensive simulation studies are conducted to evaluate the operating characteristics of the proposed design. The results demonstrate its favorable performance across various scenarios. By presenting a robust and flexible framework for dose testing and optimization, our novel multi-dose randomized basket trial design represents a valuable advancement in oncology drug development, with the potential to benefit both patients and drug developers.

# Functions

The repository includes two functions:

- batodo-2dose.R: The R code includes the function ```batodo_2dose()``` to obtain the operating characteristics of the BaTODO design for 2 doses by simulating trials.
  
  ```rscript
  batodo_2dose(ntrial, peff, ptox, n1, nsample, theta0, theta0T, weightu, delta1,cutoff_utility)
  ```

- batodo-3dose.R: The R code includes the function `batodo_3dose()` to obtain the operating characteristics of the BaTODO design for 3 doses by simulating trials.
  
  ```rscript
  batodo_2dose(ntrial, peff, ptox, n1, nsample, theta0, theta0T, weightu, delta1,cutoff_utility)
  ```

# Inputs

- `ntrial`: The total number of trials to be simulated.

- `nsample`: Maximum sample size for each dose arm.

- `n1`: Stage 1 sample size for each dose arm.

- `peff`: True toxicity rate.

- `ptox`: True efficacy rate.

- `theta0T`: Null toxicity rate.

- `theta0`: Null efficacy rate.

- `weightu`: penalty of toxicity in the utility function.

- `delta1`: maximum acceptable difference in utility.

- `cutoff_utility`: cutoff value for identifying optimal doses.

# Outputs

- `batodo_2dose()` and `batodo_3dose()`will return the operating characteristics of the BaTODOdesign as a data frame, including:
  
  ```
   (1) the number of toxicities, efficacies, and patients treated at each dose arm at interim analysis (interim_data);
   (2) the estimate of toxicity and efficacy of each dose arm at interim analysis (interim_estmate);
   (3) the number of toxicities, efficacies, and patients treated at each dose arm at finalanalysis (final_data);
   (4) the estimate of toxicity and efficacy of each dose arm at final analysis (final_estmate);
   (5) the percentage of idenfitying as optimal dose (optimal_selection);  
   (6) the percentage of idenfitying as promising dose (promising_selection).
  ```

# Example

We consider using the `batodo_2dose()` as an illustration.

- A maximum of 20 patients will be recruited to each dose arm and the interim analysis will be conducted when outcome of half of the patients has been observed. Suppose 4 tumor types are tested and 2 dose levels are considered for each tumor type. The null toxicity and efficacy rate are 40% and 20% respectively.To risk-benefit tradeoff, the penalty of toxicity is set as 0.25 and the cutoff value is 0.4. We can use the following code to simulate the scenario 1 in Table 1.

```rscript
peff <- rbind(c(0.20, 0.20), c(0.20, 0.20), c(0.20, 0.20), c(0.20, 0.20))
ptox <- rbind(c(0.40, 0.40), c(0.40, 0.40), c(0.40, 0.40), c(0.40, 0.40))

batodo_2dose(ntrial, peff, ptox, n1, nsample, theta0, theta0T, weightu, delta1,cutoff_utility)

  -----------------------output------------------------

$interim_data
             Eff-dose1 Eff-dose2 Tox-dose1 Tox-dose2 N-dose1 N-dose2
Tumor type 1         2         1         6         3      10      10
Tumor type 2         3         4         7         5      10      10
Tumor type 3         2         2         2         6      10      10
Tumor type 4         0         3         3         3      10      10

$interim_estmate
             Eff-dose1 Eff-dose2 Tox-dose1 Tox-dose2
Tumor type 1      0.20      0.17      0.51      0.36
Tumor type 2      0.28      0.33      0.60      0.50
Tumor type 3      0.20      0.21      0.31      0.51
Tumor type 4      0.09      0.25      0.36      0.36

$final_data
             Eff-dose1 Eff-dose2 Tox-dose1 Tox-dose2 N-dose1 N-dose2
Tumor type 1         2         2         6         6      10      20
Tumor type 2         3         6         7        10      10      20
Tumor type 3         7         2         5         6      20      10
Tumor type 4         0         5         3         5      10      20

$final_estmate
             Eff-dose1 Eff-dose2 Tox-dose1 Tox-dose2
Tumor type 1      0.21      0.15      0.50      0.34
Tumor type 2      0.27      0.27      0.61      0.50
Tumor type 3      0.31      0.22      0.31      0.50
Tumor type 4      0.11      0.23      0.33      0.29

$optimal_selection
             Dose1 Dose2
Tumor type 1     0     0
Tumor type 2     0     0
Tumor type 3     0     0
Tumor type 4     0     0

$promising_selection
             Dose1 Dose2
Tumor type 1     0     0
Tumor type 2     0     0
Tumor type 3     0     0
Tumor type 4     0     0
```

# Authors and Reference

* Jingyi Zhang, Fangrong Yan, and Ruitao Lin
* Zhang, J., Yan, F., and Lin, R. (2023) “BaTODO: A Novel Multi-Dose Randomized Basket Trial Design for Dose Testing and Optimization in Oncology Drug Development.”
