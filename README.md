# bayesreml

## ​​AAGI-AU-RD-BFA: Exploring a flexible Bayesian framework for ASReml-compatible inference and prediction ​ 

Authors: 
- Max Moldovan (Adelaide University),
- Emi Tanaka (Australian National University),
- Anabel Deltell (University of Valencia),
- Beata Sznajder (Adelaide University),
- Francis Hui (Australian National University)



## Background

The software package ASReml-R (The VSNi Team, 2023) offers a powerful and efficient framework for fitting complex linear mixed models, particularly in the fields of agriculture, plant breeding, and quantitative genetics. Its core strength lies in its ability to handle large unbalanced hierarchical datasets with multiple sources of random variation, often encountered in grain production practices. Unlike many general-purpose modelling tools, ASReml-R provides accurate control over variance structures, including spatial, temporal and autoregressive correlations, and supports complex model terms such as factor-level-specific variances, genotype-by-environment-by-management (GxExM) interactions, and non-linear associations. This flexibility combined with its speed and memory efficiency make it suited for small-to-moderately sized datasets and complex model structures.  

However, in the context of future industry data complexity, it is starting to become clear that there are several constraining issues with using current and future versions of ASReml-R. These are outlined below: 

- It requires a commercial license form VSNi which can be a significant cost for individuals or organizations that are intending to use it regularly. These costs further rise when it’s required to be used in conjunction with HPC. 
- ASReml-R 4.2 is now being supported solely by VSNi and this has resulted in many issues including installation problems as well as bugs that have arisen through cutting edge research undertaken in AAGI. 
- It’s a linear mixed model package that, outside approximation, is inherently constrained by linearity in its underlying model.  
- Again, outside approximate methods, the underlying optimisation method requires tractability of the objective function being optimised (likelihood) to ensure straightforward parameter estimation and evaluation of model parameters. This means, in most cases, that parameters are required to be normally distributed. 
- In some cases the software becomes unwieldy for models containing large numbers of fixed effects or complex interactions.

This project attempts to overcome these constraints through the exploration and development of a next-generation, open-source Bayesian modelling framework closely tailored to the current and future needs of the Australian grains industry. This approach overcomes the constraints through: 

- Having easily accessible open-source Bayesian software that will significantly reduce user costs. 
- Accessible support from software package owners. 
- Flexibility to specify non-linearity in the Bayesian hierarchical models. 
- Flexibility to specify non-normality of parameters and tackle intractable analyses. 
- Flexibility in the size and complexity of data structures.

The proposed solution will explore the open-source suite of R package Bayesian tools available on the market that include, greta, Stan (brms), NIMBLE, JAGS (rjags), and INLA. To ensure industry relevance these will be benchmarked against ASReml-R using relevant data from the industry. The last three dot points will also allow us to explore the potential extended modelling capabilities of each of the software packages to help prepare a guide for their future use with the Australian grains industry. Further details and uses cases can be found in the Analytics & Methodology Justification section. 
