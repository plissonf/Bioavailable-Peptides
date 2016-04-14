# Bioavailable-Peptides
Mining the features that contribute to oral bioavailability in peptide middle space
Project from Postdoctoral position at the University of Queensland (Australia)

Datasets: Final dataset of 576 observations (conformers) with 34 variables. Assembled 96 chemical structures from literature review, generated 3D chemical information in 6 replicates and 34 variables using Maestro (Schrodinger)

Features / Variables (34): Feature engineering (polarity ratios, moment of inertia), Correlation analysis

Model selection/tuning: Exploratory data analysis

R libraries: reshape2, ggplot 2, lattice, corrgram, gclus, TSP, cluster, seriation

Results: Agreed with reported observations in scientific literature

Summary: 
	
	•	Rule-of-Five Breakdown: Distributions of all peptides were represented in 6 scatterplots for 6 rules-of-5 (hydrogen bond donors, hydrogen bond acceptors, lipophilicity index logP, molecular weight, rotatable bonds, polar surface area). We observed that low polarity, low HBD and moderate lipophilicity favor medium to high F% (oral bioavailability)

	•	Dehydration energy and size: We measured the dehydration energy as the free-energy penalty associated with the desolation of peptide groups upon insertion in the membrane (simulated with CHCl3 environment) and, the radius of gyration as the indicator of size. A scatterplot between the two properties showed that a peptide of greater size (R high) and/or of higher polarity (high PSA) value would require more energy for the desolation of polar groups.

	•	Polarity and Non-Polarity Ratios: We created new features based on the ratios of the different components of solvent-accessible surface area (SASA) - polar, non-polar, neutral, others. Peptides with low polar (~20%) and large non-polar surface area (~75%) display medium-high F% in both solvent simulations - CHCl3 and H2O.

	•	Shape: From 3 principal moments of inertia, we created a ternary plot to represent shape diversity among peptide conformers datasets in both solvent simulations - CHCl3 and H2O. Each angle represent one extreme adopted geometry - rod, sphere and disk. Comparing with F% values, none of the geometries were preferred. 

	•	We evaluated the RMSD values between H2O-based and CHCl3-based conformations of a same peptide for all peptide datasets. We noted that peptides with either rigid structure (low RMSD) or highly flexible backbone (high RMSD) exhibit medium high F%.


