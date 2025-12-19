# MG_Manuscript_Model
This repository contains the code for simulating the models presented in "Biophysical odeling elucidates mechanistic principles for rational molecular glue design", _Under Review_.

## Contents

**Models:** Folder containing subfolders with code for simulating each model using the modelstructure class definition. 1433model: subfolder containing code for simulating 1433 stabilizer model. Combined Degrader Model: subfolder containing code for simulating degrader model with and without degradation terms.  FKBP model: subfolder containing code for simulating FKBP stabilizer model. Stabilizer Model: subfolder containing code for simulating stabilizer model with covalent modification.

**calculate_combined_degrader_fitness:** Function for calculating degrader model fitness to BRET, AS, iMRM, and HiBiT data.

**estimate_degrader_model_parameters:** Script for estimating degrader model parameters by fitting to data.

**fit_degrader_model_data:** Function for fitting degrader models to data.

**initialize_degrader_combined_data:** Function for initializing MG degrader data in format required for fitting.  Data taken from Sperling et al. and Yamanaka et al.

**modelstructure:** Class definition for model class.

**plot_sperling_EM_fit_dist:** Script for plotting and analyzing parameter values for MG degrader model fit to Sperling et al data.

**plot_yamanaka_EM_fit_dist:** Script for plotting and analyzing parameter values for MG degrader model fit to Yamanaka et al data.

**fit_stabilizer_model_schulze_data**: Script for estimating parameters of the CYPA stabilizer model by fitting to TR-FRET (ternary complex formation) and covalent modification (k_obs) data for Cmpd3, Cmpd4, and RMC-4998.

**plot_schulze_EM_fit**: Script for plotting the ensemble modeling results of the Schulze data, including dose-response curves for both TR-FRET and covalent kinetic data with confidence intervals.

**fit_stabilizer_model_vickery_data**: Script for fitting the 14-3-3 cooperative binding model to BRET data for various targets (ERa, CRAF, TAZ) using an ensemble modeling approach with particle swarm optimization.

**plot_vickery_EM_fit**: Script for visualizing the 14-3-3 model fits, generating relative BRET signal plots with ensemble statistics (mean and 95% confidence intervals) for ERa, CRAF, and TAZ.

**fit_stabilizer_model_deutscher_data**: Script for estimating binding parameters (KD) for the FKBP stabilizer model by fitting to AlphaScreen experimental data for compounds Cmpd5, Cmpd6, and Cmpd7.

**plot_deutscher_EM_fit**: Script for plotting the FKBP model results, showing the normalized BRET/AlphaScreen signal fits and calculating ensemble statistics for the optimized parameters.

**se:** Function for calculating standard squared error.

## References

Adam S. Sperling, Michael Burgess, Hasmik Keshishian, Jessica A. Gasser, Shruti Bhatt, Max Jan, Mikołaj Słabicki, Rob S. Sellar, Emma C. Fink, Peter G. Miller, Brian J. Liddicoat, Quinlan L. Sievers, Rohan Sharma, Dylan N. Adams, Elyse A. Olesinski, Mariateresa Fulciniti, Namrata D. Udeshi, Eric Kuhn, Anthony Letai, Nikhil C. Munshi, Steven A. Carr, Benjamin L. Ebert; Patterns of substrate affinity, competition, and degradation kinetics underlie biological activity of thalidomide analogs. Blood 2019; 134 (2): 160–170. doi: https://doi.org/10.1182/blood.2019000789

Yamanaka, S., Furihata, H., Yanagihara, Y. et al. Lenalidomide derivatives and proteolysis-targeting chimeras for controlling neosubstrate degradation. Nat Commun 14, 4683 (2023). https://doi.org/10.1038/s41467-023-40385-9

Schulze, C. J. et al. Chemical Remodeling of a Cellular Chaperone to Target the Active State of Mutant KRAS. Science 2023, 381 (6659), 794–799. https://doi.org/10.1126/science.adg9652.

Deutscher, R. C. E. et al. Discovery of Fully Synthetic FKBP12-mTOR Molecular Glues. Chem. Sci. 2025, 16 (10), 4256–4263. https://doi.org/10.1039/D4SC06917J.

Vickery, H. R. et al. Development of a NanoBRET Assay for Evaluation of 14-3-3σ Molecular Glues. SLAS Discovery 2024, 29 (5), 100165. https://doi.org/10.1016/j.slasd.2024.100165.


