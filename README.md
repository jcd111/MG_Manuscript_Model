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

**se:** Function for calculating standard squared error.

## References

Adam S. Sperling, Michael Burgess, Hasmik Keshishian, Jessica A. Gasser, Shruti Bhatt, Max Jan, Mikołaj Słabicki, Rob S. Sellar, Emma C. Fink, Peter G. Miller, Brian J. Liddicoat, Quinlan L. Sievers, Rohan Sharma, Dylan N. Adams, Elyse A. Olesinski, Mariateresa Fulciniti, Namrata D. Udeshi, Eric Kuhn, Anthony Letai, Nikhil C. Munshi, Steven A. Carr, Benjamin L. Ebert; Patterns of substrate affinity, competition, and degradation kinetics underlie biological activity of thalidomide analogs. Blood 2019; 134 (2): 160–170. doi: https://doi.org/10.1182/blood.2019000789

Yamanaka, S., Furihata, H., Yanagihara, Y. et al. Lenalidomide derivatives and proteolysis-targeting chimeras for controlling neosubstrate degradation. Nat Commun 14, 4683 (2023). https://doi.org/10.1038/s41467-023-40385-9

