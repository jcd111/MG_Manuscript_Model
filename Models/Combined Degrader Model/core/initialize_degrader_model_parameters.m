function [parameters,bounds,parameter_names] = initialize_degrader_model_parameters()
% Function for initializing simple reversibe binding and protein
% degradation model parameters


parameters = struct();
%% Association rates and dissociation constants.

% Association rate constant for E + P <-> EP
% Units: /nM/hr
% Source: Schreiber et al., assumed to be identical for all species.
parameters.k1 = 0.36;
bounds.k1_bounds = [0.036 3.6];

% Dissociation constant of E + P <-> EP for IKZF1
% Units: nM
% Source: estimated
parameters.KD1_IKZF1 = 1000;
bounds.KD1_IKZF1_bounds = [1e2 1e7];

% Dissociation constant of E + P <-> EP for IKZF3
% Units: nM
% Source: estimated
parameters.KD1_IKZF3 = 1000;
bounds.KD1_IKZF3_bounds = [1e2 1e7];

% Dissociation constant of E + P <-> EP for GSPT1
% Units: nM
% Source: estimated
parameters.KD1_GSPT1 = 1000;
bounds.KD1_GSPT1_bounds = [1e2 1e7];

% Dissociation constant of E + P <-> EP for ZFP91
% Units: nM
% Source: estimated
parameters.KD1_ZFP91 = 1000;
bounds.KD1_ZFP91_bounds = [1e2 1e7];

% Dissociation constant of E + P <-> EP for RNF166
% Units: nM
% Source: estimated
parameters.KD1_RNF166 = 1000;
bounds.KD1_RNF166_bounds = [1e2 1e7];

% Dissociation constant of E + P <-> EP for ZNF692
% Units: nM
% Source: estimated
parameters.KD1_ZNF692 = 1000;
bounds.KD1_ZNF692_bounds = [1e2 1e7];

% Dissociation constant of E + P <-> EP for SALL4
% Units: nM
% Source: estimated
parameters.KD1_SALL4 = 1000;
bounds.KD1_SALL4_bounds = [1e2 1e7];

% Dissociation constant of E + P <-> EP for PLZF
% Units: nM
% Source: estimated
parameters.KD1_PLZF = 1000;
bounds.KD1_PLZF_bounds = [1e2 1e7];

% Association Rate for EP + G <-> EP
% Units: /nM/hr
% Source: Schreiber et al., assumed to be identical for all species.
parameters.k2 = 0.36;
bounds.k2_bounds = [1e-2 1e1]*0.36;

% Association rate constant for E + G <-> EG for thalidomide
% Units: /nM/hr
% Source: Schreiber et al., assumed to be identical for all species.
parameters.k3 = 0.36;
bounds.k3_bounds = [1e-2 1e1]*0.36;

% Dissociation constant of E + G <-> EG for pomalidomide
% Units: nM
% Source: Akuffo et al., J Biol Chem, 2018.
parameters.KD3_pom = 156.6;
bounds.KD3_pom_bounds = [156.6 156.6];

% Dissociation constant of E + G <-> EG for lenalidomide
% Units: nM
% Source: Akuffo et al., J Biol Chem, 2018.
parameters.KD3_len = 177.8;
bounds.KD3_len_bounds = [177.8 177.8];

% Dissociation constant of E + G <-> EG for avadomide
% Units: nM
% Source: Akuffo et al., J Biol Chem, 2018.
parameters.KD3_ava = 6600;
bounds.KD3_ava_bounds = [6600 6600];

% Dissociation constant of E + G <-> EG for cc-885
% Units: nM
% Source: Akuffo et al., J Biol Chem, 2018.
parameters.KD3_cc885 = 200;
bounds.KD3_cc885_bounds = [200 200];

% Association rate for EG + P <-> EGP
% Units: /nM/hr
% Source: Schreiber et al., assumed to be identical for all species.
parameters.k4 = 0.36;
bounds.k4_bounds = [1e-2 1e1]*0.36;

% Cooperativity for IKZF1
% Units: nM
% Source: estimated
parameters.alpha_len_IKZF1 = 1;
bounds.alpha_len_IKZF1_bounds = [1e0 1e6];
parameters.alpha_pom_IKZF1 = 1;
bounds.alpha_pom_IKZF1_bounds = [1e0 1e6];
parameters.alpha_ava_IKZF1 = 1;
bounds.alpha_ava_IKZF1_bounds =  [1e0 1e6];
parameters.alpha_cc885_IKZF1 = 1;
bounds.alpha_cc885_IKZF1_bounds =  [1e0 1e6];

% Cooperativity for IKZF3
% Units: nM
% Source: estimated
parameters.alpha_len_IKZF3 = 1;
bounds.alpha_len_IKZF3_bounds = [1e0 1e6];
parameters.alpha_pom_IKZF3 = 1;
bounds.alpha_pom_IKZF3_bounds = [1e0 1e6];
parameters.alpha_ava_IKZF3 = 1;
bounds.alpha_ava_IKZF3_bounds = [1e0 1e6];
parameters.alpha_cc885_IKZF3 = 1;
bounds.alpha_cc885_IKZF3_bounds = [1e0 1e6];

% Cooperativity for GSPT1
% Units: nM
% Source: estimated
parameters.alpha_len_GSPT1 = 1;
bounds.alpha_len_GSPT1_bounds = [1e0 1e6];
parameters.alpha_pom_GSPT1 = 1;
bounds.alpha_pom_GSPT1_bounds = [1e0 1e6];
parameters.alpha_ava_GSPT1 = 1;
bounds.alpha_ava_GSPT1_bounds =  [1e0 1e6];
parameters.alpha_cc885_GSPT1 = 1;
bounds.alpha_cc885_GSPT1_bounds = [1e0 1e6];


% Cooperativity for ZFP91
% Units: nM
% Source: estimated
parameters.alpha_len_ZFP91 = 1;
bounds.alpha_len_ZFP91_bounds = [1e0 1e6];
parameters.alpha_pom_ZFP91 = 1;
bounds.alpha_pom_ZFP91_bounds = [1e0 1e6];
parameters.alpha_ava_ZFP91 = 1;
bounds.alpha_ava_ZFP91_bounds = [1e0 1e6];
parameters.alpha_cc885_ZFP91 = 1;
bounds.alpha_cc885_ZFP91_bounds = [1e0 1e6];

% Cooperativity for RNF166
% Units: nM
% Source: estimated
parameters.alpha_len_RNF166 = 1;
bounds.alpha_len_RNF166_bounds = [1e0 1e6];
parameters.alpha_pom_RNF166 = 1;
bounds.alpha_pom_RNF166_bounds = [1e0 1e6];
parameters.alpha_ava_RNF166 = 1;
bounds.alpha_ava_RNF166_bounds = [1e0 1e6];
parameters.alpha_cc885_RNF166 = 1;
bounds.alpha_cc885_RNF166_bounds = [1e0 1e6];

% Cooperativity for ZNF692
% Units: nM
% Source: estimated
parameters.alpha_len_ZNF692 = 1;
bounds.alpha_len_ZNF692_bounds = [1e0 1e6];
parameters.alpha_pom_ZNF692 = 1;
bounds.alpha_pom_ZNF692_bounds = [1e0 1e6];
parameters.alpha_ava_ZNF692 = 1;
bounds.alpha_ava_ZNF692_bounds = [1e0 1e6];
parameters.alpha_cc885_ZNF692 = 1;
bounds.alpha_cc885_ZNF692_bounds = [1e0 1e6];

% Cooperativity for SALL4
parameters.alpha_pom_SALL4 = 1;
bounds.alpha_pom_SALL4_bounds = [1e0 1e6];
parameters.alpha_len_SALL4 = 1;
bounds.alpha_len_SALL4_bounds = [1e0 1e6];

% Cooperativity for PLZF
parameters.alpha_pom_PLZF = 1;
bounds.alpha_pom_PLZF_bounds = [1e0 1e6];
parameters.alpha_len_PLZF = 1;
bounds.alpha_len_PLZF_bounds = [1e0 1e6];

% Scaling for in vitro binding.
% Units: Unitless
% Source: Estiamted
parameters.dK = 1;
bounds.dK_bounds = [1e0 1e3];

%% Degradation Parameters

% Free P degradation rate
parameters.g1 = 0.01;
bounds.g1_bounds = [1e-3 1e0];

% Ubiqutination rate constant for lenalidomide + IKZF1
% Units: unitless
% Source: estimated
parameters.kcat_len_IKZF1 = 1;
bounds.kcat_len_IKZF1_bounds = [0.1 1000];

% Ubiqutination rate constant for lenalidomide + CK1a
% Units: unitless
% Source: estimated
parameters.kcat_len_CK1a = 1;
bounds.kcat_len_CK1a_bounds = [0.1 1000];

% Ubiqutination rate constant for lenalidomide + RNF166
% Units: unitless
% Source: estimated
parameters.kcat_len_RNF166 = 1;
bounds.kcat_len_RNF166_bounds = [0.1 1000];

% Ubiqutination rate constant for lenalidomide + ZFP91
% Units: unitless
% Source: estimated
parameters.kcat_len_ZFP91 = 1;
bounds.kcat_len_ZFP91_bounds = [0.1 1000];

% Ubiqutination rate constant for lenalidomide + ZNF692
% Units: unitless
% Source: estimated
parameters.kcat_len_ZNF692 = 1;
bounds.kcat_len_ZNF692_bounds = [0.1 1000];

% Ubiqutination rate constant for lenalidomide + GSPT1
% Units: unitless
% Source: estimated
parameters.kcat_len_GSPT1 = 1;
bounds.kcat_len_GSPT1_bounds = [0.1 1000];

% Ubiqutination rate constant for lenalidomide + IKZF3
% Units: unitless
% Source: estimated
parameters.kcat_len_IKZF3 = 1;
bounds.kcat_len_IKZF3_bounds = [0.1 1000];

% Ubiqutination rate constant for pomalidomide + IKZF1
% Units: unitless
% Source: estimated
parameters.kcat_pom_IKZF1 = 1;
bounds.kcat_pom_IKZF1_bounds = [0.1 1000];


% Ubiqutination rate constant for pomalidomide + CK1a
% Units: unitless
% Source: estimated
parameters.kcat_pom_CK1a = 1;
bounds.kcat_pom_CK1a_bounds = [0.1 1000];

% Ubiqutination rate constant for pomalidomide + RNF166
% Units: unitless
% Source: estimated
parameters.kcat_pom_RNF166 = 1;
bounds.kcat_pom_RNF166_bounds = [0.1 1000];

% Ubiqutination rate constant for pomalidomide + ZFP91
% Units: unitless
% Source: estimated
parameters.kcat_pom_ZFP91 = 1;
bounds.kcat_pom_ZFP91_bounds = [0.1 1000];

% Ubiqutination rate constant for pomalidomide + ZNF692
% Units: unitless
% Source: estimated
parameters.kcat_pom_ZNF692 = 1;
bounds.kcat_pom_ZNF692_bounds = [0.1 1000];

% Ubiqutination rate constant for pomalidomide + GSPT1
% Units: unitless
% Source: estimated
parameters.kcat_pom_GSPT1 = 1;
bounds.kcat_pom_GSPT1_bounds = [0.1 1000];

% Ubiqutination rate constant for pomalidomide + IKZF3
% Units: unitless
% Source: estimated
parameters.kcat_pom_IKZF3 = 1;
bounds.kcat_pom_IKZF3_bounds = [0.1 1000];

% Ubiqutination rate constant for avadomide + IKZF1
% Units: unitless
% Source: estimated
parameters.kcat_ava_IKZF1 = 1;
bounds.kcat_ava_IKZF1_bounds = [0.1 1000];


% Ubiqutination rate constant for avadomide + CK1a
% Units: unitless
% Source: estimated
parameters.kcat_ava_CK1a = 1;
bounds.kcat_ava_CK1a_bounds = [0.1 1000];

% Ubiqutination rate constant for avadomide + RNF166
% Units: unitless
% Source: estimated
parameters.kcat_ava_RNF166 = 1;
bounds.kcat_ava_RNF166_bounds = [0.1 1000];


% Ubiqutination rate constant for avadomide + ZFP91
% Units: unitless
% Source: estimated
parameters.kcat_ava_ZFP91 = 1;
bounds.kcat_ava_ZFP91_bounds = [0.1 1000];

% Ubiqutination rate constant for avadomide + ZNF692
% Units: unitless
% Source: estimated
parameters.kcat_ava_ZNF692 = 1;
bounds.kcat_ava_ZNF692_bounds = [0.1 1000];

% Ubiqutination rate constant for avadomide + GSPT1
% Units: unitless
% Source: estimated
parameters.kcat_ava_GSPT1 = 1;
bounds.kcat_ava_GSPT1_bounds = [0.1 1000];

% Ubiqutination rate constant for avadomide + IKZF3
% Units: unitless
% Source: estimated
parameters.kcat_ava_IKZF3 = 1;
bounds.kcat_ava_IKZF3_bounds = [0.1 1000];

% Ubiqutination rate constant for cc8-85 + IKZF1
% Units: unitless
% Source: estimated
parameters.kcat_cc885_IKZF1 = 1;
bounds.kcat_cc885_IKZF1_bounds = [0.1 1000];


% Ubiqutination rate constant for cc-885 + CK1a
% Units: unitless
% Source: estimated
parameters.kcat_cc885_CK1a = 1;
bounds.kcat_cc885_CK1a_bounds = [0.1 1000];

% Ubiqutination rate constant for cc-885 + RNF166
% Units: unitless
% Source: estimated
parameters.kcat_cc885_RNF166 = 1;
bounds.kcat_cc885_RNF166_bounds = [0.1 1000];


% Ubiqutination rate constant for cc-885 + ZFP91
% Units: unitless
% Source: estimated
parameters.kcat_cc885_ZFP91 = 1;
bounds.kcat_cc885_ZFP91_bounds = [0.1 1000];

% Ubiqutination rate constant for cc-885 + ZNF692
% Units: unitless
% Source: estimated
parameters.kcat_cc885_ZNF692 = 1;
bounds.kcat_cc885_ZNF692_bounds = [0.1 1000];

% Ubiqutination rate constant for cc-885 + GSPT1
% Units: unitless
% Source: estimated
parameters.kcat_cc885_GSPT1 = 1;
bounds.kcat_cc885_GSPT1_bounds = [0.1 1000];

% Ubiqutination rate constant for cc-885 + IKZF3
% Units: unitless
% Source: estimated
parameters.kcat_cc885_IKZF3 = 1;
bounds.kcat_cc885_IKZF3_bounds = [0.1 1000];


% Ubiqutination rate constant for lenalidomide + SALL4
% Units: unitless
% Source: estimated
parameters.kcat_len_SALL4 = 1;
bounds.kcat_len_SALL4_bounds = [0.1 1000];

% Ubiqutination rate constant for pomalidomide + SALL4
% Units: unitless
% Source: estimated
parameters.kcat_pom_SALL4 = 1;
bounds.kcat_pom_SALL4_bounds = [0.1 1000];

% Ubiqutination rate constant for lenalidomide + PLZF
% Units: unitless
% Source: estimated
parameters.kcat_len_PLZF = 1;
bounds.kcat_len_PLZF_bounds = [0.1 1000];

% Ubiqutination rate constant for lenalidomide + PLZF
% Units: unitless
% Source: estimated
parameters.kcat_pom_PLZF = 1;
bounds.kcat_pom_PLZF_bounds = [0.1 1000];

%% Species abundances

% Total P for IKZF1 in MRM experiments
parameters.Pt_IKZF1_MRM = 40;
bounds.Pt_IKZF1_MRM_bounds = [1e-2 100];

% Total P for CK1a in MRM experiments
parameters.Pt_CK1a_MRM = 40;
bounds.Pt_CK1a_MRM_bounds = [1e-2 100];

% Total P for RNF166 in MRM experiments
parameters.Pt_RNF166_MRM = 40;
bounds.Pt_RNF166_MRM_bounds = [1e-2 100];

% Total P for ZFP91 in MRM experiments
parameters.Pt_ZFP91_MRM = 40;
bounds.Pt_ZFP91_MRM_bounds = [1e-2 100];

% Total P for ZNF692 in MRM experiments
parameters.Pt_ZNF692_MRM = 40;
bounds.Pt_ZNF692_MRM_bounds = [1e-2 100];

% Total P for GSPT1 in MRM experiments
parameters.Pt_GSPT1_MRM = 40;
bounds.Pt_GSPT1_MRM_bounds = [1e-2 100];

% Total P for IKZF3 in MRM experiments
parameters.Pt_IKZF3_MRM = 40;
bounds.Pt_IKZF3_MRM_bounds = [1e-2 100];

% Total P in BRET experiments
parameters.Pt_BRET = 10;
bounds.Pt_BRET_bounds = [1e0 1e2];

% Total P in AlphaScreen experiments.
parameters.Pt_AS = 100;
bounds.Pt_AS_bounds = [1e1 1e3];

% Total P in HIBIT Experiment
parameters.Pt_HIBIT = 1;
bounds.Pt_HIBIT_bounds = [1e-1 1e3];

% Total E in MRM experiments
parameters.Et_MRM = 1;
bounds.Et_MRM_bounds = [1e0 1e2];

% Total E in BRET experiments
parameters.Et_BRET = 1;
bounds.Et_BRET_bounds = [1e-2 1e1];

% Total E in AlphaScreen Experiments
parameters.Et_AS = 100;
bounds.Et_AS_bounds = [1e-1 1e3];

% Total E in HiBit Experiments
parameters.Et_HIBIT = 1;
bounds.Et_HIBIT_bounds = [1e-1 1e3];

% Total G
% Given by dose
parameters.Gt = 0;
bounds.Gt_bounds = nan(1,2);

%% Parameters governing drug dosed and target being measured.

% Drug being dosed.
parameters.G_name = "thal";

% Target being measured.
parameters.P_name = "IKZF1";

%% Experimental Setup

% Type of experiment
parameters.experiment_type = "MRM";


parameter_names = fieldnames(parameters);

end