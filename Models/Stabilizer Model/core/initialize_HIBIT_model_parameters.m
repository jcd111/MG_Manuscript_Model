function [parameters,bounds,parameter_names] = initialize_HIBIT_model_parameters()
% Function for initializing simple reversibe binding and protein
% degradation model parameters


    parameters = struct();


    %% Association rates and dissociation constants.

    % Association rate constant for E + P <-> EP
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k1 = 0.36;
    bounds.k1_bounds = [0.036 3.6];

    % Dissociation constant of E + P <-> EP for SALL4
    % Units: nM
    % Source: Mean and SD taken from AlphaScreen Ensemble Fitting
    parameters.KD1_SALL4 = 1.0713e5;
    bounds.KD1_SALL4_bounds = 1.0713e5 + [-1 1]*6.2948e4;

    % Dissociation constant of E + P <-> EP for IKZF1
    % Units: nM
    % Source: estimated
    parameters.KD1_IKZF1 = 1.8542e5;
    bounds.KD1_IKZF1_bounds = 1.8542e5 + [-1 1]*1.1087e5;

    % Dissociation constant of E + P <-> EP for PLZF
    % Units: nM
    % Source: estimated
    parameters.KD1_PLZF = 1.5069e5;
    bounds.KD1_PLZF_bounds = 1.5069e5 + [-1 1]*9.0810e4;
    
    % Association rate for EP + G <-> EGP for thalidomide and IKZF1
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k2_thal_IKZF1 = 0.36;
    bounds.k2_thal_IKZF1_bounds = [1e-2 1e1]*0.36;

    % Association rate for EP + G <-> EGP for thalidomide and SALL4
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k2_thal_SALL4 = 0.36;
    bounds.k2_thal_SALL4_bounds = [1e-2 1e1]*0.36;

    % Association rate for EP + G <-> EGP for thalidomide and PLZF
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k2_thal_PLZF = 0.36;
    bounds.k2_thal_PLZF_bounds = [1e-2 1e1]*0.36;

    % Association rate for EP + G <-> EGP for pomalidomide and IKZF1
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k2_pom_IKZF1 = 0.36;
    bounds.k2_pom_IKZF1_bounds = [1e-2 1e1]*0.36;

    % Association rate for EP + G <-> EGP for pomalidomide and SALL4
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k2_pom_SALL4 = 0.36;
    bounds.k2_pom_SALL4_bounds = [1e-2 1e1]*0.36;

    % Association rate for EP + G <-> EGP for pomalidomide and PLZF
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k2_pom_PLZF = 0.36;
    bounds.k2_pom_PLZF_bounds = [1e-2 1e1]*0.36;


    % Association rate for EP + G <-> EGP for lenalidomide and IKZF1
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k2_len_IKZF1 = 0.36;
    bounds.k2_len_IKZF1_bounds = [1e-2 1e1]*0.36;

    % Association rate for EP + G <-> EGP for lenalidomide and SALL4
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k2_len_SALL4 = 0.36;
    bounds.k2_len_SALL4_bounds = [1e-2 1e1]*0.36;

    % Association rate for EP + G <-> EGP for lenalidomide and PLZF
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k2_len_PLZF = 0.36;
    bounds.k2_len_PLZF_bounds = [1e-2 1e1]*0.36;

    % Dissociation constant of EP + G <-> EGP for thalidomide + SALL4
    % Units: nM
    % Source: estimated
    parameters.KD2_thal_SALL4 = 77.2184;
    bounds.KD2_thal_SALL4_bounds = 77.2184 + [-1 1]*8.2998;

    % Dissociation constant of EP + G <-> EGP for thalidomide + IKZF1
    % Units: nM
    % Source: estimated
    parameters.KD2_thal_IKZF1 = 316.1941;
    bounds.KD2_thal_IKZF1_bounds = 316.1941 + [-1 1]*13.2314;

    % Dissociation constant of EP + G <-> EGP for thalidomide + PLZF
    % Units: nM
    % Source: estimated
    parameters.KD2_thal_PLZF = 540.3116;
    bounds.KD2_thal_PLZF_bounds = 540.3116 + [-1 1]*17.2884;

    % Dissociation constant of EP + G <-> EGP for pomalidomide + SALL4
    % Units: nM
    % Source: estimated
    parameters.KD2_pom_SALL4 = 22.1191;
    bounds.KD2_pom_SALL4_bounds = 22.1191 + [-1 1]*3.2682;

    % Dissociation constant of EP + G <-> EGP for pomalidomide + IKZF1
    % Units: nM
    % Source: estimated
    parameters.KD2_pom_IKZF1 = 39.0702;
    bounds.KD2_pom_IKZF1_bounds = 29.0702 + [-1 1]*3.1205;

    % Dissociation constant of EP + G <-> EGP for pomalidomide + PLZF
    % Units: nM
    % Source: estimated
    parameters.KD2_pom_PLZF = 75.2382;
    bounds.KD2_pom_PLZF_bounds = 75.2382 + [-1 1]*4.5452;

    % Dissociation constant of EP + G <-> EGP for lenalidomide + SALL4
    % Units: nM
    % Source: estimated
    parameters.KD2_len_SALL4 = 46.3994;
    bounds.KD2_len_SALL4_bounds = 46.3994 + [-1 1]*4.9076;

    % Dissociation constant of EP + G <-> EGP for lenalidomide + IKZF1
    % Units: nM
    % Source: estimated
    parameters.KD2_len_IKZF1 = 40.1617;
    bounds.KD2_len_IKZF1_bounds = 40.1617 + [-1 1]*2.9664;

    % Dissociation constant of EP + G <-> EGP for lenalidomide + PLZF
    % Units: nM
    % Source: estimated
    parameters.KD2_len_PLZF = 96.8389;
    bounds.KD2_len_PLZF_bounds = 96.8389 + [-1 1]*4.7083;

    % Association rate constant for E + G <-> EG for thalidomide
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k3_thal = 0.36;
    bounds.k3_thal_bounds = [1e-2 1e1]*0.36;

    % Association rate constant for E + G <-> EG for pomalidomide
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k3_pom = 0.36;
    bounds.k3_pom_bounds = [1e-2 1e1]*0.36;

    % Association rate constant for E + G <-> EG for lenalidomide
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k3_len = 0.36;
    bounds.k3_len_bounds = [1e-2 1e1]*0.36;

    % Dissociation constant of E + G <-> EG for thalidomide
    % Units: nM
    % Source: Akuffo et al., J Biol Chem, 2018.
    parameters.KD3_thal = mean([65 38]*1000);
    bounds.KD3_thal_bounds = [26 105]*1000;
    
    % Dissociation constant of E + G <-> EG for pomalidomide
    % Units: nM
    % Source: Akuffo et al., J Biol Chem, 2018.
    parameters.KD3_pom = mean([16 16]*1000);
    bounds.KD3_pom_bounds = [12 20]*1000;

    % Dissociation constant of E + G <-> EG for lenalidomide
    % Units: nM
    % Source: Akuffo et al., J Biol Chem, 2018.
    parameters.KD3_len = mean([11 13]*1000);
    bounds.KD3_len_bounds = [8 15]*1000;

    % Association rate for EG + P <-> EGP for thalidomide and IKZF1
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k4_thal_IKZF1 = 0.36;
    bounds.k4_thal_IKZF1_bounds = [1e-2 1e1]*0.36;

    % Association rate for EG + P <-> EGP for thalidomide and SALL4
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k4_thal_SALL4 = 0.36;
    bounds.k4_thal_SALL4_bounds = [1e-2 1e1]*0.36;

    % Association rate for EG + P <-> EGP for thalidomide and PLZF
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k4_thal_PLZF = 0.36;
    bounds.k4_thal_PLZF_bounds = [1e-2 1e1]*0.36;

    % Association rate for EG + P <-> EGP for pomalidomide and IKZF1
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k4_pom_IKZF1 = 0.36;
    bounds.k4_pom_IKZF1_bounds = [1e-2 1e1]*0.36;

    % Association rate for EG + P <-> EGP for pomalidomide and SALL4
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k4_pom_SALL4 = 0.36;
    bounds.k4_pom_SALL4_bounds = [1e-2 1e1]*0.36;

    % Association rate for EG + P <-> EGP for pomalidomide and PLZF
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k4_pom_PLZF = 0.36;
    bounds.k4_pom_PLZF_bounds = [1e-2 1e1]*0.36;


    % Association rate for EG + P <-> EGP for lenalidomide and IKZF1
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k4_len_IKZF1 = 0.36;
    bounds.k4_len_IKZF1_bounds = [1e-2 1e1]*0.36;

    % Association rate for EG + P <-> EGP for lenalidomide and SALL4
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k4_len_SALL4 = 0.36;
    bounds.k4_len_SALL4_bounds = [1e-2 1e1]*0.36;

    % Association rate for EG + P <-> EGP for lenalidomide and PLZF
    % Units: /nM/hr
    % Source: Schreiber et al., assumed to be identical for all species.
    parameters.k4_len_PLZF = 0.36;
    bounds.k4_len_PLZF_bounds = [1e-2 1e1]*0.36;

    %% Degradation Parameters

    % Free P degradation rate
    parameters.g1 = 0.01;
    bounds.g1_bounds = [1e-3 1e0];
    % 
    % % Ubiquitination rate
    % parameters.g_EP = 10;
    % bounds.g_EP_bounds = [1e-2 1e2];

    % Ubiquitination rate
    parameters.kcat_EGP_thal_IKZF1 = 10;
    bounds.g_EP_thal_IKZF1_bounds = [1e-2 1e2];

    % Ubiquitination rate
    parameters.kcat_EGP_thal_SALL4 = 10;
    bounds.g_EP_thal_SALL4_bounds = [1e-2 1e2];

    % Ubiquitination rate
    parameters.kcat_EGP_thal_PLZF = 10;
    bounds.g_EP_thal_PLZF_bounds = [1e-2 1e2];

    % Ubiquitination rate
    parameters.kcat_EGP_pom_IKZF1 = 10;
    bounds.kcat_EGP_pom_IKZF1_bounds = [1e-2 1e2];

    % Ubiquitination rate
    parameters.kcat_EGP_pom_SALL4 = 10;
    bounds.kcat_EGP_pom_SALL4_bounds = [1e-2 1e2];

    % Ubiquitination rate
    parameters.kcat_EGP_pom_PLZF = 10;
    bounds.kcat_EGP_pom_PLZF_bounds = [1e-2 1e2];

    % Ubiquitination rate
    parameters.kcat_EGP_len_IKZF1 = 10;
    bounds.kcat_EGP_len_IKZF1_bounds = [1e-2 1e2];

    % Ubiquitination rate
    parameters.kcat_EGP_len_SALL4 = 10;
    bounds.kcat_EGP_len_SALL4_bounds = [1e-2 1e2];

    % Ubiquitination rate
    parameters.kcat_EGP_len_PLZF = 10;
    bounds.kcat_EGP_len_PLZF_bounds = [1e-2 1e2];


    % Ubiquitination rate
    parameters.kcat_EP_IKZF1 = 10;
    bounds.kcat_EP_IKZF1_bounds = [1e-2 1e2];

    % Ubiquitination rate
    parameters.kcat_EP_SALL4 = 10;
    bounds.kcat_EP_SALL4_bounds = [1e-2 1e2];

    % Ubiquitination rate
    parameters.kcat_EP_PLZF = 10;
    bounds.kcat_EP_PLZF_bounds = [1e-2 1e2];


    % Degradation due to ubiqutination.
    parameters.g2 = 1;
    bounds.g2_bounds = [1e-2 1e2];

    %% Species abundances

    % Total P for IKZF1
    parameters.Pt_IKZF1 = 40;
    bounds.Pt_IKZF1_bounds = [10 100];

    % Total P for SALL4
    parameters.Pt_SALL4 = 40;
    bounds.Pt_SALL4_bounds = [10 100];

    % Total P for PLZF
    parameters.Pt_PLZF = 40;
    bounds.Pt_PLZF_bounds = [10 100];

    % Total E
    parameters.Et = 40;
    bounds.Et_bounds = [10 100];

    % Total G
    % Given by dose
    parameters.Gt = 0;
    bounds.Gt_bounds = nan(1,2);

    %% Parameters governing drug dosed and target being measured.

    % Drug being dosed.
    parameters.G_name = "thal";

    % Target being measured.
    parameters.P_name = "IKZF1";


    parameter_names = fieldnames(parameters);

end