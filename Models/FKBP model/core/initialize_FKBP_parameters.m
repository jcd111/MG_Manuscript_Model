function [parameters, bounds, parameter_names] = initialize_FKBP_parameters()
% Function for setting parameters and bounds for simple equilibrium model
    %% Parameters Governing Protein  Abundance

    % Total amount of MG
    % Units: nM
    % Source: Set by dose when simulating
    parameters.Gt = 5000;
    bounds.Gt_bounds = nan(1,2);
    
    % Total amount of effector protein.
    % Units: nM
    % Source: estimated
    parameters.Et = 20;
    bounds.Et_bounds = [5 2000];
    
    % Total amount of Target Protein.
    % Units: nM
    % Source: Estiamted
    parameters.Pt = 40;
    bounds.Pt_bounds = [5 2000];
   
   
    %% Dissociation Constants
       

    % Dissociation constant of E + G <-> EG for Compound 5
    % Units: nM
    % Source: estimated
    parameters.KD3_Cmpd5 = 1000;
    bounds.KD3_Cmpd5_bounds = [1e2 1e7];

    % Dissociation constant of E + G <-> EG for Compound 6
    % Units: nM
    % Source: estimated
    parameters.KD3_Cmpd6= 1000;
    bounds.KD3_Cmpd6_bounds = [1e2 1e7];

    % Dissociation constant of E + G <-> EG for Compound 7
    % Units: nM
    % Source: estimated
    parameters.KD3_Cmpd7 = 1000;
    bounds.KD3_Cmpd7_bounds = [1e2 1e7];


    % Dissociation constant of EG + P <-> EGP for Cmpd and FRB
    % Units: nM
    % Source: estimated
    parameters.KD4_Cmpd5_FRB = 1;
    bounds.KD4_Cmpd5_FRB_bounds = [1e-4 1e7];
    parameters.KD4_Cmpd6_FRB = 1;
    bounds.KD4_Cmpd6_FRB_bounds = [1e-4 1e7];
    parameters.KD4_Cmpd7_FRB = 1;
    bounds.KD4_Cmpd7_FRB_bounds =  [1e-4 1e7];




    %% Parameters governing drug dosed and target being measured.

    % Drug being dosed.
    parameters.G_name = "Cmpd5";

    % Target being measured.
    parameters.P_name = "FRB";

    
    %% Storing parameter names
    
    parameter_names = fieldnames(parameters);


    
end