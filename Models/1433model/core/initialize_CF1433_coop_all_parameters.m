function [parameters, bounds, parameter_names] = initialize_CF1433_coop_all_parameters()
% Function for setting parameters and bounds for simple equilibrium model
    %% Parameters Governing Protein  Abundance

    % Total amount of MG
    % Units: nM
    % Source: Set by dose when simulating
    parameters.Gt = 0;
    bounds.Gt_bounds = nan(1,2);
    
    % Total amount of effector protein.
    % Units: nM
    % Source: estimated
    parameters.Et = 500;
    bounds.Et_bounds = [5 2000];
    
    % Total amount of Target Protein.
    % Units: nM
    % Source: Estiamted
    parameters.Pt_CRAF = 40;
    bounds.Pt_CRAF_bounds = [5 2000];
    parameters.Pt_TAZ = 40;
    bounds.Pt_TAZ_bounds = [5 2000];
    parameters.Pt_ERa = 40;
    bounds.Pt_ERa_bounds = [5 2000];
   
    %% Dissociation Constants
       
  
    % Dissociation constant of E + P <-> EP for IKZF1
    % Units: nM
    % Source: estimated
    parameters.KD1_CRAF = 1000;
    bounds.KD1_CRAF_bounds = [1e2 1e7];

    % Dissociation constant of E + P <-> EP for IKZF3
    % Units: nM
    % Source: estimated
    parameters.KD1_TAZ = 1000;
    bounds.KD1_TAZ_bounds = [1e2 1e7];


    % Dissociation constant of E + P <-> EP for GSPT1
    % Units: nM
    % Source: estimated
    parameters.KD1_ERa= 1000;
    bounds.KD1_ERa_bounds = [1e2 1e7];

    parameters.alpha_CRAF01=1000;
    bounds.alpha_CRAF01_bounds = [1e-1 1e4];
    parameters.alpha_CRAF02=1;
    bounds.alpha_CRAF02_bounds = [1e-1 1e4];
    parameters.alpha_CRAF03=1;
    bounds.alpha_CRAF03_bounds = [1e-1 1e4];
    parameters.alpha_CRAF04=1;
    bounds.alpha_CRAF04_bounds = [1e-1 1e4];


    parameters.alpha_TAZ01=1;
    bounds.alpha_TAZ01_bounds = [1e-1 1e4];
    parameters.alpha_TAZ02=1;
    bounds.alpha_TAZ02_bounds = [1e-1 1e4];
    parameters.alpha_TAZ03=1;
    bounds.alpha_TAZ03_bounds = [1e-1 1e4];


    parameters.alpha_ERa01=1;
    bounds.alpha_ERa01_bounds = [1e-1 1e4];
    parameters.alpha_ERa02=1;
    bounds.alpha_ERa02_bounds = [1e-1 1e4];



 
    % % Dissociation constant of EP + G <-> EGP for IKZF1
    % % Units: nM
    % % Source: estimated
    % parameters.KD2_CRAF01_CRAF = 1;
    % bounds.KD2_CRAF01_CRAF_bounds = [1e-4 1e7];
    % parameters.KD2_CRAF02_CRAF = 1;
    % bounds.KD2_CRAF02_CRAF_bounds = [1e-4 1e7];
    % parameters.KD2_CRAF03_CRAF = 1;
    % bounds.KD2_CRAF03_CRAF_bounds =  [1e-4 1e7];
    % parameters.KD2_CRAF04_CRAF = 1;
    % bounds.KD2_CRAF04_CRAF_bounds =  [1e-4 1e7];
    % 
    % 
    % % Dissociation constant of EP + G <-> EGP for GSPT1
    % % Units: nM
    % % Source: estimated
    % parameters.KD2_TAZ01_TAZ = 1;
    % bounds.KD2_TAZ01_TAZ_bounds = [1e-4 1e7];
    % parameters.KD2_TAZ02_TAZ = 1;
    % bounds.KD2_TAZ02_TAZ_bounds = [1e-4 1e7];
    % parameters.KD2_TAZ03_TAZ = 1;
    % bounds.KD2_TAZ03_TAZ_bounds =  [1e-4 1e7];
    % parameters.KD2_TAZ04_TAZ = 1;
    % bounds.KD2_TAZ04_TAZ_bounds =  [1e-4 1e7];
    % 
    % 
    % % Dissociation constant of EP + G <-> EGP for ZFP91
    % % Units: nM
    % % Source: estimated
    % parameters.KD2_ERa01_ERa = 1;
    % bounds.KD2_ERa01_ERa_bounds = [1e-4 1e7];
    % parameters.KD2_ERa02_ERa = 1;
    % bounds.KD2_ERa02_ERa_bounds = [1e-4 1e7];
    % parameters.KD2_ERa03_ERa = 1;
    % bounds.KD2_ERa03_ERa_bounds =  [1e-4 1e7];
    % parameters.KD2_ERa04_ERa = 1;
    % bounds.KD2_ERa04_ERa_bounds =  [1e-4 1e7];


    
    % Dissociation constant of E + G <-> EG for Avadomide
    % Units: nM
    parameters.KD3_CRAF01= 1000;
    bounds.KD3_CRAF01_bounds = [100,10000];
    parameters.KD3_CRAF02= 1000;
    bounds.KD3_CRAF02_bounds = [100,10000];
    parameters.KD3_CRAF03= 1000;
    bounds.KD3_CRAF03_bounds = [100,10000];
    parameters.KD3_CRAF04= 1000;
    bounds.KD3_CRAF04_bounds = [100,10000];

    parameters.KD3_TAZ01= 1000;
    bounds.KD3_TAZ01_bounds = [100,10000];
    parameters.KD3_TAZ02= 1000;
    bounds.KD3_TAZ02_bounds = [100,10000];
    parameters.KD3_TAZ03= 1000;
    bounds.KD3_TAZ03_bounds = [100,10000];
    parameters.KD3_TAZ04= 1000;
    bounds.KD3_TAZ04_bounds = [100,10000];

    parameters.KD3_ERa01= 1000;
    bounds.KD3_ERa01_bounds = [100,10000];
    parameters.KD3_ERa02= 1000;
    bounds.KD3_ERa02_bounds = [100,10000];
    parameters.KD3_ERa03= 1000;
    bounds.KD3_ERa03_bounds = [100,10000];
    parameters.KD3_ERa04= 1000;
    bounds.KD3_ERa04_bounds = [100,10000];





    %% Parameters governing drug dosed and target being measured.

    % Drug being dosed.
    parameters.G_name = "CRAF01";

    % Target being measured.
    parameters.P_name = "CRAF";

    
    %% Storing parameter names
    
    parameter_names = fieldnames(parameters);


    
end