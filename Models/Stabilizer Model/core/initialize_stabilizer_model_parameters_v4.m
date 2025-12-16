function [parameters,bounds,parameter_names] = initialize_stabilizer_model_parameters_v4()
% Function for initializing simple reversibe binding and protein
% degradation model parameters


parameters = struct();


%% Association/dissociation rates

% E + P -> EP
parameters.kf1 = 3;
bounds.kf1_bounds = [1e-4 1];

parameters.KD1 = 1;
bounds.KD1_bounds = [1 1e4];


parameters.kf2 = 2;
bounds.kf2_bounds = [1e-4 1];


parameters.KD2 = 1;
bounds.KD2_bounds = [1 1e4];



parameters.kcat = 0.2;
bounds.kcat_bounds = [1e-4 1e1];


% Total P for PLZF
parameters.Pt = 40;
bounds.Pt_bounds = [10 100];

% Total E
parameters.Et = 40;
bounds.Et_bounds = [10 100];

% Total G
% Given by dose
parameters.Gt = 0;
bounds.Gt_bounds = nan(1,2);

parameters.exp_name= "TRFRET";


parameter_names = fieldnames(parameters);

end
