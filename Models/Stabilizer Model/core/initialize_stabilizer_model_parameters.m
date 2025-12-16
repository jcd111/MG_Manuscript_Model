function [parameters,bounds,parameter_names] = initialize_stabilizer_model_parameters()
% Function for initializing simple reversibe binding and protein
% degradation model parameters


parameters = struct();


%% Association/dissociation rates

% E + P -> EP
parameters.kf1 = 3;
bounds.kf1_bounds = [1e-4 1];

parameters.kb1 = 4;
bounds.kb1_bounds = [0.001 1e2];


parameters.kf2 = 2;
bounds.kf2_bounds = [1e-4 1];


parameters.kb2 = 1;
bounds.kb2_bounds = [0.001 1e2];


parameters.kcat = 2;
bounds.kcat_bounds = [1e-2 1e4];


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

parameter_names = fieldnames(parameters);

end
