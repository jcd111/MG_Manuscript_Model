%% Initializing Parameters
clear;
addpath(genpath('Models/1433model'))
addpath(genpath('Modeling Functions'))

[parameters, bounds, parameter_names] = initialize_CF1433_coop_all_parameters();

%% Setting  eval function

eval_function = @(obj,tspan) CF1433_coop_all_eval_function(obj);

%% Setting up experimental data

experimental_data = {};

%% Constructing model.

model = modelstructure('parameters',parameters,'parameter_names',parameter_names,...
    'bounds',bounds,'experimental_data',experimental_data,'species_names',...
    ["EGP","EP","EG","G","E","P"], 'eval_function',eval_function);

save('Models/1433model/1433_coop_all_model_v1','model')

%%
load('1433_coop_all_model_v1.mat')
model.parameters.G_name = "CRAF01";
model.parameters.P_name = "CRAF";
BaselinePredictions = model.evaluate;
EP_Baseline = BaselinePredictions(2);
EGP_Baseline = BaselinePredictions(1);
y_baseline = EP_Baseline + EGP_Baseline;
G_range = logspace(-4,4);
y_predicted = zeros(length(G_range),1);

for ii = 1:length(G_range)
    model.parameters.Gt = G_range(ii);
    tmp1 = model.evaluate;
    y_predicted(ii) = (tmp1(2) + tmp1(1))/y_baseline;
    
end

figure(1)
clf
hold on
plot(G_range,y_predicted,'b','LineWidth',1.5)
set(gca,'Xscale','log')
