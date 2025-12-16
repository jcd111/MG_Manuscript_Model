%% Loading Models

clear;
close all;
clc;
addpath(genpath("Modeling Functions"))
addpath(genpath("Models/Combined Degrader Model"))
addpath(genpath("Combined Degrader Fitting"))

ConditionsToOpt = 1:90;
ConditionsToOpt([3 8 11 12 13 14 18 20 23 28 39 40 41 46]) = [];
% ConditionsToOpt = [1,2,4,5,31,32,33,34];

filenames = ["2025_09_23_yamanaka_fit_n50_PSW_EM_SEw100_V1";
    "2025_09_23_yamanaka_fit_n50_PSW_EM_SEw100_V2";
    "2025_09_23_yamanaka_fit_n50_PSW_EM_SEw100_V3"];

save_plots = true;
%% Selecting Parameter Sets within 5% of maximum value.

% combining ensmeble together
ensemble = [];
fitness = [];
for ii = 1:length(filenames)
    load(filenames(ii));
    ensemble = [ensemble, outstruct.OptimizedEnsemble];
    fitness = [fitness; outstruct.fitness];
end

% Sorting and selecting best sets

[fitness,ind] = sort(fitness);
ensemble = ensemble(:,ind);

fitness = fitness(fitness <= 1.05*min(fitness));
ensemble = ensemble(:,fitness <= 1.05*min(fitness));

model = outstruct.model;
EnsembleParameters = outstruct.fitted_parameters;

% plotting histogram of fitness values

figure(13)
histogram(fitness,10)
xlabel("SSE")
ylabel("Frequency")
title("Distribution of FItness Values")
if save_plots
    saveas(gcf,"Combined Degrader Fitting/Finalized Fitness Figures/Yamanaka et al/Fitness Histogram.svg")
end


%% Calculating Ensemble Statistics

sd = std(ensemble,0,2);
avg = mean(ensemble,2);
CV = sd./avg;

stat_summary = table(EnsembleParameters',avg,sd,CV);

%% Simulating HIBIT Experiments

HIBIT_data = model.experimental_data(ConditionsToOpt(65:70));
tspan = [0 16];
G_range = logspace(-4,7,100);
P_sim = cell(length(HIBIT_data),1);
EP_ind = find(strcmp(model.species_names,"EP"),1);
EGP_ind = find(strcmp(model.species_names,"EGP"),1);
P_ind = find(strcmp(model.species_names,"P"),1);

for jj = 1:size(ensemble,2)
    model = set_model_parameters(model, ensemble(:,jj),EnsembleParameters);
    for ii = 1:length(P_sim)
        condition = HIBIT_data{ii};
        model.parameters.P_name = condition{6};
        model.parameters.G_name = condition{5};
        model.parameters.experiment_type = condition{8};
        for kk = 1:length(G_range)
            model.parameters.Gt = G_range(kk);
            ModelPredictions = model.evaluate(tspan);
            y_baseline = model.parameters.(append("Pt_HIBIT"));
            EP = ModelPredictions(end,EP_ind);
            EGP = ModelPredictions(end,EGP_ind);
            P = ModelPredictions(end,P_ind);
            P_sim{ii}(kk,jj) = (EP + EGP + P)./y_baseline;
        end
    end
end

%% Plotting HIBIT Fits

c = [0.8500 0.3250 0.0980; 0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];

% IKZF1 Titration (1:2)
figure(20)
clf
hold on
for ii = 1:2
    P_cond = P_sim{ii};
    P_mean = mean(P_cond,2);
    plot(G_range,P_mean,"LineWidth",0.5, "Color",c(ii,:));
    fillspecs = {"Color",c(ii,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range,P_cond, [0.05 0.95],fillspecs,'off');
end
for ii = 1:2
    data_cond = HIBIT_data{ii};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii,:));
end

xlabel('Dose (nM)')
ylabel('Normalized HiBiT Signal')
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Location","Best");
title(l,'Treatment')
title("IKZF1 HiBiT Response")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
ylim([0 1])
yticks(0:0.5:1)
xticks([10^-4 10^-2 10^0 10^2 10^4 10^6])
xlim([10^-4 10^6])
legend box off
if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Yamanaka et al/HIBIT/IKZF1.pdf","ContentType","vector")
end


% SALL4 Titration (3:4)
figure(21)
clf
hold on
for ii = 3:4
    P_cond = P_sim{ii};
    P_mean = mean(P_cond,2);
    plot(G_range,P_mean,"LineWidth",0.5, "Color",c(ii-2,:));
    fillspecs = {"Color",c(ii-2,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range,P_cond, [0.05 0.95],fillspecs,'off');
end
for ii = 3:4
    data_cond = HIBIT_data{ii};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii-2,:));
end

xlabel('Dose (nM)')
ylabel('Normalized HiBiT Signal')
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Location","Best");
title(l,'Treatment')
title("SALL4 HiBiT Response")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
ylim([0 1])
yticks(0:0.5:1)
xticks([10^-4 10^-2 10^0 10^2 10^4 10^6])
xlim([10^-4 10^6])
legend box off
if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Yamanaka et al/HIBIT/SALL4.pdf","ContentType","vector")
end
% PLZF Titration (5:6)
figure(22)
clf
hold on
for ii = 5:6
    P_cond = P_sim{ii};
    P_mean = mean(P_cond,2);
    plot(G_range,P_mean,"LineWidth",0.5, "Color",c(ii-4,:));
    fillspecs = {"Color",c(ii-4,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range,P_cond, [0.05 0.95],fillspecs,'off');
end
for ii = 5:6
    data_cond = HIBIT_data{ii};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii-4,:));
end

xlabel('Dose (nM)')
ylabel('Normalized HiBiT Signal')
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Location","Best");
title(l,'Treatment')
title("PLZF HiBiT Response")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
ylim([0 1])
yticks(0:0.5:1)
xticks([10^-4 10^-2 10^0 10^2 10^4 10^6])
xlim([10^-4 10^6])
legend box off
legend box off
if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Yamanaka et al/HIBIT/PLZF.pdf","ContentType","vector")
end
%% Simulating AS Experiments.

AS_data = model.experimental_data(ConditionsToOpt(71:end));
tspan = [0 1];
G_range = logspace(0,6,100);
P_sim = cell(length(AS_data),1);
EP_ind = find(strcmp(model.species_names,"EP"),1);
EGP_ind = find(strcmp(model.species_names,"EGP"),1);
for jj = 1:size(ensemble,2)
    model = set_model_parameters(model, ensemble(:,jj),EnsembleParameters);
    for ii = 1:length(P_sim)
        condition = AS_data{ii};
        model.parameters.P_name = condition{6};
        model.parameters.G_name = condition{5};
        model.parameters.experiment_type = condition{8};
        model.parameters.Gt = 0;
        BaselinePredictions = model.evaluate(tspan);
        y_baseline = BaselinePredictions(end,EP_ind) + BaselinePredictions(end,EGP_ind);
        for kk = 1:length(G_range)
            model.parameters.Gt = G_range(kk);
            ModelPredictions = model.evaluate(tspan);
            EP = ModelPredictions(end,EP_ind);
            EGP = ModelPredictions(end,EGP_ind);
            P_sim{ii}(kk,jj) = (EP + EGP)./y_baseline;
        end
    end
end

%% Plotting Results

c = [0.8500 0.3250 0.0980; 0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];

% IKZF1 Titration (1:2)
figure(23)
clf
hold on
for ii = 1:2
    P_cond = P_sim{ii};
    P_mean = mean(P_cond,2);
    plot(G_range,P_mean,"LineWidth",0.5, "Color",c(ii,:));
    fillspecs = {"Color",c(ii,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range,P_cond, [0.05 0.95],fillspecs,'off');
end
for ii = 1:2
    data_cond = AS_data{ii};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii,:));
end

xlabel('Dose (nM)')
ylabel('Normalized AlphaScreen Signal')
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Location","Best");
title(l,'Treatment')
title("IKZF1 AlphaScreen Response")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
ylim([0 120])
yticks(0:20:120)
xticks([10^0 10^2 10^4 10^6])
xlim([10^0 10^6])
legend box off
if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Yamanaka et al/AS/IKZF1.pdf","ContentType","vector")
end
% SALL4 Titration (3:4)
figure(24)
clf
hold on
for ii = 3:4
    P_cond = P_sim{ii};
    P_mean = mean(P_cond,2);
    plot(G_range,P_mean,"LineWidth",0.5, "Color",c(ii-2,:));
    fillspecs = {"Color",c(ii-2,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range,P_cond, [0.05 0.95],fillspecs,'off');
end
for ii = 3:4
    data_cond = AS_data{ii};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii-2,:));
end

xlabel('Dose (nM)')
ylabel('Normalized AlphaScreen Signal')
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Location","Best");
title(l,'Treatment')
title("SALL4 AlphaScreen Response")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
ylim([0 120])
yticks(0:20:120)
xticks([10^0 10^2 10^4 10^6])
xlim([10^0 10^6])
legend box off
if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Yamanaka et al/AS/SALL4.pdf","ContentType","vector")
end
% PLZF Titration (5:6)
figure(25)
clf
hold on
for ii = 5:6
    P_cond = P_sim{ii};
    P_mean = mean(P_cond,2);
    plot(G_range,P_mean,"LineWidth",0.5, "Color",c(ii-4,:));
    fillspecs = {"Color",c(ii-4,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range,P_cond, [0.05 0.95],fillspecs,'off');
end
for ii = 5:6
    data_cond = AS_data{ii};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii-4,:));
end

xlabel('Dose (nM)')
ylabel('Normalized AlphaScreen Signal')
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Location","Best");
title(l,'Treatment')
title("PLZF AlphaScreen Response")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
ylim([0 120])
yticks(0:20:120)
xticks([10^0 10^2 10^4 10^6])
xlim([10^0 10^6])
legend box off
if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Yamanaka et al/AS/PLZF.pdf","ContentType","vector")
end