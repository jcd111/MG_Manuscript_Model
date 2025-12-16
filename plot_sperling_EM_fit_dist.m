%% Loading Models

clear;
% close all;
clc;
addpath(genpath("Modeling Functions"))
addpath(genpath("Models/Combined Degrader Model"))
addpath(genpath("Combined Degrader Fitting"))

ConditionsToOpt = 1:90;
ConditionsToOpt([3 8 11 12 13 14 18 20 23 28 39 40 41 46]) = [];
% ConditionsToOpt = [1,2,4,5,31,32,33,34];

filenames = ["2025_09_22_sperling_fit_n25_PSW_EM_SE_V1";
    "2025_09_22_sperling_fit_n25_PSW_EM_SE_V2";
    "2025_09_22_sperling_fit_n25_PSW_EM_SE_V3"];

save_plots = false;

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

fitness = fitness(fitness <= 1.00*min(fitness));
ensemble = ensemble(:,fitness <= 1.00*min(fitness));

model = outstruct.model;
EnsembleParameters = outstruct.fitted_parameters;

% plotting histogram of fitness values

figure(20)
histogram(fitness,10)
xlabel("SSE")
ylabel("Frequency")
xlim([30 35])
title("Distribution of FItness Values")

if save_plots
    saveas(gcf,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/Fitness Histogram","epsc")
end


%% Calculating Ensemble Statistics

sd = std(ensemble,0,2);
avg = mean(ensemble,2);
CV = sd./avg;

stat_summary = table(EnsembleParameters',avg,sd,CV);

%% Simulating Timecourse Data

timecourse_data = model.experimental_data(ConditionsToOpt(1:20));
tspan = linspace(0,30,1000);
P_sim = cell(length(timecourse_data),1);
EP_ind = find(strcmp(model.species_names,"EP"),1);
EGP_ind = find(strcmp(model.species_names,"EGP"),1);
P_ind = find(strcmp(model.species_names,"P"),1);

for jj = 1:size(ensemble,2)
    model = set_model_parameters(model, ensemble(:,jj), EnsembleParameters);
    for ii = 1:length(P_sim)
        condition = timecourse_data{ii};
        model.parameters.P_name = condition{6};
        model.parameters.G_name = condition{5};
        model.parameters.Gt = condition{4};
        model.parameters.experiment_type = condition{8};
        ModelPredictions = model.evaluate(tspan);
        y_baseline = model.parameters.(append("Pt_",condition{6},"_MRM"));
        EP = ModelPredictions(:,EP_ind);
        EGP = ModelPredictions(:,EGP_ind);
        P = ModelPredictions(:,P_ind);
        P_sim{ii}(:,jj) = (EP + EGP + P)./y_baseline;

    end
end


%% Plotting Results


c = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
% IKZF1 Fit
figure(1)
clf
hold on
for ii = 1:4
    P_cond = P_sim{ii};
    P_mean = mean(P_cond,2);
    plot(tspan,P_mean,"LineWidth",0.5, "Color",c(ii,:));
    fillspecs = {"Color",c(ii,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(tspan,P_cond, [0.05 0.95],fillspecs,'off');
end

for ii = 1:4
    data_cond = timecourse_data{ii};
    t_exp = data_cond{1};
    c_exp = data_cond{2};
    c_err = data_cond{3};

    errorbar(t_exp,c_exp,c_err,'.','MarkerSize',10,'Color',c(ii,:))
end

xlabel("Time (hours)")
ylabel("Normalized iMRM Signal")
title("IKZF1 Degradation Timecourse")
l = legend("Lenalidomide (1{\mu}M)","Pomalidomide (1{\mu}M)","Avadomide (1{\mu}M)","cc885 (0.1{\mu}M)","Location","Best","FontSize",6);
title(l,'Treatment')
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
ylim([0 1])
yticks(0:0.5:1)
xticks(0:10:30)
legend box off

if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/Timecourse/IKZF1.pdf","ContentType","vector")
end

% IKZF3 Fit
figure(2)
clf
hold on
for ii = 1:4
    P_cond = P_sim{ii+4};
    P_mean = mean(P_cond,2);
    plot(tspan,P_mean,"LineWidth",0.5, "Color",c(ii,:));
    fillspecs = {"Color",c(ii,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(tspan,P_cond, [0.05 0.95],fillspecs,'off');
end

for ii = 1:4
    data_cond = timecourse_data{ii+4};
    t_exp = data_cond{1};
    c_exp = data_cond{2};
    c_err = data_cond{3};

    errorbar(t_exp,c_exp,c_err,'.','MarkerSize',10,'Color',c(ii,:))
end

xlabel("Time (hours)")
ylabel("Normalized iMRM Signal")
l = legend("Lenalidomide (1{\mu}M)","Pomalidomide (1{\mu}M)","Avadomide (1{\mu}M)","cc885 (0.1{\mu}M)","Location","Best");
title("IKZF3 Degradation Timecourse")
title(l,'Treatment')
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
ylim([0 1])
yticks(0:0.5:1)
xticks(0:10:30)
legend box off
if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/Timecourse/IKZF3.pdf","ContentType","vector")
end

% GSPT1 Fit
figure(3)
clf
hold on
P_cond = P_sim{9};
P_mean = mean(P_cond,2);
plot(tspan,P_mean,"LineWidth",0.5, "Color",c(4,:));
fillspecs = {"Color",c(4,:), "FaceAlpha",0.3,"EdgeAlpha",0};
plot_ensemble_CI(tspan,P_cond, [0.05 0.95],fillspecs,'off');
data_cond = timecourse_data{9};
t_exp = data_cond{1};
c_exp = data_cond{2};
c_err = data_cond{3};
errorbar(t_exp,c_exp,c_err,'.','MarkerSize',10,'Color',c(4,:))

xlabel("Time (hours)")
ylabel("Normalized iMRM Signal")
l = legend("cc885 (0.1{\mu}M)","Location","Best");
title(l,'Treatment')
title("GSPT1 Degradation Timecourse")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
ylim([0 1])
yticks(0:0.5:1)
xticks(0:10:30)
legend box off

if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/Timecourse/GSPT1.pdf","ContentType","vector")
end

% ZFP91 Fit
figure(4)
clf
hold on
for ii = 1:3
    P_cond = P_sim{ii+9};
    P_mean = mean(P_cond,2);
    plot(tspan,P_mean,"LineWidth",0.5, "Color",c(ii,:));
    fillspecs = {"Color",c(ii,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(tspan,P_cond, [0.05 0.95],fillspecs,'off');
end

for ii = 1:3
    data_cond = timecourse_data{ii+9};
    t_exp = data_cond{1};
    c_exp = data_cond{2};
    c_err = data_cond{3};

    errorbar(t_exp,c_exp,c_err,'.','MarkerSize',10,'Color',c(ii,:))
end

xlabel("Time (hours)")
ylabel("Normalized iMRM Signal")
l = legend("Lenalidomide (1{\mu}M)","Pomalidomide (1{\mu}M)","Avadomide (1{\mu}M)","Location","Best");
title(l,'Treatment')
title("ZFP91 Degradation Timecourse")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
ylim([0 1])
yticks(0:0.5:1)
xticks(0:10:30)
legend box off
if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/Timecourse/ZFP91.pdf","ContentType","vector")
end

% RNF166 Fit
figure(5)
clf
hold on
for ii = 1:4
    P_cond = P_sim{ii+13};
    P_mean = mean(P_cond,2);
    plot(tspan,P_mean,"LineWidth",0.5, "Color",c(ii,:));
    fillspecs = {"Color",c(ii,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(tspan,P_cond, [0.05 0.95],fillspecs,'off');
end

for ii = 1:4
    data_cond = timecourse_data{ii+13};
    t_exp = data_cond{1};
    c_exp = data_cond{2};
    c_err = data_cond{3};

    errorbar(t_exp,c_exp,c_err,'.','MarkerSize',10,'Color',c(ii,:))
end

xlabel("Time (hours)")
ylabel("Normalized iMRM Signal")
l = legend("Lenalidomide (1{\mu}M)","Pomalidomide (1{\mu}M)","Avadomide (1{\mu}M)","cc885 (0.1{\mu}M)","Location","Best");
title(l,'Treatment')
title("RNF166 Degradation Timecourse")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
ylim([0 1])
yticks(0:0.5:1)
xticks(0:10:30)
legend box off
if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/Timecourse/RNF166.pdf","ContentType","vector")
end


% ZNF692 Fit
figure(6)
clf
hold on
for ii = 1:4
    P_cond = P_sim{ii+16};
    P_mean = mean(P_cond,2);
    plot(tspan,P_mean,"LineWidth",0.5, "Color",c(ii,:));
    fillspecs = {"Color",c(ii,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(tspan,P_cond, [0.05 0.95],fillspecs,'off');
end

for ii = 1:4
    data_cond = timecourse_data{ii+16};
    t_exp = data_cond{1};
    c_exp = data_cond{2};
    c_err = data_cond{3};

    errorbar(t_exp,c_exp,c_err,'.','MarkerSize',10,'Color',c(ii,:))
end

xlabel("Time (hours)")
ylabel("Normalized iMRM Signal")
l = legend("Lenalidomide (1{\mu}M)","Pomalidomide (1{\mu}M)","Avadomide (1{\mu}M)","cc885 (0.1{\mu}M)","Location","Best");
title(l,'Treatment')
title("ZNF692 Degradation Timecourse")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
ylim([0 1])
yticks(0:0.5:1)
xticks(0:10:30)
legend box off
if save_plots
   exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/Timecourse/ZNF692.pdf","ContentType","vector")
end


%% Simulating Dose Response Data


titration_data = model.experimental_data(ConditionsToOpt(21:40));
tspan = [0 6];
G_range = logspace(-4, 7,100);
P_sim = cell(length(titration_data),1);
EP_ind = find(strcmp(model.species_names,"EP"),1);
EGP_ind = find(strcmp(model.species_names,"EGP"),1);
P_ind = find(strcmp(model.species_names,"P"),1);

for jj = 1:size(ensemble,2)
    model = set_model_parameters(model, ensemble(:,jj),EnsembleParameters);
    for ii = 1:length(P_sim)
        condition = titration_data{ii};
        model.parameters.P_name = condition{6};
        model.parameters.G_name = condition{5};
        model.parameters.experiment_type = condition{8};
        for kk = 1:length(G_range)
            model.parameters.Gt = G_range(kk);
            ModelPredictions = model.evaluate(tspan);
            y_baseline = model.parameters.(append("Pt_",condition{6},"_MRM"));
            EP = ModelPredictions(end,EP_ind);
            EGP = ModelPredictions(end,EGP_ind);
            P = ModelPredictions(end,P_ind);
            P_sim{ii}(kk,jj) = (EP + EGP + P)./y_baseline;
        end
    end
end

%% Plotting Results

c = [ 0.8500 0.3250 0.0980; 0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
% IKZF1 Titration (1:4)
figure(7)
clf
hold on
for ii = 1:4
    P_cond = P_sim{ii};
    P_mean = mean(P_cond,2);
    plot(G_range,P_mean,"LineWidth",0.5, "Color",c(ii,:));
    fillspecs = {"Color",c(ii,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range,P_cond, [0.05 0.95],fillspecs,'off');
end
for ii = 1:4
    data_cond = titration_data{ii};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii,:));
end

xlabel('Dose (nM)')
ylabel('Normalized iMRM Signal')
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Avadomide","cc885","Location","Best");
title(l,'Treatment')
title("IKZF1 Dose Response")
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
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/Dose Response/IKZF1.pdf","ContentType","vector")
end


% IKZF3 Titration (5:8)
figure(8)
clf
hold on
for ii = 1:4
    P_cond = P_sim{ii + 4};
    P_mean = mean(P_cond,2);
    plot(G_range,P_mean,"LineWidth",0.5, "Color",c(ii,:));
    fillspecs = {"Color",c(ii,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range,P_cond, [0.05 0.95],fillspecs,'off');
end
for ii = 1:4
    data_cond = titration_data{ii + 4};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii,:));
end

xlabel('Dose (nM)')
ylabel('Normalized iMRM Signal')
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Avadomide","cc885","Location","Best");
title(l,'Treatment')
title("IKZF3 Dose Response")
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
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/Dose Response/IKZF3.pdf","ContentType","vector")
end
% GSPT1 Titration (9)
figure(9)
clf
hold on
P_cond = P_sim{9};
P_mean = mean(P_cond,2);
plot(G_range,P_mean,"LineWidth",0.5, "Color",c(4,:));
fillspecs = {"Color",c(4,:), "FaceAlpha",0.3,"EdgeAlpha",0};
plot_ensemble_CI(G_range,P_cond, [0.05 0.95],fillspecs,'off');
data_cond = titration_data{9};
G_exp = data_cond{4};
y_exp = data_cond{2};
y_err = data_cond{3};
errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(4,:));

xlabel('Dose (nM)')
ylabel('Normalized iMRM Signal')
set(gca,'Xscale','log')
l = legend("cc885","Location","Best");
title(l,'Treatment')
title("GSPT1 Dose Response")
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
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/Dose Response/GSPT1.pdf","ContentType","vector")
end

% ZFP91 Titration (10:13)
figure(10)
clf
hold on
for ii = 1:3
    P_cond = P_sim{ii + 9};
    P_mean = mean(P_cond,2);
    plot(G_range,P_mean,"LineWidth",0.5, "Color",c(ii,:));
    fillspecs = {"Color",c(ii,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range,P_cond, [0.05 0.95],fillspecs,'off');
end
for ii = 1:3
    data_cond = titration_data{ii + 9};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii,:));
end

xlabel('Dose (nM)')
ylabel('Normalized iMRM Signal')
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Avadomide","Location","Best");
title(l,'Treatment')
title("ZFP91 Dose Response")
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
   exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/Dose Response/ZFP91.pdf","ContentType","vector")
end
% RNF166 titration (14:17)
figure(11)
clf
hold on
for ii = 1:4
    P_cond = P_sim{ii + 12};
    P_mean = mean(P_cond,2);
    plot(G_range,P_mean,"LineWidth",0.5, "Color",c(ii,:));
    fillspecs = {"Color",c(ii,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range,P_cond, [0.05 0.95],fillspecs,'off');
end
for ii = 1:4
    data_cond = titration_data{ii + 12};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii,:));
end

xlabel('Dose (nM)')
ylabel('Normalized iMRM Signal')
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Avadomide","cc885","Location","Best");
title(l,'Treatment')
title("RNF166 Dose Response")
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
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/Dose Response/RNF166.pdf","ContentType","vector")
end

% ZNF692 titration (18:21)
figure(12)
clf
hold on
for ii = 1:4
    P_cond = P_sim{ii+16};
    P_mean = mean(P_cond,2);
    plot(G_range,P_mean,"LineWidth",0.5, "Color",c(ii,:));
    fillspecs = {"Color",c(ii,:), "FaceAlpha",0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range,P_cond, [0.05 0.95],fillspecs,'off');
end
for ii = 1:4
    data_cond = titration_data{ii+16};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii,:));
end

xlabel('Dose (nM)')
ylabel('Normalized iMRM Signal')
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Avadomide","cc885","Location","Best");
title(l,'Treatment')
title("ZNF692 Dose Response")
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
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/Dose Response/ZNF692.pdf","ContentType","vector")
end


%% Simulating BRET Experiments.


BRET_data = model.experimental_data(ConditionsToOpt(41:64));
tspan = [0 6];
G_range = logspace(-3,5,100);
P_sim = cell(length(BRET_data),1);
EP_ind = find(strcmp(model.species_names,"EP"),1);
EGP_ind = find(strcmp(model.species_names,"EGP"),1);
for jj = 1:size(ensemble,2)
    model = set_model_parameters(model, ensemble(:,jj),EnsembleParameters);
    for ii = 1:length(P_sim)
        condition = BRET_data{ii};
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

%% Plotting BRET Experiments
c = [ 0.8500 0.3250 0.0980; 0 0.4470 0.7410; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560];
% IKZF1 Titration (1:4)
figure(13)
clf
hold on
for ii = 1:4
    P_cond = P_sim{ii};
    P_mean = mean(P_cond,2);
    plot(G_range, P_mean, "LineWidth",0.5,"Color",c(ii,:));
    fillspecs = {"Color", c(ii,:), "FaceAlpha", 0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range, P_cond, [0.05 0.95], fillspecs,"off");
end
for ii = 1:4
    data_cond = BRET_data{ii};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii,:));
end

xlabel('Dose (nM)')
ylabel("Normalized BRET Signal")
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Avadomide","cc885","Location","Best");
legend box off
title(l,'Treatment')
title("IKZF1 BRET Response")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
xticks([10^-4 10^-2 10^0 10^2 10^4 10^6])
xlim([10^-4 10^6])
ylim([0 10])
yticks(0:5:15)
if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/BRET/IKZF1.pdf","ContentType","vector")
end
% IKZF3 Titration (5:8)
figure(14)
clf
hold on
for ii = 5:8
    P_cond = P_sim{ii};
    P_mean = mean(P_cond,2);
    plot(G_range, P_mean, "LineWidth",0.5,"Color",c(ii-4,:));
    fillspecs = {"Color", c(ii-4,:), "FaceAlpha", 0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range, P_cond, [0.05 0.95], fillspecs,"off");
end
for ii = 5:8
    data_cond = BRET_data{ii};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii-4,:));
end

xlabel('Dose (nM)')
ylabel("Normalized BRET Signal")
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Avadomide","cc885","Location","Best");
legend box off
title(l,'Treatment')
title("IKZF3 BRET Response")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
xticks([10^-4 10^-2 10^0 10^2 10^4 10^6])
xlim([10^-4 10^6])
ylim([0 10])
yticks(0:5:15)
if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/BRET/IKZF3.pdf","ContentType","vector")
end
% GSPT1 Titration (9:12)
figure(15)
clf
hold on
for ii = 9:12
    P_cond = P_sim{ii};
    P_mean = mean(P_cond,2);
    plot(G_range, P_mean, "LineWidth",0.5,"Color",c(ii-8,:));
    fillspecs = {"Color", c(ii-8,:), "FaceAlpha", 0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range, P_cond, [0.05 0.95], fillspecs,"off");
end
for ii = 9:12
    data_cond = BRET_data{ii};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii-8,:));
end

xlabel('Dose (nM)')
ylabel("Normalized BRET Signal")
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Avadomide","cc885","Location","Best");
legend box off
title(l,'Treatment')
title("GSPT1 BRET Response")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
xticks([10^-4 10^-2 10^0 10^2 10^4 10^6])
xlim([10^-4 10^6])
ylim([0 10])
yticks(0:5:15)
if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/BRET/GSPT1.pdf","ContentType","vector")
end
% ZFP91 Titration (13:16)
figure(16)
clf
hold on
for ii = 13:16
    P_cond = P_sim{ii};
    P_mean = mean(P_cond,2);
    plot(G_range, P_mean, "LineWidth",0.5,"Color",c(ii-12,:));
    fillspecs = {"Color", c(ii-12,:), "FaceAlpha", 0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range, P_cond, [0.05 0.95], fillspecs,"off");
end
for ii = 13:16
    data_cond = BRET_data{ii};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii-12,:));
end

xlabel('Dose (nM)')
ylabel("Normalized BRET Signal")
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Avadomide","cc885","Location","Best");
legend box off
title(l,'Treatment')
title("ZFP91 BRET Response")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
xticks([10^-4 10^-2 10^0 10^2 10^4 10^6])
xlim([10^-4 10^6])
ylim([0 10])
yticks(0:5:15)
if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/BRET/ZFP91.pdf","ContentType","vector")
end
%  RNF166 Titration (17:20)
figure(18)
clf
hold on
for ii = 17:20
    P_cond = P_sim{ii};
    P_mean = mean(P_cond,2);
    plot(G_range, P_mean, "LineWidth",0.5,"Color",c(ii-16,:));
    fillspecs = {"Color", c(ii-16,:), "FaceAlpha", 0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range, P_cond, [0.05 0.95], fillspecs,"off");
end
for ii = 17:20
    data_cond = BRET_data{ii};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii-16,:));
end

xlabel('Dose (nM)')
ylabel("Normalized BRET Signal")
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Avadomide","cc885","Location","Best");
legend box off
title(l,'Treatment')
title("RNF166 BRET Response")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
xticks([10^-4 10^-2 10^0 10^2 10^4 10^6])
xlim([10^-4 10^6])
ylim([0 10])
yticks(0:5:15)
if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/BRET/RNF166.pdf","ContentType","vector")
end

%  ZNF692 Titration (21:24)
figure(19)
clf
hold on
for ii = 21:24
    P_cond = P_sim{ii};
    P_mean = mean(P_cond,2);
    plot(G_range, P_mean, "LineWidth",0.5,"Color",c(ii-20,:));
    fillspecs = {"Color", c(ii-20,:), "FaceAlpha", 0.3,"EdgeAlpha",0};
    plot_ensemble_CI(G_range, P_cond, [0.05 0.95], fillspecs,"off");
end
for ii = 21:24
    data_cond = BRET_data{ii};
    G_exp = data_cond{4};
    y_exp = data_cond{2};
    y_err = data_cond{3};
    errorbar(G_exp,y_exp,y_err,'.','MarkerSize',10, "Color",c(ii-20,:));
end

xlabel('Dose (nM)')
ylabel("Normalized BRET Signal")
set(gca,'Xscale','log')
l = legend("Pomalidomide","Lenalidomide","Avadomide","cc885","Location","Best");
legend box off
title(l,'Treatment')
title("ZNF692 BRET Response")
l.IconColumnWidth = 10;
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10)
xticks([10^-4 10^-2 10^0 10^2 10^4 10^6])
xlim([10^-4 10^6])
ylim([0 15])
yticks(0:5:15)
if save_plots
    exportgraphics(gca,"Combined Degrader Fitting/Finalized Fitness Figures/Sperling et al/BRET/ZNF692.pdf","ContentType","vector")
end

