% Loading Model
clear 
addpath(genpath("../Modeling Functions"))
addpath(genpath("../Models/FKBP model"))
addpath(genpath("../FKBP fitting_v2"))
load("2025_11_04_FKBP_cmpd7_KD3_fixed.mat")
model = outstruct.model;
experimental_data = model.experimental_data;

%%

OptimizedEnsemble = outstruct.OptimizedEnsemble;

Parameters = ["KD2_Cmpd5_FRB","KD2_Cmpd6_FRB","KD2_Cmpd7_FRB"]';

%% Calculating Stats
OptimizedEnsemble=OptimizedEnsemble(:,1:27);
sd = std(OptimizedEnsemble,0,2);
avg = mean(OptimizedEnsemble,2);
prc_025=prctile(OptimizedEnsemble,2.5,2);
prc_975=prctile(OptimizedEnsemble,97.5,2);
CV = sd./avg;

stat_summary = table(Parameters,avg,sd,CV,prc_025,prc_975);
stat_summary
%% Calculating some ratios.

numParams = numel(Parameters);

% Dynamically assign each row of OptimizedEnsemble to a variable named after the corresponding parameter
for i = 1:numParams
    eval([Parameters{i} ' = OptimizedEnsemble(i,:);']);
end


%%  Plotting IKZF1 degradation data

G_range = logspace(-3,6.3,1000);


c = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4941 0.1843 0.5569];
drugs = ["Cmpd5","Cmpd6","Cmpd7"];
targets = ["FRB"];

model.parameters.Et=20;
model.parameters.Gt=5000;
% Plotting IKZF1 fits iterating through each parameter set in ensemble.
for target_i=1:length(targets)
    figure(target_i)
    clf
    hold on
    for drug_i=1:length(drugs)
        y_ensemble = zeros(length(G_range),size(OptimizedEnsemble,2));
        for kk = 1:size(OptimizedEnsemble,2)
            for param = 1:length(Parameters)
                model.parameters.(Parameters(param)) = OptimizedEnsemble(param,kk);
            end


            % model.parameters.Gt = 0;
            data_tmp = experimental_data{drug_i};
            % G_exp = data_tmp{1};
            model.parameters.G_name = data_tmp{4};
            model.parameters.P_name = data_tmp{5};
            % BaselinePredictions = model.evaluate;
            % EP_Baseline = BaselinePredictions(2);
            % EGP_Baseline = BaselinePredictions(1);
            % y_baseline = EP_Baseline + EGP_Baseline;
            y_baseline=model.parameters.Et;
            y_predicted = zeros(size(G_range));
            for jj = 1:length(G_range)
                model.parameters.Pt = G_range(jj);
                ModelPredictions = model.evaluate;
                % EP = ModelPredictions(2);
                EGP = ModelPredictions(1);
                y_predicted(jj) = (EGP)./y_baseline;
            end
            y_ensemble(:,kk) = y_predicted;
        end
        fillspecs = {"Color",c(drug_i,:),"FaceAlpha",0.3,"EdgeAlpha",0};
        plot_ensemble_CI(G_range,y_ensemble,[0.05 0.95],fillspecs,"off")
        plot(G_range,mean(y_ensemble,2),"LineWidth",1.5,"Color",c(drug_i,:))
    end

    for drug_i=1:3
        data_tmp = experimental_data{drug_i};
        G_exp = data_tmp{1};
        y_exp = data_tmp{2};
        y_err = data_tmp{3};
        errorbar(G_exp, y_exp,y_err,'.','MarkerSize',10,"Color",c(drug_i,:),"HandleVisibility","off")
    end
    set(gca,'Xscale',"log")
    title(targets(target_i))
    l = legend(drugs,"Location","NorthWest");
    legend box off
    title(l,'Drug')
    xlabel("Dose (nM)")
    ylabel("Normalized BRET Signal = EGP/E_T")

end

%% Estimating the CV values of the data 
ylim([0,1])
yticks([0.:0.5:1])
xlim([10^2,10^6])
set(gca,'FontSize',10)
set(gca,'TickDir','out')
set(gca,'LineWidth',1)
set(gcf,'Position',[0,0,250,187.5])
set(gca,'TickLength',[0.05, 0.01])

% set(gcf,'FontSize',12)
set(gca,"TickDir","out")
set(gca,"TickLength",[0.025 0.01])
set(gcf,"Position",[0,0,250,187.5])
set(gca,"FontSize",10) 
xticks([10.^(0:1:5)])

% exportgraphics(gcf, 'FKBP.pdf', 'ContentType', 'vector');