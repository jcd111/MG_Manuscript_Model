%% Loading Model
clear
close all
addpath(genpath("Models/FKBP model"))
load("FKBPmodel_v2.mat")

savename = "2025_11_04_FKBP_cmpd7_KD3_fixed";

%% Setting Up Experimental Data



% Lenalidomide and IKZF1
% Glue Dose
G = 1000.*[1.002345428
2.011243132
3.98589108
8.034069737
16.17673696
32.11750262
63.97069708
128.0909823
];
% AlphaScreen Measurements
y = [191.2343755
195.0801815
192.2176286
194.0111753
193.320064
194.4770033
192.0766448
199.2195827];

% Converting to mean and sd.
y_mean = (y(:,1)-190)/80;
% y_sd = y(:,2);
y_sd=zeros(size(y_mean));
model.experimental_data{1} = {G,y_mean,y_sd,"Cmpd5","FRB"};

% Len and IKZF3
% Glue Dose
G = 1000.*[1.002345428
2.011243132
3.98589108
8.034069737
16.17673696
32.11750262
63.97069708
128.0909823
];

% AlphaScreen Measurements
y = [191.0519686
194.6092665
193.2219567
195.7771067
193.6107522
196.94204
196.0576209
202.2768967];
% Converting to mean and sd.
y_mean = (y(:,1)-190)/80;
% y_sd = y(:,2);
y_sd=zeros(size(y_mean));

model.experimental_data{2} = {G,y_mean,y_sd,"Cmpd6","FRB"};


% Lenalidomide and GSPT1
% Glue Dose
G = 1000.*[1.002345428
2.011243132
3.98589108
8.034069737
16.17673696
32.11750262
63.97069708
128.0909823
];

% AlphaScreen Measurements
% AlphaScreen Measurements
y = [190.9509544
194.8723394
197.3061271
202.1867834
206.8072737
217.0496754
221.4572365
237.544653];
% Converting to mean and sd.
y_mean = (y(:,1)-190)/80;
% y_sd = y(:,2);
y_sd=zeros(size(y_mean));

model.experimental_data{3} =  {G,y_mean,y_sd,"Cmpd7","FRB"};





%% Running optimization.

model.parameters.Et=20;
model.parameters.Gt=5000;
model.parameters.KD1_Cmpd5=5.5;
model.parameters.KD1_Cmpd6=86;
model.parameters.KD1_Cmpd7=6.2;

ParamsToOpt = ["KD2_Cmpd5_FRB","KD2_Cmpd6_FRB","KD2_Cmpd7_FRB"];
% ParamsToOpt=["KD1_Cmpd5","KD1_Cmpd6","KD2_Cmpd5_FRB","KD2_Cmpd6_FRB","KD2_Cmpd7_FRB"];


FitFunction = @calculate_FKBP_fitness;
ErrorFunction = @se;

numberOfVariables=length(ParamsToOpt);


% OptOptions = optimoptions('patternsearch','UseParallel',true,'MaxFunctionEvaluations',5000*length(ParamsToOpt),'MaxIterations',500*length(ParamsToOpt));
% OptOptions = optimoptions('patternsearch','UseParallel',true,MaxIterations=1000*numberOfVariables,MaxFunctionEvaluations=10000*numberOfVariables);

tic
OptOptions = optimoptions("particleswarm",'UseParallel',true,'MaxIterations',3000*length(ParamsToOpt));

outstruct = fit_BRET_data(model,"ParamsToOpt",ParamsToOpt,"FitFunction",FitFunction,...
    "ErrorFunction",@se,"OptScale","log","OptMethod","particleswarm_ensemble", ...
    "NSets",50, "Npools",1,'OptOptions',OptOptions,"savename",savename,"UseParallel",false);


toc
