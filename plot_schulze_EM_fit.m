%% 파일 2: 통합 결과 플로팅 (TR-FRET & Covalent Data)
clear; close all; clc;

% 컬러 설정 (Cmpd3: Blue, Cmpd4: Orange, RMC: Yellow)
c = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
G_range_log = logspace(0, 3.4, 100);
G_range_lin = linspace(0, 3500, 120);

% 데이터 로드 리스트
files = {
    'cmpd3_EMfit_se_250530', ...
    'cmpd4_EMfit_se_250530', ...
    'rmc_EMfit_se_250530'
};
labels = {"Cmpd3", "Cmpd4", "RMC4998"};

%% Figure 1: TR-FRET Dose-Response
figure(1); clf; hold on;

for i = 1:length(files)
    % 메인 결과 로드
    load(fullfile("CYPA stabilizer fitting", files{i} + "_v2.mat"));
    
    % 시뮬레이션 평균선 및 신뢰구간
    plot(G_range, mean(P_sim, 2), "LineWidth", 1.5, "Color", c(i,:));
    fillspecs = {"Color", c(i,:), "FaceAlpha", 0.1, "EdgeAlpha", 0};
    plot_ensemble_CI(G_range, P_sim, [0.05 0.95], fillspecs, 'off');
    
    % 실험 데이터 (Errorbar)
    errorbar(model.experimental_data{1}{1}, model.experimental_data{1}{2}, 0, ...
        '.', 'MarkerSize', 12, 'Color', c(i,:));
end

set(gca, "Xscale", 'log', "FontSize", 10, "TickDir", "out", "TickLength", [0.025 0.01]);
xlabel("Dose (nM)"); ylabel("KRAS:CYPA (%)"); ylim([0, 110]);
xticks([10.^(0:1:5)]);
legend(labels, "Location", "SouthEast", "Title", "Drug");
title("TR-FRET Combined Data");
set(gcf, "Position", [0, 0, 350, 260]);

%% Figure 2: Covalent k_obs Data
figure(2); clf; hold on;

for i = 1:length(files)
    % Covalent 결과 로드
    load(fullfile("CYPA stabilizer fitting", files{i} + "_cov_v2.mat"));
    
    % 시뮬레이션 평균선 및 신뢰구간
    plot(G_range, mean(rates, 2), "LineWidth", 1.5, "Color", c(i,:));
    fillspecs = {"Color", c(i,:), "FaceAlpha", 0.1, "EdgeAlpha", 0};
    plot_ensemble_CI(G_range, rates, [0.05 0.95], fillspecs, 'off');
    
    % 실험 데이터 (Scatter)
    scatter(model.experimental_data{2}{1}, model.experimental_data{2}{2}, 15, ...
        'MarkerFaceColor', c(i,:), 'MarkerEdgeAlpha', 0);
end

set(gca, "FontSize", 10, "TickDir", "out", "TickLength", [0.025 0.01]);
xlabel("Dose (nM)"); ylabel("k_{obs} (s^{-1})"); ylim([0, 0.15]); xlim([0, 3500]);
legend(labels, "Location", "NorthWest", "Title", "Drug");
title("Covalent Combined Data");
set(gcf, "Position", [400, 0, 350, 260]);
