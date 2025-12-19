%% 14-3-3 Cooperative Binding Model Plotting (Combined: ERa, CRAF, TAZ)
clear 
close all
addpath(genpath("Models/1433model"))
addpath(genpath("1433 fitting"))

% 타겟별 설정 정보 정의
targets_info = {
    struct('name', "ERa", ...
           'file', "2025_11_04_1433_coop_ERa_EMfitting_PSn10000k50_SE_v1.mat", ...
           'drugs', ["ERa01", "ERa02"], ...
           'params', ["Et","Pt","KD1_ERa", "alpha_ERa01","alpha_ERa02","KD3_ERa01","KD3_ERa02"]', ...
           'ylim', [0.9, 2.7], 'yticks', 1:0.5:2.5, 'slice', []), ...
    struct('name', "CRAF", ...
           'file', "2025_11_04_1433_coop_CRAF_EMfitting_PSn10000k50_SE_v1.mat", ...
           'drugs', ["CRAF01", "CRAF02", "CRAF03", "CRAF04"], ...
           'params', ["Et","Pt","KD1_CRAF","alpha_CRAF01","alpha_CRAF02","alpha_CRAF03","alpha_CRAF04","KD3_CRAF01","KD3_CRAF02","KD3_CRAF03","KD3_CRAF04"]', ...
           'ylim', [0.8, 1.7], 'yticks', 0.8:0.2:1.6, 'slice', 23:50), ...
    struct('name', "TAZ", ...
           'file', "2025_11_04_1433_coop_TAZ_EMfitting_PSn10000k50_SE_v1.mat", ...
           'drugs', ["TAZ01", "TAZ02", "TAZ03"], ...
           'params', ["Et","Pt","KD1_TAZ","alpha_TAZ01","alpha_TAZ02","alpha_TAZ03","KD3_TAZ01","KD3_TAZ02","KD3_TAZ03"]', ...
           'ylim', [0.9, 2.7], 'yticks', 1:0.5:2.5, 'slice', [])
};

G_range = logspace(0, 5, 1000);
colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4941 0.1843 0.5569];

% 각 타겟에 대해 루프 실행
for t_idx = 1:length(targets_info)
    info = targets_info{t_idx};
    fprintf('Processing: %s\n', info.name);
    
    % 데이터 로드
    if ~exist(info.file, 'file')
        warning('File %s not found. Skipping %s.', info.file, info.name);
        continue;
    end
    load(info.file)
    
    model = outstruct.model;
    experimental_data = model.experimental_data;
    OptimizedEnsemble = outstruct.OptimizedEnsemble;
    
    % CRAF와 같이 특정 앙상블 범위만 사용하는 경우 처리
    if ~isempty(info.slice)
        OptimizedEnsemble = OptimizedEnsemble(:, info.slice);
    end
    
    Parameters = info.params;

    %% 통계 계산 (Calculating Stats)
    sd = std(OptimizedEnsemble, 0, 2);
    avg = mean(OptimizedEnsemble, 2);
    prc_025 = prctile(OptimizedEnsemble, 2.5, 2);
    prc_975 = prctile(OptimizedEnsemble, 97.5, 2);
    CV = sd ./ avg;
    stat_summary = table(Parameters, avg, sd, CV, prc_025, prc_975);
    
    % 개별 변수 할당 (필요 시)
    for i = 1:numel(Parameters)
        eval([Parameters{i} ' = OptimizedEnsemble(i,:);']);
    end

    %% 시뮬레이션 및 플로팅
    figure(t_idx)
    clf
    hold on
    drugs = info.drugs;
    
    for drug_i = 1:length(drugs)
        y_ensemble = zeros(length(G_range), size(OptimizedEnsemble, 2));
        
        for kk = 1:size(OptimizedEnsemble, 2)
            % 모델 파라미터 업데이트
            for p_idx = 1:length(Parameters)
                model.parameters.(Parameters(p_idx)) = OptimizedEnsemble(p_idx, kk);
            end

            % Baseline (Gt = 0) 계산
            model.parameters.Gt = 0;
            data_tmp = experimental_data{drug_i};
            model.parameters.G_name = data_tmp{4};
            model.parameters.P_name = data_tmp{5};
            
            BaselinePredictions = model.evaluate;
            y_baseline = BaselinePredictions(1) + BaselinePredictions(2); % EGP + EP
            
            % Dose-response 계산
            y_predicted = zeros(size(G_range));
            for jj = 1:length(G_range)
                model.parameters.Gt = G_range(jj);
                ModelPredictions = model.evaluate;
                y_predicted(jj) = (ModelPredictions(1) + ModelPredictions(2)) ./ y_baseline;
            end
            y_ensemble(:, kk) = y_predicted;
        end
        
        % 앙상블 신뢰구간 및 평균 플로팅
        fillspecs = {"Color", colors(drug_i, :), "FaceAlpha", 0.3, "EdgeAlpha", 0};
        plot_ensemble_CI(G_range, y_ensemble, [0.05 0.95], fillspecs, "off")
        plot(G_range, mean(y_ensemble, 2), "LineWidth", 1.5, "Color", colors(drug_i, :))
    end

    % 실험 데이터 (Errorbar) 추가
    for drug_i = 1:length(drugs)
        data_tmp = experimental_data{drug_i};
        markerSize = 15;
        if info.name == "CRAF", markerSize = 10; end % CRAF 파일 기준
        errorbar(data_tmp{1}, data_tmp{2}, data_tmp{3}, '.', 'MarkerSize', markerSize, ...
                 "Color", colors(drug_i, :), "HandleVisibility", "off")
    end

    %% 그래프 서식 설정
    set(gca, 'Xscale', "log")
    title(info.name)
    l = legend(drugs, "Location", "NorthWest");
    legend box off
    title(l, 'Drug')
    xlabel("Dose (nM)")
    ylabel("Relative BRET Signal")
    ylim(info.ylim)
    yticks(info.yticks)
    
    % 공통 서식
    set(gca, 'FontSize', 14, 'TickDir', 'out', 'LineWidth', 1.2)
    
    % CRAF 전용 추가 서식 및 파일 저장
    if info.name == "CRAF"
        set(gca, "TickLength", [0.025 0.01])
        set(gcf, "Position", [0, 0, 250, 187.5])
        set(gca, "FontSize", 10) 
        xticks(10.^(0:1:5))
    end
end