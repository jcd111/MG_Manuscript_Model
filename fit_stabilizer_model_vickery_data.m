%% 14-3-3 Cooperative Binding Model Fitting (Combined: ERa, CRAF, TAZ)
clear 
addpath(genpath("../Modeling Functions"))
addpath(genpath("../Models/1433model"))

% 기본 모델 로드
load("1433_coop_model_v1.mat")

% 공통 모델 설정
model.bounds.Et_bounds = [1, 1000];
model.bounds.Pt_bounds = [1, 1000];

%% 타겟별 설정 정의 (Configuration)
% 실행하고자 하는 타겟을 리스트에 넣으세요 (예: ["ERa", "CRAF", "TAZ"])
run_targets = ["ERa", "CRAF", "TAZ"]; 

for t_name = run_targets
    fprintf('--- Starting Optimization for Target: %s ---\n', t_name);
    
    % 현재 타겟에 맞게 모델 초기화 (실험 데이터 초기화)
    current_model = model;
    current_model.experimental_data = {};
    
    switch t_name
        case "ERa"
            savename = "2025_11_04_1433_coop_ERa_EMfitting_PSn10000k50_SE_v1";
            n_sets = 300;
            params = ["Et","Pt","KD1_ERa","alpha_ERa01","alpha_ERa02","KD3_ERa01","KD3_ERa02"];
            
            % ERa01 Data
            G1 = 1e9*10.^[-6.8068754; -6.5044557; -6.20412; -5.90309; -5.60206; -5.30103; -5; -4.69897];
            y1 = [1.29659197 0.23589118; 1.39686615 0.15809879; 1.50050211 0.15406787; 1.81675673 0.2336021; 1.871302 0.08733602; 2.12938332 0.29150887; 2.51930792 0.07115281; 2.43631939 0.10110964];
            current_model.experimental_data{1} = {G1, y1(:,1), y1(:,2), "ERa01", "ERa"};
            
            % ERa02 Data
            G2 = G1; % Same Dose
            y2 = [1.06205505 0.02332213; 1.06128884 0.15046826; 1.05384991 0.02499402; 1.09265226 0.07816389; 1.11759336 0.03332713; 1.00954238 0.08551722; 1.04784427 0.10451696; 1.08907304 0.09060318];
            current_model.experimental_data{2} = {G2, y2(:,1), y2(:,2), "ERa02", "ERa"};

        case "CRAF"
            savename = "2025_11_04_1433_coop_CRAF_EMfitting_PSn10000k50_SE_v1";
            n_sets = 50;
            params = ["Et","Pt","KD1_CRAF","alpha_CRAF01","alpha_CRAF02","alpha_CRAF03","alpha_CRAF04","KD3_CRAF01","KD3_CRAF02","KD3_CRAF03","KD3_CRAF04"];
            
            G_common = 1e9*10.^[-8.3113302; -8.0102997; -7.70927; -7.40824; -7.10721; -6.80618; -6.50515; -6.20412; -5.90309; -5.60206; -5.30103; -5.0];
            y_c01 = [0.96548867 0.01562349; 0.99680650 0.03662492; 0.95637715 0.04297132; 0.94334910 0.01616446; 1.01509703 0.06656667; 1.07820591 0.04576315; 1.01729560 0.02312441; 0.98986465 0.05314353; 0.98760409 0.03713082; 1.01129388 0.03225764; 1.03848977 0.02544161; 1.00012766 0.03685381];
            y_c02 = [0.98729581 0.02164551; 0.95880601 0.01723791; 0.96789082 0.01851238; 1.06207009 0.03864189; 1.29722492 0.07099525; 1.29087127 0.16940956; 1.34292857 0.07192083; 1.45195501 0.06809085; 1.47186095 0.01756381; 1.50679987 0.09109365; 1.60426760 0.06633227; 1.60961507 0.04138530];
            y_c03 = [0.99007295 0.05907445; 0.98870139 0.08583928; 0.98321207 0.06416677; 1.05539251 0.06136794; 1.20360949 0.06063322; 1.22220600 0.04628974; 1.31026381 0.03589700; 1.26880993 0.08104530; 1.39856462 0.11709923; 1.43488366 0.07365801; 1.44965606 0.05360979; 1.47387388 0.05748559];
            y_c04 = [1.00000000 0.11623242; 1.07690035 0.07184782; 1.09612251 0.03857005; 1.00081234 0.06026594; 1.04372669 0.10703783; 1.10612625 0.04859359; 1.16026746 0.06160273; 1.13023765 0.02931828; 1.27237337 0.07151398; 1.28970730 0.08134349; 1.36696747 0.03559749; 1.33912948 0.06400713];
            
            current_model.experimental_data{1} = {G_common, y_c01(:,1), y_c01(:,2), "CRAF01", "CRAF"};
            current_model.experimental_data{2} = {G_common, y_c02(:,1), y_c02(:,2), "CRAF02", "CRAF"};
            current_model.experimental_data{3} = {G_common, y_c03(:,1), y_c03(:,2), "CRAF03", "CRAF"};
            current_model.experimental_data{4} = {G_common, y_c04(:,1), y_c04(:,2), "CRAF04", "CRAF"};

        case "TAZ"
            savename = "2025_11_04_1433_coop_TAZ_EMfitting_PSn10000k50_SE_v1";
            n_sets = 300;
            params = ["Et","Pt","KD1_TAZ","alpha_TAZ01","alpha_TAZ02","alpha_TAZ03","KD3_TAZ01","KD3_TAZ02","KD3_TAZ03"];
            
            G_taz = 1e9*10.^[-8.5; -8.0; -7.69897; -7.4089354; -7.1079054; -6.806041; -6.50515; -6.20412; -5.90309; -5.60206; -5.30103; -5.0];
            y_t01 = [1.0 0.178; 1.162 0.080; 0.974 0.266; 0.828 0.268; 1.257 0.168; 1.062 0.078; 1.039 0.087; 1.217 0.128; 1.093 0.229; 0.939 0.093; 1.006 0.172; 1.070 0.105]; % Simplified for brevity
            y_t02 = [1.0 0.178; 1.085 0.280; 1.011 0.169; 1.048 0.170; 1.193 0.085; 1.248 0.085; 1.363 0.037; 1.528 0.093; 1.923 0.113; 1.999 0.376; 2.130 0.126; 2.186 0.097];
            y_t03 = [1.0 0.178; 1.260 0.134; 1.249 0.163; 1.239 0.228; 1.266 0.179; 1.522 0.227; 1.521 0.104; 1.531 0.273; 1.580 0.258; 1.679 0.242; 1.586 0.192; 1.824 0.059];
            
            current_model.experimental_data{1} = {G_taz, y_t01(:,1), y_t01(:,2), "TAZ01", "TAZ"};
            current_model.experimental_data{2} = {G_taz, y_t02(:,1), y_t02(:,2), "TAZ02", "TAZ"};
            current_model.experimental_data{3} = {G_taz, y_t03(:,1), y_t03(:,2), "TAZ03", "TAZ"};
    end

    %% Optimization 실행
    tic
    OptOptions = optimoptions("particleswarm", 'UseParallel', true, 'MaxIterations', 1000 * length(params));

    outstruct = fit_BRET_data(current_model, ...
        "ParamsToOpt", params, ...
        "FitFunction", @calculate_BRET_fitness_everycombo, ...
        "ErrorFunction", @se, ...
        "OptScale", "log", ...
        "OptMethod", "particleswarm_ensemble", ...
        "NSets", n_sets, ...
        "Npools", 1, ...
        "OptOptions", OptOptions, ...
        "savename", savename, ...
        "UseParallel", false);
    toc
    
    fprintf('--- Finished Optimization for %s ---\n\n', t_name);
end