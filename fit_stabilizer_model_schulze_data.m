%% Unified Stabilizer EM Fitting (Cmpd3, Cmpd4, RMC-4998)
clear; clc;

% 경로 추가
addpath(genpath("Models/Stabilizer Model"))

% 기본 모델 로드
load("Stabilizer_model_v4")

% 화합물 리스트 정의
compounds = {'Cmpd3', 'Cmpd4', 'RMC'};

for i = 1:length(compounds)
    cpd = compounds{i};
    fprintf('--- Running Optimization for: %s ---\n', cpd);
    
    % 모델 초기화 및 공통 바운드 설정
    model.bounds.kf1_bounds = [1e-7, 1];
    model.bounds.kf2_bounds = [1e-7, 1];
    
    % 화합물별 고유 설정
    switch cpd
        case 'Cmpd3'
            savename = "CYPA stabilizer fitting/2025_06_06_stabilizer_EMfitting_1_10000_kf_fixed_se_cmpd3";
            model.bounds.kcat_bounds = [1e-6, 2e-1];
            model.parameters.kf1 = 1e-4;
            model.parameters.kf2 = 1e-4;
            ParamsToOpt = ["kb1", "kb2", "kcat"];
            n_sets = 10000;
            
            % Experimental Data (TR-FRET)
            G_fret = 10.^[0.0516519691853996; 0.5299672648244; 1.0093481754358; 1.49142830352522; ...
                           1.9643326809657; 2.43682000462306; 2.9180438957618; 3.39294842910096];
            y_fret = [0.0584034459858239; 3.14307878480376; 12.8296096953979; 30.1193034871952; ...
                      67.5455320641382; 96.1731065971635; 105.825219018626; 99.9454901170799];
            % Experimental Data (Covalent)
            G_cov = 10^3*[0.107337438413987; 0.207531798509672; 0.398201973051095; 0.794816640436933; ...
                          1.59663924973271; 3.20579724215275];
            y_cov = [0.000271754910324043; 0.000257813893969818; 0.0000571349671763937; ...
                     0.000275804343213621; 0.000559264645484102; 0.000317648483072597];

        case 'Cmpd4'
            savename = "CYPA stabilizer fitting/2025_06_06_stabilizer_EMfitting_1_10000_kf_fixed_se_cmpd4";
            model.bounds.kcat_bounds = [1e-6, 2e-1];
            model.parameters.kf1 = 1e-4;
            model.parameters.kf2 = 1e-4;
            ParamsToOpt = ["kb1", "kb2", "kcat"];
            n_sets = 10000;
            
            % Experimental Data (TR-FRET)
            G_fret = 10.^[0.0527817026305809; 0.535634089493926; 1.00356644819337; 1.49014167841876; ...
                           1.96244373940919; 2.4363411020476; 2.91561420283751; 3.39109296533113];
            y_fret = [5.45891230121821; 11.7448113139034; 21.8997713216114; 53.288482832847; ...
                      96.4654736393993; 102.781049663538; 105.86409381236; 99.8478810765915];
            % Experimental Data (Covalent)
            G_cov = 10^3*[0.0949097639918919; 0.20407120443088; 0.394674852489017; 0.7970730646033; ...
                          1.59721750279323; 3.2019979309182];
            y_cov = [0.000919432173897271; 0.00140137796014599; 0.00233508860313828; ...
                     0.00409153011899287; 0.0057638193576109; 0.0077221630664443];

        case 'RMC'
            savename = "CYPA stabilizer fitting/2025_06_06_stabilizer_EMfitting_1_10000_kf_fixed_se_rmc";
            model.parameters.kf1 = 1e-4;
            model.parameters.kf2 = 1e-4;
            ParamsToOpt = ["kf1", "kf2", "kcat"]; % rmc2 파일에서 수정된 파라미터 기준
            n_sets = 5000; % rmc2 파일 기준
            
            % Experimental Data (TR-FRET)
            G_fret = 10.^[0.0568078576786148; 0.543218288805079; 1.01818104987376; 1.5026344563116; ...
                           1.98770521028111; 2.45395907943249; 2.93197556492251; 3.41742861217634];
            y_fret = [0.0784571124218646; 4.83103470611712; 20.078939381468; 62.707704509654; ...
                      81.9955140493035; 91.3619461719859; 95.3871071805477; 98.3390339373152];
            % Experimental Data (Covalent)
            G_cov = 10^3*[0.0507567808819845; 0.0968486292819558; 0.19531893188245; 0.400115382378208; ...
                          0.789948701542757; 1.58944790958227; 3.19352527562015];
            y_cov = [0.0111393251534454; 0.0173697962611312; 0.0363905528910886; 0.0670696012043494; ...
                     0.0861140797800228; 0.121225453952439; 0.134412064633379];
    end
    
    % 공통 단백질 농도 설정
    Pt_FRET = 50;
    Et_FRET = 20000;
    
    % 모델에 데이터 삽입
    model.experimental_data{1} = {G_fret, y_fret, 0, cpd, "KRAS", "TRFRET", Pt_FRET, Et_FRET};
    model.experimental_data{2} = {G_cov, y_cov, 0, cpd, "KRAS", "covalent", Pt_FRET, 25000};
    
    % 최적화 옵션 설정 및 실행
    FitFunction = @calculate_stabilizer_fitness_v2;
    ErrorFunction = @se;
    
    tic
    numberOfVariables = length(ParamsToOpt);
    OptOptions = optimoptions('patternsearch', 'UseParallel', true);
    
    % rmc의 경우 MaxIterations 옵션 추가 (rmc2 파일 기준)
    if strcmp(cpd, 'RMC')
        OptOptions.MaxIterations = 500 * numberOfVariables;
        OptOptions.MaxFunctionEvaluations = 5000 * numberOfVariables;
    end
    
    outstruct = fit_dose_response_data(model, "ParamsToOpt", ParamsToOpt, "FitFunction", FitFunction, ...
        "ErrorFunction", ErrorFunction, "OptScale", "log", "OptMethod", "patternsearch", 'tspan', [0, 4*3600], ...
        "NSets", n_sets, "BestNSets", 50, 'OptOptions', OptOptions, 'SaveName', savename);
    toc
    
    clear G_fret y_fret G_cov y_cov savename ParamsToOpt n_sets;
end