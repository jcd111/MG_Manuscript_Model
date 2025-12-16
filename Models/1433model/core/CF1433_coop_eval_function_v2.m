function y = CF1433_coop_eval_function_v2(obj)
% Function for evaluating complex formation model for multiple datasets.
% Solves equilibrium of binding using fsolve and reports all species values
% for different KD values.

    % Getting drug and target info.
    P_name = obj.parameters.P_name;
    G_name = obj.parameters.G_name;
    
    % Setting KD3 based on P name
    obj.parameters.KD1 = obj.parameters.(append("KD1_",P_name));
    % Setting KD2 based on G name and P name
    obj.parameters.alpha = obj.parameters.(append("alpha_",G_name));
    % Setting KD3 based on G name
    obj.parameters.KD3 = obj.parameters.(append("KD3_",G_name));
    % solving conservation equations
    [EGP, EP, EG,G, E, P] = solve_conservation_equations(obj.parameters);
    % storing outputs.
    y = [EGP, EP, EG, G, E, P];
end


function sol = calculate_conservation_equations(E, parameters)
% Function calculating the equations that need to be solved for the model.

    % Intializing parameters.
    KD1 = parameters.KD1;
    KD2 = parameters.KD3/parameters.alpha;
    KD3 = parameters.KD3;
    Gt = parameters.Gt;
    Pt = parameters.Pt;
    Et = parameters.EP_ratio*parameters.Pt;

    % Calculating constants
    A1 = (1 + E/KD3)*(E/(KD1*KD2));
    A2 = (1 + E/KD1)*(E/(KD1*KD2));
    B1 = (1 + E/KD1)*(1 + E/KD3) + (Pt - Gt)*(E/(KD1*KD2));
    B2 = (1 + E/KD1)*(1 + E/KD3) + (Gt - Pt)*(E/(KD1*KD2));
    C1 = -(1 + E/KD1)*Gt;
    C2 = -(1 + E/KD3)*Pt;

    % Calculating G and P
    G = (-B1 + sqrt(B1^2 - 4*A1*C1))/(2*A1);
    P = (-B2 + sqrt(B2^2 - 4*A2*C2))/(2*A2);

    % Calculating conservation equation.
    sol = E*(1 + P/KD1 + G/KD3 + (G*P)/(KD1*KD2)) - Et;
end

function  [EGP, EP, EG, G, E, P] = solve_conservation_equations(parameters)
% Function for simulating equilibrium model of MG binding kinetics.
    
    % Initializing other parameters.
    KD1 = parameters.KD1;
    KD2 = parameters.KD3/parameters.alpha;
    KD3 = parameters.KD3;
    Gt = parameters.Gt;
    Et = parameters.EP_ratio*parameters.Pt;
    Pt = parameters.Pt;
    
    % Solving for G and P based on conservation equations.
    func = @(y) calculate_conservation_equations(y,parameters);
    y0 = [eps Et-eps];
    options = optimset('Display','off');
    E = fzero(func,y0,options);
    % Calculating other species
    A1 = (1 + E/KD3)*(E/(KD1*KD2));
    A2 = (1 + E/KD1)*(E/(KD1*KD2));
    B1 = (1 + E/KD1)*(1 + E/KD3) + (Pt - Gt)*(E/(KD1*KD2));
    B2 = (1 + E/KD1)*(1 + E/KD3) + (Gt-Pt)*(E/(KD1*KD2));
    C1 = -(1 + E/KD1)*Gt;
    C2 = -(1 + E/KD3)*Pt;

    % Calculating G and P
    G = (-B1 + sqrt(B1^2 - 4*A1*C1))/(2*A1);
    P = (-B2 + sqrt(B2^2 - 4*A2*C2))/(2*A2);

    % calculating EGP, EP, and EG
    EP = E*P/KD1;
    EGP = E*P*G/(KD1*KD2);
    EG = E*G/KD3;
end
