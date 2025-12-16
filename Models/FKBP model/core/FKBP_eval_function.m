function y = FKBP_eval_function(obj)
% Function for evaluating complex formation model for multiple datasets.
% Solves equilibrium of binding using fsolve and reports all species values
% for different KD values.

    % Getting drug and target info.
    P_name = obj.parameters.P_name;
    G_name = obj.parameters.G_name;
    
    % Setting KD3 based on P name
    % obj.parameters.KD1 = obj.parameters.(append("KD1_",P_name));
    % Setting KD2 based on G name and P name
    obj.parameters.KD4 = obj.parameters.(append("KD4_",G_name,"_",P_name));
    % Setting KD3 based on G name
    obj.parameters.KD3 = obj.parameters.(append("KD3_",G_name));
    % solving conservation equations
    [EGP, EG,G, E, P] = solve_conservation_equations(obj.parameters);
    % storing outputs.
    y = [EGP, EG, G, E, P];
end


function sol = calculate_conservation_equations(E, parameters)
% Function calculating the equations that need to be solved for the model.

    % Intializing parameters.

    KD3 = parameters.KD3;
    KD4 = parameters.KD4;
    Gt = parameters.Gt;
    Pt = parameters.Pt;
    Et = parameters.Et;

    % Calculating constants
    % A1 = (1 + E/KD3)*(E/(KD1*KD2));
    % A2 = (1 + E/KD1)*(E/(KD1*KD2));
    % B1 = (1 + E/KD1)*(1 + E/KD3) + (Pt - Gt)*(E/(KD1*KD2));
    % B2 = (1 + E/KD1)*(1 + E/KD3) + (Gt - Pt)*(E/(KD1*KD2));
    % C1 = -(1 + E/KD1)*Gt;
    % C2 = -(1 + E/KD3)*Pt;

    % Calculating G and P
    % G = (-B1 + sqrt(B1^2 - 4*A1*C1))/(2*A1);
    % P = (-B2 + sqrt(B2^2 - 4*A2*C2))/(2*A2);
    G = Gt-Et+E;
    P = Pt/(1+E*G/(KD3*KD4));
    % Calculating conservation equation.
    sol = E*(1 + G/KD3 + (G*P)/(KD3*KD4)) - Et;
end

function  [EGP, EG, G, E, P] = solve_conservation_equations(parameters)
% Function for simulating equilibrium model of MG binding kinetics.
    
    % Initializing other parameters.
    KD3 = parameters.KD3;
    KD4 = parameters.KD4;
    Gt = parameters.Gt;
    Pt = parameters.Pt;
    Et = parameters.Et;
    
    % Solving for G and P based on conservation equations.
    func = @(y) calculate_conservation_equations(y,parameters);
    y0 = [eps Et-eps];
    options = optimset('Display','off');
    E = fzero(func,y0,options);

    % Calculating other species
    G = Gt-Et+E;
    P = Pt/(1+E*G/(KD3*KD4));
    % Calculating conservation equation.
   

    % calculating EGP, EP, and EG
    EGP = E*P*G/(KD3*KD4);
    EG = E*G/KD3;
end
