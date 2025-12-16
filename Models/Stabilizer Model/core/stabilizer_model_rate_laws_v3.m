function dydt = stabilizer_model_rate_laws_v3(t,y,p)
% Function for calcualting rates laws for simple reversible binding +
% degradation model.


    % Initalizing species
    E = y(1);               % Free target protein, nM
    G = y(2);               % Free effector protein, nM
    P = y(3);               % Free molecular glue, nM
    EG = y(4);              % E:G complex, nM
    EGP = y(5);             % E:G:P ternary complex, nM
    EGP_cov=y(6);
    % Pub = y(7);
    % Initializing Parameters
    % Association/dissociation parameters
    % k1 = p.k1;          % E + P -> EP
    % KD1 = p.KD1;
    % k2 = p.k2;          % EP + G -> EGP
    % KD2 = p.KD2;
    % k3 = p.k3;          % E + G -> EG
    % KD3 = p.KD3;
    % k4 = p.k4;          % EG + P -> EGP
    % 
    % Degradation parameters
    
    kf1=p.kf1;
    kb1=kf1*p.KD1;

    kf2=p.kf2;
    kb2=kf2*p.KD2;

    kcat=p.kcat;


    % Calculating rates of change.

    dydt = zeros(size(y));

    % E

    dydt(1) = -kf1*E*G+kb1*EG;
    % G 
    dydt(2) = -kf1*E*G+kb1*EG;
    % P
    dydt(3) = -kf2*EG*P+kb2*EGP;

    % EG
    dydt(4) = +kf1*E*G-kb1*EG-kf2*EG*P+kb2*EGP;

    % EGP
    dydt(5) = +kf2*EG*P-kb2*EGP-kcat*EGP;

    % EGP_cov
    dydt(6) = +kcat*EGP;


end