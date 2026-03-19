function dydt = degrader_model_rate_laws(t,y,p)
% Function for calcualting rates laws for simple reversible binding +
% degradation model.


% Initalizing species
P = y(1);               % Free target protein, nM
E = y(2);               % Free effector protein, nM
G = y(3);               % Free molecular glue, nM
EP = y(4);              % E:P complex, nM
EG = y(5);              % E:G complex, nM
EGP = y(6);             % E:G:P ternary complex, nM

% Initializing Parameters
% Association/dissociation parameters
k3 = p.k3;          % E + P -> EP
KD3 = p.KD3;
k4 = p.k4;          % EP + G -> EGP
k1 = p.k1;          % E + G -> EG
KD1 = p.KD1;
k2 = p.k2;          % EG + P -> EGP
alpha = p.alpha;          % Cooperativity for enhanced binding interactions
dK = p.dK;                % Scaling of binding affinities from in vitro to in vivo (AS/HiBiT Combo Only).

% Degradation parameters
g1 = p.g1;            % Free P degradation rate.
kcat = p.kcat;
v_P = p.v_P;          % Basal P production rate.
% calculating KD2 and KD4 based on cooperativity.
KD4 = KD1/alpha;
KD2 = KD3/alpha;

% Calculating rates of change.
dydt = zeros(size(y));

% P
dydt(1) = v_P - g1*P - k3*E*P + k3*KD1*EP - k2*EG*P + k2*KD2*EGP;

% E
dydt(2) = kcat*EGP - k3*E*P + k3*KD3*EP - k1*E*G + k1*KD1*EG;

% G
dydt(3) = kcat*EGP - k1*E*G + k1*KD3*EG - k4*EP*G + k4*KD4*EGP;

% EP
dydt(4) = k3*E*P - k3*KD3*EP + k4*KD4*EGP - k4*EP*G;

% EG
dydt(5) = k1*E*G - k1*KD1*EG + k2*KD2*EGP - k2*EG*P;

% EGP
dydt(6) = -kcat*EGP + k2*EG*P - k2*KD2*EGP + k4*EP*G - k4*KD4*EGP;


end