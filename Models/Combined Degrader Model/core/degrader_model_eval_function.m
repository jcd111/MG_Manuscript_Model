function y = degrader_model_eval_function(obj,tspan)
% Function for evaluating simple reversible binding model with P
% degrdataion.
% Setting binding affinites based on the drug and targets being used.
% Getting drug and target info.
P_name = obj.parameters.P_name;
G_name = obj.parameters.G_name;

% Setting KD3 based on P name
obj.parameters.KD1 = obj.parameters.(append("KD1_",P_name));
% Setting KD2 based on G name and P name
obj.parameters.alpha = obj.parameters.(append("alpha_",G_name,"_",P_name));
% Setting KD3 based on G name
obj.parameters.KD3 = obj.parameters.(append("KD3_",G_name));

% Setting Et and Pt values based on experiment type
if strcmp(obj.parameters.experiment_type,"MRM")
    obj.parameters.Pt = obj.parameters.(append("Pt_",P_name,"_",obj.parameters.experiment_type));
    obj.parameters.Et = obj.parameters.(append("Et_",obj.parameters.experiment_type));
elseif strcmp(obj.parameters.experiment_type,"BRET")
    obj.parameters.kcat = 0;
    obj.parameters.g1 = 0;
    obj.parameters.Pt = obj.parameters.(append("Pt_",obj.parameters.experiment_type));
    obj.parameters.Et = obj.parameters.(append("Et_",obj.parameters.experiment_type));
elseif strcmp(obj.parameters.experiment_type,"AS")
    obj.parameters.kcat = 0;
    obj.parameters.g1 = 0;
    obj.parameters.Pt = obj.parameters.(append("Pt_",obj.parameters.experiment_type));
    obj.parameters.Et = obj.parameters.(append("Et_",obj.parameters.experiment_type));
    % Scaling binding affinities by dK
    obj.parameters.KD1 = obj.parameters.dK*obj.parameters.KD1;
    obj.parameters.KD3 = obj.parameters.dK*obj.parameters.KD3;
elseif strcmp(obj.parameters.experiment_type,"HIBIT")
    obj.parameters.Pt = obj.parameters.(append("Pt_",obj.parameters.experiment_type));
    obj.parameters.Et = obj.parameters.(append("Et_",obj.parameters.experiment_type));
end


% Setting kcat based on G and P
obj.parameters.kcat = obj.parameters.(append("kcat_",G_name,"_",P_name));


% Determining initial amount of EP binding based on equilibrium.
Pt = obj.parameters.Pt; Et = obj.parameters.Et; KD1 = obj.parameters.KD1; k1 = obj.parameters.k1;
kcat_EP = obj.parameters.kcat_EP; g1 = obj.parameters.g1;

Kd_eff = (kcat_EP + k1*KD1)/k1;

P0 = 0.5*(Pt - Kd_eff - Et + sqrt((Pt - Kd_eff - Et).^2 + 4*Pt*Kd_eff));
EP0 = Pt - P0;
E0 = Et - EP0;

obj.parameters.v_P = g1*P0 + k1*E0*P0 - k1*KD1*EP0;
% Setting initial conditions.
obj.y0 = zeros(1,length(obj.species_names));
obj.y0(strcmp(obj.species_names,"P")) = P0;
obj.y0(strcmp(obj.species_names,"EP")) = EP0;
obj.y0(strcmp(obj.species_names,"E")) = E0;
obj.y0(strcmp(obj.species_names,"G")) = obj.parameters.Gt;



% Simulating model, checking whether t = 0 is already included as a
% timepoint.
options = odeset('NonNegative',1,'AbsTol',1e-8);
dydt = @(t,y) obj.rate_laws(t,y,obj.parameters);
if tspan(1) == 0
    [~,y] = ode15s(dydt,tspan,obj.y0,options);
else
    % If 0 is not included, adding to beginning of tspan then
    % removing in output.
    % if tspan is a column vector, changing to row
    if size(tspan,2) == 1
        tspan = tspan';
    end
    [~,y_tmp] = ode15s(dydt,[0,tspan],obj.y0,options);
    y = y_tmp(2:end,:);
end

end



