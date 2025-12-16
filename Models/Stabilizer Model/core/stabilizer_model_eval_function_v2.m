function y = stabilizer_model_eval_function_v2(obj,tspan)
% Function for evaluating simple reversible binding model with P
% degrdataion.
      % Setting binding affinites based on the drug and targets being used.
      % Getting drug and target info.

      % 
      % obj.parameters.a_P = g1*P0 + k1*E0*P0 -k1*KD1*EP0;
      % Setting initial conditions.
      obj.y0 = zeros(1,length(obj.species_names));
      obj.y0(strcmp(obj.species_names,"P")) = obj.parameters.Pt;
      obj.y0(strcmp(obj.species_names,"EP")) = 0;
      obj.y0(strcmp(obj.species_names,"E")) = obj.parameters.Et;
      obj.y0(strcmp(obj.species_names,"G")) = obj.parameters.Gt;


      % Simulating with ODE15s
      options = odeset('NonNegative',1,'AbsTol',1e-12);
      dydt = @(t,y) obj.rate_laws(t,y,obj.parameters);
      [t,y] = ode15s(dydt,tspan,obj.y0,options);

end



