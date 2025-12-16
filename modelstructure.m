classdef modelstructure
   properties
       parameters {mustBeA(parameters,'struct')} = struct()
       parameter_names {mustBeA(parameter_names,["cell","string"])} = {}
       bounds {mustBeA(bounds,'struct')} = struct()
       experimental_data {mustBeA(experimental_data,'cell')} = {}
       species_names{mustBeA(species_names,["cell","string"])} = {};
       eval_function {mustBeA(eval_function,'function_handle')} = @(obj,tspan) normal_eval(obj,tspan)
       rate_laws {mustBeA(rate_laws,'function_handle')} = @(x,y,k) x 
       y0 double = 0
   end
   
   
   methods
       function obj = modelstructure(varargin)
           p = inputParser;
           checkName = @(x) ischar(x) || isstring(x) || iscell(x);
           checkNameOrCell = @(x) ischar(x) || isstring(x) || ...
               (iscell(x) && ( ischar(x{1}) || isstring(x{1}) ) );
           checkBool = @(x) islogical(x) || isnumeric(x);
           checkFunction = @(x)isa(x,'function_handle');
           checkNumOrCell = @(x) isnumeric(x) || iscell(x);
           checkStruct = @(x) isa(x,'struct');
           
           addParameter(p,'parameters',struct(),checkStruct);
           addParameter(p,'parameter_names',{},checkName);
           addParameter(p,'bounds',struct(),checkStruct);
           addParameter(p,'experimental_data',{},@iscell);
           addParameter(p,'eval_function',@(obj,tspan) normal_eval(obj,tspan), checkFunction);
           addParameter(p,'species_names',{},checkNameOrCell);
           addParameter(p,'rate_laws',@(x,y,k) x,checkFunction);
           p.KeepUnmatched = false;
           p.CaseSensitive = false;
           
           parse(p,varargin{:});
           
           obj.parameter_names = p.Results.parameter_names;
           obj.parameters = p.Results.parameters;
           obj.bounds = p.Results.bounds;
           obj.experimental_data = p.Results.experimental_data;
           obj.eval_function = p.Results.eval_function;
           obj.species_names = p.Results.species_names;
           obj.rate_laws = p.Results.rate_laws;
           
       end
       
       function [y] = evaluate(obj,tspan)   
           % If tspan exists, input to eval function
           if exist('tspan','var')
                y = obj.eval_function(obj,tspan);
           % If not, no time dependence     
           else
               y = obj.eval_function(obj);
           end
       end
       
   end

end