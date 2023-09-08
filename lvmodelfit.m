function [psol,sumsq] = lvmodelfit(tspan,data,lb,ub,pguess,scale)
% Sarah Sundius, February 6, 2023
%
% Function to fit a Lotka-Volterra model to experimental data
% Inputs:   tspan = time points
%           data = true densities (experimental or simulated), oriented
%           (r,c) = (time,species)
%           psize = number of fit parameters
%           lb/ub = lower/upper bounds on fit parameters, if no bound = []
%           pguess = initial parameter guess
%           scale = boolean to use log scale or not (1 = log scaled)
%
% Outputs:  psol = fit parameter values
%           sumsq = sum of squares

% Get the dimensions of input data (rows = time, cols = species)
numsp = size(data,2);

% Prepare the optimization problem
N0 = data(1,:);
p = optimvar('p',numsp,numsp+1,'LowerBound',lb,'UpperBound',ub);
myfcn = fcn2optimexpr(@solvelvode,p,numsp,tspan,N0,scale);
obj = sum(sum((myfcn' - data).^2));
prob = optimproblem('Objective',obj);
opts = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxFunctionEvaluations',1000,'StepTolerance',1e-12,'FunctionTolerance',1e-12);

% Solve the problem
p0.p = pguess;
[psol,sumsq] = solve(prob,p0,Options=opts);


end