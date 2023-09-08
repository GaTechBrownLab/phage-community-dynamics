function solpts = solvelvode(p,numsp,tspan,N0,scale)
% Sarah Sundius, February 6, 2023
%
% Function to solve the Lotka-Volterra ODE --> using ode45? write own
% solver?
% Inputs:   p = set of parameters
%           numsp = number of species
%           tspan = time points
%           N0 = initial density
%           scale = boolean to use log scale or not (1 = log scaled)
%
% Output:   solpts = time series solution to model

sol = ode45(@(t,N)lvgrowthode(t,N,p,numsp,scale),tspan,N0);
solpts = deval(sol,tspan);


end