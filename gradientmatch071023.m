function [r,alpha,beta,psol] = gradientmatch071023(t,N,lb,ub,p0,logderiv)
% Sarah Sundius, July 10, 2023
%
% Function to fit a Lotka-Volterra model by gradient matching 
% Inputs:   t = times
%           N = bacterial densities - (r,c) = (time,numsp)
%           lb/ub = lower/upper bounds on fit parameters - [r; A]
%           p0 = initial parameter guess - [r; A]
%           logderiv = how to do derivative est ('gradient' or
%           'splinefit' - Liao et al.)
% Outputs:  r = growth rates
%           alpha = interspecific interaction coeffs (scaled by Nmax)
%           beta = all interaction coeffs (non-scaled)
%           psol = raw solution to optimization prob

Nmax = max(N);
N = N./Nmax;

% Calculate log derivative - gradients
if strcmp(logderiv,'gradient')
    nt = size(t,2);
    ns = size(N,2);
    lnN = log(N)';
    dlnNdt = zeros(ns,nt);
    for i = 1:ns
        dlnNdt(i,:) = gradient(lnN(i,:))./gradient(t);
    end

% Calculate log derivative - spline fitting
elseif strcmp(logderiv,'splinefit')
    nt = size(t,2);
    ns = size(N,2);
    dlnNdt = zeros(ns,nt);   %   [# of species] x [# of time]
    for i=1:ns
        % using spline to find the correct derivative at each N data
        pp = spline(t,log(N(:,i)));
        pder = fnder(pp, 1);
        dlnNdt(i,:) = ppval(pder, t);
    end
    
end

% Define function whose sum of squares is minimized for lsqnonlin, Cx-d
Nmat = ones(length(t),ns+1);
for i = 1:ns
    Nmat(:,i+1) = N(:,i);
end
dlnNdtfcn = @(p) sum(sum((Nmat*p - dlnNdt').^2));   %form 1: p = [r1 r2; b11 b21; b12 b22...]

% Solve lsqnonlin
pinitial = [p0(1,:); (p0(2:end,:))'];
lb = [lb(1,:); (lb(2:end,:))'];
ub = [ub(1,:); (ub(2:end,:))'];
opts = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxIterations',1000,'Display','off');
psol = lsqnonlin(dlnNdtfcn,pinitial,lb,ub,opts);

% Rearrange parameters
% scaled by max version                    
r = psol(1,:)';                     %growth rates
alpha = psol(2:end,:)';             %interspecific interaction coeff
beta = alpha./Nmax;                %all interaction coeffs 


end