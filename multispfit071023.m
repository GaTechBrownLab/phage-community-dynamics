function [pred,beta,alpha,R2,method] = multispfit071023(data,time,lb,ub,pguess,scale)
% Sarah Sundius, February 13, 2023
%
% Function to fit a Lotka-Volterra model for multi-species data, compares
% methods and selects the best based off average R-squared
% Inputs:   data = bacterial densities, (r,c) = (time,species)
%           time = time points
%           pguess = initial guess for parameters - [r; A] 
%           lb/ub = lower/upper parameter bound - [r; A] 
%           scale = boolean to use log scaled or not (1 = log scaled)
% Outputs:  pred = predicted bacterial densities
%           beta = raw model coefficients (i.e. const, lin, quad)
%           A = interaction coefficients
%           k = carrying capacities
%           R2 = R-squared (measure of model accuracy)
%           method = method producing the best results

% Get number of species
[nt,ns] = size(data);

time2 = time;   % unsmoothed data
data2 = data;   % unsmoothed data

% Method 1: Gradient matching - check gradient and spline fitting est for log deriv
[r2g,a2g,b2g,psol2g] = gradientmatch071023(time2,data2,lb,ub,pguess,'gradient');
pred_2g = glv_simulation(time,data(1,:),[b2g r2g]);
R2_2g = Rsquare(data,pred_2g);
%R2_2g = RMSE2(data,pred_2g);
avg_2g = mean(R2_2g);

[r2s,a2s,b2s,psol2s] = gradientmatch071023(time2,data2,lb,ub,pguess,'splinefit');
pred_2s = glv_simulation(time,data(1,:),[b2s r2s]);
R2_2s = Rsquare(data,pred_2s);
%R2_2s = RMSE2(data,pred_2s);
avg_2s = mean(R2_2s);

if avg_2g > avg_2s 
    r2 = r2g;
    b2 = b2g;
    a2 = a2g;
    pred_2 = pred_2g;
    R2_2 = R2_2g;
    method2 = 'gradient match, gradient';
else
    r2 = r2s;
    b2 = b2s;
    a2 = a2s;
    pred_2 = pred_2s;
    R2_2 = R2_2s;
    method2 = 'gradient match, spline';
end

% Method 2: Liao et al. 2020 - check input deriv (gradient) vs. their method (spline fitting) 
Nmax = max(data2);
Nscale = data2./Nmax;

lb = [lb(2:ns+1,:) lb(1,:)'];
ub = [ub(2:ns+1,:) ub(1,:)'];
dLguess = zeros(ns,length(time2));
for i = 1:ns
    dLguess(i,:) = gradient(log(Nscale(:,i)))./gradient(time2');
end

[optBeta3g,initialGuessDL3g] = glv_linreg(time2,Nscale,lb,ub,'logderiv',dLguess);
Beta3g = [optBeta3g(1:ns,1:ns)./Nmax optBeta3g(:,ns+1)];
pred_3g = glv_simulation(time,data2(1,:),Beta3g);
R2_3g = Rsquare(data,pred_3g);
%R2_3g = RMSE2(data,pred_3g);
avg_3g = mean(R2_3g);

[optBeta3s,initialGuessDL3s] = glv_linreg(time2,Nscale,lb,ub);
Beta3s = [optBeta3s(1:ns,1:ns)./Nmax optBeta3s(:,ns+1)];
pred_3s = glv_simulation(time,data2(1,:),Beta3s);
R2_3s = Rsquare(data,pred_3s);
%R2_3s = RMSE2(data,pred_3s);
avg_3s = mean(R2_3s);

if avg_3g > avg_3s 
    b3 = Beta3g;
    a3 = optBeta3g(1:ns,1:ns);
    pred_3 = pred_3g;
    R2_3 = R2_3g;
    method3 = 'Liao, gradient';
else
    b3 = Beta3s;
    a3 = optBeta3s(1:ns,1:ns);
    pred_3 = pred_3s;
    R2_3 = R2_3s;
    method3 = 'Liao, spline';
end

% Method selection and assignment
ranked = sort([mean(R2_2),mean(R2_3)],'descend');

if ranked(1,1) == mean(R2_2)
    pred = pred_2; 
    R2 = R2_2;
    beta = [b2 r2];
    alpha = a2;
    method = method2;

else
    pred = pred_3;
    R2 = R2_3;
    beta = b3;
    alpha = a3;
    method = method3;
    
end

end