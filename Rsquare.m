function R2 = Rsquare(data,pred)
% function to calculate the R^2 values for predicted data
% input: observational data and fitted/predicted values (time is cols for
% both)
% output: R^2 value

% calculate residuals
y = data;
f = pred;
e = y - f;

% calculate the mean of the observed data
ybar = mean(y,1);

% calculate the sum of squares of the residuals
SSres = sum(e.^2,1);
    
% calculate the total sum of squares
SStot = sum((y-ybar).^2,1);

% calculate R2
R2 = 1 - (SSres./SStot);


end