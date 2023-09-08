function err = RMSE2(data,pred)
% function to calculate the rmse values for predicted data
% input: observational data (time is cols), fitted/predicted values (time is cols)
% output: rmse or normalized rmse value

numdata = size(data,1);
err = sqrt(sum((data-pred).^2)./numdata);


end