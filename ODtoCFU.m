function cfu_ml = ODtoCFU(OD600,species)
% function to convert OD600 to cfu/ml for Ellinor's MC experiment data
% input: OD600 = time series of OD600 measurements, species = 'PA14',
% 'SA', 'AB', 'BC' - type of bacteria
% output: cfu_ml = corresponding cfu/ml measurements, converted via
% standard curves


% Read data table
bac_data = readtable('CFUvsOD600 Ellinor.csv','VariableNamingRule','preserve');

% Parse data + convert to doubles
if strcmp(species,'PA14')
    x = bac_data{1:5,5};  %OD600(nm)
    y = bac_data{1:5,4};  %cfu/ml
elseif strcmp(species,'SA')
    x = bac_data{6:10,5};  %OD600(nm)
    y = bac_data{6:10,4};  %cfu/ml
elseif strcmp(species,'AB')
    x = bac_data{11:14,5};  %OD600(nm)
    y = bac_data{11:14,4};  %cfu/ml
elseif strcmp(species,'BC')
    x = bac_data{15:19,5};  %OD600(nm)
    y = bac_data{15:19,4};  %cfu/ml
end

% Fit model
mdl = fitlm(x,y);

% Predict cfu/ml for input OD600 values
cfu_ml = predict(mdl,OD600);

% checks
%figure
%hold on
%plot(mdl)
%xlabel('OD600 (nm)')
%ylabel('cfu/ml')
%title(species)




end
