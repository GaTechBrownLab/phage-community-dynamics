function coeff_table = makecoefftable(existtable,ns,betaMat,alphaMat,R2,treatment,phage)
% Sarah Sundius, Sept 7, 2023
% Function to make a table that saves the coefficients from model fitting
% Inputs:   existtable = existing table of coefficients, 0 if no existing
%           ns = number of species
%           betaMat = interaction/growth rate matrix, ns x ns+1, [b11 b12
%           ... r1; b21 b22 ... r2; etc]
%           alphaMat = rescaled interaction coeff matrix (bij/bii), ns x ns
%           [a11 a12 ...; a21 a22 ...]
%           R2 = R-squared values, ns x 1 vector, [r21; r22; ...]
%           treatment = treatment name, string
%           phage = 'Y' or 'N' 
% Output:   coeff_table = table of model fitting parameters for each
%           treatment


% Set up the table (if no existtable)
if ~istable(existtable)
    if ns == 2
        varNames_coeff = {'Treatment','Phage','b11','b12','b21','b22','r1','r2','a11','a12','a21','a22','R21','R22','avgR2'};
        varTypes_coeff = ['string','string',repelem({'double'},length(varNames_coeff)-2)];
        coeff_table = table('Size',[15,length(varNames_coeff)],'VariableTypes',varTypes_coeff,'VariableNames',varNames_coeff);
    elseif ns == 3
        varNames_coeff = {'Treatment','Phage','b11','b12','b13','b21','b22','b23','b31','b32','b33','r1','r2','r3','a11','a12','a13','a21','a22','a23','a31','a32','a33','R21','R22','R23','avgR2'};
        varTypes_coeff = ['string','string',repelem({'double'},length(varNames_coeff)-2)];
        coeff_table = table('Size',[13,length(varNames_coeff)],'VariableTypes',varTypes_coeff,'VariableNames',varNames_coeff);
    elseif ns == 4
        varNames_coeff = {'Treatment','Phage','b11','b12','b13','b14','b21','b22','b23','b24','b31','b32','b33','b34','b41','b42','b43','b44','r1','r2','r3','r4','a11','a12','a13','a14','a21','a22','a23','a24','a31','a32','a33','a34','a41','a42','a43','a44','R21','R22','R23','R24','avgR2'};
        varTypes_coeff = ['string','string',repelem({'double'},length(varNames_coeff)-2)];
        coeff_table = table('Size',[4,length(varNames_coeff)],'VariableTypes',varTypes_coeff,'VariableNames',varNames_coeff);   
    end

else
    coeff_table = existtable;

end

% Prepare variables for input
betaline = reshape(betaMat(1:ns,1:ns)',1,ns.^2);
alphaline = reshape(alphaMat',1,ns.^2);
rline = reshape(betaMat(:,ns+1),1,ns);
treatmentline = [betaline rline alphaline R2' mean(R2)];

% Add data
if ns == 2
    for i = 1:15
        if ismissing(coeff_table(i,'Treatment'))
            coeff_table(i,'Treatment') = {treatment};
            coeff_table(i,'Phage') = {phage};
            coeff_table(i,3:end) = array2table(treatmentline);
            break
        end
    end

elseif ns == 3
    for i = 1:13
        if ismissing(coeff_table(i,'Treatment'))
            coeff_table(i,'Treatment') = {treatment};
            coeff_table(i,'Phage') = {phage};
            coeff_table(i,3:end) = array2table(treatmentline);
            break
        end
    end
    
else
    for i = 1:4
        if ismissing(coeff_table(i,'Treatment'))
            coeff_table(i,'Treatment') = {treatment};
            coeff_table(i,'Phage') = {phage};
            coeff_table(i,3:end) = array2table(treatmentline);
            break
        end
    end
end

end