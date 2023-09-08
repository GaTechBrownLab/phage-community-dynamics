function dens_table = makedenstable2(time,PA14N,PA14Y,CRKON,CRKOY,Treatment,Target)
% Sarah Sundius, February 16, 2023
% Function to make a table that saves predicted/simulated densities
% Treatments: +AB, +BC, +SA, +AB+SA, +BC+SA, +BC+AB, +MC
% Targets: AB, BC, SA

% Set up the table
varNames_dens = {'Treatment','Target','Time','Quantity','Phage'};
varTypes_dens = {'string','string','double','double','string'};
treatment1 = cell2mat(strcat('PA14',Treatment));
treatment2 = cell2mat(strcat('CRKO',Treatment));

ns = size(PA14N,2);
nt = size(time,2);
rows = (ns*nt)*4;

dens_table = table('Size',[rows,5],'VariableTypes',varTypes_dens,'VariableNames',varNames_dens);
dens_table.Time = repmat(time',[ns*4,1]);
dens_table.Phage = char([repelem({'N'},rows/2) repelem({'Y'},rows/2)]');
dens_table.Treatment = categorical([repelem(strcat('PA14',Treatment),rows/4) repelem(strcat('CRKO',Treatment),rows/4) repelem(strcat('PA14',Treatment),rows/4) repelem(strcat('CRKO',Treatment),rows/4)]');
dens_table.Target = categorical(repmat([repelem({'PA14'},nt) repelem(Target,nt) repelem({'CRKO'},nt) repelem(Target,nt)]',[2,1]));

% Add bacterial densities
% strain: PA14
dens_table(dens_table.Phage=='N'&dens_table.Target=='PA14',:).Quantity = PA14N(:,1);
dens_table(dens_table.Phage=='Y'&dens_table.Target=='PA14',:).Quantity = PA14Y(:,1);
for i = 1:ns-1
    dens_table(dens_table.Treatment==treatment1&dens_table.Phage=='N'&dens_table.Target==Target{i},:).Quantity = PA14N(:,i+1);
    dens_table(dens_table.Treatment==treatment1&dens_table.Phage=='Y'&dens_table.Target==Target{i},:).Quantity = PA14Y(:,i+1);
end   

% strain: CRKO
dens_table(dens_table.Phage=='N'&dens_table.Target=='CRKO',:).Quantity = CRKON(:,1);
dens_table(dens_table.Phage=='Y'&dens_table.Target=='CRKO',:).Quantity = CRKOY(:,1);
for i = 1:ns-1
    dens_table(dens_table.Treatment==treatment2&dens_table.Phage=='N'&dens_table.Target==Target{i},:).Quantity = CRKON(:,i+1);
    dens_table(dens_table.Treatment==treatment2&dens_table.Phage=='Y'&dens_table.Target==Target{i},:).Quantity = CRKOY(:,i+1);
end



end