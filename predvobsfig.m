function predvobsfig(PA14_N,PA14_Y,CRKO_N,CRKO_Y,species,titletxt)
% Sarah Sundius, February 8, 2023
%
% Function to make plot of predicted vs. observed values as compared to 1:1
% ratio line 
% Inputs:   PA14_N, CRKO_N = 2D cell array, {1} observed data, {2} predicted/sim data, no phage
%           PA14_Y, CRKO_Y = 2D cell array, {1} observed data, {2} predicted/sim data, with phage
%           all data oriented (r,c) = (time,numsp)
%           species = cell array of species names ('PA14','CRKO','AB','BC','SA')
%           titletxt = cell array with table title (ex. '+ S. aureus')
% Outputs:  doesn't return anything just makes the plot

% Set colors for species 
numsp = size(species,2);
colors = cell(1,numsp);
for i = 1:numsp
    if strcmp(species{i},'PA14') || strcmp(species{i},'CRKO')
        colors{i} = '#06141F';
    elseif strcmp(species{i},'AB')
        colors{i} = '#CD4F38';
    elseif strcmp(species{i},'BC')
        colors{i} = '#3D4F7D';
    elseif strcmp(species{i},'SA')
        colors{i} = '#E48C2A';
    end
end

% Make 1:1 line
line1to1 = [10^0,10^5,10^10];

% Make figure
figure
hold on
box on
plot(log10(line1to1),log10(line1to1),'--','LineWidth',1.5,'Color',[0.7 0.7 0.7],'DisplayName','pred = obs')

for i = 1:numsp
    plot(log10(PA14_N{1}(:,i)),log10(PA14_N{2}(:,i)),'.','MarkerSize',15,'Color',colors{i},'DisplayName',strcat(species{i},', no phage'))
    plot(log10(PA14_Y{1}(:,i)),log10(PA14_Y{2}(:,i)),'o','LineWidth',1,'MarkerSize',5,'Color',colors{i},'DisplayName',strcat(species{i},', w/ phage'))
    %if i == 1
    %    plot(log10(CRKO_N{1}(:,i)),log10(CRKO_N{2}(:,i)),'.','MarkerSize',15,'Color',colors{i},'DisplayName',strcat('CRKO, no phage'))
    %    plot(log10(CRKO_Y{1}(:,i)),log10(CRKO_Y{2}(:,i)),'o','LineWidth',1,'MarkerSize',5,'Color',colors{i},'DisplayName',strcat(species{i},'CRKO, w/ phage'))
    %else
    %    plot(log10(CRKO_N{1}(:,i)),log10(CRKO_N{2}(:,i)),'.','MarkerSize',15,'Color',colors{i},'DisplayName',strcat(species{i},', no phage'))
    %    plot(log10(CRKO_Y{1}(:,i)),log10(CRKO_Y{2}(:,i)),'o','LineWidth',1,'MarkerSize',5,'Color',colors{i},'DisplayName',strcat(species{i},', w/ phage'))
    %end
end

xlim([0 10])
xticks([2,4,6,8,10])
xticklabels({'10^2','10^4','10^6','10^8'})
ylim([0 10])
yticks([2,4,6,8,10])
yticklabels({'10^2','10^4','10^6','10^8'})
xlabel('Observed bacterial density, log10(cfu/ml)')
ylabel('Predicted bacterial density, log10(cfu/ml)')
legend('Location','northeastoutside')
title(titletxt,'FontWeight','normal','FontAngle','italic','FontSize',12)

