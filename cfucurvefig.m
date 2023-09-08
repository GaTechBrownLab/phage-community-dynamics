function cfucurvefig(time,data,pred,species,titletxt)
% Sarah Sundius, February 8, 2023
%
% Function to plot bacterial density (cfu/ml) for observed and predicted trajectories vs. time for a community of
% up to four species
% Inputs:   time = time points, (r,c) = (time,1)
%           data = observed data, (r,c) = (time,numsp)
%           pred = model predicted/simulated data, (r,c) = (time,numsp)
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

% Make figure
figure
hold on
grid on
box on
plot([0,11*24],[2,2],'--','LineWidth',1,'Color',[0.7 0.7 0.7],'HandleVisibility','off')

for i = 1:numsp
    plot(time,log10(data(:,i)),'.--','LineWidth',1.5,'MarkerSize',15,'Color',colors{i},'DisplayName',strcat(species{i},' obs'))
    plot(time,log10(pred(:,i)),'-','LineWidth',1.5,'Color',colors{i},'DisplayName',strcat(species{i},' pred'))
    %plot(time,(data(:,i)),'.--','LineWidth',1.5,'MarkerSize',15,'Color',colors{i},'DisplayName',strcat(species{i},' obs'))
    %plot(time,(pred(:,i)),'-','LineWidth',1.5,'Color',colors{i},'DisplayName',strcat(species{i},' pred'))
end

xlim([0,11*24])
xticks([0,1,3,7,10].*24)
xticklabels({'0','3','5','7','10'})
ylim([0,11])
yticks([2,4,6,8,10])
yticklabels({'10^2','10^4','10^6','10^8'})
%ylim([0 1.5])
xlabel('Days post infection')
ylabel('Bacterial density, cfu/ml')
legend('Location','northeastoutside')
title(titletxt,'FontWeight','normal','FontAngle','italic','FontSize',12)


end