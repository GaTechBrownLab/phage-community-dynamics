function cfucurvefig_panels(time,data,pred,species,titletxt,avgR2,pos)
% Sarah Sundius, February 15, 2023
%
% Function to make subplots of bacterial density (cfu/ml) for observed and predicted trajectories vs. time for a community of
% up to four species
% Inputs:   time = time points, (r,c) = (time,1)
%           data = observed data, (r,c) = (time,numsp)
%           pred = model predicted/simulated data, (r,c) = (time,numsp)
%           species = cell array of species names ('PA14','CRKO','AB','BC','SA')
%           titletxt = cell array with table title (ex. '+ S. aureus')
%           avgR2 = average R-squared of model prediction
%           pos = tile number/position for tiled layout
% Outputs:  doesn't return anything just makes the plot

% Set colors for species 
numsp = size(species,2);
colors = cell(1,numsp);
for i = 1:numsp
    if strcmp(species{i},'PA14') || strcmp(species{i},'CRKO')
        colors{i} = '#06141F';
        leg{i} = 'PA';
    elseif strcmp(species{i},'AB')
        colors{i} = '#CD4F38';
        leg{i} = species{i};
    elseif strcmp(species{i},'BC')
        colors{i} = '#3D4F7D';
        leg{i} = species{i};
    elseif strcmp(species{i},'SA')
        colors{i} = '#E48C2A';
        leg{i} = species{i};
    end
end

% Make figure
nexttile(pos)
hold on
grid on
box on
plot([0,11*24],[2,2],'--','LineWidth',1,'Color',[0.7 0.7 0.7],'HandleVisibility','off')

for i = 1:numsp
    plot(time,log10(data(:,i)),'.--','LineWidth',1.5,'MarkerSize',15,'Color',colors{i},'HandleVisibility','off')
    plot(time,log10(pred(:,i)),'-','LineWidth',1.5,'Color',colors{i},'DisplayName',leg{i})
end

xlim([0,11*24])
xticks([0,1,3,7,10].*24)
ax = get(gca,'XTickLabel');
set(gca,'XTickLabel',ax,'fontsize',8,'FontWeight','Normal')
set(gca,'XTickLabelMode','auto')
xticklabels({'0','3','5','7','10'})

ylim([0,11])
yticks([2,4,6,8,10])
ay = get(gca,'YTickLabel');
set(gca,'YTickLabel',ay,'fontsize',8,'FontWeight','Normal')
set(gca,'YTickLabelMode','auto')
yticklabels({'10^2','10^4','10^6','10^8'})

title(titletxt{1},'FontWeight','normal','FontSize',9)

if abs(avgR2) < 10
    text(140,10,['R^2 = ',num2str(avgR2,'%.2f')],'FontWeight','normal','FontSize',8)
elseif abs(avgR2) < 100
    text(125,10,['R^2 = ',num2str(avgR2,'%.2f')],'FontWeight','normal','FontSize',8)
else
    text(115,10,['R^2 = ',num2str(avgR2,'%.2f')],'FontWeight','normal','FontSize',8)
end

end