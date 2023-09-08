%% Multi-Species Model Fitting
% Prediction of full community dynamics from 2- and 3-species data - final modeling

% Author: Sarah Sundius
% Date: September 3, 2023

clear; clc;

%% Load data 
% Using data curves from Rafael 08/18/2022

% Input data
data_ss = readtable('Bacterial species SA PA AB BC growth in LB 24h at 37C_18_08_2022.xlsx');

% Convert to matrix
data_ss = data_ss{:,:};

% Separate variables
time = data_ss(1:end,1);
PA14 = data_ss(1:end,2:13);
AB = data_ss(1:end,14:25);
BC = data_ss(1:end,26:37);
SA = data_ss(1:end,38:49);

% Average replicates and smooth
PA14avg = smooth(mean(PA14,2));
ABavg = smooth(mean(AB,2));
BCavg = smooth(mean(BC,2));
SAavg = smooth(mean(SA,2));

% Convert OD600 to cfu/ml
cfuPA14 = ODtoCFU(PA14avg,'PA14');
cfuAB = ODtoCFU(ABavg,'AB');
cfuBC = ODtoCFU(BCavg,'BC');
cfuSA = ODtoCFU(SAavg,'SA');


%% Fit single species curves
% Uses ode45 to solve for ln(N) from dln(N)/dt = r(1-N/k), least squares
% regression

% Set up to save parameters
rsave1 = zeros(1,4);
ksave1 = zeros(1,4);
r21 = zeros(1,4);

% Omit negative densities
maxPA = max(cfuPA14);
maxAB = max(cfuAB);
maxBC = max(cfuBC);
maxSA = max(cfuSA);

cfuPA14s = log(cfuPA14(2:end,1)./maxPA);
cfuABs = log(cfuAB(2:end,1)./maxAB);
cfuBCs = log(cfuBC(3:end,1)./maxBC);
cfuSAs = log(cfuSA(2:end,1)./maxSA);

cfuPA142 = log(cfuPA14(2:end,1));
cfuAB2 = log(cfuAB(2:end,1));
cfuBC2 = log(cfuBC(3:end,1));
cfuSA2 = log(cfuSA(2:end,1));

time = time(2:end,1);
time2 = time(2:end,1);
scale = 1; %set methods to use log scaled data

% Fit to logistic growth model and save parameters
PA140 = cfuPA142(1,1);
PA14guess = [1 1]; 
[psolPA14,sumsqPA14] = lvmodelfit(time,cfuPA14s,[],[],PA14guess,scale);
rsave1(1,1) = psolPA14.p(1);
ksave1(1,1) = psolPA14.p(2).*maxPA;

AB0 = cfuAB2(1,1);
ABguess = [1 1];
[psolAB,sumsqAB] = lvmodelfit(time,cfuABs,[],[],ABguess,scale);
rsave1(1,2) = psolAB.p(1);
ksave1(1,2) = psolAB.p(2).*maxAB;

BC0 = cfuBC2(1,1);
BCguess = [1 1];
[psolBC,sumsqBC] = lvmodelfit(time2,cfuBCs,[],[],BCguess,scale);
rsave1(1,3) = psolBC.p(1);
ksave1(1,3) = psolBC.p(2).*maxBC;

SA0 = cfuSA2(1,1);
SAguess = [1 1];
[psolSA,sumsqSA] = lvmodelfit(time,cfuSAs,[],[],SAguess,scale);
rsave1(1,4) = psolSA.p(1);
ksave1(1,4) = psolSA.p(2).*maxSA;

% Simulate logistic growth using fit parameters 
PA14pred_cfu = solvelvode([rsave1(1,1),ksave1(1,1)],1,time,PA140,scale);
ABpred_cfu = solvelvode([rsave1(1,2),ksave1(1,2)],1,time,AB0,scale);
BCpred_cfu = solvelvode([rsave1(1,3),ksave1(1,3)],1,time2,BC0,scale);
SApred_cfu = solvelvode([rsave1(1,4),ksave1(1,4)],1,time,SA0,scale);

% Calculate R^2
r21(1,1) = Rsquare(cfuPA142,PA14pred_cfu');
r21(1,2) = Rsquare(cfuAB2,ABpred_cfu');
r21(1,3) = Rsquare(cfuBC2,BCpred_cfu');
r21(1,4) = Rsquare(cfuSA2,SApred_cfu');

% Calculate doubling time
DT1 = 60.*(log(2)./rsave1);

% Plot, compare actual and predicted
%{
figure
subplot(2,2,1)
hold on
plot(time,log10(exp(cfuPA142)),'.-','Color','#06141F','LineWidth',1,'MarkerSize',10)
plot(time,log10(exp(PA14pred_cfu)),'-','Color','#06141F','LineWidth',1)
ylim([0,11])
yticks([2,4,6,8,10])
yticklabels({'10^2','10^4','10^6','10^8'})
xlim([0 11])
xlabel('time')
ylabel('log10(cfu/ml)')
title('PA14')
legend('Data','Predicted','Location','southeast')

subplot(2,2,2)
hold on
plot(time,log10(exp(cfuAB2)),'.-','Color','#CD4F38','LineWidth',1,'MarkerSize',10)
plot(time,log10(exp(ABpred_cfu)),'-','Color','#CD4F38','LineWidth',1)
ylim([0,11])
yticks([2,4,6,8,10])
yticklabels({'10^2','10^4','10^6','10^8'})
xlim([0 11])
xlabel('time')
ylabel('log10(cfu/ml)')
title('AB')
legend('Data','Predicted','Location','southeast')

subplot(2,2,3)
hold on
plot(time2,log10(exp(cfuBC2)),'.-','Color','#3D4F7D','LineWidth',1,'MarkerSize',10)
plot(time2,log10(exp(BCpred_cfu)),'-','Color','#3D4F7D','LineWidth',1)
ylim([0,11])
yticks([2,4,6,8,10])
yticklabels({'10^2','10^4','10^6','10^8'})
xlim([0 11])
xlabel('time')
ylabel('log10(cfu/ml)')
title('BC')
legend('Data','Predicted','Location','southeast')

subplot(2,2,4)
hold on
plot(time,log10(exp(cfuSA2)),'.-','Color','#E48C2A','LineWidth',1,'MarkerSize',10)
plot(time,log10(exp(SApred_cfu)),'-','Color','#E48C2A','LineWidth',1)
ylim([0,11])
yticks([2,4,6,8,10])
yticklabels({'10^2','10^4','10^6','10^8'})
xlim([0 11])
xlabel('time')
ylabel('log10(cfu/ml)')
title('SA')
legend('Data','Predicted','Location','southeast')
%}


%% Read in multi-species data - exp w/o PA

% 10-day multi species growth curves - w/o PA

% Input data 
sheet_name = sheetnames('10-day exp. community data wo PA sarah.xlsx');
for k=1:numel(sheet_name)
    data_pw10{k} = readtable('10-day exp. community data wo PA sarah.xlsx','Sheet',sheet_name(k),'VariableNamingRule','preserve');
end

impcolspw = [3,4,6,5];

% naming: treatment_targetpw (pairwise)
% Monoculture 
SA_SApw = parsedata(data_pw10{1}(1:25,impcolspw));
AB_ABpw = parsedata(data_pw10{2}(1:25,impcolspw));
BC_BCpw = parsedata(data_pw10{3}(1:25,impcolspw));

% 2-species experiment
ABSA_ABpw = parsedata(data_pw10{4}(1:25,impcolspw));
ABSA_SApw = parsedata(data_pw10{4}(26:50,impcolspw));

BCSA_BCpw = parsedata(data_pw10{5}(1:25,impcolspw));
BCSA_SApw = parsedata(data_pw10{5}(26:50,impcolspw));

BCAB_BCpw = parsedata(data_pw10{6}(1:25,impcolspw));
BCAB_ABpw = parsedata(data_pw10{6}(26:50,impcolspw));

% 3-species experiment
ABBCSA_ABpw = parsedata(data_pw10{7}(1:25,impcolspw));
ABBCSA_BCpw = parsedata(data_pw10{7}(26:50,impcolspw));
ABBCSA_SApw = parsedata(data_pw10{7}(51:75,impcolspw));

% Average replicates
SA_SApwq = mean(SA_SApw,2);
AB_ABpwq = mean(AB_ABpw,2);
BC_BCpwq = mean(BC_BCpw,2);

ABSA_ABpwq = mean(ABSA_ABpw,2);
ABSA_SApwq = mean(ABSA_SApw,2);

BCSA_BCpwq = mean(BCSA_BCpw,2);
BCSA_SApwq = mean(BCSA_SApw,2);

BCAB_BCpwq = mean(BCAB_BCpw,2);
BCAB_ABpwq = mean(BCAB_ABpw,2);

ABBCSA_ABpwq = mean(ABBCSA_ABpw,2);
ABBCSA_BCpwq = mean(ABBCSA_BCpw,2);
ABBCSA_SApwq = mean(ABBCSA_SApw,2);


%% Read in multi-species data - no PA

% 10-day multi species growth curves - no PA

% Input data 
% import data from .xlsx file - multiple sheets
sheet_name = sheetnames('10-day exp. community data no PA sarah.xlsx');
for k=1:numel(sheet_name)
    data_ms10{k} = readtable('10-day exp. community data no PA sarah.xlsx','Sheet',sheet_name(k),'VariableNamingRule','preserve');
end

% Separate by treatment, phage/no phage
impcolsN = [3,4,7,9];
impcolsY = [13,14,17,19];

% naming: treatment_strainN (no phage), treatment_strainY (yes phage) -
% target is treatment if single species treatment, no PA
% Strain: PA14
% 2-species experiments 
SA_PA14Nqf = parsedata(data_ms10{1}(1:25,impcolsN));
SA_PA14Yqf = parsedata(data_ms10{1}(:,impcolsY));

BC_PA14Nqf = parsedata(data_ms10{4}(1:25,impcolsN));
BC_PA14Yqf = parsedata(data_ms10{4}(:,impcolsY));

AB_PA14Nqf = parsedata(data_ms10{5}(1:25,impcolsN));
AB_PA14Yqf = parsedata(data_ms10{5}(:,impcolsY));

% Average replicates (quantities final)
SA_PA14Nq = mean(SA_PA14Nqf,2);
AB_PA14Nq = mean(AB_PA14Nqf,2);
BC_PA14Nq = mean(BC_PA14Nqf,2);

SA_PA14Yq = mean(SA_PA14Yqf,2);
AB_PA14Yq = mean(AB_PA14Yqf,2);
BC_PA14Yq = mean(BC_PA14Yqf,2);

% 3-species experiments (target_treatment_PAstrain, N=no phage)
AB_ABSA_PA14Nqf = parsedata(data_ms10{6}(1:25,impcolsN));
SA_ABSA_PA14Nqf = parsedata(data_ms10{6}(26:50,impcolsN));
AB_ABSA_PA14Yqf = parsedata(data_ms10{6}(1:49,impcolsY));
SA_ABSA_PA14Yqf = parsedata(data_ms10{6}(50:98,impcolsY));

BC_BCSA_PA14Nqf = parsedata(data_ms10{8}(1:25,impcolsN));
SA_BCSA_PA14Nqf = parsedata(data_ms10{8}(26:50,impcolsN));
BC_BCSA_PA14Yqf = parsedata(data_ms10{8}(1:49,impcolsY));
SA_BCSA_PA14Yqf = parsedata(data_ms10{8}(50:98,impcolsY));

AB_BCAB_PA14Nqf = parsedata(data_ms10{9}(1:25,impcolsN));
BC_BCAB_PA14Nqf = parsedata(data_ms10{9}(26:50,impcolsN));
AB_BCAB_PA14Yqf = parsedata(data_ms10{9}(1:49,impcolsY));
BC_BCAB_PA14Yqf = parsedata(data_ms10{9}(50:98,impcolsY));

% Average replicates (quantities final)
AB_ABSA_PA14Nq = mean(AB_ABSA_PA14Nqf,2);
SA_ABSA_PA14Nq = mean(SA_ABSA_PA14Nqf,2);
BC_BCSA_PA14Nq = mean(BC_BCSA_PA14Nqf,2);
SA_BCSA_PA14Nq = mean(SA_BCSA_PA14Nqf,2);
AB_BCAB_PA14Nq = mean(AB_BCAB_PA14Nqf,2);
BC_BCAB_PA14Nq = mean(BC_BCAB_PA14Nqf,2);

AB_ABSA_PA14Yq = mean(AB_ABSA_PA14Yqf,2);
SA_ABSA_PA14Yq = mean(SA_ABSA_PA14Yqf,2);
BC_BCSA_PA14Yq = mean(BC_BCSA_PA14Yqf,2);
SA_BCSA_PA14Yq = mean(SA_BCSA_PA14Yqf,2);
AB_BCAB_PA14Yq = mean(AB_BCAB_PA14Yqf,2);
BC_BCAB_PA14Yq = mean(BC_BCAB_PA14Yqf,2);

% Full multi-species community (target_treatment_PAstrain, N=no phage)
AB_MC_PA14Nqf = parsedata(data_ms10{10}(1:25,impcolsN));
BC_MC_PA14Nqf = parsedata(data_ms10{10}(26:50,impcolsN));
SA_MC_PA14Nqf = parsedata(data_ms10{10}(51:75,impcolsN));

AB_MC_PA14Yqf = parsedata(data_ms10{10}(1:49,impcolsY));
BC_MC_PA14Yqf = parsedata(data_ms10{10}(50:98,impcolsY));
SA_MC_PA14Yqf = parsedata(data_ms10{10}(99:147,impcolsY));

% Average replicates (quantities final)
AB_MC_PA14Nq = mean(AB_MC_PA14Nqf,2);
BC_MC_PA14Nq = mean(BC_MC_PA14Nqf,2);
SA_MC_PA14Nq = mean(SA_MC_PA14Nqf,2);

AB_MC_PA14Yq = mean(AB_MC_PA14Yqf,2);
BC_MC_PA14Yq = mean(BC_MC_PA14Yqf,2);
SA_MC_PA14Yq = mean(SA_MC_PA14Yqf,2);


% Strain: CR-KO
% 2-species experiments 
SA_CRKONqf = parsedata(data_ms10{3}(:,impcolsN));
SA_CRKOYqf = parsedata(data_ms10{3}(:,impcolsY));

BC_CRKONqf = parsedata(data_ms10{12}(:,impcolsN));
BC_CRKOYqf = parsedata(data_ms10{12}(:,impcolsY));

AB_CRKONqf = parsedata(data_ms10{14}(:,impcolsN));
AB_CRKOYqf = parsedata(data_ms10{14}(:,impcolsY));

% Average replicates (quantities final)
SA_CRKONq = mean(SA_CRKONqf,2);
AB_CRKONq = mean(AB_CRKONqf,2);
BC_CRKONq = mean(BC_CRKONqf,2);

SA_CRKOYq = mean(SA_CRKOYqf,2);
AB_CRKOYq = mean(AB_CRKOYqf,2);
BC_CRKOYq = mean(BC_CRKOYqf,2);

% 3-species experiments (target_treatment_PAstrain, N=no phage)
AB_ABSA_CRKONqf = parsedata(data_ms10{2}(1:25,impcolsN));
SA_ABSA_CRKONqf = parsedata(data_ms10{2}(26:end,impcolsN));
AB_ABSA_CRKOYqf = parsedata(data_ms10{2}(1:25,impcolsY));
SA_ABSA_CRKOYqf = parsedata(data_ms10{2}(26:end,impcolsY));

BC_BCSA_CRKONqf = parsedata(data_ms10{7}(1:25,impcolsN));
SA_BCSA_CRKONqf = parsedata(data_ms10{7}(26:end,impcolsN));
BC_BCSA_CRKOYqf = parsedata(data_ms10{7}(1:25,impcolsY));
SA_BCSA_CRKOYqf = parsedata(data_ms10{7}(26:end,impcolsY));

AB_BCAB_CRKONqf = parsedata(data_ms10{13}(1:25,impcolsN));
BC_BCAB_CRKONqf = parsedata(data_ms10{13}(26:end,impcolsN));
AB_BCAB_CRKOYqf = parsedata(data_ms10{13}(1:25,impcolsY));
BC_BCAB_CRKOYqf = parsedata(data_ms10{13}(26:end,impcolsY));

% Average replicates (quantities final)
AB_ABSA_CRKONq = mean(AB_ABSA_CRKONqf,2);
SA_ABSA_CRKONq = mean(SA_ABSA_CRKONqf,2);
BC_BCSA_CRKONq = mean(BC_BCSA_CRKONqf,2);
SA_BCSA_CRKONq = mean(SA_BCSA_CRKONqf,2);
AB_BCAB_CRKONq = mean(AB_BCAB_CRKONqf,2);
BC_BCAB_CRKONq = mean(BC_BCAB_CRKONqf,2);

AB_ABSA_CRKOYq = mean(AB_ABSA_CRKOYqf,2);
SA_ABSA_CRKOYq = mean(SA_ABSA_CRKOYqf,2);
BC_BCSA_CRKOYq = mean(BC_BCSA_CRKOYqf,2);
SA_BCSA_CRKOYq = mean(SA_BCSA_CRKOYqf,2);
AB_BCAB_CRKOYq = mean(AB_BCAB_CRKOYqf,2);
BC_BCAB_CRKOYq = mean(BC_BCAB_CRKOYqf,2);

% Full multi-species community (target_treatment_PAstrain, N=no phage)
AB_MC_CRKONqf = parsedata(data_ms10{11}(1:25,impcolsN));
BC_MC_CRKONqf = parsedata(data_ms10{11}(26:50,impcolsN));
SA_MC_CRKONqf = parsedata(data_ms10{11}(51:end,impcolsN));

AB_MC_CRKOYqf = parsedata(data_ms10{11}(1:25,impcolsY));
BC_MC_CRKOYqf = parsedata(data_ms10{11}(26:50,impcolsY));
SA_MC_CRKOYqf = parsedata(data_ms10{11}(51:end,impcolsY));

% Average replicates (quantities final)
AB_MC_CRKONq = mean(AB_MC_CRKONqf,2);
BC_MC_CRKONq = mean(BC_MC_CRKONqf,2);
SA_MC_CRKONq = mean(SA_MC_CRKONqf,2);

AB_MC_CRKOYq = mean(AB_MC_CRKOYqf,2);
BC_MC_CRKOYq = mean(BC_MC_CRKOYqf,2);
SA_MC_CRKOYq = mean(SA_MC_CRKOYqf,2);


%% Read in multi-species data - PA only 

% 10-day multi species growth curves - PA, no phage

% Input data 
% import data from .xlsx file - multiple sheets
sheet_name = sheetnames('10-day exp. PA densities without phage sarah.xlsx');
for k=1:numel(sheet_name)
    data_pa10N{k} = readtable('10-day exp. PA densities without phage sarah.xlsx','Sheet',sheet_name(k),'VariableNamingRule','preserve');
end

% Separate by treatment, target (strain)
impcolsPA14 = 3:6;
impcolsCRKO = 10:13;

% naming: treatmentp_targetN (no phage) - only PA, target is strain 
% Strain: PA14
% 2-species experiments
PA14p_PA14Nqf = parsedata(data_pa10N{1}(:,impcolsPA14));
SAp_PA14Nqf = parsedata(data_pa10N{2}(:,impcolsPA14));
ABp_PA14Nqf = parsedata(data_pa10N{3}(:,impcolsPA14));
BCp_PA14Nqf = parsedata(data_pa10N{4}(:,impcolsPA14));

% Average replicates (quantities final)
PA14p_PA14Nq = mean(PA14p_PA14Nqf,2);
SAp_PA14Nq = mean(SAp_PA14Nqf,2);
ABp_PA14Nq = mean(ABp_PA14Nqf,2);
BCp_PA14Nq = mean(BCp_PA14Nqf,2);

% 3-species experiments
ABSAp_PA14Nqf = parsedata(data_pa10N{5}(:,impcolsPA14));
BCSAp_PA14Nqf = parsedata(data_pa10N{6}(:,impcolsPA14));
BCABp_PA14Nqf = parsedata(data_pa10N{7}(:,impcolsPA14));

% Average replicates (quantities final)
ABSAp_PA14Nq = mean(ABSAp_PA14Nqf,2);
BCSAp_PA14Nq = mean(BCSAp_PA14Nqf,2);
BCABp_PA14Nq = mean(BCABp_PA14Nqf,2);

% Full multi-species community experiment
MCp_PA14Nqf = parsedata(data_pa10N{8}(:,impcolsPA14));

% Average replicates (quantities final)
MCp_PA14Nq = mean(MCp_PA14Nqf,2);


% Strain: CR-KO
% 2-species experiments
PA14p_CRKONqf = parsedata(data_pa10N{1}(:,impcolsCRKO));
SAp_CRKONqf = parsedata(data_pa10N{2}(:,impcolsCRKO));
ABp_CRKONqf = parsedata(data_pa10N{3}(:,impcolsCRKO));
BCp_CRKONqf = parsedata(data_pa10N{4}(:,impcolsCRKO));

% Average replicates (quantities final)
PA14p_CRKONq = mean(PA14p_CRKONqf,2);
SAp_CRKONq = mean(SAp_CRKONqf,2);
ABp_CRKONq = mean(ABp_CRKONqf,2);
BCp_CRKONq = mean(BCp_CRKONqf,2);

% 3-species experiments
ABSAp_CRKONqf = parsedata(data_pa10N{5}(:,impcolsCRKO));
BCSAp_CRKONqf = parsedata(data_pa10N{6}(:,impcolsCRKO));
BCABp_CRKONqf = parsedata(data_pa10N{7}(:,impcolsCRKO));

% Average replicates (quantities final)
ABSAp_CRKONq = mean(ABSAp_CRKONqf,2);
BCSAp_CRKONq = mean(BCSAp_CRKONqf,2);
BCABp_CRKONq = mean(BCABp_CRKONqf,2);

% Full multi-species community experiment
MCp_CRKONqf = parsedata(data_pa10N{8}(:,impcolsCRKO));

% Average replicates (quantities final)
MCp_CRKONq = mean(MCp_CRKONqf,2);


% 10-day multi species growth curves - PA, with phage

% Input data 
% import data from .xlsx file - multiple sheets
sheet_name = sheetnames('10-day exp. PA densities with phage sarah.xlsx');
for k=1:numel(sheet_name)
    data_pa10Y{k} = readtable('10-day exp. PA densities with phage sarah.xlsx','Sheet',sheet_name(k),'VariableNamingRule','preserve');
end

% Separate by treatment, target (strain)
impcolsPA14 = 3:6;
impcolsCRKO = 10:13;

% naming: treatmentp_targetY (with phage) - only PA, target is strain 
% Strain: PA14
PA14p_PA14Yqf = parsedata(data_pa10Y{1}(:,impcolsPA14));
SAp_PA14Yqf = parsedata(data_pa10Y{2}(:,impcolsPA14));
ABp_PA14Yqf = parsedata(data_pa10Y{3}(:,impcolsPA14));
BCp_PA14Yqf = parsedata(data_pa10Y{4}(:,impcolsPA14));

% Average replicates (quantities final)
PA14p_PA14Yq = mean(PA14p_PA14Yqf,2);
SAp_PA14Yq = mean(SAp_PA14Yqf,2);
ABp_PA14Yq = mean(ABp_PA14Yqf,2);
BCp_PA14Yq = mean(BCp_PA14Yqf,2);

% 3-species experiments
ABSAp_PA14Yqf = parsedata(data_pa10Y{5}(:,impcolsPA14));
BCSAp_PA14Yqf = parsedata(data_pa10Y{6}(:,impcolsPA14));
BCABp_PA14Yqf = parsedata(data_pa10Y{7}(:,impcolsPA14));

% Average replicates (quantities final)
ABSAp_PA14Yq = mean(ABSAp_PA14Yqf,2);
BCSAp_PA14Yq = mean(BCSAp_PA14Yqf,2);
BCABp_PA14Yq = mean(BCABp_PA14Yqf,2);

% Full multi-species community experiment
MCp_PA14Yqf = parsedata(data_pa10Y{8}(:,impcolsPA14));

% Average replicates (quantities final)
MCp_PA14Yq = mean(MCp_PA14Yqf,2);


% Strain: CR-KO
PA14p_CRKOYqf = parsedata(data_pa10Y{1}(1:25,impcolsCRKO));
SAp_CRKOYqf = parsedata(data_pa10Y{2}(1:25,impcolsCRKO));
ABp_CRKOYqf = parsedata(data_pa10Y{3}(1:25,impcolsCRKO));
BCp_CRKOYqf = parsedata(data_pa10Y{4}(1:25,impcolsCRKO));

% Average replicates (quantities final)
PA14p_CRKOYq = mean(PA14p_CRKOYqf,2);
SAp_CRKOYq = mean(SAp_CRKOYqf,2);
ABp_CRKOYq = mean(ABp_CRKOYqf,2);
BCp_CRKOYq = mean(BCp_CRKOYqf,2);

% 3-species experiments
ABSAp_CRKOYqf = parsedata(data_pa10Y{5}(1:25,impcolsCRKO));
BCSAp_CRKOYqf = parsedata(data_pa10Y{6}(1:25,impcolsCRKO));
BCABp_CRKOYqf = parsedata(data_pa10Y{7}(1:25,impcolsCRKO));

% Average replicates (quantities final)
ABSAp_CRKOYq = mean(ABSAp_CRKOYqf,2);
BCSAp_CRKOYq = mean(BCSAp_CRKOYqf,2);
BCABp_CRKOYq = mean(BCABp_CRKOYqf,2);

% Full multi-species community experiment
MCp_CRKOYqf = parsedata(data_pa10Y{8}(1:25,impcolsCRKO));

% Average replicates (quantities final)
MCp_CRKOYq = mean(MCp_CRKOYqf,2);


%% Fit two species curves - no PA

% Set up to fit and keep track of parameters
time = [0,1,3,7,10];
tspan = time*24;

r22pw = zeros(2,3);
Asave2pw = zeros(2,3);

% Fit interaction parameters (assuming growth rates from single species dynamics)
% AB + SA 
ABSA_pw = [ABSA_ABpwq ABSA_SApwq];  
pguess = [rsave1(1,2) rsave1(1,4); 0 0; 0 0];        
lb = [rsave1(1,2) rsave1(1,4); -inf -inf; -inf -inf]; 
ub = [rsave1(1,2) rsave1(1,4); 0 0; 0 0];            
[ABSA_pwpred,BetaABSA_pw,AlphaABSA_pw,r221,method1] = multispfit071023(ABSA_pw,tspan,lb,ub,pguess,1);
Asave2pw(:,1) = [AlphaABSA_pw(1,2); AlphaABSA_pw(2,1)];
r22pw(:,1) = r221';

% Save coefficients
twosp_coeffs = makecoefftable(0,2,BetaABSA_pw,AlphaABSA_pw,r221','AB+SA','N');

% Make density vs. time figure
%cfucurvefig(tspan,ABSA_pw,ABSA_pwpred,{'AB','SA'},{'A. baumannii + S. aureus'})


% BC + SA  
BCSA_pw = [BCSA_BCpwq BCSA_SApwq];  
pguess = [rsave1(1,3) rsave1(1,4); 0 0; 0 0];        
lb = [rsave1(1,3) rsave1(1,4); -inf -inf; -inf -inf]; 
ub = [rsave1(1,3) rsave1(1,4); 0 0; 0 0];            
[BCSA_pwpred,BetaBCSA_pw,AlphaBCSA_pw,r222,method2] = multispfit071023(BCSA_pw,tspan,lb,ub,pguess,1); 
Asave2pw(:,2) = [AlphaBCSA_pw(1,2); AlphaBCSA_pw(2,1)];
r22pw(:,2) = r222';

% Save coefficients
twosp_coeffs = makecoefftable(twosp_coeffs,2,BetaBCSA_pw,AlphaBCSA_pw,r222','BC+SA','N');

% Make density vs. time figure
%cfucurvefig(tspan,BCSA_pw,BCSA_pwpred,{'BC','SA'},{'B. cenocepacia + S. aureus'})


% BC + AB
BCAB_pw = [BCAB_BCpwq([1,2,4,5],:) BCAB_ABpwq([1,2,4,5],:)];   
pguess = [rsave1(1,3) rsave1(1,2); 0 0; 0 0];        
lb = [rsave1(1,3) rsave1(1,2); -inf -inf; -inf -inf]; 
ub = [rsave1(1,3) rsave1(1,2); 0 0; 0 0];            
[BCAB_pwpred,BetaBCAB_pw,AlphaBCAB_pw,r223,method3] = multispfit071023(BCAB_pw,tspan(:,[1,2,4,5]),lb,ub,pguess,1); 
Asave2pw(:,3) = [AlphaBCAB_pw(1,2); AlphaBCAB_pw(2,1)];
r22pw(:,3) = r223';

% Save coefficients
twosp_coeffs = makecoefftable(twosp_coeffs,2,BetaBCAB_pw,AlphaBCAB_pw,r223','BC+AB','N');

% Make density vs. time figure
%cfucurvefig(tspan(:,[1,2,4,5]),BCAB_pw,BCAB_pwpred,{'BC','AB'},{'B. cenocepacia + A. baumannii'})


% All 3 - AB, BC, SA
a11AB = mean([AlphaABSA_pw(1,1),AlphaBCAB_pw(2,2)]);
a11BC = mean([AlphaBCSA_pw(1,1),AlphaBCAB_pw(1,1)]);
a11SA = mean([AlphaABSA_pw(2,2),AlphaBCSA_pw(2,2)]);

ABBCSA_pw = [ABBCSA_ABpwq ABBCSA_BCpwq ABBCSA_SApwq];

% Fit model
pguess = [rsave1(1,2) rsave1(1,3) rsave1(1,4); a11AB AlphaBCAB_pw(2,1) AlphaABSA_pw(1,2); AlphaBCAB_pw(1,2) a11BC AlphaBCSA_pw(1,2); AlphaABSA_pw(2,1) AlphaBCSA_pw(2,1) a11SA];
lb = [rsave1(1,2) rsave1(1,3) rsave1(1,4); -inf -inf -inf; -inf -inf -inf; -inf -inf -inf];
ub = [rsave1(1,2) rsave1(1,3) rsave1(1,4); 0 0 0; 0 0 0; 0 0 0];
[ABBCSA_pwpred,BetaABBCSA_pw,AlphaABBCSA_pw,r224,method4] = multispfit071023(ABBCSA_pw,tspan,lb,ub,pguess,1); 

% Save coefficients
threesp_coeffs = makecoefftable(0,3,BetaABBCSA_pw,AlphaABBCSA_pw,r224','AB+BC+SA','N');

% Make density vs. time figure
%cfucurvefig(tspan,ABBCSA_pw,ABBCSA_pwpred,{'AB','BC','SA'},{'A. baumannii + B. cenocepacia + S. aureus'})


%% Fit two species curves - PA14 

% Set up to fit and keep track of parameters
time = [0,1,3,7,10];
tspan = time*24;

r22 = zeros(2,3);
r22P = zeros(2,3);
Asave2 = zeros(2,3);
Asave2P = Asave2;


% Fit interaction parameters (assuming growth rates from single species dynamics)
% PA14 + AB w/o phage
PA14_AB_N = [ABp_PA14Nq AB_PA14Nq];    
pguess = [rsave1(1,1) rsave1(1,2); 0 0; 0 0];        
lb = [rsave1(1,1) rsave1(1,2); -inf -inf; -inf -inf]; 
ub = [rsave1(1,1) rsave1(1,2); 0 0; 0 0];            
[PA14_AB_Npred,BetaAB,AlphaAB,r221,method1] = multispfit071023(PA14_AB_N,tspan,lb,ub,pguess,1); 
Asave2(:,1) = [AlphaAB(1,2); AlphaAB(2,1)];
r22(:,1) = r221';

% Save coefficients
twosp_coeffs = makecoefftable(twosp_coeffs,2,BetaAB,AlphaAB,r221','PA14+AB','N');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_AB_N,PA14_AB_Npred,{'PA14','AB'},{'+ A. baumannii'})

% PA14 + AB w/ phage
PA14_AB_Y = [ABp_PA14Yq AB_PA14Yq];
[PA14_AB_Ypred,BetaABP,AlphaABP,r221P,method1P] = multispfit071023(PA14_AB_Y,tspan,lb,ub,pguess,1); 
Asave2P(:,1) = [AlphaABP(1,2); AlphaABP(2,1)];
r22P(:,1) = r221P';

% Save coefficients
twosp_coeffs = makecoefftable(twosp_coeffs,2,BetaABP,AlphaABP,r221P','PA14+AB','Y');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_AB_Y,PA14_AB_Ypred,{'PA14','AB'},{'+ A. baumannii'})
%cfucurvefig(tspan,PA14_AB_Y(:,1:2),PA14_AB_Ypred(:,1:2),{'PA14','AB'},{'+ A. baumannii'})

% Make observed vs. predicted figure
%PA14_AB_Nc = {PA14_AB_N,PA14_AB_Npred};
%PA14_AB_Yc = {PA14_AB_Y,PA14_AB_Ypred};
%predvobsfig(PA14_AB_Nc,PA14_AB_Yc,[],[],{'PA14','AB'},{'+ A. baumannii'})


% PA14 + BC w/o phage
PA14_BC_N = [BCp_PA14Nq BC_PA14Nq];      
pguess = [rsave1(1,1) rsave1(1,3); 0 0; 0 0];            
lb = [rsave1(1,1) rsave1(1,3); -inf -inf; -inf -inf];    
ub = [rsave1(1,1) rsave1(1,3); 0 0; 0 0];               
[PA14_BC_Npred,BetaBC,AlphaBC,r222,method2] = multispfit071023(PA14_BC_N,tspan,lb,ub,pguess,1); 
Asave2(:,2) = [AlphaBC(1,2); AlphaBC(2,1)];
r22(:,2) = r222';

% Save coefficients
twosp_coeffs = makecoefftable(twosp_coeffs,2,BetaBC,AlphaBC,r222','PA14+BC','N');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_BC_N,PA14_BC_Npred,{'PA14','BC'},{'+ B. cenocepacia'})

% PA + BC w/ phage
PA14_BC_Y = [BCp_PA14Yq BC_PA14Yq];     
[PA14_BC_Ypred,BetaBCP,AlphaBCP,r222P,method2P] = multispfit071023(PA14_BC_Y,tspan,lb,ub,pguess,1);
Asave2P(:,2) = [AlphaBCP(1,2); AlphaBCP(2,1)];
r22P(:,2) = r222P';

% Save coefficients
twosp_coeffs = makecoefftable(twosp_coeffs,2,BetaBCP,AlphaBCP,r222P','PA14+BC','Y');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_BC_Y,PA14_BC_Ypred,{'PA14','BC'},{'+ B. cenocepacia'})
%cfucurvefig(tspan,PA14_BC_Y(:,1:2),PA14_BC_Ypred(:,1:2),{'PA14','BC'},{'+ B. cenocepacia'})

% Make observed vs. predicted figure
%PA14_BC_Nc = {PA14_BC_N,PA14_BC_Npred};
%PA14_BC_Yc = {PA14_BC_Y,PA14_BC_Ypred};
%predvobsfig(PA14_BC_Nc,PA14_BC_Yc,[],[],{'PA14','BC'},{'+ B. cenocepacia'})


% PA14 + SA w/o phage
PA14_SA_N = [SAp_PA14Nq SA_PA14Nq];    
pguess = [rsave1(1,1) rsave1(1,4); 0 0; 0 0];            
lb = [rsave1(1,1) rsave1(1,4); -inf -inf; -inf -inf];   
ub = [rsave1(1,1) rsave1(1,4); 0 0; 0 0];               
[PA14_SA_Npred,BetaSA,AlphaSA,r223,method3] = multispfit071023(PA14_SA_N,tspan,lb,ub,pguess,1); 
Asave2(:,3) = [AlphaSA(1,2); AlphaSA(2,1)];
r22(:,3) = r223';

% Save coefficients
twosp_coeffs = makecoefftable(twosp_coeffs,2,BetaSA,AlphaSA,r223','PA14+SA','N');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_SA_N,PA14_SA_Npred,{'PA14','SA'},{'+ S. aureus'})

% PA14 + SA w/ phage
PA14_SA_Y = [SAp_PA14Yq SA_PA14Yq]; 
[PA14_SA_Ypred,BetaSAP,AlphaSAP,r223P,method3P] = multispfit071023(PA14_SA_Y,tspan,lb,ub,pguess,1);  
Asave2P(:,3) = [AlphaSAP(1,2); AlphaSAP(2,1)];
r22P(:,3) = r223P';

% Save coefficients
twosp_coeffs = makecoefftable(twosp_coeffs,2,BetaSAP,AlphaSAP,r223P','PA14+SA','Y');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_SA_Y,PA14_SA_Ypred,{'PA14','SA'},{'+ S. aureus'})
%cfucurvefig(tspan,PA14_SA_Y(:,1:2),PA14_SA_Ypred(:,1:2),{'PA14','SA'},{'+ S. aureus'})

% Make observed vs. predicted figure
%PA14_SA_Nc = {PA14_SA_N,PA14_SA_Npred};
%PA14_SA_Yc = {PA14_SA_Y,PA14_SA_Ypred};
%predvobsfig(PA14_SA_Nc,PA14_SA_Yc,[],[],{'PA14','SA'},{'+ S. aureus'})


%% Fit two species curves - CRKO 

% Set up to fit and keep track of parameters
time = [0,1,3,7,10];
tspan = time*24;

r22C = zeros(2,3);
r22PC = zeros(2,3);
Asave2C = zeros(2,3);
Asave2PC = Asave2C;


% Fit interaction parameters (assuming growth rates single species dynamics)
% PA14 + AB w/o phage
PA14_AB_NC = [ABp_CRKONq AB_CRKONq];    
pguess = [rsave1(1,1) rsave1(1,2); 0 0; 0 0];       
lb = [rsave1(1,1) rsave1(1,2); -inf -inf; -inf -inf]; 
ub = [rsave1(1,1) rsave1(1,2); 0 0; 0 0];            
[PA14_AB_NpredC,BetaABC,AlphaABC,r221,method1C] = multispfit071023(PA14_AB_NC,tspan,lb,ub,pguess,1); 
Asave2C(:,1) = [AlphaABC(1,2); AlphaABC(2,1)];
r22C(:,1) = r221';

% Save coefficients
twosp_coeffs = makecoefftable(twosp_coeffs,2,BetaABC,AlphaABC,r221','CRKO+AB','N');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_AB_NC,PA14_AB_NpredC,{'CRKO','AB'},{'+ A. baumannii'})
  
% PA14 + AB w/ phage
PA14_AB_YC = [ABp_CRKOYq AB_CRKOYq];   
[PA14_AB_YpredC,BetaABPC,AlphaABPC,r221P,method1PC] = multispfit071023(PA14_AB_YC,tspan,lb,ub,pguess,1);  
Asave2PC(:,1) = [AlphaABPC(1,2); AlphaABPC(2,1)];
r22PC(:,1) = r221P';

% Save coefficients
twosp_coeffs = makecoefftable(twosp_coeffs,2,BetaABPC,AlphaABPC,r221P','CRKO+AB','Y');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_AB_YC,PA14_AB_YpredC,{'CRKO','AB'},{'+ A. baumannii'})
%cfucurvefig(tspan,PA14_AB_YC(:,1:2),PA14_AB_YpredC(:,1:2),{'CRKO','AB'},{'+ A. baumannii'})

% Make observed vs. predicted figure
%PA14_AB_NcC = {PA14_AB_NC,PA14_AB_NpredC};
%PA14_AB_YcC = {PA14_AB_YC,PA14_AB_YpredC};
%predvobsfig(PA14_AB_NcC,PA14_AB_YcC,[],[],{'CRKO','AB'},{'+ A. baumannii'})


% PA14 + BC w/o phage
PA14_BC_NC = [BCp_CRKONq BC_CRKONq];      
pguess = [rsave1(1,1) rsave1(1,3); 0 0; 0 0];            
lb = [rsave1(1,1) rsave1(1,3); -inf -inf; -inf -inf];    
ub = [rsave1(1,1) rsave1(1,3); 0 0; 0 0];               
[PA14_BC_NpredC,BetaBCC,AlphaBCC,r222,method2C] = multispfit071023(PA14_BC_NC,tspan,lb,ub,pguess,1); 
Asave2C(:,2) = [AlphaBCC(1,2); AlphaBCC(2,1)];
r22C(:,2) = r222';

% Save coefficients
twosp_coeffs = makecoefftable(twosp_coeffs,2,BetaBCC,AlphaBCC,r222','CRKO+BC','N');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_BC_NC,PA14_BC_NpredC,{'CRKO','BC'},{'+ B. cenocepacia'})

% PA + BC w/ phage
PA14_BC_YC = [BCp_CRKOYq BC_CRKOYq];  
[PA14_BC_YpredC,BetaBCPC,AlphaBCPC,r222P,method2PC] = multispfit071023(PA14_BC_YC,tspan,lb,ub,pguess,1);  
Asave2PC(:,2) = [AlphaBCPC(1,2); AlphaBCPC(2,1)];
r22PC(:,2) = r222P';

% Save coefficients
twosp_coeffs = makecoefftable(twosp_coeffs,2,BetaBCPC,AlphaBCPC,r222P','CRKO+BC','Y');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_BC_YC,PA14_BC_YpredC,{'CRKO','BC'},{'+ B. cenocepacia'})
%cfucurvefig(tspan,PA14_BC_YC(:,1:2),PA14_BC_YpredC(:,1:2),{'CRKO','BC'},{'+ B. cenocepacia'})

% Make observed vs. predicted figure
%PA14_BC_NcC = {PA14_BC_NC,PA14_BC_NpredC};
%PA14_BC_YcC = {PA14_BC_YC,PA14_BC_YpredC};
%predvobsfig(PA14_BC_NcC,PA14_BC_YcC,[],[],{'CRKO','BC'},{'+ B. cenocepacia'})


% PA14 + SA w/o phage
PA14_SA_NC = [SAp_CRKONq SA_CRKONq];    
pguess = [rsave1(1,1) rsave1(1,4); 0 0; 0 0];            
lb = [rsave1(1,1) rsave1(1,4); -inf -inf; -inf -inf];   
ub = [rsave1(1,1) rsave1(1,4); 0 0; 0 0];               
[PA14_SA_NpredC,BetaSAC,AlphaSAC,r223,method3C] = multispfit071023(PA14_SA_NC,tspan,lb,ub,pguess,1); 
Asave2C(:,3) = [AlphaSAC(1,2); AlphaSAC(2,1)];
r22C(:,3) = r223';

% Save coefficients
twosp_coeffs = makecoefftable(twosp_coeffs,2,BetaSAC,AlphaSAC,r223','CRKO+SA','N');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_SA_NC,PA14_SA_NpredC,{'CRKO','SA'},{'+ S. aureus'})

% PA14 + SA w/ phage
PA14_SA_YC = [SAp_CRKOYq SA_CRKOYq];    
[PA14_SA_YpredC,BetaSAPC,AlphaSAPC,r223P,method3PC] = multispfit071023(PA14_SA_YC,tspan,lb,ub,pguess,1);
Asave2PC(:,3) = [AlphaSAPC(1,2); AlphaSAPC(2,1)];
r22PC(:,3) = r223P';

% Save coefficients
twosp_coeffs = makecoefftable(twosp_coeffs,2,BetaSAPC,AlphaSAPC,r223P','CRKO+SA','Y');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_SA_YC,PA14_SA_YpredC,{'CRKO','SA'},{'+ S. aureus'})
%cfucurvefig(tspan,PA14_SA_YC(:,1:2),PA14_SA_YpredC(:,1:2),{'CRKO','SA'},{'+ S. aureus'})

% Make observed vs. predicted figure
%PA14_SA_NcC = {PA14_SA_NC,PA14_SA_NpredC};
%PA14_SA_YcC = {PA14_SA_YC,PA14_SA_YpredC};
%predvobsfig(PA14_SA_NcC,PA14_SA_YcC,[],[],{'CRKO','SA'},{'+ S. aureus'})


%% Write file with 2 species predicted densities

% Make density tables 
densAB_table = makedenstable2(time,PA14_AB_Npred,PA14_AB_Ypred,PA14_AB_NpredC,PA14_AB_YpredC,{'+AB'},{'AB'});
densBC_table = makedenstable2(time,PA14_BC_Npred,PA14_BC_Ypred,PA14_BC_NpredC,PA14_BC_YpredC,{'+BC'},{'SA'});
densSA_table = makedenstable2(time,PA14_SA_Npred,PA14_SA_Ypred,PA14_SA_NpredC,PA14_SA_YpredC,{'+SA'},{'BC'});

%filename_dens = 'preddensities_090323.xlsx';
%writetable(densAB_table,filename_dens,'Sheet','Two Species','Range','A1')
%writetable(densBC_table,filename_dens,'Sheet','Two Species','WriteMode','Append')
%writetable(densSA_table,filename_dens,'Sheet','Two Species','WriteMode','Append')


%% Fit three species curves - PA14 

% set a11 guess (PA-PA interaction coeff)
a11_ABSA = AlphaSA(1,1);
a11_BCSA = AlphaSA(1,1);
a11_BCAB = AlphaBC(1,1);
a11_ABSAP = AlphaSAP(1,1);
a11_BCSAP = AlphaSAP(1,1);
a11_BCABP = AlphaBCP(1,1);

% Fit interaction parameters (assuming growth rates from single species dynamics)
% PA14 + AB + SA (no phage)
PA14_ABSA_N = [ABSAp_PA14Nq AB_ABSA_PA14Nq SA_ABSA_PA14Nq]; 
pguess = [rsave1(1,1) rsave1(1,2) rsave1(1,4); a11_ABSA AlphaAB(1,2) AlphaSA(1,2); AlphaAB(2,1) AlphaAB(2,2) AlphaABSA_pw(1,2); AlphaSA(2,1) AlphaABSA_pw(2,1) AlphaSA(2,2)];
lb = [rsave1(1,1) rsave1(1,2) rsave1(1,4); -inf -inf -inf; -inf -inf -inf; -inf -inf -inf];    
ub = [rsave1(1,1) rsave1(1,2) rsave1(1,4); 0 0 0; 0 0 0; 0 0 0]; 
[PA14_ABSA_Npred,BetaABSA,AlphaABSA,r2_ABSA,methodABSA] = multispfit071023(PA14_ABSA_N,tspan,lb,ub,pguess,1);  

% Save coefficients
threesp_coeffs = makecoefftable(threesp_coeffs,3,BetaABSA,AlphaABSA,r2_ABSA','PA14+AB+SA','N');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_ABSA_N,PA14_ABSA_Npred,{'PA14','AB','SA'},{'+ A. baumannii and S. aureus'})

% PA14 + AB + SA (with phage)
PA14_ABSA_Y = [ABSAp_PA14Yq AB_ABSA_PA14Yq SA_ABSA_PA14Yq]; 
pguess = [rsave1(1,1) rsave1(1,2) rsave1(1,4); a11_ABSAP AlphaABP(1,2) AlphaSAP(1,2); AlphaABP(2,1) AlphaABP(2,2) AlphaABSA_pw(1,2); AlphaSAP(2,1) AlphaABSA_pw(2,1) AlphaSAP(2,2)]; 
lb = [rsave1(1,1) rsave1(1,2) rsave1(1,4); -inf -inf -inf; -inf -inf -inf; -inf -inf -inf];    
ub = [rsave1(1,1) rsave1(1,2) rsave1(1,4); 0 0 0; 0 0 0; 0 0 0]; 
[PA14_ABSA_Ypred,BetaABSAP,AlphaABSAP,r2_ABSAP,methodABSAP] = multispfit071023(PA14_ABSA_Y,tspan,lb,ub,pguess,1);  

% Save coefficients
threesp_coeffs = makecoefftable(threesp_coeffs,3,BetaABSAP,AlphaABSAP,r2_ABSAP','PA14+AB+SA','Y');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_ABSA_Y,PA14_ABSA_Ypred,{'PA14','AB','SA'},{'+ A. baumannii and S. aureus'})

% Make observed vs. predicted figure
%PA14_ABSA_Nc = {PA14_ABSA_N,PA14_ABSA_Npred};
%PA14_ABSA_Yc = {PA14_ABSA_Y,PA14_ABSA_Ypred};
%predvobsfig(PA14_ABSA_Nc,PA14_ABSA_Yc,[],[],{'PA14','AB','SA'},{'+ A. baumannii and S. aureus'})


% Fit interaction parameters (assuming growth rates from single species dynamics)
% PA14 + BC + SA (no phage)
PA14_BCSA_N = [BCSAp_PA14Nq BC_BCSA_PA14Nq SA_BCSA_PA14Nq]; 
pguess = [rsave1(1,1) rsave1(1,3) rsave1(1,4); a11_BCSA AlphaBC(1,2) AlphaSA(1,2); AlphaBC(2,1) AlphaBC(2,2) AlphaBCSA_pw(1,2); AlphaSA(2,1) AlphaBCSA_pw(2,1) AlphaSA(2,2)];
lb = [rsave1(1,1) rsave1(1,3) rsave1(1,4); -inf -inf -inf; -inf -inf -inf; -inf -inf -inf];    
ub = [rsave1(1,1) rsave1(1,3) rsave1(1,4); 0 0 0; 0 0 0; 0 0 0]; 
[PA14_BCSA_Npred,BetaBCSA,AlphaBCSA,r2_BCSA,methodBCSA] = multispfit071023(PA14_BCSA_N,tspan,lb,ub,pguess,1); 

% Save coefficients
threesp_coeffs = makecoefftable(threesp_coeffs,3,BetaBCSA,AlphaBCSA,r2_BCSA','PA14+BC+SA','N');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_BCSA_N,PA14_BCSA_Npred,{'PA14','BC','SA'},{'+ B. cenocepacia and S. aureus'})

% PA14 + BC + SA (with phage)
PA14_BCSA_Y = [BCSAp_PA14Yq BC_BCSA_PA14Yq SA_BCSA_PA14Yq]; 
pguess = [rsave1(1,1) rsave1(1,3) rsave1(1,4); a11_BCSAP AlphaBCP(1,2) AlphaSAP(1,2); AlphaBCP(2,1) AlphaBCP(2,2) AlphaBCSA_pw(1,2); AlphaSAP(2,1) AlphaBCSA_pw(2,1) AlphaSAP(2,2)];
lb = [rsave1(1,1) rsave1(1,3) rsave1(1,4); -inf -inf -inf; -inf -inf -inf; -inf -inf -inf];   
ub = [rsave1(1,1) rsave1(1,3) rsave1(1,4); 0 0 0; 0 0 0; 0 0 0];  
[PA14_BCSA_Ypred,BetaBCSAP,AlphaBCSAP,r2_BCSAP,methodBCSAP] = multispfit071023(PA14_BCSA_Y,tspan,lb,ub,pguess,1); 

% Save coefficients
threesp_coeffs = makecoefftable(threesp_coeffs,3,BetaBCSAP,AlphaBCSAP,r2_BCSAP','PA14+BC+SA','Y');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_BCSA_Y,PA14_BCSA_Ypred,{'PA14','BC','SA'},{'+ B. cenocepacia and S. aureus'})

% Make observed vs. predicted figure
%PA14_BCSA_Nc = {PA14_BCSA_N,PA14_BCSA_Npred};
%PA14_BCSA_Yc = {PA14_BCSA_Y,PA14_BCSA_Ypred};
%predvobsfig(PA14_BCSA_Nc,PA14_BCSA_Yc,[],[],{'PA14','BC','SA'},{'+ B. cenocepacia and S. aureus'})


% Fit interaction parameters (assuming growth rates from single species dynamics)
% PA14 + BC + AB (no phage)
PA14_BCAB_N = [BCABp_PA14Nq BC_BCAB_PA14Nq AB_BCAB_PA14Nq];  
pguess = [rsave1(1,1) rsave1(1,3) rsave1(1,2); a11_BCAB AlphaBC(1,2) AlphaAB(1,2); AlphaBC(2,1) AlphaBC(2,2) AlphaBCAB_pw(1,2); AlphaAB(2,1) AlphaBCAB_pw(2,1) AlphaAB(2,2)];
lb = [rsave1(1,1) rsave1(1,3) rsave1(1,2); -inf -inf -inf; -inf -inf -inf; -inf -inf -inf];    
ub = [rsave1(1,1) rsave1(1,3) rsave1(1,2); 0 0 0; 0 0 0; 0 0 0]; 
[PA14_BCAB_Npred,BetaBCAB,AlphaBCAB,r2_BCAB,methodBCAB] = multispfit071023(PA14_BCAB_N,tspan,lb,ub,pguess,1); 

% Save coefficients
threesp_coeffs = makecoefftable(threesp_coeffs,3,BetaBCAB,AlphaBCAB,r2_BCAB','PA14+BC+AB','N');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_BCAB_N,PA14_BCAB_Npred,{'PA14','BC','AB'},{'+ B. cenocepacia and A. baumannii'})

% PA14 + BC + AB (with phage)
PA14_BCAB_Y = [(BCABp_PA14Yq) (BC_BCAB_PA14Yq) (AB_BCAB_PA14Yq)];  
pguess = [rsave1(1,1) rsave1(1,3) rsave1(1,2); a11_BCABP AlphaBCP(1,2) AlphaABP(1,2); AlphaBCP(2,1) AlphaBCP(2,2) AlphaBCAB_pw(1,2); AlphaABP(2,1) AlphaBCAB_pw(2,1) AlphaABP(2,2)];
lb = [rsave1(1,1) rsave1(1,3) rsave1(1,2); -inf -inf -inf; -inf -inf -inf; -inf -inf -inf];    
ub = [rsave1(1,1) rsave1(1,3) rsave1(1,2); 0 0 0; 0 0 0; 0 0 0]; 
[PA14_BCAB_Ypred,BetaBCABP,AlphaBCABP,r2_BCABP,methodBCABP] = multispfit071023(PA14_BCAB_Y,tspan,lb,ub,pguess,1); 

% Save coefficients
threesp_coeffs = makecoefftable(threesp_coeffs,3,BetaBCABP,AlphaBCABP,r2_BCABP','PA14+BC+AB','Y');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_BCAB_Y,PA14_BCAB_Ypred,{'PA14','BC','AB'},{'+ B. cenocepacia and A. baumannii'})

% Make observed vs. predicted figure
%PA14_BCAB_Nc = {PA14_BCAB_N,PA14_BCAB_Npred};
%PA14_BCAB_Yc = {PA14_BCAB_Y,PA14_BCAB_Ypred};
%predvobsfig(PA14_BCAB_Nc,PA14_BCAB_Yc,[],[],{'PA14','BC','AB'},{'+ B. cenocepacia and A. baumannii'})


%% Fit three species curves - CRKO

% Set up to fit and keep track of parameters
time = [0,1,3,7,10];
tspan = time*24;

% set a11 guess (PA-PA interaction coeff)
a11_ABSAC = AlphaSAC(1,1);
a11_BCSAC = AlphaSAC(1,1);
a11_BCABC = AlphaBCC(1,1);
a11_ABSAPC = AlphaSAPC(1,1);
a11_BCSAPC = AlphaSAPC(1,1);
a11_BCABPC = AlphaBCPC(1,1);

% Fit interaction parameters (assuming growth rates from single species dynamics)
% PA14 + AB + SA (no phage)
PA14_ABSA_NC = [ABSAp_CRKONq AB_ABSA_CRKONq SA_ABSA_CRKONq]; 
pguess = [rsave1(1,1) rsave1(1,2) rsave1(1,4); a11_ABSAC AlphaABC(1,2) AlphaSAC(1,2); AlphaABC(2,1) AlphaABC(2,2) AlphaABSA_pw(1,2); AlphaSAC(2,1) AlphaABSA_pw(2,1) AlphaSAC(2,2)]; 
lb = [rsave1(1,1) rsave1(1,2) rsave1(1,4); -inf -inf -inf; -inf -inf -inf; -inf -inf -inf];    
ub = [rsave1(1,1) rsave1(1,2) rsave1(1,4); 0 0 0; 0 0 0; 0 0 0]; 
[PA14_ABSA_NpredC,BetaABSAC,AlphaABSAC,r2_ABSAC,methodABSAC] = multispfit071023(PA14_ABSA_NC,tspan,lb,ub,pguess,1);  

% Save coefficients
threesp_coeffs = makecoefftable(threesp_coeffs,3,BetaABSAC,AlphaABSAC,r2_ABSAC','CRKO+AB+SA','N');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_ABSA_NC,PA14_ABSA_NpredC,{'CRKO','AB','SA'},{'+ A. baumannii and S. aureus'})

% PA14 + AB + SA (with phage)
PA14_ABSA_YC = [ABSAp_CRKOYq AB_ABSA_CRKOYq SA_ABSA_CRKOYq]; 
pguess = [rsave1(1,1) rsave1(1,2) rsave1(1,4); a11_ABSAPC AlphaABPC(1,2) AlphaSAPC(1,2); AlphaABPC(2,1) AlphaABPC(2,2) AlphaABSA_pw(1,2); AlphaSAPC(2,1) AlphaABSA_pw(2,1) AlphaSAPC(2,2)];
lb = [rsave1(1,1) rsave1(1,2) rsave1(1,4); -inf -inf -inf; -inf -inf -inf; -inf -inf -inf];    
ub = [rsave1(1,1) rsave1(1,2) rsave1(1,4); 0 0 0; 0 0 0; 0 0 0];
[PA14_ABSA_YpredC,BetaABSAPC,AlphaABSAPC,r2_ABSAPC,methodABSAPC] = multispfit071023(PA14_ABSA_YC,tspan,lb,ub,pguess,1);

% Save coefficients
threesp_coeffs = makecoefftable(threesp_coeffs,3,BetaABSAPC,AlphaABSAPC,r2_ABSAPC','CRKO+AB+SA','Y');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_ABSA_YC,PA14_ABSA_YpredC,{'CRKO','AB','SA'},{'+ A. baumannii and S. aureus'})

% Make observed vs. predicted figure
%PA14_ABSA_NcC = {PA14_ABSA_NC,PA14_ABSA_NpredC};
%PA14_ABSA_YcC = {PA14_ABSA_YC,PA14_ABSA_YpredC};
%predvobsfig(PA14_ABSA_NcC,PA14_ABSA_YcC,[],[],{'CRKO','AB','SA'},{'+ A. baumannii and S. aureus'})


% Fit interaction parameters (assuming growth rates from single species dynamics)
% PA14 + BC + SA (no phage)
PA14_BCSA_NC = [BCSAp_CRKONq BC_BCSA_CRKONq SA_BCSA_CRKONq]; 
pguess = [rsave1(1,1) rsave1(1,3) rsave1(1,4); a11_BCSAC AlphaBCC(1,2) AlphaSAC(1,2); AlphaBCC(2,1) AlphaBCC(2,2) AlphaBCSA_pw(1,2); AlphaSAC(2,1) AlphaBCSA_pw(2,1) AlphaSAC(2,2)];
lb = [rsave1(1,1) rsave1(1,3) rsave1(1,4); -inf -inf -inf; -inf -inf -inf; -inf -inf -inf];    
ub = [rsave1(1,1) rsave1(1,3) rsave1(1,4); 0 0 0; 0 0 0; 0 0 0]; 
[PA14_BCSA_NpredC,BetaBCSAC,AlphaBCSAC,r2_BCSAC,methodBCSAC] = multispfit071023(PA14_BCSA_NC,tspan,lb,ub,pguess,1); 

% Save coefficients
threesp_coeffs = makecoefftable(threesp_coeffs,3,BetaBCSAC,AlphaBCSAC,r2_BCSAC','CRKO+BC+SA','N');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_BCSA_NC,PA14_BCSA_NpredC,{'CRKO','BC','SA'},{'+ B. cenocepacia and S. aureus'})

% PA14 + BC + SA (with phage)
PA14_BCSA_YC = [BCSAp_CRKOYq BC_BCSA_CRKOYq SA_BCSA_CRKOYq]; 
pguess = [rsave1(1,1) rsave1(1,3) rsave1(1,4); a11_BCSAPC AlphaBCPC(1,2) AlphaSAPC(1,2); AlphaBCPC(2,1) AlphaBCPC(2,2) AlphaBCSA_pw(1,2); AlphaSAPC(2,1) AlphaBCSA_pw(2,1) AlphaSAPC(2,2)];
lb = [rsave1(1,1) rsave1(1,3) rsave1(1,4); -inf -inf -inf; -inf -inf -inf; -inf -inf -inf];    
ub = [rsave1(1,1) rsave1(1,3) rsave1(1,4); 0 0 0; 0 0 0; 0 0 0]; 
[PA14_BCSA_YpredC,BetaBCSAPC,AlphaBCSAPC,r2_BCSAPC,methodBCSAPC] = multispfit071023(PA14_BCSA_YC,tspan,lb,ub,pguess,1); 

% Save coefficients
threesp_coeffs = makecoefftable(threesp_coeffs,3,BetaBCSAPC,AlphaBCSAPC,r2_BCSAPC','CRKO+BC+SA','Y');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_BCSA_YC,PA14_BCSA_YpredC,{'CRKO','BC','SA'},{'+ B. cenocepacia and S. aureus'})

% Make observed vs. predicted figure
%PA14_BCSA_NcC = {PA14_BCSA_NC,PA14_BCSA_NpredC};
%PA14_BCSA_YcC = {PA14_BCSA_YC,PA14_BCSA_YpredC};
%predvobsfig(PA14_BCSA_NcC,PA14_BCSA_YcC,[],[],{'CRKO','BC','SA'},{'+ B. cenocepacia and S. aureus'})


% Fit interaction parameters (assuming growth rates from single species dynamics)
% PA14 + BC + AB (no phage)
PA14_BCAB_NC = [BCABp_CRKONq BC_BCAB_CRKONq AB_BCAB_CRKONq];  
pguess = [rsave1(1,1) rsave1(1,3) rsave1(1,2); a11_BCABC AlphaBCC(1,2) AlphaABC(1,2); AlphaBCC(2,1) AlphaBCC(2,2) AlphaBCAB_pw(1,2); AlphaABC(2,1) AlphaBCAB_pw(2,1) AlphaABC(2,2)];
lb = [rsave1(1,1) rsave1(1,3) rsave1(1,2); -inf -inf -inf; -inf -inf -inf; -inf -inf -inf];    
ub = [rsave1(1,1) rsave1(1,3) rsave1(1,2); 0 0 0; 0 0 0; 0 0 0]; 
[PA14_BCAB_NpredC,BetaBCABC,AlphaBCABC,r2_BCABC,methodBCABC] = multispfit071023(PA14_BCAB_NC,tspan,lb,ub,pguess,1); 

% Save coefficients
threesp_coeffs = makecoefftable(threesp_coeffs,3,BetaBCABC,AlphaBCABC,r2_BCABC','CRKO+BC+AB','N');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_BCAB_NC,PA14_BCAB_NpredC,{'CRKO','BC','AB'},{'+ B. cenocepacia and A. baumannii'})

% PA14 + BC + AB (with phage)
PA14_BCAB_YC = [(BCABp_CRKOYq) (BC_BCAB_CRKOYq) (AB_BCAB_CRKOYq)];  
pguess = [rsave1(1,1) rsave1(1,3) rsave1(1,2); a11_BCABPC AlphaBCPC(1,2) AlphaABPC(1,2); AlphaBCPC(2,1) AlphaBCPC(2,2) AlphaBCAB_pw(1,2); AlphaABPC(2,1) AlphaBCAB_pw(2,1) AlphaABPC(2,2)];
lb = [rsave1(1,1) rsave1(1,3) rsave1(1,2); -inf -inf -inf; -inf -inf -inf; -inf -inf -inf];    
ub = [rsave1(1,1) rsave1(1,3) rsave1(1,2); 0 0 0; 0 0 0; 0 0 0]; 
[PA14_BCAB_YpredC,BetaBCABPC,AlphaBCABPC,r2_BCABPC,methodBCABPC] = multispfit071023(PA14_BCAB_YC,tspan,lb,ub,pguess,1); 

% Save coefficients
threesp_coeffs = makecoefftable(threesp_coeffs,3,BetaBCABPC,AlphaBCABPC,r2_BCABPC','CRKO+BC+AB','Y');

% Make density vs. time figure
%cfucurvefig(tspan,PA14_BCAB_YC,PA14_BCAB_YpredC,{'CRKO','BC','AB'},{'+ B. cenocepacia and A. baumannii'})

% Make observed vs. predicted figure
%PA14_BCAB_NcC = {PA14_BCAB_NC,PA14_BCAB_NpredC};
%PA14_BCAB_YcC = {PA14_BCAB_YC,PA14_BCAB_YpredC};
%predvobsfig(PA14_BCAB_NcC,PA14_BCAB_YcC,[],[],{'CRKO','BC','AB'},{'+ B. cenocepacia and A. baumannii'})


%% Add 3 species predicted densities to new sheet of existing file

% Make density tables 
densABSA_table = makedenstable2(time,PA14_ABSA_Npred,PA14_ABSA_Ypred,PA14_ABSA_NpredC,PA14_ABSA_YpredC,{'+AB+SA'},{'AB','SA'});
densBCSA_table = makedenstable2(time,PA14_BCSA_Npred,PA14_BCSA_Ypred,PA14_BCSA_NpredC,PA14_BCSA_YpredC,{'+BC+SA'},{'BC','SA'});
densBCAB_table = makedenstable2(time,PA14_BCAB_Npred,PA14_BCAB_Ypred,PA14_BCAB_NpredC,PA14_BCAB_YpredC,{'+BC+AB'},{'BC','AB'});

%writetable(densABSA_table,filename_dens,'Sheet','Three Species','Range','A1')
%writetable(densBCSA_table,filename_dens,'Sheet','Three Species','WriteMode','Append')
%writetable(densBCAB_table,filename_dens,'Sheet','Three Species','WriteMode','Append')


%% Save the best 3 species values for beta --> 4 sp prediction (PA14)

% no phage - for prediction, try different combos of beta parameters
%{
% (1) all 3 species combos with PA
b11 = mean([BetaABSA(1,1),BetaBCSA(1,1),BetaBCAB(1,1)]);
b12 = mean([BetaABSA(1,2),BetaBCAB(1,3)]);
b13 = mean([BetaBCSA(1,2),BetaBCAB(1,2)]);
b14 = mean([BetaABSA(1,3),BetaBCSA(1,3)]);

b21 = mean([BetaABSA(2,1),BetaBCAB(3,1)]);
b22 = mean([BetaABSA(2,2),BetaBCAB(3,3)]);
b23 = BetaBCAB(3,2);
b24 = BetaABSA(2,3);

b31 = mean([BetaBCAB(2,1),BetaBCSA(2,1)]);
b32 = BetaBCAB(2,3);
b33 = mean([BetaBCAB(2,2),BetaBCSA(2,2)]);
b34 = BetaBCSA(2,3);

b41 = mean([BetaABSA(3,1),BetaBCSA(3,1)]);
b42 = BetaABSA(3,2);
b43 = BetaBCSA(3,2);
b44 = mean([BetaABSA(3,3),BetaBCSA(3,3)]);

% (2) all 3 species combos
b11 = mean([BetaABSA(1,1),BetaBCSA(1,1),BetaBCAB(1,1)]);
b12 = mean([BetaABSA(1,2),BetaBCAB(1,3)]);
b13 = mean([BetaBCSA(1,2),BetaBCAB(1,2)]);
b14 = mean([BetaABSA(1,3),BetaBCSA(1,3)]);

b21 = mean([BetaABSA(2,1),BetaBCAB(3,1)]);
b22 = mean([BetaABSA(2,2),BetaBCAB(3,3),BetaABBCSA_pw(1,1)]);
b23 = mean([BetaBCAB(3,2),BetaABBCSA_pw(1,2)]);
b24 = mean([BetaABSA(2,3),BetaABBCSA_pw(1,3)]);

b31 = mean([BetaBCAB(2,1),BetaBCSA(2,1)]);
b32 = mean([BetaBCAB(2,3),BetaABBCSA_pw(2,1)]);
b33 = mean([BetaBCAB(2,2),BetaBCSA(2,2),BetaABBCSA_pw(2,2)]);
b34 = mean([BetaBCSA(2,3),BetaABBCSA_pw(2,3)]);

b41 = mean([BetaABSA(3,1),BetaBCSA(3,1)]);
b42 = mean([BetaABSA(3,2),BetaABBCSA_pw(3,1)]);
b43 = mean([BetaBCSA(3,2),BetaABBCSA_pw(3,2)]);
b44 = mean([BetaABSA(3,3),BetaBCSA(3,3),BetaABBCSA_pw(3,3)]);
%}

% (3) All 2 species combos, all 3 species combos with PA
b11 = mean([BetaAB(1,1),BetaBC(1,1),BetaSA(1,1),BetaABSA(1,1),BetaBCSA(1,1),BetaBCAB(1,1)]);
b12 = mean([BetaAB(1,2),BetaABSA(1,2),BetaBCAB(1,3)]);
b13 = mean([BetaBC(1,2),BetaBCSA(1,2),BetaBCAB(1,2)]);
b14 = mean([BetaSA(1,2),BetaABSA(1,3),BetaBCSA(1,3)]);

b21 = mean([BetaAB(2,1),BetaABSA(2,1),BetaBCAB(3,1)]);
b22 = mean([BetaAB(2,2),BetaABSA(2,2),BetaBCAB(3,3),BetaABSA_pw(1,1),BetaBCAB_pw(2,2)]);
b23 = mean([BetaBCAB_pw(2,1),BetaBCAB(3,2)]);
b24 = mean([BetaABSA_pw(1,2),BetaABSA(2,3)]);

b31 = mean([BetaBC(2,1),BetaBCAB(2,1),BetaBCSA(2,1)]);
b32 = mean([BetaBCAB_pw(1,2),BetaBCAB(2,3)]);
b33 = mean([BetaBC(2,2),BetaBCAB(2,2),BetaBCSA(2,2),BetaBCSA_pw(1,1),BetaBCAB_pw(1,1)]);
b34 = mean([BetaBCSA_pw(1,2),BetaBCSA(2,3)]);

b41 = mean([BetaSA(2,1),BetaABSA(3,1),BetaBCSA(3,1)]);
b42 = mean([BetaABSA_pw(2,1),BetaABSA(3,2)]);
b43 = mean([BetaBCSA_pw(2,1),BetaBCSA(3,2)]);
b44 = mean([BetaSA(2,2),BetaABSA(3,3),BetaBCSA(3,3),BetaABSA_pw(2,2),BetaBCSA_pw(2,2)]);

%{
% (4) All 2 species combos, all 3 species combos
b11 = mean([BetaAB(1,1),BetaBC(1,1),BetaSA(1,1),BetaABSA(1,1),BetaBCSA(1,1),BetaBCAB(1,1)]);
b12 = mean([BetaAB(1,2),BetaABSA(1,2),BetaBCAB(1,3)]);
b13 = mean([BetaBC(1,2),BetaBCSA(1,2),BetaBCAB(1,2)]);
b14 = mean([BetaSA(1,2),BetaABSA(1,3),BetaBCSA(1,3)]);

b21 = mean([BetaAB(2,1),BetaABSA(2,1),BetaBCAB(3,1)]);
b22 = mean([BetaAB(2,2),BetaABSA(2,2),BetaBCAB(3,3),BetaABSA_pw(1,1),BetaBCAB_pw(2,2),BetaABBCSA_pw(1,1)]);
b23 = mean([BetaBCAB_pw(2,1),BetaBCAB(3,2),BetaABBCSA_pw(1,2)]);
b24 = mean([BetaABSA_pw(1,2),BetaABSA(2,3),BetaABBCSA_pw(1,3)]);

b31 = mean([BetaBC(2,1),BetaBCAB(2,1),BetaBCSA(2,1)]);
b32 = mean([BetaBCAB_pw(1,2),BetaBCAB(2,3),BetaABBCSA_pw(2,1)]);
b33 = mean([BetaBC(2,2),BetaBCAB(2,2),BetaBCSA(2,2),BetaBCSA_pw(1,1),BetaBCAB_pw(1,1),BetaABBCSA_pw(2,2)]);
b34 = mean([BetaBCSA_pw(1,2),BetaBCSA(2,3),BetaABBCSA_pw(2,3)]);

b41 = mean([BetaSA(2,1),BetaABSA(3,1),BetaBCSA(3,1)]);
b42 = mean([BetaABSA_pw(2,1),BetaABSA(3,2),BetaABBCSA_pw(3,1)]);
b43 = mean([BetaBCSA_pw(2,1),BetaBCSA(3,2),BetaABBCSA_pw(3,2)]);
b44 = mean([BetaSA(2,2),BetaABSA(3,3),BetaBCSA(3,3),BetaABSA_pw(2,2),BetaBCSA_pw(2,2),BetaABBCSA_pw(3,3)]);

% (5) All 2 species combos with PA, all 3 species combos with PA
b11 = mean([BetaAB(1,1),BetaBC(1,1),BetaSA(1,1),BetaABSA(1,1),BetaBCSA(1,1),BetaBCAB(1,1)]);
b12 = mean([BetaAB(1,2),BetaABSA(1,2),BetaBCAB(1,3)]);
b13 = mean([BetaBC(1,2),BetaBCSA(1,2),BetaBCAB(1,2)]);
b14 = mean([BetaSA(1,2),BetaABSA(1,3),BetaBCSA(1,3)]);

b21 = mean([BetaAB(2,1),BetaABSA(2,1),BetaBCAB(3,1)]);
b22 = mean([BetaAB(2,2),BetaABSA(2,2),BetaBCAB(3,3)]);
b23 = BetaBCAB(3,2);
b24 = BetaABSA(2,3);

b31 = mean([BetaBC(2,1),BetaBCAB(2,1),BetaBCSA(2,1)]);
b32 = BetaBCAB(2,3);
b33 = mean([BetaBC(2,2),BetaBCAB(2,2),BetaBCSA(2,2)]);
b34 = BetaBCSA(2,3);

b41 = mean([BetaSA(2,1),BetaABSA(3,1),BetaBCSA(3,1)]);
b42 = BetaABSA(3,2);
b43 = BetaBCSA(3,2);
b44 = mean([BetaSA(2,2),BetaABSA(3,3),BetaBCSA(3,3)]);

% (6) All 2 species combos with PA, all 3 species combos
b11 = mean([BetaAB(1,1),BetaBC(1,1),BetaSA(1,1),BetaABSA(1,1),BetaBCSA(1,1),BetaBCAB(1,1)]);
b12 = mean([BetaAB(1,2),BetaABSA(1,2),BetaBCAB(1,3)]);
b13 = mean([BetaBC(1,2),BetaBCSA(1,2),BetaBCAB(1,2)]);
b14 = mean([BetaSA(1,2),BetaABSA(1,3),BetaBCSA(1,3)]);

b21 = mean([BetaAB(2,1),BetaABSA(2,1),BetaBCAB(3,1)]);
b22 = mean([BetaAB(2,2),BetaABSA(2,2),BetaBCAB(3,3),BetaABBCSA_pw(1,1)]);
b23 = mean([BetaBCAB(3,2),BetaABBCSA_pw(1,2)]);
b24 = mean([BetaABSA(2,3),BetaABBCSA_pw(1,3)]);

b31 = mean([BetaBC(2,1),BetaBCAB(2,1),BetaBCSA(2,1)]);
b32 = mean([BetaBCAB(2,3),BetaABBCSA_pw(2,1)]);
b33 = mean([BetaBC(2,2),BetaBCAB(2,2),BetaBCSA(2,2),BetaABBCSA_pw(2,2)]);
b34 = mean([BetaBCSA(2,3),BetaABBCSA_pw(2,3)]);

b41 = mean([BetaSA(2,1),BetaABSA(3,1),BetaBCSA(3,1)]);
b42 = mean([BetaABSA(3,2),BetaABBCSA_pw(3,1)]);
b43 = mean([BetaBCSA(3,2),BetaABBCSA_pw(3,2)]);
b44 = mean([BetaSA(2,2),BetaABSA(3,3),BetaBCSA(3,3),BetaABBCSA_pw(3,3)]);
%}

% Predicted community interaction matrix, growth rates from single species growth curves
BetaMC = [b11 b12 b13 b14 rsave1(1,1); 
    b21 b22 b23 b24 rsave1(1,2);
    b31 b32 b33 b34 rsave1(1,3);
    b41 b42 b43 b44 rsave1(1,4)];

% Set up to fit and keep track of parameters
time = [0,1,3,7,10];
tspan = time*24;

% Simulate full community dynamics
PA14_MC_N = [MCp_PA14Nq AB_MC_PA14Nq BC_MC_PA14Nq SA_MC_PA14Nq];  %no phage; 1=PA, 2=AB, 3=BC, 4=SA

% Alpha prediction (interaction coeffs scaled by max single species density for each species)
AlphaMC = BetaMC(1:4,1:4).*[maxPA maxAB maxBC maxSA];

% Predict from pairwise average
PA14_MC_Npred = glv_simulation(tspan,PA14_MC_N(1,:),BetaMC);

% Calculate R-squared
R2save_MC = Rsquare(PA14_MC_N,PA14_MC_Npred);

% Save coefficients
multisp_coeffs = makecoefftable(0,4,BetaMC,AlphaMC,R2save_MC','PA14+all','N');

% Make density vs. time plot
%cfucurvefig(tspan,PA14_MC_N,PA14_MC_Npred,{'PA14','AB','BC','SA'},{'Polyculture'})


% with phage - for prediction
% (1) all 3 species combos with PA (only combos with phage implied)
b11P = mean([BetaABSAP(1,1),BetaBCSAP(1,1),BetaBCABP(1,1)]);
b12P = mean([BetaABSAP(1,2),BetaBCABP(1,3)]);
b13P = mean([BetaBCSAP(1,2),BetaBCABP(1,2)]);
b14P = mean([BetaABSAP(1,3),BetaBCSAP(1,3)]);

b21P = mean([BetaABSAP(2,1),BetaBCABP(3,1)]);
b22P = mean([BetaABSAP(2,2),BetaBCABP(3,3)]);
b23P = BetaBCABP(3,2);
b24P = BetaABSAP(2,3);

b31P = mean([BetaBCABP(2,1),BetaBCSAP(2,1)]);
b32P = BetaBCABP(2,3);
b33P = mean([BetaBCSAP(2,2),BetaBCABP(2,2)]);
b34P = BetaBCSAP(2,3);

b41P = mean([BetaABSAP(3,1),BetaBCSAP(3,1)]);
b42P = BetaABSAP(3,2);
b43P = BetaBCSAP(3,2);
b44P = mean([BetaABSAP(3,3),BetaBCSAP(3,3)]);

% Predicted community interaction matrix, growth rates from single species growth curves
BetaMCP = [b11P b12P b13P b14P rsave1(1,1); 
    b21P b22P b23P b24P rsave1(1,2);
    b31P b32P b33P b34P rsave1(1,3);
    b41P b42P b43P b44P rsave1(1,4)];

% Set up to fit and keep track of parameters
time = [0,1,3,7,10];
tspan = time*24;

% Simulate full community dynamics 
PA14_MC_Y = [MCp_PA14Yq AB_MC_PA14Yq BC_MC_PA14Yq SA_MC_PA14Yq];  %w/ phage; 1=PA, 2=AB, 3=BC, 4=SA

% Alpha prediction (interaction coeffs scaled by max single species density for each species)
AlphaMCP = BetaMCP(1:4,1:4).*[maxPA maxAB maxBC maxSA];

% Predict from pairwise average
PA14_MC_Ypred = glv_simulation(tspan,PA14_MC_Y(1,:),BetaMCP);

% Calculate R-squared
R2save_MCP = Rsquare(PA14_MC_Y,PA14_MC_Ypred);

% Save coefficients
multisp_coeffs = makecoefftable(multisp_coeffs,4,BetaMCP,AlphaMCP,R2save_MCP','PA14+all','Y');

% Make density vs. time plot
%cfucurvefig(tspan,PA14_MC_Y,PA14_MC_Ypred,{'PA14','AB','BC','SA'},{'Polyculture'})

% Make observed vs. predicted figure
%PA14_MC_Nc = {PA14_MC_N,PA14_MC_Npred};
%PA14_MC_Yc = {PA14_MC_Y,PA14_MC_Ypred};
%predvobsfig(PA14_MC_Nc,PA14_MC_Yc,[],[],{'PA14','AB','BC','SA'},{'Polyculture'})


%% Save the best 3 species values for beta --> 4 sp prediction (CRKO)

% no phage - predict, try different combos of beta parameters
% (1) all 3 species combos with PA
b11 = mean([BetaABSAC(1,1),BetaBCSAC(1,1),BetaBCABC(1,1)]);
b12 = mean([BetaABSAC(1,2),BetaBCABC(1,3)]);
b13 = mean([BetaBCSAC(1,2),BetaBCABC(1,2)]);
b14 = mean([BetaABSAC(1,3),BetaBCSAC(1,3)]);

b21 = mean([BetaABSAC(2,1),BetaBCABC(3,1)]);
b22 = mean([BetaABSAC(2,2),BetaBCABC(3,3)]);
b23 = BetaBCABC(3,2);
b24 = BetaABSAC(2,3);

b31 = mean([BetaBCABC(2,1),BetaBCSAC(2,1)]);
b32 = BetaBCABC(2,3);
b33 = mean([BetaBCABC(2,2),BetaBCSAC(2,2)]);
b34 = BetaBCSAC(2,3);

b41 = mean([BetaABSAC(3,1),BetaBCSAC(3,1)]);
b42 = BetaABSAC(3,2);
b43 = BetaBCSAC(3,2);
b44 = mean([BetaABSAC(3,3),BetaBCSAC(3,3)]);

%{
% (2) all 3 species combos
b11 = mean([BetaABSAC(1,1),BetaBCSAC(1,1),BetaBCABC(1,1)]);
b12 = mean([BetaABSAC(1,2),BetaBCABC(1,3)]);
b13 = mean([BetaBCSAC(1,2),BetaBCABC(1,2)]);
b14 = mean([BetaABSAC(1,3),BetaBCSAC(1,3)]);

b21 = mean([BetaABSAC(2,1),BetaBCABC(3,1)]);
b22 = mean([BetaABSAC(2,2),BetaBCABC(3,3),BetaABBCSA_pw(1,1)]);
b23 = mean([BetaBCABC(3,2),BetaABBCSA_pw(1,2)]);
b24 = mean([BetaABSAC(2,3),BetaABBCSA_pw(1,3)]);

b31 = mean([BetaBCABC(2,1),BetaBCSAC(2,1)]);
b32 = mean([BetaBCABC(2,3),BetaABBCSA_pw(2,1)]);
b33 = mean([BetaBCABC(2,2),BetaBCSAC(2,2),BetaABBCSA_pw(2,2)]);
b34 = mean([BetaBCSAC(2,3),BetaABBCSA_pw(2,3)]);

b41 = mean([BetaABSAC(3,1),BetaBCSAC(3,1)]);
b42 = mean([BetaABSAC(3,2),BetaABBCSA_pw(3,1)]);
b43 = mean([BetaBCSAC(3,2),BetaABBCSA_pw(3,2)]);
b44 = mean([BetaABSAC(3,3),BetaBCSAC(3,3),BetaABBCSA_pw(3,3)]);

% (3) All 2 species combos, all 3 species combos with PA
b11 = mean([BetaABC(1,1),BetaBCC(1,1),BetaSAC(1,1),BetaABSAC(1,1),BetaBCSAC(1,1),BetaBCABC(1,1)]);
b12 = mean([BetaABC(1,2),BetaABSAC(1,2),BetaBCABC(1,3)]);
b13 = mean([BetaBCC(1,2),BetaBCSAC(1,2),BetaBCABC(1,2)]);
b14 = mean([BetaSAC(1,2),BetaABSAC(1,3),BetaBCSAC(1,3)]);

b21 = mean([BetaABC(2,1),BetaABSAC(2,1),BetaBCABC(3,1)]);
b22 = mean([BetaABC(2,2),BetaABSAC(2,2),BetaBCABC(3,3),BetaABSA_pw(1,1),BetaBCAB_pw(2,2)]);
b23 = mean([BetaBCAB_pw(2,1),BetaBCABC(3,2)]);
b24 = mean([BetaABSA_pw(1,2),BetaABSAC(2,3)]);

b31 = mean([BetaBCC(2,1),BetaBCABC(2,1),BetaBCSAC(2,1)]);
b32 = mean([BetaBCAB_pw(1,2),BetaBCABC(2,3)]);
b33 = mean([BetaBCC(2,2),BetaBCABC(2,2),BetaBCSAC(2,2),BetaBCSA_pw(1,1),BetaBCAB_pw(1,1)]);
b34 = mean([BetaBCSA_pw(1,2),BetaBCSAC(2,3)]);

b41 = mean([BetaSAC(2,1),BetaABSAC(3,1),BetaBCSAC(3,1)]);
b42 = mean([BetaABSA_pw(2,1),BetaABSAC(3,2)]);
b43 = mean([BetaBCSA_pw(2,1),BetaBCSAC(3,2)]);
b44 = mean([BetaSAC(2,2),BetaABSAC(3,3),BetaBCSAC(3,3),BetaABSA_pw(2,2),BetaBCSA_pw(2,2)]);

% (4) All 2 species combos, all 3 species combos
b11 = mean([BetaABC(1,1),BetaBCC(1,1),BetaSAC(1,1),BetaABSAC(1,1),BetaBCSAC(1,1),BetaBCABC(1,1)]);
b12 = mean([BetaABC(1,2),BetaABSAC(1,2),BetaBCABC(1,3)]);
b13 = mean([BetaBCC(1,2),BetaBCSAC(1,2),BetaBCABC(1,2)]);
b14 = mean([BetaSAC(1,2),BetaABSAC(1,3),BetaBCSAC(1,3)]);

b21 = mean([BetaABC(2,1),BetaABSAC(2,1),BetaBCABC(3,1)]);
b22 = mean([BetaABC(2,2),BetaABSAC(2,2),BetaBCABC(3,3),BetaABSA_pw(1,1),BetaBCAB_pw(2,2),BetaABBCSA_pw(1,1)]);
b23 = mean([BetaBCAB_pw(2,1),BetaBCABC(3,2),BetaABBCSA_pw(1,2)]);
b24 = mean([BetaABSA_pw(1,2),BetaABSAC(2,3),BetaABBCSA_pw(1,3)]);

b31 = mean([BetaBCC(2,1),BetaBCABC(2,1),BetaBCSAC(2,1)]);
b32 = mean([BetaBCAB_pw(1,2),BetaBCABC(2,3),BetaABBCSA_pw(2,1)]);
b33 = mean([BetaBCC(2,2),BetaBCABC(2,2),BetaBCSAC(2,2),BetaBCSA_pw(1,1),BetaBCAB_pw(1,1),BetaABBCSA_pw(2,2)]);
b34 = mean([BetaBCSA_pw(1,2),BetaBCSAC(2,3),BetaABBCSA_pw(2,3)]);

b41 = mean([BetaSAC(2,1),BetaABSAC(3,1),BetaBCSAC(3,1)]);
b42 = mean([BetaABSA_pw(2,1),BetaABSAC(3,2),BetaABBCSA_pw(3,1)]);
b43 = mean([BetaBCSA_pw(2,1),BetaBCSAC(3,2),BetaABBCSA_pw(3,2)]);
b44 = mean([BetaSAC(2,2),BetaABSAC(3,3),BetaBCSAC(3,3),BetaABSA_pw(2,2),BetaBCSA_pw(2,2),BetaABBCSA_pw(3,3)]);

% (5) All 2 species combos with PA, all 3 species combos with PA
b11 = mean([BetaABC(1,1),BetaBCC(1,1),BetaSAC(1,1),BetaABSAC(1,1),BetaBCSAC(1,1),BetaBCABC(1,1)]);
b12 = mean([BetaABC(1,2),BetaABSAC(1,2),BetaBCABC(1,3)]);
b13 = mean([BetaBCC(1,2),BetaBCSAC(1,2),BetaBCABC(1,2)]);
b14 = mean([BetaSAC(1,2),BetaABSAC(1,3),BetaBCSAC(1,3)]);

b21 = mean([BetaABC(2,1),BetaABSAC(2,1),BetaBCABC(3,1)]);
b22 = mean([BetaABC(2,2),BetaABSAC(2,2),BetaBCABC(3,3)]);
b23 = BetaBCABC(3,2);
b24 = BetaABSAC(2,3);

b31 = mean([BetaBCC(2,1),BetaBCABC(2,1),BetaBCSAC(2,1)]);
b32 = BetaBCABC(2,3);
b33 = mean([BetaBCC(2,2),BetaBCABC(2,2),BetaBCSAC(2,2)]);
b34 = BetaBCSAC(2,3);

b41 = mean([BetaSAC(2,1),BetaABSAC(3,1),BetaBCSAC(3,1)]);
b42 = BetaABSAC(3,2);
b43 = BetaBCSAC(3,2);
b44 = mean([BetaSAC(2,2),BetaABSAC(3,3),BetaBCSAC(3,3)]);

% (6) All 2 species combos with PA, all 3 species combos
b11 = mean([BetaABC(1,1),BetaBCC(1,1),BetaSAC(1,1),BetaABSAC(1,1),BetaBCSAC(1,1),BetaBCABC(1,1)]);
b12 = mean([BetaABC(1,2),BetaABSAC(1,2),BetaBCABC(1,3)]);
b13 = mean([BetaBCC(1,2),BetaBCSAC(1,2),BetaBCABC(1,2)]);
b14 = mean([BetaSAC(1,2),BetaABSAC(1,3),BetaBCSAC(1,3)]);

b21 = mean([BetaABC(2,1),BetaABSAC(2,1),BetaBCABC(3,1)]);
b22 = mean([BetaABC(2,2),BetaABSAC(2,2),BetaBCABC(3,3),BetaABBCSA_pw(1,1)]);
b23 = mean([BetaBCABC(3,2),BetaABBCSA_pw(1,2)]);
b24 = mean([BetaABSAC(2,3),BetaABBCSA_pw(1,3)]);

b31 = mean([BetaBCC(2,1),BetaBCABC(2,1),BetaBCSAC(2,1)]);
b32 = mean([BetaBCABC(2,3),BetaABBCSA_pw(2,1)]);
b33 = mean([BetaBCC(2,2),BetaBCABC(2,2),BetaBCSAC(2,2),BetaABBCSA_pw(2,2)]);
b34 = mean([BetaBCSAC(2,3),BetaABBCSA_pw(2,3)]);

b41 = mean([BetaSAC(2,1),BetaABSAC(3,1),BetaBCSAC(3,1)]);
b42 = mean([BetaABSAC(3,2),BetaABBCSA_pw(3,1)]);
b43 = mean([BetaBCSAC(3,2),BetaABBCSA_pw(3,2)]);
b44 = mean([BetaSAC(2,2),BetaABSAC(3,3),BetaBCSAC(3,3),BetaABBCSA_pw(3,3)]);
%}

% Predicted community interaction matrix, growth rates from single species growth curves
BetaMCC = [b11 b12 b13 b14 rsave1(1,1); 
    b21 b22 b23 b24 rsave1(1,2);
    b31 b32 b33 b34 rsave1(1,3);
    b41 b42 b43 b44 rsave1(1,4)];

% Set up to fit and keep track of parameters
time = [0,1,3,7,10];
tspan = time*24;

% Simulate full community dynamics
PA14_MC_NC = [MCp_CRKONq AB_MC_CRKONq BC_MC_CRKONq SA_MC_CRKONq];  %no phage; 1=PA, 2=AB, 3=BC, 4=SA

% Alpha prediction (interaction coeffs scaled by max single species density for each species)
AlphaMCC = BetaMCC(1:4,1:4).*[maxPA maxAB maxBC maxSA];

% Predict from pairwise average
PA14_MC_NpredC = glv_simulation(tspan,PA14_MC_NC(1,:),BetaMCC);

% Calculate R-squared
R2save_MCC = Rsquare(PA14_MC_NC,PA14_MC_NpredC);

% Save coefficients
multisp_coeffs = makecoefftable(multisp_coeffs,4,BetaMCC,AlphaMCC,R2save_MCC','CRKO+all','N');

% Make density vs. time plot
%cfucurvefig(tspan,PA14_MC_NC,PA14_MC_NpredC,{'CRKO','AB','BC','SA'},{'Polyculture'})


% with phage - predict
% (1) all 3 species combos with PA (only combos with phage implied)
b11P = mean([BetaABSAPC(1,1),BetaBCSAPC(1,1),BetaBCABPC(1,1)]);
b12P = mean([BetaABSAPC(1,2),BetaBCABPC(1,3)]);
b13P = mean([BetaBCSAPC(1,2),BetaBCABPC(1,2)]);
b14P = mean([BetaABSAPC(1,3),BetaBCSAPC(1,3)]);

b21P = mean([BetaABSAPC(2,1),BetaBCABPC(3,1)]);
b22P = mean([BetaABSAPC(2,2),BetaBCABPC(3,3)]);
b23P = BetaBCABPC(3,2);
b24P = BetaABSAPC(2,3);

b31P = mean([BetaBCABPC(2,1),BetaBCSAPC(2,1)]);
b32P = BetaBCABPC(2,3);
b33P = mean([BetaBCSAPC(2,2),BetaBCABPC(2,2)]);
b34P = BetaBCSAPC(2,3);

b41P = mean([BetaABSAPC(3,1),BetaBCSAPC(3,1)]);
b42P = BetaABSAPC(3,2);
b43P = BetaBCSAPC(3,2);
b44P = mean([BetaABSAPC(3,3),BetaBCSAPC(3,3)]);

% Predicted community interaction matrix, growth rates from single species growth curves
BetaMCPC = [b11P b12P b13P b14P rsave1(1,1); 
    b21P b22P b23P b24P rsave1(1,2);
    b31P b32P b33P b34P rsave1(1,3);
    b41P b42P b43P b44P rsave1(1,4)];

% Set up to fit and keep track of parameters
time = [0,1,3,7,10];
tspan = time*24;

% Simulate full community dynamics
PA14_MC_YC = [MCp_CRKOYq AB_MC_CRKOYq BC_MC_CRKOYq SA_MC_CRKOYq];  %w/ phage; 1=PA, 2=AB, 3=BC, 4=SA

% Alpha prediction (interaction coeffs scaled by max single species density for each species)
AlphaMCPC = BetaMCPC(1:4,1:4).*[maxPA maxAB maxBC maxSA];

% Predict from pairwsie average
PA14_MC_YpredC = glv_simulation(tspan,PA14_MC_YC(1,:),BetaMCPC);

% Calculate R-squared
R2save_MCPC = Rsquare(PA14_MC_YC,PA14_MC_YpredC);

% Save coefficients
multisp_coeffs = makecoefftable(multisp_coeffs,4,BetaMCPC,AlphaMCPC,R2save_MCPC','CRKO+all','Y');

% Make density vs. time plot
%cfucurvefig(tspan,PA14_MC_YC,PA14_MC_YpredC,{'CRKO','AB','BC','SA'},{'Polyculture'})

% Make observed vs. predicted figure
%PA14_MC_NcC = {PA14_MC_NC,PA14_MC_NpredC};
%PA14_MC_YcC = {PA14_MC_YC,PA14_MC_YpredC};
%predvobsfig(PA14_MC_NcC,PA14_MC_YcC,[],[],{'CRKO','AB','BC','SA'},{'Polyculture'})


%% Add full community species predicted densities to new sheet of existing file

% Make density tables 
densMC_table = makedenstable2(time,PA14_MC_Npred,PA14_MC_Ypred,PA14_MC_NpredC,PA14_MC_YpredC,{'+all'},{'AB','BC','SA'});

%writetable(densMC_table,filename_dens,'Sheet','Full Community','Range','A1')


%% Write a file with all coefficients from model fitting

%filename_coeffs = 'coefficients_090323.xlsx';
%writetable(twosp_coeffs,filename_coeffs,'Sheet','Two Species','Range','A1')
%writetable(threesp_coeffs,filename_coeffs,'Sheet','Three Species','Range','A1')
%writetable(multisp_coeffs,filename_coeffs,'Sheet','Full Community','Range','A1')


%% Make panel observed vs. model figures (Fig 7,8)

% Calculate average R-squared
avgR2_SA = mean(r22(:,3));
avgR2_AB = mean(r22(:,1));
avgR2_BC = mean(r22(:,2));
avgR2_SAC = mean(r22C(:,3));
avgR2_ABC = mean(r22C(:,1));
avgR2_BCC = mean(r22C(:,2));
avgR2_SAP = mean(r22P(:,3));
avgR2_ABP = mean(r22P(:,1));
avgR2_BCP = mean(r22P(:,2));
avgR2_SAPC = mean(r22PC(:,3));
avgR2_ABPC = mean(r22PC(:,1));
avgR2_BCPC = mean(r22PC(:,2));

avgR2_ABSA = mean(r2_ABSA);
avgR2_BCSA = mean(r2_BCSA);
avgR2_BCAB = mean(r2_BCAB);
avgR2_ABSAC = mean(r2_ABSAC);
avgR2_BCSAC = mean(r2_BCSAC);
avgR2_BCABC = mean(r2_BCABC);
avgR2_ABSAP = mean(r2_ABSAP);
avgR2_BCSAP = mean(r2_BCSAP);
avgR2_BCABP = mean(r2_BCABP);
avgR2_ABSAPC = mean(r2_ABSAPC);
avgR2_BCSAPC = mean(r2_BCSAPC);
avgR2_BCABPC = mean(r2_BCABPC);

avgR2_MC = mean(R2save_MC);
avgR2_MCC = mean(R2save_MCC);
avgR2_MCP = mean(R2save_MCP);
avgR2_MCPC = mean(R2save_MCPC);

pos = [1:7 9:15];

% Figure 7 - No phage
fig1 = figure;
hold on
axis tight
T1 = tiledlayout(4,4);
ylabel(T1,'Bacterial density, cfu/ml','FontSize',10);
xlabel(T1,'Days post infection','FontSize',10);
set(gcf,'Unit','Inches','Position',[0,0,7.25,7.25]);

% PA14
cfucurvefig_panels(tspan,PA14_SA_N,PA14_SA_Npred,{'PA14','SA'},{'PA14 + SA'},avgR2_SA,pos(1))
cfucurvefig_panels(tspan,PA14_AB_N,PA14_AB_Npred,{'PA14','AB'},{'PA14 + AB'},avgR2_AB,pos(2))
cfucurvefig_panels(tspan,PA14_BC_N,PA14_BC_Npred,{'PA14','BC'},{'PA14 + BC'},avgR2_BC,pos(3))
cfucurvefig_panels(tspan,PA14_ABSA_N,PA14_ABSA_Npred,{'PA14','AB','SA'},{'PA14 + SA + AB'},avgR2_ABSA,pos(4))
cfucurvefig_panels(tspan,PA14_BCSA_N,PA14_BCSA_Npred,{'PA14','BC','SA'},{'PA14 + SA + BC'},avgR2_BCSA,pos(5))
cfucurvefig_panels(tspan,PA14_BCAB_N,PA14_BCAB_Npred,{'PA14','BC','AB'},{'PA14 + AB + BC'},avgR2_BCAB,pos(6))
cfucurvefig_panels(tspan,PA14_MC_N,PA14_MC_Npred,{'PA14','AB','BC','SA'},{'PA14 Polyculture'},avgR2_MC,pos(7))

% CRKO
cfucurvefig_panels(tspan,PA14_SA_NC,PA14_SA_NpredC,{'CRKO','SA'},{'CRISPR-KO + SA'},avgR2_SAC,pos(8))
cfucurvefig_panels(tspan,PA14_AB_NC,PA14_AB_NpredC,{'CRKO','AB'},{'CRISPR-KO + AB'},avgR2_ABC,pos(9))
cfucurvefig_panels(tspan,PA14_BC_NC,PA14_BC_NpredC,{'CRKO','BC'},{'CRISPR-KO + BC'},avgR2_BCC,pos(10))
cfucurvefig_panels(tspan,PA14_ABSA_NC,PA14_ABSA_NpredC,{'CRKO','AB','SA'},{'CRISPR-KO + SA + AB'},avgR2_ABSAC,pos(11))
cfucurvefig_panels(tspan,PA14_BCSA_NC,PA14_BCSA_NpredC,{'CRKO','BC','SA'},{'CRISPR-KO + SA + BC'},avgR2_BCSAC,pos(12))
cfucurvefig_panels(tspan,PA14_BCAB_NC,PA14_BCAB_NpredC,{'CRKO','BC','AB'},{'CRISPR-KO + BC + AB'},avgR2_BCABC,pos(13))
cfucurvefig_panels(tspan,PA14_MC_NC,PA14_MC_NpredC,{'CRKO','AB','BC','SA'},{'CRISPR-KO Polyculture'},avgR2_MCC,pos(14))

% Add legend
leg1 = legend('Orientation','Vertical');
leg1.Layout.Tile = 16;
neworder = [1,4,2,3];
leg1.AutoUpdate = 'off';
leg1.PlotChildren = leg1.PlotChildren(neworder);

% Save figure
%exportgraphics(fig1,'Fig7_MatlabFinal2.eps','Resolution',300,'ContentType','vector')
%exportgraphics(fig1,'Fig7_MatlabFinal.tif','Resolution',300,'ContentType','vector')


% Figure 8 - Phage
fig2 = figure;
hold on
T2 = tiledlayout(4,4);
ylabel(T2,'Bacterial density, cfu/ml','FontSize',10);
xlabel(T2,'Days post infection','FontSize',10);
set(gcf,'Unit','Inches','Position',[0,0,7.25,7.25]);

% PA14
cfucurvefig_panels(tspan,PA14_SA_Y,PA14_SA_Ypred,{'PA14','SA'},{'PA14 + SA'},avgR2_SAP,pos(1))
cfucurvefig_panels(tspan,PA14_AB_Y,PA14_AB_Ypred,{'PA14','AB'},{'PA14 + AB'},avgR2_ABP,pos(2))
cfucurvefig_panels(tspan,PA14_BC_Y,PA14_BC_Ypred,{'PA14','BC'},{'PA14 + BC'},avgR2_BCP,pos(3))
cfucurvefig_panels(tspan,PA14_ABSA_Y,PA14_ABSA_Ypred,{'PA14','AB','SA'},{'PA14 + SA + AB'},avgR2_ABSAP,pos(4))
cfucurvefig_panels(tspan,PA14_BCSA_Y,PA14_BCSA_Ypred,{'PA14','BC','SA'},{'PA14 + SA + BC'},avgR2_BCSAP,pos(5))
cfucurvefig_panels(tspan,PA14_BCAB_Y,PA14_BCAB_Ypred,{'PA14','BC','AB'},{'PA14 + AB + BC'},avgR2_BCABP,pos(6))
cfucurvefig_panels(tspan,PA14_MC_Y,PA14_MC_Ypred,{'PA14','AB','BC','SA'},{'PA14 Polyculture'},avgR2_MCP,pos(7))

% CRKO
cfucurvefig_panels(tspan,PA14_SA_YC,PA14_SA_YpredC,{'CRKO','SA'},{'CRISPR-KO + SA'},avgR2_SAPC,pos(8))
cfucurvefig_panels(tspan,PA14_AB_YC,PA14_AB_YpredC,{'CRKO','AB'},{'CRISPR-KO + AB'},avgR2_ABPC,pos(9))
cfucurvefig_panels(tspan,PA14_BC_YC,PA14_BC_YpredC,{'CRKO','BC'},{'CRISPR-KO + BC'},avgR2_BCPC,pos(10))
cfucurvefig_panels(tspan,PA14_ABSA_YC,PA14_ABSA_YpredC,{'CRKO','AB','SA'},{'CRISPR-KO + SA + AB'},avgR2_ABSAPC,pos(11))
cfucurvefig_panels(tspan,PA14_BCSA_YC,PA14_BCSA_YpredC,{'CRKO','BC','SA'},{'CRISPR-KO + SA + BC'},avgR2_BCSAPC,pos(12))
cfucurvefig_panels(tspan,PA14_BCAB_YC,PA14_BCAB_YpredC,{'CRKO','BC','AB'},{'CRISPR-KO + BC + AB'},avgR2_BCABPC,pos(13))
cfucurvefig_panels(tspan,PA14_MC_YC,PA14_MC_YpredC,{'CRKO','AB','BC','SA'},{'CRISPR-KO Polyculture'},avgR2_MCPC,pos(14))

% Add legend
leg2 = legend('Orientation','Vertical');
leg2.Layout.Tile = 16;
leg2.AutoUpdate = 'off';
leg2.PlotChildren = leg2.PlotChildren(neworder);

% Save figure
%exportgraphics(fig2,'Fig8_MatlabFinal2.eps','Resolution',300,'ContentType','vector')
%exportgraphics(fig2,'Fig8_MatlabFinal.tif','Resolution',300,'ContentType','vector')


%% Long time predictions (Fig S9)

tlong = (0:1:1000)';
%tlong = (0:1:2000)';

% IC: PA common 
PA14_MC_Nlong = glv_simulation(tlong,PA14_MC_N(1,:),BetaMC);
PA14_MC_NClong = glv_simulation(tlong,PA14_MC_NC(1,:),BetaMCC);
PA14_MC_Ylong = glv_simulation(tlong,PA14_MC_Y(1,:),BetaMCP);
PA14_MC_YClong = glv_simulation(tlong,PA14_MC_YC(1,:),BetaMCPC);

% IC: PA rare
%PA14_MC_Nlong = glv_simulation(tlong,[1e3 PA14_MC_N(1,2:end)],BetaMC);
%PA14_MC_NClong = glv_simulation(tlong,[1e3 PA14_MC_NC(1,2:end)],BetaMCC);
%PA14_MC_Ylong = glv_simulation(tlong,[1e3 PA14_MC_Y(1,2:end)],BetaMCP);
%PA14_MC_YClong = glv_simulation(tlong,[1e3 PA14_MC_YC(1,2:end)],BetaMCPC);

simcolors = {'#06141F','#CD4F38','#3D4F7D','#E48C2A'};

fig3 = figure;
hold on
T3 = tiledlayout(2,2);
ylabel(T3,'Bacterial density, cfu/ml','FontSize',10);
xlabel(T3,'Days post infection','FontSize',10);
set(gcf,'Unit','Inches','Position',[0,0,6,5]);

% PA14
nexttile
hold on
grid on
box on
for i = 1:4
    plot(tlong,log10(PA14_MC_Nlong(:,i)),'LineWidth',1.5,'Color',simcolors{i})
end
plot([0,1000],[2,2],'--','LineWidth',1,'Color',[0.7 0.7 0.7],'HandleVisibility','off')
%plot([0,2000],[2,2],'--','LineWidth',1,'Color',[0.7 0.7 0.7],'HandleVisibility','off')
ylim([0,11])
yticks([2,4,6,8,10])
yticklabels({'10^2','10^4','10^6','10^8'})
xticks([0,10,20,30,40].*24)
xticklabels({'0','10','20','30','40'})
%xticks([0,20,40,60,80].*24)
%xticklabels({'0','20','40','60','80'})
title('PA14 Polyculture','FontWeight','normal','FontSize',9)

% CRKO
nexttile
hold on
grid on
box on
for i = 1:4
    plot(tlong,log10(PA14_MC_NClong(:,i)),'LineWidth',1.5,'Color',simcolors{i})
end
plot([0,1000],[2,2],'--','LineWidth',1,'Color',[0.7 0.7 0.7],'HandleVisibility','off')
%plot([0,2000],[2,2],'--','LineWidth',1,'Color',[0.7 0.7 0.7],'HandleVisibility','off')
ylim([0,11])
yticks([2,4,6,8,10])
yticklabels({'10^2','10^4','10^6','10^8'})
xticks([0,10,20,30,40].*24)
xticklabels({'0','10','20','30','40'})
%xticks([0,20,40,60,80].*24)
%xticklabels({'0','20','40','60','80'})
title('CRISPR-KO Polyculture','FontWeight','normal','FontSize',9)

% PA14 + Phage
nexttile
hold on
grid on
box on
for i = 1:4
    plot(tlong,log10(PA14_MC_Ylong(:,i)),'LineWidth',1.5,'Color',simcolors{i})
end
plot([0,1000],[2,2],'--','LineWidth',1,'Color',[0.7 0.7 0.7],'HandleVisibility','off')
%plot([0,2000],[2,2],'--','LineWidth',1,'Color',[0.7 0.7 0.7],'HandleVisibility','off')
ylim([0,11])
yticks([2,4,6,8,10])
yticklabels({'10^2','10^4','10^6','10^8'})
xticks([0,10,20,30,40].*24)
xticklabels({'0','10','20','30','40'})
%xticks([0,20,40,60,80].*24)
%xticklabels({'0','20','40','60','80'})
title('PA14 Polyculture + phage','FontWeight','normal','FontSize',9)

% CRKO + Phage
nexttile
hold on
grid on
box on
for i = 1:4
    plot(tlong,log10(PA14_MC_YClong(:,i)),'LineWidth',1.5,'Color',simcolors{i})
end
plot([0,1000],[2,2],'--','LineWidth',1,'Color',[0.7 0.7 0.7],'HandleVisibility','off')
%plot([0,2000],[2,2],'--','LineWidth',1,'Color',[0.7 0.7 0.7],'HandleVisibility','off')
ylim([0,11])
yticks([2,4,6,8,10])
yticklabels({'10^2','10^4','10^6','10^8'})
xticks([0,10,20,30,40].*24)
xticklabels({'0','10','20','30','40'})
%xticks([0,20,40,60,80].*24)
%xticklabels({'0','20','40','60','80'})
title('CRISPR-KO Polyculture + phage','FontWeight','normal','FontSize',9)

% Add legend
leg3 = legend({'PA','AB','BC','SA'},'Orientation','Vertical');
leg3.Layout.Tile = 'East';
leg3.AutoUpdate = 'off';
leg3.PlotChildren = leg3.PlotChildren(neworder);

% Save figures
%exportgraphics(fig3,'FigS9_MatlabFinal2.eps','Resolution',300,'ContentType','vector')
%exportgraphics(fig3,'FigS9_MatlabFinal.tif','Resolution',300,'ContentType','vector')
