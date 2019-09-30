%% Pick the file
[filename, pathaname]=uigetfile('*xlsx','Pick the data in the weighted means file');


%% cdc42 total density and error bars acquisition
% bem1 deleted data
%total density at 24 hours
db0=xlsread(filename,'Error bars cdc42','O11');% 0% Gal total density
db1=xlsread(filename,'Error bars cdc42','P11');% 0.06% Gal
db2=xlsread(filename,'Error bars cdc42','Q11');% 0.1% Gal
%error bars
stb0=xlsread(filename,'Error bars cdc42','O13');% 0% Gal
stb1=xlsread(filename,'Error bars cdc42','P13');% 0.06% Gal
stb2=xlsread(filename,'Error bars cdc42','Q13');% 0.1% Gal


% bem1bem3 deleted data
db130=xlsread(filename,'Error bars cdc42','O7');% 0% Gal total density
db131=xlsread(filename,'Error bars cdc42','P7');% 0.06% Gal
db132=xlsread(filename,'Error bars cdc42','Q7');% 0.1% Gal
%error bars
stb130=xlsread(filename,'Error bars cdc42','O9');% 0% Gal
stb131=xlsread(filename,'Error bars cdc42','P9');% 0.06% Gal
stb132=xlsread(filename,'Error bars cdc42','Q9');% 0.1% Gal

% WT+cdc42 data

wt0=xlsread(filename,'Error bars cdc42','O3');% 0% Gal total density
wt1=xlsread(filename,'Error bars cdc42','P3');% 0.06% Gal
wt2=xlsread(filename,'Error bars cdc42','Q3');% 0.1% Gal
%error bars
stwt0=xlsread(filename,'Error bars cdc42','O5');% 0% Gal
stwt1=xlsread(filename,'Error bars cdc42','P5');% 0.06% Gal
stwt2=xlsread(filename,'Error bars cdc42','Q5');% 0.1% Gal

%% Plot of the cdc42 per genotypes per gal at 24h

db=[db0,db1,db2];
stdb=[stb0,stb1,stb2];

db13=[db130,db131,db132];
stdb13=[stb130,stb131,stb132];

wt=[wt0,wt1,wt2];
stwt=[stwt0,stwt1,stwt2];
figure;
errorbar(db/48,stdb/48,'Color',[0 0 0],'LineWidth',1)
hold on
scatter([1 2 3],db/48,200,'filled','DisplayName','dBem1');
hold on 
line([1 2 3],[1 1 1],'Color',[1 0 0],'LineWidth',1)
set(gca, 'XScale','linear', 'YScale','log')
grid('on')

hold on
errorbar(db13/48,stdb13/48,'Color',[0 0 0],'LineWidth',1)
hold on
scatter([1,2, 3],db13/48,200,'filled','DisplayName','dBem1dBem3');
hold on 
line([1,2,3],[1 1 1],'Color',[1 0 0],'LineWidth',1)
set(gca, 'XScale','linear', 'YScale','log')
grid('on')


hold on 

errorbar(wt/48,stwt/48,'Color',[0 0 0],'LineWidth',1)
hold on
scatter([1,2,3],wt/48,200,'filled','DisplayName','WT+cdc42');
hold on 
line([1,2,3],[1 1 1],'Color',[1 0 0],'LineWidth',1)
set(gca, 'XScale','linear', 'YScale','log')
grid('on')


