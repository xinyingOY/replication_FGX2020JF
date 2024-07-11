%# matrix_generation

%# input the type of matrix you would like to generate


%# the following code is taken from the file named main.m

clear; clc;
close all;

addpath('../data')

% factor
allfactors = csvread('factors.csv',1,0);
date = allfactors(:,1);
rf = allfactors(:,2);
factors = allfactors(:,3:end);
allf = factors;

L = length(date);
P = size(factors,2);


% other information
summary = readtable('summary.csv');
factorname = summary.Row;
factorname_full = summary.Descpription;

%% form a smaller set of portfolios for bivariate sorted porfolios

%# rewritten for Variable: Ri
%# port_3x2
% Ri = data_extract1(rf,P);
%# port_5x5
% Ri = data_extract2(rf,P);
%# port202
% Ri = data_extract3(rf,P);


Ri_1 = data_extract1(rf,P);

Ri_2 = data_extract2(rf,P);

Ri_3 = data_extract3(rf,P);


%% form a smaller time span
%# add a function for matrix interception
%# make a count for which fits the time
%# all matrix is 498

%# create a variable timestamp to store 
%# the time you would like to depart
time_stamp = 20021231;


date1 = date(date <= time_stamp);
n1 = size(date1,1);





%% continue to generate matrix
generate_matrix(Ri_1, allf, "port_3x2","1976_2017",factorname);
generate_matrix(Ri_2, allf, "port_5x5","1976_2017",factorname);
generate_matrix(Ri_3, allf, "port_202","1976_2017",factorname);

generate_matrix(Ri_1(1:n1,:), allf(1:n1,:), "port_3x2","1976_2002",factorname);
generate_matrix(Ri_2(1:n1,:), allf(1:n1,:), "port_5x5","1976_2002",factorname);
generate_matrix(Ri_3(1:n1,:), allf(1:n1,:), "port_202","1976_2002",factorname);

generate_matrix(Ri_1(n1+1:end,:), allf(n1+1:end,:), "port_3x2","2003_2017",factorname);
generate_matrix(Ri_2(n1+1:end,:), allf(n1+1:end,:), "port_5x5","2003_2017",factorname);
generate_matrix(Ri_3(n1+1:end,:), allf(n1+1:end,:), "port_202","2003_2017",factorname);



% generate_matrix(Ri=Ri_2,allf=allf,name1 = "port_5x5",name2 = "1976_2017",factorname = factorname);
% generate_matrix(Ri=Ri_3,allf=allf,name1 = "port_202",name2 = "1976_2017",factorname = factorname);

% generate_matrix(Ri=Ri_1(1:n1,:),allf=allf(1:n1,:),name1 = "port_3x2",name2 = "1976_2002",factorname = factorname);
% generate_matrix(Ri=Ri_2(1:n1,:),allf=allf(1:n1,:),name1 = "port_5x5",name2 = "1976_2002",factorname = factorname);
% generate_matrix(Ri=Ri_3(1:n1,:),allf=allf(1:n1,:),name1 = "port_202",name2 = "1976_2002",factorname = factorname);
% 
% generate_matrix(Ri=Ri_1(n1+1:end,:),allf=allf(n1+1:end,:),name1 = "port_3x2",name2 = "2003_2017",factorname = factorname);
% generate_matrix(Ri=Ri_2(n1+1:end,:),allf=allf(n1+1:end,:),name1 = "port_5x5",name2 = "2003_2017",factorname = factorname);
% generate_matrix(Ri=Ri_3(n1+1:end,:),allf=allf(n1+1:end,:),name1 = "port_202",name2 = "2003_2017",factorname = factorname);

% Ri=Ri_1;
% allf = allf;
% name1 = "port_3x2";
% name2 = "1976_2017";
% 
% f = size(allf',1); % 150
% 
% tmpf = nancov([allf,Ri]);
% cov_f = tmpf((f+1):end,1:f);
% 
% ER    = nanmean(Ri',2); % 750 1
% 
% %# merge factors and ER
% v = [cov_f ER]; %# 750 151
% 
% factorname_v = [factorname' "ER"];
% 
% v_new = array2table(v,'VariableNames',factorname_v);
% 
% writetable(v_new, strcat(name1,"_",name2,".csv"));

