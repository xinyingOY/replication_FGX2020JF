%% 这个得到的是得到图1的相关参数


clear; clc;
close all;

addpath('../glmnet_matlab')
addpath('../functions')
addpath('../../data')

% fix the random seed for any additional CV
seed_num = 100;

% factor
allfactors = csvread('factors.csv',1,0);
date = allfactors(:,1);
rf = allfactors(:,2);
factors = allfactors(:,3:end);

L = length(date);
P = size(factors,2);

% test portfolios
port_3x2 = csvread('port_3x2.csv',0,1);
port_3x2 = port_3x2 - rf*ones(1,size(port_3x2,2)); % excess return

% other information
summary = readtable('summary.csv');
factorname = summary.Row;
factorname_full = summary.Descpription;
year_pub = summary.Year;
year_end = summary.Year_end;

port_3x2_id = readtable('port_3x2_id.csv');

mkt_ind = find(ismember(factorname,'MktRf'));
smb_ind = find(ismember(factorname,'SMB'));
hml_ind = find(ismember(factorname,'HML'));

%% form a smaller set of portfolios for bivariate sorted porfolios
kk = 10; % minimun number of stocks in a portfolio

include_3x2 = find(port_3x2_id.min_stk6>=kk)';
port_3x2b = [];

for i = 1:P
    if ismember(i,include_3x2)
        port_3x2b = [port_3x2b port_3x2(:,(i*6-5):(i*6))];
    end
end

Ri = port_3x2b; % test asset

% load tune_center for 200 randome seeds selected by cross-validations
% For each randome seed in 1:200, we run a cross-validation and find the
% best tuning parameter.
% Then we take the average for the 200 selected tuning parameters at the
% log scale as tune center, because we plot heat maps on log(lambda).

% choose control factors before 2012
ContrlList = find(year_pub < 2012);
ControlFactor = factors(:,ContrlList);


% test factors since 2012
TestList = find(year_pub >= 2012);
TestFactor = factors(:,TestList);


%% 
% 设置参数
Kfld = 10;  % 折数

alpha = 1;  % LASSO 惩罚项（1 代表 LASSO）
NumLambda = 100;  % 测试的正则化参数数量
seednum = 100;
% 初始化结果存储
num_run=200;
optimal_parameters = NaN(num_run, 2);  % 保存每次模拟的最优参数
for k = 1:200
    
    disp(k)
    
    gt = TestFactor(:,1)'; % test factor
    ht = ControlFactor'; % control factor
    
    result = DS_TSCV(Ri', gt, ht, alpha, NumLambda, seednum+k, Kfld)
    optimal_parameters(k, 1) = result.tune1;
    optimal_parameters(k, 2) = result.tune2;
    
    
end
cd ../main/new_tune
save('tune_figure1_new.mat', 'optimal_parameters');


