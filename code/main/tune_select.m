%% 这里得到的是待测因子正文中使用的参数。
%% 插入数据

clear; clc;
close all;

addpath('../glmnet_matlab')
addpath('../functions')
addpath('../../data')

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

% choose control factors before 2012
ContrlList = find(year_pub < 2012);
ControlFactor = factors(:,ContrlList);
FF3 = factors(:,[mkt_ind,smb_ind,hml_ind])';

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
all_tune_results = []; % 调整维度
num_run=200;

% test factor individually
for j = 1:length(TestList)
% for j = 1:1

    disp(j)

% depart1
% j = 1;

gt = TestFactor(:,j)'; % test factor
ht = ControlFactor'; % control factor

optimal_parameters = NaN(num_run, 3);  % 保存每次模拟的最优参数

    for i = 1:num_run
        result = DS_TSCV(Ri', gt, ht, alpha, NumLambda, seednum+i, Kfld)
        
        optimal_parameters(i, 1) = result.tune1;
        optimal_parameters(i, 2) = result.tune2;
        optimal_parameters(i, 3) = j;
    end
    all_tune_results = [all_tune_results; optimal_parameters];
    end
    %% 

% % 使用 splitapply 对每个分组计算平均值
index_column = all_tune_results(:, end);  % 提取索引列
data_to_average = all_tune_results(:, 1:end-1);  % 剔除索引列
average_tune = splitapply(@mean, data_to_average, index_column);

%%%与原论文参数对比
load('tune_main.mat');
log_tunecenter=log(tune_center);

%%%这里是先平均再取对数的结果
log_tune=log(average_tune);
cd ../main/new_tune

save('all_tune_results.mat', 'all_tune_results');
save('average_tune.mat', 'average_tune');
save('log_tune.mat', 'log_tune');

%%%这里是先取对数再平均的结果
log_alltune=log(all_tune_results);
data_to_average = log_alltune(:, 1:end-1);  % 剔除索引列
log_average_tune = splitapply(@mean, data_to_average, index_column);
save('log_average_tune.mat', 'log_average_tune');


