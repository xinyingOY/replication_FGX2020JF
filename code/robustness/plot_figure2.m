%  Robustness to the Choice of Tuning Parameters
% 即待测因子在200次模拟的参数选择表现
% 提供了paper中图2的参数表现，但没有绘制热力图，只是提供了200次模拟的参数分布情况



clear; clc;
close all;

addpath('../glmnet_matlab')
addpath('../functions')
addpath('../../data')
addpath('../main/new_tune')

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

Ri = port_3x2; % test asset

% load tune_center for 200 randome seeds selected by cross-validations
% For each randome seed in 1:200, we run a cross-validation and find the
% best tuning parameter.
% Then we take the average for the 200 selected tuning parameters at the
% log scale as tune center, because we plot heat maps on log(lambda).

load('all_tune_results.mat');
load('log_average_tune.mat');


% choose control factors before 2012
ContrlList = find(year_pub < 2012);
ControlFactor = factors(:,ContrlList);


% test factors since 2012
TestList = find(year_pub >= 2012);
TestFactor = factors(:,TestList);


%% 画参数选择图
% 创建一个新的图形窗口
fig = figure('Position', [100, 100, 800, 800]);  % [left, bottom, width, height]
numPlots = 15;

log_alltune=log(all_tune_results);

for i = 1:numPlots
    % 生成一个逻辑向量，指示第三列是否为i
    logical_index = all_tune_results(:, 3) == i;

    % 使用逻辑索引筛选矩阵的特定行
    x = log_alltune(logical_index, 1);
    y = log_alltune(logical_index, 2);
    
%     logical_index1= tune_tstat(:, 2) == i;
%     heatmap_tstast = tune_tstat(logical_index1,1);
    
    % 创建子图
    subplot(5, 3, i);
    
    % 绘制热力图
%     获取唯一的 x 和 y 坐标值
%     unique_x = unique(x);
%     unique_y = unique(y);
%     
%     sorted_x = sort(unique_x);
%     sorted_y = sort(unique_y);
%     [X, Y] = meshgrid(sorted_x, sorted_y);
%     
% 
%     
%     初始化颜色矩阵
%     color_matrix = zeros(size(X));
% 
%     将数据点的值大小填充到对应的坐标位置上
%     for j = 1:length(x)
%         查找数据点在坐标网格中的位置
%         idx_x = find(sorted_x == x(j));
%         idx_y = find(sorted_y == y(j));
% 
%         填充颜色矩阵
%         color_matrix(idx_y, idx_x) = heatmap_tstast(j);
%     end
% 
%     pcolor(X, Y, color_matrix);
%     shading interp;  % 设置颜色插值方式
%     colorbar;
%     
%     hold on;  % 保持当前图形以添加更多元素
    % 绘制点图
    scatter(x, y,'k','filled','SizeData', 5);
    
    hold on;
    % 添加特殊点，红色叉叉标记
    select_tune1 = log_average_tune(i,1); 
    select_tune2 = log_average_tune(i,2); 
    scatter(select_tune1, select_tune2, 'rx');  % 'rx' 表示红色叉叉
    
    title(['Scatter Plot ', num2str(i)]);
    hold off;  % 不再保持当前图形，以便下一个子图
end

% 设置整体标题
sgtitle('Robustness to tuning parameters');

% print('scater_plot.pdf', '-dpdf', '-bestfit','-r300');  % '-r300' 指定分辨率为 300 dpi

cd ../../output/output_new/robustness
print('tune_plot.pdf', '-dpdf', '-bestfit','-r300');  % '-r300' 指定分辨率为 300 dpi