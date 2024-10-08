% This is the Matlab code for paper
%**************************************************************************
% This code is constructed in Matlab R2020a
clear all
clc

%**************************************************************************
% 实证分析
% 新输出的复刻结果在'output/output_new'文件夹中
% 原文输出结果在'output/output_original'文件夹中
%**************************************************************************

%%% 得到原文表1结果 table 1
% Note: 这里使用的是原论文参数。复刻模型参数，可以运行'main/tune_select.m'
run code/main/main.m

%%% 得到原文图1结果 figure1
% Note: 这里使用的是原论文参数。复刻模型参数，可以运行'main/tune_select_figure1.m'
cd ../../../code
run main/plot_figure1.m

%%% 得到图2结果（未作热力图版本） figure2
% Note: 这里使用的是复刻的参数。复刻模型参数，可以运行'main/tune_select.m'
cd ../../../code
run robustness/plot_figure2.m
%**************************************************************************
% % 下面为原始论文输出的结果，且尚未进行参数复刻的部分
cd ../../../code
run robustness/robustness.m

% 这个运行时间稍微有点长
cd ../../../code
run robustness/robustness_stepwise.m
