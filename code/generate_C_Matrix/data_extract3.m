function[Ri] = data_extract3(rf,~)
    %# different portfolios
    %# robustness
    %# the following code is taken from the file named robustness.m
    
    addpath('../data')
    
    port_202 = csvread('port202.csv',0,1)/100;
    port_202 = port_202 - rf*ones(1,size(port_202,2));
    
    Ri = port_202;
    
    %# 498 202