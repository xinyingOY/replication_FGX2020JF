function[Ri] = data_extract1(rf,P)
    %# different portfolios
    %# main
    %# the following code is taken from the file named main.m
    
    addpath('../data')
    
    

    port_3x2 = csvread('port_3x2.csv',0,1);
    port_3x2 = port_3x2 - rf*ones(1,size(port_3x2,2)); % excess return
    %# 498 780
    port_3x2_id = readtable('port_3x2_id.csv');
    
    kk = 10; % minimun number of stocks in a portfolio
    
    include_3x2 = find(port_3x2_id.min_stk6>=kk)';
    port_3x2b = [];
    
    for i = 1:P
        if ismember(i,include_3x2)
            port_3x2b = [port_3x2b port_3x2(:,(i*6-5):(i*6))];
        end
    end
    
    Ri = port_3x2b; % test asset #% size(Ri) 498 750
