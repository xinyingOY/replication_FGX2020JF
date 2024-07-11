function[Ri] = data_extract2(rf,P)
    %# different portfolios
    %# robustness
    %# the following code is taken from the file named robustness.m
    
    addpath('../data')
    
    kk = 10;

    port_5x5 = csvread('port_5x5.csv',0,1);
    port_5x5 = port_5x5 - rf*ones(1,size(port_5x5,2)); % excess return
    port_5x5_id = readtable('port_5x5_id.csv');
    
    include_5x5 = find(port_5x5_id.min_stk>=kk)';
    port_5x5b = [];
    
    for i = 1:P
        if ismember(i,include_5x5)
            port_5x5b = [port_5x5b port_5x5(:,(i*25-24):(i*25))];
        end
    end
    
    Ri = port_5x5b;
    
    %# 498 1825