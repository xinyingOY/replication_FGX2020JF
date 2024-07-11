%# a function
function[test_num] = generate_matrix(Ri,allf,name1,name2,factorname)
    f = size(allf',1); % 150
    
    test_num = size(Ri);
   
    

    tmpf = nancov([allf,Ri]);
    cov_f = tmpf((f+1):end,1:f);
    
    ER    = nanmean(Ri',2); % 750 1
    
    %# merge factors and ER
    v = [cov_f ER]; %# 750 151
    
    factorname_v = [factorname' "ER"];
    
    v_new = array2table(v,'VariableNames',factorname_v);
    
    writetable(v_new, strcat("../output_new/C_matrix/",name1,"_",name2,".csv"));
   
