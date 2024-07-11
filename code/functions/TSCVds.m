%%
% The function for cross-validation over time
% only for 3rd selection in the DS function
%

function output = TSCVds(Ri, gt, ht, lambda, Kfld, Jrep,alpha,seednum)

if isempty(seednum)
    seednum = 101;
end

% data information
[p,T] = size(ht);

L = length(lambda);

cvm3 = NaN(L,Kfld,Jrep);
cvm33 = [];

nomissing = (sum(isnan([ht;gt]),1)==0)';

for j = 1:Jrep

    rng(seednum+j)
    indices = crossvalind('Kfold',T,Kfld);

    for k = 1:Kfld

        % divide the train and test samples
        test = (indices == k);
        train = (indices ~= k);

        ht_train = ht(:,train & nomissing);
        gt_train = gt(:,train & nomissing);

        ht_test = ht(:,test);
        gt_test = gt(:,test);

        opts3 = struct('intr',false,'standardize',true,'lambda',lambda,'alpha',alpha);
        model3 = glmnet(ht_train', gt_train','gaussian',opts3);

        gt_pred = ht_test'*model3.beta;

        LL3 = length(model3.lambda);

        cvm3(1:LL3,k,j) = nanmean((repmat(gt_test',1,LL3) - gt_pred).^2,1)';
    end

    cvm33 = [cvm33, cvm3(:,:,j)];
end

cv_sd3 = std(cvm33')/sqrt(Kfld*Jrep);
cvm333 = mean(cvm33,2);
[~,l_sel3] = min(cvm333);

cvm33ub = cvm333(l_sel3) + cv_sd3(l_sel3);
l3_1se = find(cvm333(1:l_sel3) >= cvm33ub, 1,'last');
if isempty(l3_1se)
    l3_1se = l_sel3;
end

% to reestimate the model with all data
% refit the model

opts33 = struct('intr',false,'standardize',true,'lambda',lambda([l3_1se l_sel3]),'alpha',alpha);
model3 = glmnet(ht(:,nomissing)', gt(nomissing)','gaussian',opts33);

sel3 = find(model3.beta(:,2) ~= 0);
output.sel3 = sel3;
output.lambda3 = lambda(l_sel3);


sel3_1se = find(model3.beta(:,1) ~= 0);
output.sel3_1se = sel3_1se;
output.lambda3_1se = lambda(l3_1se);

end
