% This Version: Jan 4, 2019. 
% @copyright Guanhao Feng, Stefano Giglio and Dacheng Xiu

% We implement our program via the matlab version of glmnet 
% <http://web.stanford.edu/~hastie/glmnet_matlab/> on Matlab 2014b. 

%%
% The main function for Double-Selection implementation
% no cross-validation for 1st and 2nd selections
%
function result = DS(Ri, gt, ht, tune1, tune2, alpha, seednum)

% dim of g is 1!

% depart
% Ri = Ri';
% tune1 = -log(tune_center(1,1));
% tune2 = -log(tune_center(1,2));
% alpha = 1;
% seednum = 100;


if isempty(alpha)
    alpha = 1; % default is lasso
end

% data information
n = size(Ri,1);
p = size(ht,1);
d = size(gt,1);

tmp1  = nancov([gt',Ri']);
cov_g = tmp1((d+1):end,1:d);
tmp2  = nancov([ht',Ri']);
cov_h = tmp2((p+1):end,1:p);

ER    = nanmean(Ri,2);

beta = NaN(n,p);
for i = 1:p
  beta(:,i) = cov_h(:,i)/nanvar(ht(i,:));
end
penalty = mean(beta.^2,1);
penalty = penalty./mean(penalty); % normalize the level

lambda0 = exp(linspace(0,-35,100));

% 1st selection in cross-sectional regression
opts1 = struct('standardize',false,'lambda',exp(-tune1),'alpha',alpha);
model1 = glmnet(cov_h*(diag(penalty)), ER,'gaussian',opts1);
model1_est = glmnetCoef(model1);
sel1 = find(model1_est(2:(p+1))~=0)';
err1 = mean((ER - [ones(n,1),cov_h*(diag(penalty))]*model1_est).^2);

% 2nd selection
sel2 = [];
err2 = NaN(d,1);
opts2 = struct('standardize',false,'lambda',exp(-tune2),'alpha',alpha);
for i = 1:d
    model2 = glmnet(cov_h*(diag(penalty)),cov_g(:,i),'gaussian',opts2);
    model2_est = glmnetCoef(model2);
    sel2 = [sel2; find(model2_est(2:end)~=0)];
    err2(i) = mean((cov_g(:,i) - [ones(n,1),cov_h*(diag(penalty))]*model2_est).^2);
end
sel2 = unique(sel2)';

% 3rd selection for avar zt
sel3 = [];
for i = 1:d
    TSCVout =  TSCVds(Ri, gt(i,:), ht, lambda0, 10, 1, alpha, seednum);
    sel3 = [sel3; TSCVout.sel3_1se];
end
sel3 = unique(sel3);

% post-selection estimation and inference
dsout = inferds(Ri, gt, ht, sel1, sel2, sel3);
ssout = inferds(Ri, gt, ht, sel1, [], sel3);

% output for Double Selection
result.lambdag_ds = dsout.lambdag;
result.se_ds = dsout.se;
result.gamma_ds = dsout.gamma;

% output for Single Selection
result.lambdag_ss = ssout.lambdag;
result.se_ss = ssout.se;
result.gamma_ss = ssout.gamma;

% % selection results
result.sel1 = sel1;
result.sel2 = sel2;
result.sel3 = sel3;
select1 = union(sel1,sel2);
result.select = select1;
result.err1 = err1;
result.err2 = err2;

end




