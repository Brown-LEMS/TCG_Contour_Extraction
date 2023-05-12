clear all; close all;
addpath (genpath('../util/'));

prefix = 'SE_SEL_';

%% train the beta for merging classifier

load([prefix 'Features.mat']);
load([prefix 'Y.mat']);

Features_1 = [ones(size(Features,2),1) Features'];
Y = Y';

Features_test = Features_1(1:length(Y)/2, :);
Features_train = Features_1(length(Y)/2+1:end, :);
Y_test = Y(1:length(Y)/2);
Y_train = Y(length(Y)/2+1:end);

cnt_pos_test = sum(find(Y_test));
cnt_neg_test = sum(find(Y_test==0));

% train beta
fstd = std(Features_train);
fstd = fstd + (fstd==0);
fmean = mean(Features_train);
fmean(1) = 0;

Features_train = (Features_train - repmat(fmean,size(Features_train,1),1))./ repmat(fstd,size(Features_train,1),1);

% fit the model
fprintf(2,'Training model Random Forest...\n');

model = classRF_train(Features_train,Y_train,500);

save_name = [prefix 'RF_model.mat'];
save(save_name,'model','fmean','fstd');

fprintf(2,'Testing Random Forest...\n');
% normalize test data
Features_test = (Features_test - repmat(fmean,size(Features_test,1),1))./ repmat(fstd,size(Features_test,1),1);

Y_hat = classRF_predict(Features_test,model);

% estimates = zeros(size(p_test));
% estimates(find(p_test>0.5)) = 1;
% 
% 
error = sum(abs(Y_hat - Y_test))/length(Y_hat)


