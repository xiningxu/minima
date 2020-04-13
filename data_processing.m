%% Initialization 
clear; clc;
folder_name = './output/';

file_name = 'minima_n5.mat'; %  n=5->'minima_n5.mat'; n=7->'minima_n7.mat'
N = 3500; % top N results: n = 7 -> N = 8000; n = 5 -> N = 3500
m = 2; % number of variables: n = 7 -> m = 4; n = 5 -> m = 2 

n = 2; % rounded rounded up to n decimal places

data = load([folder_name, file_name]);
Result =  data.Result; 
% Result = data.Result_bdry;

theta = round( Result(1:N, 2:(m+1)), n );  
theta_uni = unique( theta, 'rows'); % unique solution

%% Match the original solution
optim_loc = zeros( length(theta_uni), 1 );
optim_fval = zeros( length(theta_uni), 1);
optim_theta = zeros( length(theta_uni), m);

for i = 1:length(theta_uni)
    [~, optim_loc(i)] = ismember( theta_uni(i,:), theta, 'rows');
    optim_theta(i,:) = Result(optim_loc(i), 2:end);
    optim_fval(i) = Result(optim_loc(i), 1);
end
optim = sortrows([optim_fval, optim_theta, optim_loc], 1);
