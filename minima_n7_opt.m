%% Initialization
clear; clc;
folder_name = './output/';
N = 5e+4; % number of trails

% solution
Theta = zeros(4, N); 
Theta_bdry = zeros(4,N);

% objective function value at the solution
Fval = zeros(1, N); 
Fval_bdry = zeros(1,N);

%% Optimization

% bound constraints
lb = - pi * ones(4,1); % lower bound
ub = pi * ones(4,1); % upper bound

% options
options = optimset('display','notify','TolFun',1e-6,'MaxIter',1e+3,...
    'Algorithm', 'interior-point','UseParallel', true);

% initial points
x0 = - pi + 2 * pi * rand(4, N); 

% minimize function fval_n7 subject to nonliear constraints
parfor i = 1:N
    % subject to inequality constraint -> nonlin_con_n7
    [Theta(:, i), Fval(i),~,~] = fmincon(@fval_n7,x0(:,i),...
        [],[],[],[],lb,ub,@nonlin_con_n7,options);
    % subject to equality constraint -> nonlin_con_eq_n7
    [Theta_bdry(:, i), Fval_bdry(i),~,~] = fmincon(@fval_n7,x0(:,i),...
        [],[],[],[],lb,ub,@nonlin_con_eq_n7,options);
end
% sort results
Result = sortrows([Fval; Theta]', 1);
Result_bdry = sortrows([Fval_bdry; Theta_bdry]', 1);
% save([folder_name, 'minima7.mat'], 'Result','Result_bdry');

%% ========== Objective Function for n = 7 ========== %%
function val = fval_n7(x)
    x2 = x(1); x5 = x(2); x6 = x(3); x7 = x(4);

    X1 = [x2; x5; x6; x7; x2-x5; x2-x6; ...
        x2-x7; x5-x6; x5-x7; x6-x7];
    X2 = [x2, x5, x6, x7; ...
        x2, x2-x5, x2-x6, x2-x7;...
        x5, x5-x2, x5-x6, x5-x7;...
        x6, x6-x2, x6-x5, x6-x7;...
        x7, x7-x2, x7-x5, x7-x6];
    r2 = 5 + 2 * sum( cos(X1) );  r = sqrt( 4 ./ r2 - 1);

    val1 = sqrt( 2 * (1 - cos(X1)) ); val1 = sum(val1);
    
    val3 = abs(- 1 - 2 * sum(cos(X1))); val3 = sqrt(val3);

    val21 = 3 + sum(cos(X2), 2);
    val22 = r * sum( sin(X2), 2);
    val2 = sqrt( val21 + val22 ) + sqrt( val21 - val22 );
    val2 = sum( val2);

    val = val1 + val2 + val3;
end

%% ========== Non-linear Contraints for n = 7 ========== %%
function val = nonlin_fun_n7(x)
    theta2 = x(1); theta5 = x(2); 
    theta6 = x(3); theta7 = x(4);

    Y = [theta2; theta5; theta6; theta7; ...
    theta2-theta5; theta2-theta6; theta2-theta7;...
    theta5-theta6; theta5-theta7; theta6-theta7];

    val = sum( cos(Y) ) + 1/2;
end

%% ======== Inequlity Contraints for n = 7 ============= %%
function [val, ceq] = nonlin_con_n7(x)
    val = nonlin_fun_n7(x);
    ceq = [];
end
%% ========= Equlity Contraints for n = 7 ============= %%
function [val, ceq] = nonlin_con_eq_n7(x)
    val = [];
    ceq = nonlin_fun_n7(x);
end
