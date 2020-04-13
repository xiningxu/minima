%% Initialization
clear; clc;
result_folder =  './output/'; % folder name for saving
N = 1e+4; % number of trials

% solution 
Theta = zeros(2,N); 
Theta_bdry = zeros(2,N);

% objective function value at the solution 
Fval = zeros(1,N); 
Fval_bdry = zeros(1,N);

%% OPIMIZATION

% bound constraints
lb = zeros(2,1); % lower bound
ub = pi * ones(2,1); % upper bound

% options
options = optimset( 'display', 'notify', 'TolFun', 1e-6, 'MaxIter', 1e+3,...
    'Algorithm', 'interior-point', 'UseParallel', true); 
% initial points
x0 = pi * rand(2,N); 

% minimize function fval_n5 subject to nonliear constraints
parfor i = 1:N
    % subject to inequality constraint -> nonlin_con_n5
    [Theta(:,i),Fval(i),~,~] = fmincon( @fval_n5, x0(:,i),....
        [],[],[],[], lb, ub, @nonlin_con_n5, options); 
    % subject to equality constraint -> nonlin_con_eq_n5
    [Theta_bdry(:,i),Fval_bdry(i),~,~] = fmincon( @fval_n5, x0(:,i),...
        [],[],[],[], lb, ub, @nonlin_con_eq_n5, options); 
end

Result = sortrows( [Fval;Theta]', 1);
Result_bdry = sortrows( [Fval_bdry; Theta_bdry]', 1);
% save([result_folder, 'minima5.mat'], 'Result', 'Result_bdry' );

%% ========== Objective Function for n = 5 ========== %%
function val = fval_n5(x) 
theta2 = x(1); theta5 = x(2); % x2: theta2;  x5:theta5

X1 = [theta2; theta5; theta2 - theta5];
X2 = [theta2, theta5; theta2, theta2-theta5; theta5, theta5 - theta2 ];

c3 = cos(X1); 
r = 4 / ( 3 + 2 * sum(c3) ) - 1; r = sqrt(r);

val1 = sqrt( 2 * ( 1 - c3 ) ); val1 = sum( val1 );
val3 = abs( 1 - 2 * sum( c3 ) ); val3 = sqrt( val3 );

val21 = 3 + sum( cos(X2), 2);
val22 = r * sum( sin(X2), 2);
val2 = sqrt( val21 + val22 ) + sqrt( val21 - val22 );
val2 = sum( val2 );

val = val1 + val2 + val3;
end
%% ========== Non-linear Contraints for n = 5 ========== %%
function fval = non_con_fun(theta2,theta5)
    fval = cos(theta2) + cos(theta5) + cos( theta2 - theta5 ) - 1/2;
end

%% ========== Inequlity Contraints for n = 5 ========== %%
function [fval,ceq] = nonlin_con_n5(x)
    fval = non_con_fun(x(1), x(2));
    ceq = [];
end
%% =========== Equlity Contraints for n = 5 =========== %%
function [fval,ceq] = nonlin_con_eq_n5(x)
    fval = [];
    ceq = non_con_fun(x(1), x(2));
end