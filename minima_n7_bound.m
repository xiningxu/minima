clear; clc;
%% Initialization
N = 5e+4; % number of trails

Fval = zeros(1, N); % value of object function
Theta = zeros(4, N); % output of optimization
options = optimset('display','notify',...
    'TolFun',1e-6,'MaxIter',1e+3,...
    'UseParallel', true);
Result_plus = cell(3,1);
Result_minus = cell(3,1);
for k = 1:3
    % linear constraints
    lb = [pi*ones(k,1); - pi * ones(4-k,1)];
    ub = pi * ones(4,1);
    x0 = [pi*ones(k,N); - pi + 2 * pi * rand(4-k, N)]; % initial point
    
    parfor i = 1:N
        [Theta(:, i), Fval(i),~,~] = ...
            fmincon(@fval_n7,x0(:,i),[],[],[],[],lb,ub,@nonlin_con_n7,options);
    end
    
    Result_plus{k} = sortrows([Fval; Theta]', 1);  
end
% - pi bound
for k = 1:3
    % linear constraints
    lb = - pi * ones(4,1);
    ub = [- pi * ones(k,1); pi * ones(4-k,1)];
    x0 = [- pi*ones(k,N); - pi + 2 * pi * rand(4-k, N)]; % initial point
    
    parfor i = 1:N
        [Theta(:, i), Fval(i),~,~] = ...
            fmincon(@fval_n7,x0(:,i),[],[],[],[],lb,ub,@nonlin_con_n7,options);
    end
    
    Result_minus{k} = sortrows([Fval; Theta]', 1);  
end
% save( 'minima_bdry.mat', 'Result_plus', 'Result_minus');