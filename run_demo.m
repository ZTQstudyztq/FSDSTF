clear;
addpath('CM');
addpath('dataset');
addpath('measures');

%% Run dataset
filepath1 = ['.\dataset\'];
namelist = dir(strcat([filepath1],['*.mat']));
len = length(namelist);

% I2 is the row in the output table
i2 = 6;

% Choose a dataset
i = 1;
while i <= len
    % ！！！！debug point！！！！
    namelist(i).name
    load(strcat([filepath1],[namelist(i).name]));

    % Print the dataset name
    writematrix([string(namelist(i).name)],'Results.xlsx','Sheet',1,'Range',strcat(['A',num2str(i2)],[':Q',num2str(i2)]));
    i2 = i2 + 1;

    % Put the view in the form d*n
    if size(X{1},1) == size(Y,1)
        nV = length(X);
        for I = 1:nV
            X{I} = double(X{I}');
        end
    else
        nV = length(X);
        for I = 1:nV
            X{I} = double(X{I});
        end
    end
    Y = double(Y);

    %% parameter initialization
    %                   1               3            5
    % param1 = [lambda1, lambda2,| lambda3, p,| alpha, rho1,  rho2,    rho3,   rho4];
    param1 = [0.00001	0.0001	  0.1	  0.5   1.6	 0.0001	 0.0001	 0.0001  0.0001];

    L = 100; % Dimension of feature reduction
    op_param1 = zeros(1,size(param1,2));
    op_result = zeros(1,7);

    % Range of parameter values
    lambda_0 = [1e-5 1e-4 1e-3 1e-2 0.1 1 10 100 1000];
    p_0 = [0.3 0.5 0.7 0.8 0.9 1];
    % alpha_0 = [1.6];
    rho_0 = [1e-6 1e-5 1e-4 1e-3];

    %% Tune parameters for a single dataset
    i1 = 1;
    while i1 <= 3 % 1 3 5
        % ！！！！debug point！！！！
        fprintf('%s dataset on the %d tune\n',namelist(i).name, i1);
        if i1 == 1
            a = lambda_0;
            b = lambda_0;
        elseif i1 == 3
            a = lambda_0;
            b = p_0;
            % elseif i1 == 5
            %     a = alpha_0;
            %       b = rho_0;
        end
        %% Tune two parameters at a time
        j = 1;
        while j <= size(a,2)
            k = 1;
            while k <= size(b,2)
                % ！！！！debug point！！！！
                param1(i1) = a(j);
                param1(i1+1) = b(k);
                fprintf('%s dataset on the %d tune, j = %d，k = %d\n', namelist(i).name, i1, j, k);
                result = FStest(X,Y,param1(1),param1(2),param1(3),param1(4), L, param1(5),param1(6),param1(6),param1(6),param1(6));

                if result(1) > op_result(1)
                    op_param1 = param1;
                    op_result = result;
                    writematrix([i1,op_param1,op_result],'Results.xlsx','Sheet',1,'Range',strcat(['A',num2str(i2)],[':Q',num2str(i2)]));

                end

                if op_result(1) == 1
                    break;
                end
                k = k + 1;
            end
            if op_result(1) == 1
                break;
            end
            j = j + 1;
        end
        i2 = i2 + 1;
        % The results of the previous set are used for the next set
        param1 = op_param1;
        if op_result(1) == 1
            break;
        end
        i1 = i1 + 2;
    end
    i = i + 1;
    clearvars -except filepath1 i i2 len namelist;
end
