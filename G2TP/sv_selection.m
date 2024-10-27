% Selection of the singular values coming from HOSVD of the tensor S_G
%
% Output arguments
% 1. U_weigh_fun - weight function (a.k.a. convex hulls)
% 2. S_weigh_fun - reduced qLPV model (elements are to be put in the correct order)
% 3. Par_vec - vector for HOSVD of the set of LTI models (in which dimension of the state-space tensor HOSVD is to be executed)
%
% The idea of obtaining the first two outputs is based on RHOSVD which is implemented through discarding of singular values for each
% scheduling parameter and generation of convex hulls of the specified type.
% First, SVD for the weight matrix in the specified dimension is executed which allows to discard reasonably small singular values.
% Next, generation of convex hulls is done for each truncated weight matrix.
% Finally, reduced grid-based qLPV model is obtained following HOSVD-based canonical form.
%
% Since model reduction is based on the selection of singular values, it was decided to let the user choose between the auto mode of the selection and manual mode.
%
% In case of the auto mode, reduction starts from keeping 6 singular values for each scheduling parameter and continues with the iteration step of 1 until either
% the minimum possible number of singular values (same for each scheduling parameter) is reached or the maximum Vinnicombe Gap Metric for some pair is greater
% than the predefined upper bound TOL_NORM. In the latter case, a step back is produced automatically to fulfill the tolerance condition. By the way, IRNO type of weight functions
% is prescribed during ordinary iteration in sake for short running time. HOSVD-based canonical form as well as Vinnicombe Gap Metrics are obtained from a second-level subfunction SSM_red.
%
% It may happen that the resulting number of sungular values fulfilling the tolerance condition and defining successfully reduced model makes it infeasible to actually construct such a model.
% By default, LPVTools Toolbox increases the number of the kept singular values in this case automatically, however after a proper testing an extra condition was added for executing a step back
% within the iterative process, since the Vinnicombe Gap Metric seems to be smaller in this case.

function [U_weigh_fun,S_weigh_fun,Par_vec,ParNum] = sv_selection(ParNum,ParNum_keep,mode,TOL_NORM,S_G,Par_range_pr,LPV_0,wtype,Pause_T)

% Vector for High Order Singular Value Decomposition of the set of LTI models
Par_vec = [0 0 ones(1,ParNum)];

% Check for scalar shceduling parameters (SVD/HOSVD cannot be executed for scalars, so we remove the corresponding scheduling parameters)
for i = 1:ParNum
    if length(LPV_0.Domain.IVData{i,1}) == 1
        ParNum = ParNum - 1; % In this case, ParNum is not only the dimension of the state vector, but also a reduced dimension keeping non-scalar parameters only
    end
end

Par_vec = [0 0 ones(1,ParNum)]; % Reassign Par_vec variable once more in case of reduction of ParNum

if numel(mode) == numel('auto') && min(mode == 'auto') == 1 % Auto selection of the singular values
    % Using 'numel()' - word length - and 'min()' - order of letters - commands is enough to distinguish between words 'auto' and 'manual'

    disp('Warning: Auto mode of selection of the singular variables was selected')
    disp('Remark: The selection is based on the nu-gap metric')
    disp('-----------------------------------------------------------------------------')

    % Starting minimum and maximum possible numbers of the singular values for
    % each scheduling parameter
    sv_num_min_0 = 2;
    sv_num_max_0 = 6;

    % Predefine loop variables
    nugap = 0;
    sv_num_max = 0;

    % Auto correction algorithm
    while max(nugap) <= TOL_NORM

        % Define a type of the weighting functions for fast calculations
        wtype_0 = 'irno';

        [~,~,~,sv,nugap,~,~] = SSM_red(S_G,Par_vec,ParNum,sv_num_max_0,LPV_0,Par_range_pr,wtype_0,mode);

        disp(strcat('Kept number of singular values:'," ", string(sv_num_max_0)))
        disp(strcat('Maximum value of nu-gap metric:'," ", string(max(nugap))))
        disp('-----------------------------------------------------------------------------')
        pause(Pause_T)

        % Maximum number of singular values among all scheduling parameters
        for i = 1:ParNum
            sv_num(i) = numel(sv{i+2});
        end

        sv_num_max = max(sv_num);

        if max(nugap) > TOL_NORM && sv_num_max >= sv_num_min_0  % Maximum error is greater than the required tolerance whether the minimum number of the singular
            % values is reached or not ---> fail of model reduction
            disp(strcat('Warning: Tolerance requirement was not fulfilled (number of singular values = ', " ", string(sv_num_max), ')'))
            disp(strcat('Warning: Total number of kept singular values = ', " ", string(sv_num_max+1)))
            pause(Pause_T)

            % Step back to the previous number of singular values
            % Modified inputs for the step back
            wtype_0 = wtype;
            sv_num_max_0 = sv_num_max_0 + 1;

            [U_0,S_0,~,~,~,~,~] = SSM_red(S_G,Par_vec,ParNum,sv_num_max_0,LPV_0,Par_range_pr,wtype_0,mode);

            % Reassign the necessary outputs for the consistency in the following commands
            U_weigh_fun = U_0;
            S_weigh_fun = S_0;

            break % Stop the iteration process and continue with the resulted configuration

        elseif max(nugap) < TOL_NORM && sv_num_max == sv_num_min_0 % Successful model reduction

            % Reassign hte type of the weighting functions
            wtype_0 = wtype;

            [U_0,S_0,~,~,~,~,~] = SSM_red(S_G,Par_vec,ParNum,sv_num_max_0,LPV_0,Par_range_pr,wtype_0,mode);

            % Reassign the necessary outputs for the consistency in the following commands
            U_weigh_fun = U_0;
            S_weigh_fun = S_0;

            % Handling self-adjustment of the nubmer of singular values (number is increased by one; NOT in the uncertainty), because
            % implementation of the algorithm having current number is not feasible ---> step back
            % Remark: it comes from the fact that going with the increased number without self-adjustment is more accurate than with it
            for i = 1:size(U_weigh_fun,2)
                if size(U_weigh_fun{i},2) == sv_num_max_0 + 1 && i ~= size(U_weigh_fun,2) || size(U_weigh_fun{i},2) == sv_num_max_0 + 1 && ParNum == ParNum_keep % With uncertainty || without uncertainty

                    % Modified inputs for the step back
                    wtype_0 = wtype;
                    sv_num_max_0 = sv_num_max_0 + 1;

                    [U_0,S_0,~,~,~,~,~] = SSM_red(S_G,Par_vec,ParNum,sv_num_max_0,LPV_0,Par_range_pr,wtype_0,mode);

                    % Reassign the necessary outputs for the consistency in the following commands
                    U_weigh_fun = U_0;
                    S_weigh_fun = S_0;
                    disp('Extra step back was implemented')
                    disp(strcat('Warning: Total number of kept singular values = ', " ", string(sv_num_max+1)))
                    break % Stop futher iterating, as at least one self-adjusted scheduling patameter has been found
                end
            end

            disp('Warning: Optimal number of singular values for each scheduling parameter was reached')
            % -----
            break % Stop the iteration process and continue with the resulted configuration
        else
            sv_num_max_0 = sv_num_max_0 - 1; % Iterarion step
        end
    end

elseif numel(mode) == numel('manual') && min(mode == 'manual') == 1 % Manual mode

    disp('Warning: Manual mode of selection of the singular variables was selected')
    disp('Remark: The selection is based on the nu-gap metric')
    disp('-----------------------------------------------------------------------------')

    wtype_0 = wtype;
    sv_num_max_0 = []; % Predefined variable to start the SSM_red function

    [U_0,S_0,~,~,nugap,ParNum_HOSVD,ParNum_HOSVD_vec] = SSM_red(S_G,Par_vec,ParNum,sv_num_max_0,LPV_0,Par_range_pr,wtype_0,mode);

    %Kept number of singular values after the manual input
    sv_num_max_0 = ParNum_HOSVD;

    % Reassign the necessary outputs for the consistency in the following commands
    U_weigh_fun = U_0;
    S_weigh_fun = S_0;

    % Check for the self-adjustment
    for i = 1:size(U_weigh_fun,2)
        if size(U_weigh_fun{i},2) == sv_num_max_0 + 1 && i ~= size(U_weigh_fun,2) && ParNum_HOSVD_vec(i) ~= size(U_weigh_fun{i},2) || size(U_weigh_fun{i},2) == sv_num_max_0 + 1 && ParNum == ParNum_keep && ParNum_HOSVD_vec(i) ~= size(U_weigh_fun{i},2) % With uncertainty || without uncertainty

            % Modified inputs for the step back
            wtype_0 = wtype;
            sv_num_max_0 = sv_num_max_0 + 1;
            mode = 'auto'; % Switch to the auto mode in order not to choose the number of kept singular values over again

            % Polytopic representation and specified weight functions
            % Remark: variable 'ParNum_HOSVD' is not used here as an output of the SSM_red function, because the 'auto' mode does not use it,
            % and it should not be overwritten in this case
            [U_0,S_0,~,~,~,~] = SSM_red(S_G,Par_vec,ParNum,sv_num_max_0,LPV_0,Par_range_pr,wtype_0,mode);

            % Reassign the necessary outputs for the consistency in the following commands
            U_weigh_fun = U_0;
            S_weigh_fun = S_0;
            disp('Extra step back was implemented')
            disp(strcat('Warning: Total number of kept singular values = ', " ", string(ParNum_HOSVD+1)))
            break % Stop futher iterating, as at least one self-adjusted scheduling patameter has been found
        end
    end
    disp('-----------------------------------------------------------------------------')
    disp(strcat('Maximum value of nu-gap metric:'," ", string(max(nugap))))
    disp('-----------------------------------------------------------------------------')
end
end