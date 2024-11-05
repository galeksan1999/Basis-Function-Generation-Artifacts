% Polytopic representation and specified weightions
%
% Output arguments
% 1. U_0 := U_weigh_fun
% 2. S_0 = S_weigh_fun
% 3. U - matrices for each dimension (see hosvd)
% 4. sv - n-mode singular values (see hosvd)
% 5. nugap - vector of the Vinnicombe Gap Metrics
% 6. ParNum_HOSVD - kept number of singular values in the manual mode

function [U_0,S_0,U,sv,nugap,ParNum_HOSVD,ParNum_HOSVD_vec] = SSM_red(S_G,Par_vec,ParNum,sv_num_max_0,LPV_0,Par_range_pr,wtype_0,mode)

% HOSVD with no tolerance
try
   [T U sv] = hosvd(S_G, Par_vec);
catch out_mem
   if out_mem.message == 'Out of memory.'
   [T U sv] = hosvd_tall(S_G, Par_vec);
   end 
end 

tol = [0 0]; % Define an initial row-vector containing tolerances for the first two dimensions of sv cell array

if numel(mode) == numel('manual') && min(mode == 'manual') == 1 % Manual mode

    % Tolerances for each scheduling parameter
    ParNum_HOSVD_vec = [];
    for i = 1:ParNum
        disp(strcat('Singular values of:',"      ",LPV_0.Domain.IVName{i,1}))
        disp(sv{1,i+2})
        ParNum_HOSVD = input('Enter # of elements to keep:      ');
        ParNum_HOSVD_vec = [ParNum_HOSVD_vec ParNum_HOSVD];


        if ParNum_HOSVD ~= numel(sv{1,i+2})
            tol_0 = sv{1,i+2}(ParNum_HOSVD+1,1);
        else
            tol_0 = 0; %  Toleracne in case if all the elements of sv{1,i+2} are kept
        end
        tol = [tol, tol_0];
    end
    ParNum_HOSVD_vec = [0 0 ParNum_HOSVD_vec];
else % Auto mode

    ParNum_HOSVD = []; % Excessive variable in the auto mode
    ParNum_HOSVD_vec =[]; % Excessive variable in the auto mode

    % Tolerances for each scheduling parameter
    for i = 1:ParNum
        if sv_num_max_0 < numel(sv{1,i+2})
            tol_0 = sv{1,i+2}(sv_num_max_0+1,1);
        else
            tol_0 = 0; %  Toleracne in case if all the elements of sv{1,i+2} are kept
        end
        tol = [tol, tol_0];
    end
end

% Resulting HOSVD with the tolerance vector calculated from user's choices
try
   [T U sv] = hosvd(S_G, Par_vec, tol);
catch out_mem
   if out_mem.message == 'Out of memory.'
   [T U sv] = hosvd_tall(S_G, Par_vec, tol);
   end 
end 

% Generate convex representation
U_0 = genhull(U,wtype_0);

for ii = 3:numel(Par_vec)
    U_0{1,ii} = abs(U_0{1,ii});
end

% Compute the vertex systems of the convex representation
for ii = 3:numel(Par_vec)
    U_0_Inv{1,ii} = pinv(U_0{1,ii});
end

S_0 = tprod(S_G,U_0_Inv);
Stest = tprod(S_0,U_0);  % Reconstruction of the tensor S_G

% Reconstruction of LTI grid from Stest and models comparison through the Vinnicombe Gap Metric
for ii = 1:Par_range_pr
    LPV_grid_Stest(:,:,ii) = ss(Stest(1:size(LPV_0.Data.a,1),1:size(LPV_0.Data.a,2),ii),...
        Stest(1:size(LPV_0.Data.b,1),size(LPV_0.Data.a,2)+1:end,ii),...
        Stest(size(LPV_0.Data.a,1)+1:end,1:size(LPV_0.Data.a,2),ii),...
        Stest(size(LPV_0.Data.a,1)+1:end,size(LPV_0.Data.a,2)+1:end,ii));

    [gap(ii), nugap(ii)] = gapmetric(LPV_0.Data(:,:,ii),LPV_grid_Stest(:,:,ii));
end
end 