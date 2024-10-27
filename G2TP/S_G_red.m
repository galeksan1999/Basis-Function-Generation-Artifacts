% Ordered reduced S_G tensor based on the weight functions
%
% Output arguments
% 1. S_perm - reduced aLPV midel with permuted elements
% 2. Par_range_vec_red - row vector containing dimensions of the reduced model
% 3. U - weighting matrix with permuted elements in each dimension

function [S_perm,Par_range_vec_red,U]  = S_G_red(S_weigh_fun,I_perm_ordered,I_perm,ParNum,Par_vec,U_weigh_fun,LPV_0)

S_perm = S_weigh_fun(:,:,I_perm_ordered);

% Size vector storing the reduced number of the grid points of each scheduling parameter
Par_range_vec_red = [];
for i = 1:ParNum
    Par_range_vec_red(1,i) = numel(I_perm{1,i});
end

% Set the number of dimensions for the tensor Scno_perm and reshape it
dim_2 = [size(LPV_0.Data.StateUnit,1)+size(LPV_0.Data.OutputUnit,1), size(LPV_0.Data.StateUnit,1)+size(LPV_0.Data.InputUnit,1), Par_range_vec_red];

S_perm = reshape(S_perm, dim_2);

% Sort the Ucno array according to the permutation (i.e. I_perm)
for iii = 3:numel(Par_vec)
    for ii = 1:size(U_weigh_fun{1,iii},2)
        U{1,iii-2}(:,ii) = U_weigh_fun{1,iii}(:,I_perm{1,iii-2}(ii));
    end
end
end 