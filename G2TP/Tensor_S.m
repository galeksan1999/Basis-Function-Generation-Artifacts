% Discretized state-space tensor composed of the grid of LTI models

function S_G = Tensor_S(LPV_0,Par_range_pr,ParNum)

for ii = 1:Par_range_pr
    S_G(:,:,ii) = [LPV_0.Data(:,:,ii).a LPV_0.Data(:,:,ii).b; LPV_0.Data(:,:,ii).c LPV_0.Data(:,:,ii).d];
end

% Size vector storing number of the grid points of each scheduling parameter
Par_range_vec = [];
for i = 1:ParNum
    Par_range_vec(1,i) = numel(LPV_0.Domain.IVData{i,1});
end

% Set number of dimensions for S_G and reshape it
dim = [size(LPV_0.Data.StateUnit,1)+size(LPV_0.Data.OutputUnit,1), size(LPV_0.Data.StateUnit,1)+size(LPV_0.Data.InputUnit,1), Par_range_vec];
S_G = reshape(S_G, dim);
end