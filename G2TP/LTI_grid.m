% LTI grid from the vertices

function LPV_TP = LTI_grid(Par_range_pr_hull,S_perm,LPV_0,Par_range_vec_red)

% Array of the SSMs
for ii = 1:Par_range_pr_hull
    LPV_TP_undone(:,:,ii) = ss(S_perm(1:size(LPV_0.Data.a,1),1:size(LPV_0.Data.a,1),ii),...
        S_perm(1:size(LPV_0.Data.a,1),size(LPV_0.Data.a,1)+1:end,ii),...
        S_perm(size(LPV_0.Data.a,1)+1:end,1:size(LPV_0.Data.a,1),ii),...
        S_perm(size(LPV_0.Data.a,1)+1:end,size(LPV_0.Data.a,1)+1:end,ii));
end

% Reshape of LPV_TP according to dimensions of the required array of models
if size(Par_range_vec_red,2) == 1
    LPV_TP = LPV_TP_undone;
else
    LPV_TP = reshape(LPV_TP_undone, Par_range_vec_red);
end
end
