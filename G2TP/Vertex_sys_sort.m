% Sort the vertex systems according to the scheduling parameters
% 
% The goal is to sort properly the data of the reduced LPV model since after obtaining a polytopic form, 
% the order of the data connecting weight matrices in each dimension with the reduced qLPV model is lost
%
% Output arguments
% 1. idx_u_asc - cell array containing sequences of vertices put in an ascending order in each cell
% 2. I_perm - cell array containing an arrangement of the elements of
% idx_u_asc in each cell (see sort)

function [idx_u_asc, I_perm] = Vertex_sys_sort(ParNum,Par_vec,U_weigh_fun)

for iii = 3:numel(Par_vec)
    for ii = 1:size(U_weigh_fun{1,iii},2)
        [u_max{1,iii-2}(ii) index_u{1,iii-2}(ii)] = max(U_weigh_fun{1,iii}(:,ii)); % Thus, the data related to different
        % scheduling parameters is separated by means of cells. cell{1,1} - data
        % related to the 1st parameter, cell{1,2} - 2nd pararameter and so on and so forth. Cells are
        % used to account for different sizes of row vectors
    end
end

for  iii = 1:ParNum
    [idx_u_asc{1,iii} I_perm{1,iii}] = sort(index_u{1,iii}(:),'ascend');
end
end