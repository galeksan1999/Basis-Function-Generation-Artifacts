% Grid-based LPV model with the uncertainty included
%
% It may happen that the data read from some sensors contains uncertainty, so we consider the associated scheduling parameters as irrelevant 
% and want to combine them in terms of the stored data to reduce the dimension complexity.
%
% Output arguments 
% 1. Par_range_pr - product of all grid points of the corresponding scheduling parameters in each dimension
% 
% Reshape is done when the condition ParNum_keep ~= ParNum is fulfilled.	

function [LPV_0,Par_range_pr,ParNum] = sys_reshape(LPV_0,ParNum_keep,ParNum)

% All ranges of the scheduling parameters and their product
Par_range_mat = zeros(1,ParNum);
for t = 1:ParNum
    Par_range_mat(1,t) = size(LPV_0.Domain.IVData{t,1},1);
end

Par_range_pr = prod(Par_range_mat , 'all');

% Reshape of the initial system due to the uncertainty
if ParNum_keep ~= ParNum

    dim_LPV_0_GB= [Par_range_mat(1:ParNum_keep) prod(Par_range_mat(ParNum_keep+1:end))];
    LPV_0_GB_res = reshape(LPV_0.Data, dim_LPV_0_GB);
    Par_Name_res = cell(1,ParNum_keep+1); % cell array made up of names of the kept parameters and merged ones
    % Filling the cell array with the names of the kept parameters separately
    for i = 1:ParNum_keep
        Par_Name_res{1,i} = LPV_0.Domain.IVName{i,1};
    end
    % Filling the last cell of the array with the name of the merged uncertain parameters
    for i = 1:ParNum-ParNum_keep
        Par_Name_res{1,ParNum_keep+1}(1,i) = convertCharsToStrings(regexprep(LPV_0.Domain.IVName{ParNum_keep+i,1},'_','')); 
    end

    % Grid points of the parameters
    for i = 1:ParNum_keep+1
        % Restriction: range of the uncertain parameters must be the
        % same to keep dimension consistency
        Par_grid_uncert{1,i} = pgrid(convertStringsToChars(Par_Name_res{1,i}),LPV_0.Domain.IVData{i,1}');
    end
    Par_domain_uncert = rgrid(Par_grid_uncert{:,:});
    if numel(class(LPV_0)) == numel('pss')
        for i = 1:ParNum_keep
            %Remark: rate bounds for uncertain parameters are defined as infinity, since these parameters are not essential
            Par_domain_uncert.IVRateBounds(i,:) = [LPV_0.Domain.IVRateBounds(i,1)   LPV_0.Domain.IVRateBounds(i,2)];
        end
    end
    LPV_0 = pss(LPV_0_GB_res,Par_domain_uncert);
    ParNum = ParNum_keep + 1;
end
end 