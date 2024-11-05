% Dimension of the state vector containing scheduling parameters.
% The result depends on whether there is any uncertainty in the parameter data

% Output arguments:
% 1. ParNum - Dimension of the state vector  (Parnum is also a reduced dimension keeping non-scalar parameters only)
% 2. Tick_ss - Boolean variable of the function Typo_ss checking for typos in optional input arguments in case if the input system is a ss object.
% It is noted that since ParNum_keep <= ParNum, some immediate manual entering of the input data
% which cannot be stopped until completeness is required, and the typo for checking the mentioned inequality is produced separately.

function [ParNum, LPV_0,Tick_ss] = SchParams(LPV_0,ParNum_keep,wtype,mode,Basis)
if numel(class(LPV_0)) == numel('pss') % 'pss' object
    Tick_ss = nan; % Excessive variable
    ParNum = LPV_0.Domain.NumIV;
else % 'ss' object
    % Define inner structure of the class object SS_Array
    Sch_params = 'IVName';
    Range = 'IVData';
    LPV_00 = SS_Array;
    LPV_00.Data = LPV_0;
    LPV_00.Domain = struct(Sch_params,[],Range,[]);
    LPV_00.samplinggrid = LPV_0.samplinggrid;
    LPV_0 = LPV_00; % Assign 'LPV_0' variable name to the whole class object as it is in case of 'pss' object.
    % Thus, a structure of the class object is built up attaching initial array of SS models (LPV_0) to LPV_00.Data, then LPV_0 is
    % reassigned as a variable name for the whole class
    % Check if there is any typo in the avaliable input arguments
    Tick_ss = Typo_ss(LPV_0,wtype,mode,Basis);
    % Check for the typos
    if Tick_ss == 1
        ParNum = nan; % Excessive variable
        return
    end
    if ~isempty(fieldnames(LPV_0.samplinggrid)) == 1 % Auto filling the necessary data to form an 'pss' object when the 'samplinggrid' property is not empty
        LPV_0_sh_par = LPV_0.samplinggrid;
        LPV_0_sh_par_names = fieldnames(LPV_0_sh_par);
        LPV_0_sh_par_ranges = struct2cell(LPV_0_sh_par);
        ParNum = length(LPV_0_sh_par_names);
        % The rest of the Typo_ss function
        % -----
        % Incorrect writing of number of the scheduling parameters
        if ParNum_keep <= ParNum
            Tick_ss = 0;
            % Proceed
        else
            disp('Warning: incorrect number of kept scheduling parameters')
            Tick_ss = 1;
            ParNum = []; % Excessive variable needed for the output of the function
            return
        end
        % -----
        if  length(size(LPV_0_sh_par_ranges{1, 1})) == 2 % Case when there is only one grid dimension (no ohter dimensions at all OR other dimensions are actually scalars)
            % Check for the 1 parameter in the cell array is enough, because each scheduling parameter has the same dimension using 'samplinggrid' property
            for i = 1:ParNum
                LPV_0.Domain.IVName{i,1} = LPV_0_sh_par_names{i};
                LPV_0.Domain.IVData{i,1} = unique(LPV_0_sh_par_ranges{i, 1}(:,ones(1,length(size(LPV_0_sh_par_ranges{1}))-1))); % Select the data changing along rows for each scheduling parameter
            end
        else % Case when there is at least 2D grid i.e. two scheduling parameters gridded in two separate dimensions)
            for i = 1:ParNum
                LPV_0.Domain.IVName{i,1} = LPV_0_sh_par_names{i};
                dim_params_left = num2cell(ones(1,i-1)); % Left dimensional vector for each scheduling parameter
                dim_params_right = num2cell(ones(1,length(size(LPV_0_sh_par_ranges{1}))-i)); % Right dimensional vector
                if isempty(dim_params_left) == 1
                    LPV_0.Domain.IVData{i,1} = reshape(unique(LPV_0_sh_par_ranges{i, 1}(:,cell2mat(dim_params_right))),[],1);
                elseif isempty(dim_params_right) == 1
                    LPV_0.Domain.IVData{i,1} = reshape(unique(LPV_0_sh_par_ranges{i, 1}(cell2mat(dim_params_left),:)),[],1);
                else
                    LPV_0.Domain.IVData{i,1} = reshape(unique(LPV_0_sh_par_ranges{i, 1}(cell2mat(dim_params_left),:,cell2mat(dim_params_right))),[],1); % ------- ADDRESS THIS
                end
            end
        end
    else % Auto filling in case of empty 'samplinggrid' property
        % Fill the inner structure with the essential data
        ParNum = length(size(LPV_0.Data)) -2; % '-2' subtracts the dimensions of outputs and inputs successively

        % The rest of the Typo_ss function
        % -----
        % Incorrect writing of number of the scheduling parameters
        if ParNum_keep <= ParNum
            Tick_ss = 0;
            % Proceed
        else
            disp('Warning: incorrect number of kept scheduling parameters')
            Tick_ss = 1;
            ParNum = []; % Excessive variable needed for the output of the function
            return
        end
        % -----
        for i = 1:ParNum
            LPV_0.Domain.IVName{i,1} = int2str(i); 
            Sch_params_dims = size(LPV_0.Data);
            LPV_0.Domain.IVData{i,1} = linspace(-1,1,Sch_params_dims(2+i)); % Number of grid points for each scheduling parameter ranging over [-1 1]
            LPV_0.Domain.IVData{i,1} = LPV_0.Domain.IVData{i,1}';
        end
        ParNum = size(LPV_0.Domain.IVData,1);
    end
end
end