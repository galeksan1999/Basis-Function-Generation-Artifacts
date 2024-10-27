% Reduced grid-based LPV model
% BASIS_GRID and BASIS_TP are obtained utilizing Basis_gen function

% Output argument:
% 1. Tick_bad_precision - check for identical grid points in the reduced model for each scheduling parameter

function [LPV_TProd,BASIS_GRID,BASIS_TP,Tick_bad_precision] = LPV_TP_red(ParNum,LPV_0,idx_u_asc,LPV_TP,Basis,U)

% Define strorages for basic functions and the logical variable
BASIS_TP = {};
BASIS_GRID = {};
Tick_bad_precision = 0;

% Check for identical grid points in the reduced model for each scheduling parameter
for i = 1:ParNum
    if length(idx_u_asc{1,i}) ~= length(unique(idx_u_asc{1,i}))
        disp('Warning: Identical grid points in the reduced model are detected')
        disp('Warning: It is recommended to switch to the manual mode')
        Tick_bad_precision = 1;
        break
    end
end

if Tick_bad_precision == 1
    LPV_TProd = nan; % Excessive variable # 1
    BASIS_GRID = nan; % Excessive variable # 2
    BASIS_TP = nan; % Excessive variable # 3
    return
end

% Add unities in the dimensions associated with scalar scheduling
% parameters. It is required to assign time-varying grid points below
% Recover the initial size of ParNum if there are any scalars associated with the corresponding dimensions
if ParNum ~= length(LPV_0.Domain.IVName)
    for i = ParNum+1:length(LPV_0.Domain.IVName)
        idx_u_asc{1,i} = 1;
    end
    ParNum = length(LPV_0.Domain.IVName);
end

% Time-varying grid points
for i_check = 1:ParNum
    Par_grid_red{1,i_check} = pgrid(LPV_0.Domain.IVName{i_check,1},LPV_0.Domain.IVData{i_check,1}(idx_u_asc{1,i_check}));
    % Time-varying grid points for individual scheduling paramater in order to
    % create basic functions separately for them
    % (there is no optimal way to create basic functions when Par_grid_red is composed of the data from all scheduling parameters)
    Par_grid_red_0 = pgrid(LPV_0.Domain.IVName{i_check,1},LPV_0.Domain.IVData{i_check,1}(idx_u_asc{1,i_check}));
    if Basis == 1
        % Basic functions generation
        [BASIS_GRID,BASIS_TP] = Basis_gen(LPV_0,U,Par_grid_red_0,idx_u_asc,i_check,BASIS_GRID,BASIS_TP);
    end
end

% Domain based on the grid points
Par_domain_red = rgrid(Par_grid_red{:,:});

% Rate bounds
if numel(class(LPV_0)) == numel('pss')
    for i_check = 1:ParNum
        Par_domain_red.IVRateBounds(i_check,:) = [LPV_0.Domain.IVRateBounds(i_check,1)   LPV_0.Domain.IVRateBounds(i_check,2)];
    end
end

% Reduced grid-based qLPV model
LPV_TProd = pss(LPV_TP,Par_domain_red);

% Construct the reduced sampling grid
for i = 1:ParNum
    Sampling_grid_red{i}= Par_grid_red{i}.GridData;  
end 
syms z [1 ParNum]; 
eval([char(z), '= ndgrid(Sampling_grid_red{:});'])

for i = 1:ParNum
    eval(['Sampling_gr{i} = ','z',num2str(i),';'])
end

Sampling_grid_red_str = cell2struct(Sampling_gr, regexprep(LPV_TProd.Domain.IVName,' ','_'),2);

% Assign state names, input names, and output names of the initial grid-based qLPV model to the corresponding
% character variables of the reduced model
LPV_TProd.Data.StateName = LPV_0.Data.StateName;
LPV_TProd.Data.InputName = LPV_0.Data.InputName;
LPV_TProd.Data.OutputName = LPV_0.Data.OutputName;
LPV_TProd.Data.samplinggrid = Sampling_grid_red_str;
end