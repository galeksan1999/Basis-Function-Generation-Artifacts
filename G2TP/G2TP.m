function varargout = G2TP(LPV_0,varargin)

% G2TP executes quasi-linear parameter-varying (qLPV) grid-based model reduction through the Tensor Product (TP) model
% transformation. 
%
% Structure
% [LPV_TProd, BASIS_GRID, BASIS_TP, U] =  G2TP(LPV_0, wtype, ParNum_keep,mode,TOL_NORM,Pause_T,Basis)
%
% Input arguments
% 1. LPV_0 - grid-based qLPV model
% 2. wtype - type of the weight functions (see the genhull function)
% 3. ParNum_keep - number of the kept scheduling parameters which are assumed to be measurable with no uncertainty
%
% It is assumed that parameters to be kept locate in the beginning and there are no uncertain parameters between them. 
% E.g. let A_1xA_2xA_3x.....xA_N be a multidimensional grid domain of linear time invariant (LTI) models
% where A_1xA_2xA_3x....A_n are the grid points of the corresponding scheduling parameters to be kept (n<=N)
% 
% 4. mode - type of model reduction which is based on singular values
% selection ('auto' by default or 'manual')
%
% Auto mode is built on two assumptions:
% 4.1. Maximum and minimum possible numbers of the kept singular values are 6 and 2, respectively
% 4.2. Reduced qLPV model is compared to the unmodified one after each iteration on the ground of Vinnicombe Gap Metric metric 
% which calculates the distance between each pair of LTI systems
%
% Manual mode allows the user to specify how many singular values are to be kept for each scheduling parameter individually. 
% This is done through interaction with the Command window.
%
% 5. TOL_NORM - upper bound of the Vinnicombe Gap Metric between the initial
% model and reduced one after executing  reduced high order singular value decomposition (RHOSVD)
% 6. Pause_T - pause between specific calculations
% 7. Basis - Boolean variable enabling generation of basis functions
%
% Default settings stored in a separate class LPV_TP_opt
%
% Output arguments 
% 1. LPV_TProd - resulting reduced grid-based qLPV model
% 2. BASIS_GRID - basis functions of the input qLPV model
% 3. BASIS_TP - basis functions of the reduced qLPV model
% 4. U - weighting matrix with permuted elements in each dimension
%
% LPVTools and TP Tool toolboxes are used for system design  

% Input data function
[wtype, ParNum_keep, mode,TOL_NORM,Pause_T,Basis] = Num_in(nargin,varargin);

% Number of scheduling parameters
[ParNum, LPV_0,Tick_ss] = SchParams(LPV_0,ParNum_keep,wtype,mode,Basis);

% Check if there is any typo in the input arguments 
Tick_pss = Typo_pss(LPV_0,wtype,ParNum_keep,ParNum,mode,Basis);

% Tick for typos
if Tick_pss == 1 || Tick_ss == 1
    return
end 

% Grid-based qLPV model with the uncertainty included
[LPV_0,Par_range_pr,ParNum] = sys_reshape(LPV_0,ParNum_keep,ParNum);

% TP-type polytopic HOSVD-based qLPV model
% --------------------------------------------------------------------------------------------
% Discretized state-space tensor 
S_G = Tensor_S(LPV_0,Par_range_pr,ParNum);

% Finite element polytopic form 
[U_weigh_fun,S_weigh_fun,Par_vec,ParNum] = sv_selection(ParNum,ParNum_keep,mode,TOL_NORM,S_G,Par_range_pr,LPV_0,wtype,Pause_T);

% Sort the vertex systems according to the scheduling parameters
[idx_u_asc, I_perm] = Vertex_sys_sort(ParNum,Par_vec,U_weigh_fun);

% Sorting algorithm of the elements' indices
[I_perm_ordered,Par_range_pr_hull] = Indices_sort(ParNum,U_weigh_fun,I_perm);

% Ordered reduced state-space tensor based on the weight functions
[S_perm,Par_range_vec_red,U]  = S_G_red(S_weigh_fun,I_perm_ordered,I_perm,ParNum,Par_vec,U_weigh_fun,LPV_0);

% LTI grid made up of the vertices
LPV_TP = LTI_grid(Par_range_pr_hull,S_perm,LPV_0,Par_range_vec_red);
% --------------------------------------------------------------------------------------------

% Reduced grid-based qLPV model
[LPV_TProd,BASIS_GRID,BASIS_TP,Tick_bad_precision] = LPV_TP_red(ParNum,LPV_0,idx_u_asc,LPV_TP,Basis,U);

% Tick for identical grid points in the reduced model
if Tick_bad_precision == 1
    return
end 

% Display the outputs
Output(S_G,Par_vec,ParNum,LPV_0,U,Pause_T);

% Optional output arguments 
varargout{1} = LPV_TProd;
varargout{2} = BASIS_GRID;
varargout{3} = BASIS_TP;
varargout{4} = U;
end

% Add warnings, errors, and maybe cases 