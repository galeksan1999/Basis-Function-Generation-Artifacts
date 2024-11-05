% Input data formation
%
% First, the default settings stored in a separate class LPV_TP_opt are assigned. 
% Then, if some optional arguments are specified, the input data is redefined.

function [wtype, ParNum_keep, mode,TOL_NORM,Pause_T,Basis] = Num_in(nargin,varargin)

% Set of the default input data
varargin_opt = LPV_TP_opt;

% Preassign default inputs
ParNum_keep = varargin_opt.ParNum_keep;
wtype = varargin_opt.wtype;
mode = varargin_opt.SVmode;
TOL_NORM = varargin_opt.Tolerance;
Pause_T = varargin_opt.Pause;
Basis = varargin_opt.Tick_basis;

varargin = varargin{1};

if nargin == 1 % case: nothing is predefined in 'varargin'
    % Proceed with the predefined default settings
elseif  nargin == 2 % case: wtype is predefined by the user as well
    wtype = varargin{1};
elseif nargin == 3  % case: ParNum_keep is predefined by the user whether he also specifies wtype or not as well
    wtype = varargin{1};
    ParNum_keep = varargin{2};
elseif nargin == 4 % case: mode of auto selection of singular values is predefined by the user as well
    wtype = varargin{1};
    ParNum_keep = varargin{2};
    mode = varargin{3};
elseif nargin == 5  % case: tolerance between the original systems and the reduced one is predefined by the user as well
    wtype = varargin{1};
    ParNum_keep = varargin{2};
    mode = varargin{3};
    TOL_NORM = varargin{4};
elseif nargin == 6 % case: pause time between some calculations is predefined by the user as well
    wtype = varargin{1};
    ParNum_keep = varargin{2};
    mode = varargin{3};
    TOL_NORM = varargin{4};
    Pause_T = varargin{5};
elseif nargin == 7 % case: basic functions generation is predefined by the user as well
    wtype = varargin{1};
    ParNum_keep = varargin{2};
    mode = varargin{3};
    TOL_NORM = varargin{4};
    Pause_T = varargin{5};
    Basis = varargin{6};
end
end