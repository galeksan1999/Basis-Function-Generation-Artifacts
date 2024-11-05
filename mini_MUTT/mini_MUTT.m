% --------------------------------------------------------
% miniMUTT model
%
% For control architecture, see
% Tensor Product Type Polytopic LPV Model-based Active Flutter Suppression Design
% by Bela Takarics and Balint Vanek
%
% Budapest, 2024
%
% Note: Make sure you have the latest version of the G2TPTool added
% to the path for automatic generation of basis functions
% --------------------------------------------------------

clc, clear

%% I. Grid-Based, Rate-bounded Generalized Plant
% Load State-Space Model
load('miniMUTT_ss')

% Weights For the Generalized Plant

W1 = 1;
W2 = blkdiag(11,150,200);
Wy = 0.0001*blkdiag(3,3,6);
Wz = 0.001;
sys1 = tf([1 2*2*pi*3 (2*pi*3)^2],2*[1  2*2*pi*3/sqrt(200*2) (2*pi*3/sqrt(200*2))^2]);
sys2 = tf(200*[1 2*2*pi*6 (2*pi*6)^2],[1  2*2*pi* 6*55/2.75 (2*pi* 6*55/2.75)^2]);
Wu = sys1*sys2;

% Generate Generalized Plant
Gweighted = genplant(miniMUTT_ss,W1,W2,[], Wy, Wu, Wz);

% Equidistant Partition of Rate-Bounded Scheduling Parameter
Vinf = miniMUTT_ss.SamplingGrid.Vinf';
VGrid = pgrid('V',Vinf);
DomainGrid = rgrid(VGrid);
DomainGrid.V.IVRateBounds = [-3 3];

% Grid-Based Generalized LPV Plant
H = pss(Gweighted,DomainGrid);

%% II. Grid-Based, Rate-Bounded LPV Contoroller Synthesis
%%% Classical Approach %%%
% Basis Functions
f1 = basis(1,0);
f2 = basis(VGrid,1);
f3 = basis(VGrid^2,2*VGrid);

% Matrix Functions Used For LMIs
Xb_cl = [f1;f2;f3];
Yb_cl = Xb_cl;

% Controller Synthesis
opt = lpvsynOptions('BackOffFactor',1.15);
tic
[~,gamma_cl] = lpvsyn(H,3,1,Xb_cl,Yb_cl,opt);
t1 = toc;

%%% Automatic Approach  %%%
% Call G2TP Function (Auto Mode, CNO Type of Weighting Functions)
[~, BASIS_GRID, ~, U] = G2TP(H, 'cno', 1, 'auto'); % For Manual Mode, Change 4th Argument to 'manual'

% Matrix Functions Used for LMIs
XbGrid_auto = [BASIS_GRID{1}];
YbGrid_auto = XbGrid_auto;
tic
[~,gamma_auto] = lpvsyn(H,3,1,XbGrid_auto,YbGrid_auto,opt);
t2 = toc;

%% III. Print the Results
fprintf(['Classical Approach \n   ' ...
    'Elapsed time: %d seconds \n   ' ...
    '       Gamma: %d \n'], t1, gamma_cl)

fprintf(['Automatic Approach \n   ' ...
    'Elapsed time: %d seconds \n   ' ...
    '       Gamma: %d \n'], t2, gamma_auto)