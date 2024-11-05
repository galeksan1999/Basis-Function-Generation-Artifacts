% --------------------------------------------------------
% Toy Problem. Rate-dependent stabilization.
% LPV grid-based, rate-bounded controller synthesis
%
% The core model is taken from LPVTools Manual:
% G. Balas, A. Hjartarson, A. Packard, and P. Seiler, LPVTools:
% A Toolbox for Modeling, Analysis, and Synthesis of Parameter Varying
% Control Systems, Software and User s Manual, 2015, p.22-29
%
% Control architecture is given in 
% F. Wu, X. H. Yang, A. Packard, and G. Becker, “Induced L2-norm
%control for LPV systems with bounded parameter variation rates,"
%International Journal of Robust and Nonlinear Control, vol. 6, no.
% 9-10, pp. 983–998, 1996.
%
% Budapest, 2024
%
% Note: Make sure you have the latest version of the G2TPTool added
% to the path for automatic generation of basis functions
% --------------------------------------------------------

clc, clear

%% I. Grid-Based, Rate-Bounded Generalized Plant
% Number of grid points
Num_grid = [7 49];

for i = 1:length(Num_grid)
    % Equidistant Partition of Rate-Bounded Scheduling Parameter
    rho = pgrid('rho',linspace(-pi,pi,Num_grid(i)) );
    rho.RateBounds = [-5 5];

    % State-Space Representation of the Plant
    pcos = cos(rho);
    psin = sin(rho);
    A = [0.75 2 pcos psin;0 0.5 -psin pcos;0 0 -10 0; 0 0 0 -10];
    B = [0 0 0;3 0 0;0 10 0;0 0 10];
    C = [1 0 0 0;0 1 0 0];
    D = zeros(2,3);

    % Grid-Based LPV Plant:
    G = pss(A,B,C,D);

    % Weights For the Generalized Plant
    Wp = eye(2);
    Wn = ss(10*tf([1 10],[1 1000]))*eye(2);
    Wf = 1;
    Wu = (1/280)*eye(2);
    Wr = ss(tf(20,[1 0.2]))*eye(2);

    % Control Interconnection Structure
    systemnames = 'G Wp Wn Wf Wu Wr';
    input_to_G = '[ Wf; u ]';
    input_to_Wp = '[ G-Wr ]';
    input_to_Wn = '[ dn ]';
    input_to_Wf = '[ df ]';
    input_to_Wu = '[ u ]';
    input_to_Wr = '[ dr ]';
    inputvar = '[ df; dr(2); dn(2); u(2)]';
    outputvar = '[ Wu; Wp; G-Wr+Wn ]';

    % Generate generalized plant
    H = sysic;
    H_ary{i} = H;

    %% II. Grid-Based, Rate-Bounded LPV Contoroller Synthesis
    if length(rho.GridData) == Num_grid(1)
        %%% Classical Approach %%%
        % Basis Functions
        b1 = basis(1,0);
        bcos = basis(pcos,'rho',-psin);
        bsin = basis(psin,'rho',pcos);

        % Matrix Functions Used For LMIs
        Xb_cl = [b1;bcos;bsin];
        Yb_cl = Xb_cl;

        % Controller Synthesis
        opt = lpvsynOptions('BackOffFactor',1.02);
        tic
        [~,gamma_cl] = lpvsyn(H,2,2,Xb_cl,Yb_cl,opt);
        t1 = toc;
    end

    if length(rho.GridData) == Num_grid(2)
        %%% Automatic Approach  %%%
        % Call G2TP Function (Auto Mode, CNO Type of Weighting Functions)
        [~,BASIS, ~,U] =  G2TP(H_ary{2},'cno',1,'auto');
        % Reduced grid of basis functions' values
        BASIS_pmat = lpvsplit(BASIS{1}.BasisFunction,'rho',1:(Num_grid(2)-1)/(Num_grid(1)-1):Num_grid(2),'index');
        % Reduced partial derivatives of basis functions
        BASIS_partials = lpvsplit(BASIS{1}.Partials,'rho',1:(Num_grid(2)-1)/(Num_grid(1)-1):Num_grid(2),'index');
        % Rate bounds
        BASIS_pmat.Parameter.rho.RateBounds = rho.RateBounds;
        BASIS_partials.Parameter.rho.RateBounds = rho.RateBounds;
        % Interpolated basis functions
        BASIS = basis(BASIS_pmat, BASIS_partials);
        % Matrix Functions Used for LMIs
        XbGrid_auto = BASIS;
        YbGrid_auto = XbGrid_auto;

        % Controller Synthesis
        tic
        [~,gamma_auto] = lpvsyn(H_ary{1},2,2,XbGrid_auto,YbGrid_auto,opt);
        t2 = toc;
    end
end

%% III. Print the Results
fprintf(['Classical Approach \n   ' ...
    'Elapsed time: %d seconds \n   ' ...
    '       Gamma: %d \n'], t1, gamma_cl)

fprintf(['Automatic Approach \n   ' ...
    'Elapsed time: %d seconds \n   ' ...
    '       Gamma: %d \n'], t2, gamma_auto)