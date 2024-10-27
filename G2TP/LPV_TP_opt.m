% Class specially defined for the 'LPV_TP_R_F' function for the default input data

classdef LPV_TP_opt
    properties
        wtype = 'cno'; % ~= 'cno' ---> other type
        ParNum_keep = 1; % ~= 1 ---> other number
        SVmode = 'auto'; % ~= 'auto' ---> manual mode (type in 'manual')
        Tolerance = 10^(-2); % Max allowed distance between the original model and the reduced one
        Pause = 0; % Pause time to make calculations readable on the go
        Tick_basis = 1; % Enabling for basis functions genertaion
    end
end