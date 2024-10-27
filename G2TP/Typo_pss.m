% Typo check function having a 'pss' object as an input grid of LTI system
% Tick_pss: = Boolean variable

function Tick_pss = Typo_pss(LPV_0,wtype,ParNum_keep,ParNum,mode,Basis)

if numel(class(LPV_0)) == numel('SS_Array') % 'ss' object
    Tick_pss = nan;
    return

elseif numel(class(LPV_0)) == numel('pss') % 'pss' object
    % Incorrect writing of the type of the weighting functions
    wtype_list = {'eye', 'ortho', 'snnn', 'cno', 'irno', 'box'};

    for i = 1:numel(wtype_list)
        if i == 6 && numel(wtype) ~= numel(wtype_list{i}) || i == 6 && min(wtype == wtype_list{i}) ~= 1 % last iteration; either word length is different or there is at least one different letter in the corresponding position
            disp('Warning: incorrect type of weighting functions')
            Tick_pss = 1;
            return
        elseif numel(wtype) == numel(wtype_list{i}) && min(wtype == wtype_list{i}) == 1
            Tick_pss = 0;
            break
        else
            continue
        end
    end

    % Incorrect writing of number of the scheduling parameters
    if ParNum_keep <= ParNum
        Tick_pss = 0;
        % Proceed
    else
        disp('Warning: incorrect number of kept scheduling parameters')
        Tick_pss = 1;
        return
    end

    % Incorrect writing of the mode
    if not(numel(mode) == numel('manual') && min(mode == 'manual') == 1 || numel(mode) == numel('auto') && min(mode == 'auto') == 1) % Switch to the manual mode
        disp('Warning: incorrect type of the selection mode')
        Tick_pss = 1;
        return
    end

    % Incorrect writing of the tick for basis functions generation
    if Basis ~= 0 && Basis ~= 1
        disp(strcat('Warning: use logical', " ", string(1), " ", 'to generate basic functions' ))
        Tick_pss = 1;
        return
    end
end
end