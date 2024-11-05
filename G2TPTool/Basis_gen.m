% Basic functions generation
% Based on employing Piecewise Cubic Hermite Interpolating Polynomials (PCHIPs)
% 
% Steps of the approximation of the weight functons:
% 1. Applying PCHIPs for an initial approximation
% 2. Grid refinement that is based on the data obtained from the initial approximation 
% and performed two times following iterative process to make further calculations more accurate
% 3. Numerical differentiation on a refined grid. For the last point, a first order derivative is evaluated using 
% the analytical form of the corresponding PCHIP

function [BASIS_GRID,BASIS_TP] = Basis_gen(LPV_0,U,Par_grid_red_0,idx_u_asc,i_check,BASIS_GRID,BASIS_TP)

% Default settings
n_0 = 1; % grid step for the getting the fitting curve without grid refinement
Num_ref_0 = 2; % Number of refinements
Counter_refinement_0 = 2; % scale ratio in each refinement iteration

for i = i_check:size(U,2)
    figure
    hold on
    for ii = 1:size(U{1,i},2)
        % Restore loop variables
        Num_ref = 1; % Number of grid refinements
        n = 1; % initial grid step == grid of the system composed of weighting functions (:= integers)
        Counter_refinement = Counter_refinement_0; 

        % X and Y data
        x = 1:n_0:size(U{1,i},1);
        % -----
        % ' is needed to keep the consistency in dimensions
        x = x';
        y = U{1,i}(:,ii);
        % -----

        % PCHIP fitting curve
        f=fit(x,y,'pchipinterp');

        % Grid refinement
        while Num_ref <= Num_ref_0

            % Update the grid step
            n = n/Counter_refinement;
            % Clear storage variables before every new refinement
            x_refined = [];
            y_refined = [];

            x_refined = 1:n:size(U{1,i},1);
            x_refined = x_refined';

            for iii = 1:numel(y) - 1
                y_refined_0 = [y(iii); mean(y(iii:iii+1))];
                y_refined = [y_refined; y_refined_0];
                if iii == numel(y) - 1
                    y_refined = [y_refined; y(end)]; % add the last element to the vector
                end
            end

            % Reassign refined variables and obtain PCHIP fitting curve
            x = x_refined;
            y= y_refined;
            f=fit(x,y,'pchipinterp');

            % Iteration step
            Num_ref = Num_ref + 1;
        end

        % Results
        x_orig = 1:1:size(U{1,i},1);
        plot(x_orig,U{1,i}(:,ii),'LineStyle','-'); % Original curve
        plot(f,x,y); % Approximated curve
        axis([1 size(U{1,i},1) -0.1 1]); axis square
        %legend 'Original curve' 'Approximated curve'
        legend off

        ax = gca;
        h = findobj(gca,'Type','line');

        set(gca,'YLabel',[]); set(gca,'XLabel',[]);

        x_plot = h(1).XData; % Domain of the approximate function
        y_plot = h(1).YData; % Values of the approximate function

        Der = diff(y_plot)./diff(x_plot); % First derivatives according to the forward Euler method evaluated at each grid point except the last one
        % For the last point, the corresponding analytical PCHIP polynimial defined in the last segment is used

        f_piecewise_coeffs = f.p.coefs; % Number of rows := number of approximated segments (number of 3-order polynomials for each segment)
        % PCHIP Family (piece-wise cubic polynomials): f(x) = a_3(x-x_r)^2(x-x_l) + a_2(x-x_r)(x-x_l)^2 + a_1(x-x_l) + a_0, where x_l, x_r - jont left and right gridpoints, respectively. Entries in the cells - polynomial coefficients a_i

        % First derivative evaluated at the last grid point
        %  f'(x) =  a_3*3*x^2+a_2*2*x+a_1

        syms x_var
        % Symbolic function (cubic polynomial defined on the last segment)
        Polynomial = f_piecewise_coeffs(end,1)*(x_var-x_plot(end))^2*(x_var-x_plot(end-1)) +  f_piecewise_coeffs(end,2)*(x_var-x_plot(end))*(x_var-x_plot(end-1))^2 + f_piecewise_coeffs(end,3)*(x_var-x_plot(end-1)) + f_piecewise_coeffs(end,4);
        % First derivative of 'Polynomial'
        Polynomial_der = double(subs(diff(Polynomial,x_var,1),   x_var, x_plot(end)));

        % Complete vector comprised by first derivatives evaluated at each grid point
        Der = [Der Polynomial_der];
        % Leave only the derivatives evaluated at the initial grid poinst
        % (which are integers),  sort the y data and put everything in the
        % required form to execute 'basis' command
        for p = 1:numel(x_orig)
            Data_0(ii,:,p) = U{1,i}(p,ii);
            Der_upd(ii,:,p) = Der(find(x_plot == x_orig(p)));
        end
    end
    % -----
    % Basis functions for TP model

    % Domain of each scheduling parameter based on the grid points
    Par_domain_red = rgrid(Par_grid_red_0);

    %Basis generation for TP model
    Data_TP = pmat(Data_0(:,:,idx_u_asc{1,i}),Par_domain_red);
    Partials_TP = pmat(Der_upd(:,:,idx_u_asc{1,i}),Par_domain_red);

    BASIS_TP{i} = basis(Data_TP,Partials_TP);

    % Basis functions for the grid-based model
    % -----
    % Domain of each scheduling parameter based on the grid points
    Par_grid_full_0 = pgrid(LPV_0.Domain.IVName{i,1}, LPV_0.Domain.IVData{i,1});
    Par_domain_full = rgrid(Par_grid_full_0);
    %Basis generation for TP model
    Data_grid = pmat(Data_0,Par_domain_full);
    Partials_grid = pmat(Der_upd,Par_domain_full);

    BASIS_GRID{i} = basis(Data_grid,Partials_grid);
    % -----

    % if not the last iteration, go back to create Par_grid_red_0 for the
    % next scheduling parameter
    if i ~= size(U,2)
        break
    end
end
hold off
end

