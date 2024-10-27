% Sorting of the reduced qLPV model's elements
%
% Output arguments
% 1. I_perm_ordered - row vector containing all the indices put in the correct order
% 2. Par_range_pr_hull - product of dimensions of the reduces qLPV model (total number of grid points defined over the multidimensional 
% domain under consideration)
 
function [I_perm_ordered,Par_range_pr_hull] = Indices_sort(ParNum,U_weigh_fun,I_perm)

Par_range_mat_hull = zeros(1,ParNum);
for t = 1:ParNum
    Par_range_mat_hull(1,t) = size(U_weigh_fun{1,2+t},2);
end
Par_range_pr_hull = prod(Par_range_mat_hull , 'all');

% Initial row vector 
Product_array = nan; % Predefined loop variable
for ii = 1:ParNum - 1
    Product(1,ii) = prod(Par_range_mat_hull(1,1:ii));
    for j = 1:numel(I_perm{1,ii+1})
        Product_2 = Product(1,ii)*(I_perm{1,ii+1}(j,1)-1);
        Product_array = [Product_array Product_2];
    end
end
Product_array = [I_perm{1,1}' Product_array(2:end)]; % Elimination of the test value  + adding the values of
% I_perm{1,1}. All the indeces for the correct order of the elements can be obtained from this vector by providing
% the following summation

% Construct the vector for the proper partition
for i = 1:ParNum
    I_perm_numbers_cell_distr_vec(1,i) = numel(I_perm{1,i});
end

% Partition
I_perm_numbers_cell = mat2cell(Product_array,[1],I_perm_numbers_cell_distr_vec);

% Putting indices in the correct order

% Storages for the metadata and final data (define a cell array elements  of which is to be
% updated with each iteration)
Indices = [];

if numel(I_perm_numbers_cell) == 1
    Indices = I_perm_numbers_cell{1,end}; % Full range of the indeces in case of 1 scheduling parameter
elseif numel(I_perm_numbers_cell) == 2
    % Define a cell array to start a loop in case of multiple parameters
    Indices_0 = cell(1,numel(I_perm_numbers_cell{1,2}));

    for i = 1:numel(I_perm_numbers_cell{1,2})
        Indices_0{1,i} = I_perm_numbers_cell{1,1} + I_perm_numbers_cell{1,2}(1,i);
    end

    Indices = Indices_0;
    Indices = cell2mat(Indices); % Full range of the indeces in case of 2 scheduling parameters
else
    Indices_0 = cell(1,numel(I_perm_numbers_cell{1,2}));

    for i = 1:numel(I_perm_numbers_cell{1,2})
        Indices_0{1,i} = I_perm_numbers_cell{1,1} + I_perm_numbers_cell{1,2}(1,i);
    end
    % Summation to get the full range of indices in case of multiple
    % scheduling parameters
    for i = 1:numel(I_perm_numbers_cell) - 2
        for j = 1:numel(I_perm_numbers_cell{1,i+2})
            for jj = 1:numel(Indices_0)
                Indices_00= Indices_0{1,jj} + I_perm_numbers_cell{1,i+2}(1,j);
                Indices = [Indices Indices_00];
            end
        end
        Indices_0 = Indices;
        Indices_0 = mat2cell(Indices_0,[1],size(Indices_0,2));
        if i ~= numel(I_perm_numbers_cell) - 2
            Indices = [];
        else
            % Nothing changes, so we keep the indices from the last
            % iteration. Previous data should be eraised with each
            % iteration by means of 'Indices = []' command
        end
    end
end

I_perm_ordered = Indices;
end 