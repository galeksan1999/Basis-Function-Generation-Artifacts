% Output data

function Output(S_G,Par_vec,ParNum,LPV_0,U,Pause_T)

% List of singular values of each scheduling parameter
try
   [~, ~, sv_disp] = hosvd(S_G, Par_vec);
catch out_mem
   if out_mem.message == 'Out of memory.'
   [~, ~, sv_disp] = hosvd_tall(S_G, Par_vec);
   end 
end
pause(Pause_T)
format long
for i = 1:ParNum
    disp(strcat('Singular values of:',"      ",LPV_0.Domain.IVName{i,1}))
    disp(sv_disp{1,i+2})
end
format short
% Weighting functions of the reduced model
plothull(U);
end 