function gopt = setgaitoptfcn(s,gopt)
    %% Set Objective and Constraint functions
    gopt.obj_fcn = @(y) solvedifffmincon(y,s,gopt);
    gopt.con_fcn = @(y) nonlcon(y,s,gopt);
end