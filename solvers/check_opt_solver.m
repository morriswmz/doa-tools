function check_opt_solver(name)
%CHECK_OPT_SOLVER Checks if specified solver is present.
switch lower(name)
    case 'quadprog'
        if ~exist('quadprog.m', 'file')
            error('Cannot find quadratic programming optimizer. Please make sure the optimization toolbox is installed.');
        end
    case 'sdpt3'
        if ~exist('sqlp.m', 'file')
            error('Cannot locate SDPT3 solver. Please make sure you installed it correctly.');
        end
    case 'cvx'
        if ~exist('cvx_begin.m', 'file')
            error('Cannot locate CVX installation. Please make sure you have CVX installed.');
        end
end
end

