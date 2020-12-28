function [runtime, sol, CS, flagofGlobal, iter, outindex] = SloveQuasiProm(method, data, dim, solInit, epsilon, timeLimit, gap)

    if nargin < 7
        gap = -1;
    end

    if nargin < 6
        timeLimit = -1;
    end
    

    if strcmp(method, 'GoIA')
        tic;
        [sol, CS, iter, flagofGlobal, outindex] = GoIA(data, solInit, epsilon, gap, timeLimit); 
        runtime = toc;
    else
        error('GoIA-LR and GoIA-UN are for linear model estimation! Please select the method of GoIA for non-linear model fitting problems.');
    end

end