function [runtime, sol, CS, flagofGlobal, iter, outindex] = SloveLinearProm(method, data, dim, solInit, epsilon, timeLimit, gap)

    if nargin < 7
        gap = -1;
    end
    if nargin < 6
        timeLimit = -1;
    end
    if dim+1 ~= size(data, 1) && dim ~= size(data, 1)
        error('The input data shoule be dim*N !');
    end
    
    tic;
    if strcmp(method, 'GoIA-LR')
        if dim+1 ~= size(data, 1)
            data = [ data(1:end-1, :); ones(1, size(data, 2)); data(end, :) ]; 
        end        
        [sol, CS, iter, flagofGlobal, outindex] = GoIA_LR_Linear(data, solInit, epsilon, gap, timeLimit); 
    elseif strcmp(method, 'GoIA-UN')
         if dim+1 ~= size(data, 1)
            data = [data; ones(1, size(data, 2))];
         end
         [sol, CS, iter, flagofGlobal, outindex] = GoIA_UN_Linear(data, solInit, epsilon, gap, timeLimit);
    end
    runtime = toc;
    
end