function [result, REConsensusSet, iter, flagofGlobal, outindex] = GoIA_LR_Linear(data, solInit, epsilon, gap, timeLimit)
    
    if gap == -1
        error('gap cannot be negtive!');
    end

    cd(['src/utlis/GoIA_LR/LinearLR'])

    [result, iter, flagofGlobal]=Search(data, epsilon, gap, solInit, timeLimit);
    REConsensusSet = sum(abs(data(end,:)-result'*data(1:end-1, :))<epsilon+eps);
    outindex = abs(data(end,:)-result'*data(1:end-1, :))>epsilon+eps;
    cd('../../../..')

end