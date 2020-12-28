function [result, REConsensusSet, iter, flagofGlobal, outindex] = GoIA_UN_Linear(data, solInit, epsilon, gap, timeLimit)

    if gap == -1
        error('gap cannot be negtive!');
    end
    
    cd(['src/utlis/GoIA_UN/LinearUN'])
    
    [result, iter, flagofGlobal]=Search(data, epsilon, gap, solInit, timeLimit);
    REConsensusSet = sum(abs(result'*data(1:end,:))<epsilon+eps); 
    outindex = abs(result'*data(1:end,:))>epsilon+eps; 
    cd('../../../..')
    
end