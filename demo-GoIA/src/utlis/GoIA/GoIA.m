function [result, REConsensusSet, iter, flagofGlobal, outindex] = GoIA(data, solInit, epsilon, gap, timeLimit)
    if gap == -1
        error('gap cannot be negtive!');
    end

    cd(['src/utlis/GoIA/Quasi'])
    
    N = size(data, 2);
    A = []; b = []; c = []; d = [];
    for i = 1 : N
        A = [A, data(i).A];
        b = [b, data(i).b];
        c = [c, data(i).c];
        d = [d, data(i).d];
    end
    [result, iter, flagofGlobal] = Search(A, b, c, d, epsilon, gap, solInit, timeLimit);

    r = zeros(N, 1);
    dim = size(A, 2) / N;
    for i = 1 : N
        r(i) = ( A(1, dim*(i-1)+1 : dim*i) * result + b(1, i)) * (A(1, dim*(i-1)+1 : dim*i) * result + b(1, i)) +...
                 ( A(2, dim*(i-1)+1 : dim*i) * result + b(2, i)) * (A(2, dim*(i-1)+1 : dim*i) * result + b(2, i)) -...               
                 epsilon * epsilon * ( c(dim*(i-1)+1 : dim*i) * result + d(i) ) * ( c(dim*(i-1)+1 : dim*i) * result + d(i) );
    end    
    REConsensusSet = sum(r <= 0);       
    outindex = r > 0;
    cd('../../../..')

end