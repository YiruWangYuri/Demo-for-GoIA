function [result, iter, flagofGlobal] = Search(A, b, c, d, epsilon, gap, solInit, timeLimit)
    
    N = length(d);
    dimension = length(c)/N;
    iter = 1;
    bits_per_dim = 2;
    num_thread = bits_per_dim^dimension;   
    branchQueue = [];
    branchnum(iter) = 1;
    bestBranch= [-100 * ones(dimension ,1); 200 * ones(dimension, 1); 0; 1];
    best_U_branch_center = solInit;
    
    r = zeros(N, 1);
    for i = 1 : N
        r(i) = ( A(1, dimension*(i-1)+1 : dimension*i) * best_U_branch_center + b(1, i)) * (A(1, dimension*(i-1)+1 : dimension*i) * best_U_branch_center + b(1, i)) +...
                 ( A(2, dimension*(i-1)+1 : dimension*i) * best_U_branch_center + b(2, i)) * (A(2, dimension*(i-1)+1 : dimension*i) * best_U_branch_center + b(2, i)) -...               
                 epsilon * epsilon * ( c(dimension*(i-1)+1 : dimension*i) * best_U_branch_center + d(i) ) * ( c(dimension*(i-1)+1 : dimension*i) * best_U_branch_center + d(i) );
    end
    best_U = sum(r>0)/N;


    % dynamic updates of lower and upper bounds 
    % will significantly increase runtime
    
%     figure
%     draw_l = animatedline;
%     hold on
%     draw_u = animatedline;  

    while(1)
    
        ttt = toc;
        if ttt>=timeLimit
            flagofGlobal = false;
            break
        end
        
        [L,U] = getLUperBranch(A, b, c, d, epsilon, bestBranch, bits_per_dim, dimension, num_thread);
        [temp_best_U, temp_best_U_index]= min(U);
        if (temp_best_U <= best_U)
            best_U = temp_best_U;
            best_U_index = temp_best_U_index;
        else
            best_U_index = -9;
        end    

        if ( best_U_index>=0 )
            [best_U_branch_]=getChild(bestBranch, best_U_index, bits_per_dim, dimension);
            best_U_branch_center = [ best_U_branch_(1:dimension) + best_U_branch_(dimension+1:end)*0.5;          
                                                     L(best_U_index);
                                                     U(best_U_index)];
        end
        [branch_, num_] = addSubbranchTobranch(L, U, best_U, bestBranch, bits_per_dim, dimension, num_thread);
        branchQueue = [branchQueue, branch_(:, 1:num_)];
        iter = iter + 1;
        branchnum(iter) = branchnum(iter - 1) + num_;
        
        [best_L, index] = min(branchQueue(end-1, :));
        
        bestBranch = branchQueue(:, index);
        branchQueue(end-1, index)=1e10;
        branchQueue(end, index)=-1e10;
    
        LL(iter-1)=best_L;
        UU(iter-1)=best_U;   
 
%       addpoints(draw_l, iter, LL(iter-1));
%       addpoints(draw_u, iter, UU(iter-1));
%       drawnow
      
        if(best_U-best_L<=gap)
            flagofGlobal = true;
            break
        end       
    end
    
    result = best_U_branch_center(1:dimension);

    iter = size(LL, 2);
    L = [0, LL];
    U = [1, UU];
end