function data = Triangulation2Quasi(data_Pt, data_P)
% r =     ||A\theta + b ||_2
 %       ----------------------------------
 %                  c\theta + d
    if (size(data_Pt, 1) ~= 2 && size(data_Pt, 1) ~= 3)
        error('The input data shoule be 2*N !');
    end
    
    if (size(data_P, 1) ~= 1)
        error('The input data shoule be 1*numberofCamera !');
    end    
    
    if (size(data_Pt, 2) ~= size(data_P, 2))
        error('The number of camera and 2d point should be the same !');
    end

    N = size(data_Pt, 2); 
    for i = 1 : N
        P = data_P{i};
        u = data_Pt(1 : 2, i);
        data(i).A = [ P(1, 1:3) - u(1) * P(3, 1:3);
                           P(2, 1:3) - u(2) * P(3, 1:3)];
        data(i).b = [ P(1, 4) - u(1) * P(3, 4);
                           P(2, 4) - u(2) * P(3, 4)];
        data(i).c = P(3, 1:3);
        data(i).d = P(3, 4);        
    end 
 
end
