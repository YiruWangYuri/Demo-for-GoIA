function data = AffineReg2Quasi(data_Pt)
% r =     ||A\theta + b ||_2
%       ----------------------------------
%                  c\theta + d
    if (size(data_Pt, 1) ~= 6 && size(data_Pt, 1) ~= 4)
        error('The input data shoule be 6*N !');
    end
    
    N = size(data_Pt, 2); 
    
    if (size(data_Pt, 1) == 4)
        data_Pt = [data_Pt(1:2, :); ones(1, N); data_Pt(3:4, :); ones(1, N)];
    end    
    
    for i = 1 : N
        data(i).A = [ data_Pt(1:3, i)', 0, 0, 0;
                            0, 0, 0, data_Pt(1:3, i)' ];
        data(i).b = [ -data_Pt(4, i);
                           -data_Pt(5, i) ];
        data(i).c = zeros(1, 6);
        data(i).d = 1;        
    end 
 
end
