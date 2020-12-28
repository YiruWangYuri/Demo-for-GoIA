function data = HomoTo5D( data_Pt )
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
    
    data = zeros(6, N);

    x1 = data_Pt(1, :);
    y1 = data_Pt(2, :);
    
    x2 = data_Pt(4, :);
    y2 = data_Pt(5, :);
        
    data = [(x1.*y2); (y1.*y2); (y2); (-x1.*x2); (-x2.*y1); (-x2)];
 
end
