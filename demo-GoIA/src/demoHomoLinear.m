function [RESULT] = demoHomoLinear(Data, epsilonPixel, timeLimit, selectedMethod, methodList)
% demo program for homography estimation
% Our simple  approach converts H
% estimation from an 8-dimensional model fitting problem to two
% linear model fitting problems.
% Data: prepocessed data 
% epsilon: inlier threshold
% ICL-NUIM Database
%% prepare data
    disp(['*********** Decomposed homography experiments ***********'])
    disp(['------------------------------------------------------------------------------------'])
    temp = 0;
    dim1 = 5; dim2 = 3;
    Repeats = size(Data, 2);
    
    for i = 1: Repeats   

        data_originalPt = Data(i).originaldata;
        data_Pt = Data(i).data;
        data_T1 = Data(i).T1;
        data_T2 = Data(i).T2;
        data_imLrgb = Data(i).imLrgb;
        data_imRrgb = Data(i).imRrgb;
            
        %% excute each method in each scene
        currentMethod = methodList{selectedMethod};    
        
        %H-->5D: eq.(23)
        %solInit
        solInit1 = rand(dim1 + 1, 1);            
        %gap  
        gap = 0; %the relative gap between Upper and Lower bound           
                
        data1 = HomoTo5D(data_Pt);
        [runtime1, solution1, CS1, flagofGlobal1, iter1, outindex1] = SloveLinearProm(currentMethod, data1, dim1, solInit1, epsilonPixel*data_T1(1), timeLimit, gap);           

        %H-->3D: eq.(24)        
        solInit2 = rand(dim2 + 1, 1);                
        timeLimit0 = timeLimit - runtime1;

        h11 = solution1(1); h12 = solution1(2); h13 = solution1(3); 
        h21 = solution1(4); h22 = solution1(5); h23 = solution1(6);
        if (size(data_Pt, 1) == 4)
            N = size(data_Pt, 2);
            data_Pt = [data_Pt(1:2,:); ones(1, N); data_Pt(3:4,:); ones(1, N)];
        end
        tempdata = data_Pt(:, ~outindex1);
        data2 = ones(4, CS1);
        for ind = 1 : CS1
            if ( abs(tempdata(4, ind)) >= abs(tempdata(5, ind)) )
                data2(:, ind) = [tempdata(1, ind); tempdata(2, ind); 1; (h11*tempdata(1,ind) + h12*tempdata(2,ind) + h13)./tempdata(4, ind)];
            elseif ( abs(tempdata(4, ind)) < abs(tempdata(5, ind)) )
                data2(:, ind) = [tempdata(1, ind); tempdata(2, ind); 1; (h21*tempdata(1,ind) + h22*tempdata(2,ind) + h23)./tempdata(5, ind)];
            end           
        end    
                    
        if (timeLimit0 > 0)
            [runtime2, solution2, CS2, flagofGlobal2, iter2, outindex2] = SloveLinearProm(currentMethod, data2, dim2, solInit2, epsilonPixel*data_T1(1), timeLimit0, gap);           
        else
            runtime2 = 0; 
            solution2 = inf(dim2 + 1, 1);
            CS2 = CS1;
            flagofGlobal2 = false;
            iter2 = 0;
            outindex2 = false(1, CS1);
        end
        h31 = -solution2(1)/solution2(4); h32 = -solution2(2)/solution2(4); h33 = -solution2(3)/solution2(4);
        H_GoIA = [h11, h12, h13;
                            h21, h22, h23;
                            h31, h32, h33];
        H_GoIA = inv(data_T2)*H_GoIA*data_T1;
        H_GoIA = H_GoIA./H_GoIA(end, end);

        %H-->DLT              
        data_CS = data_Pt(:, ~outindex1);
        data_CS = data_CS(:, ~outindex2);
        if (CS2 >= 4 )
            epsilon_DLT = 4 * data_T1(1);
            H_DLT = homography2d( data_CS(1:3, :), data_CS(4:6, :) );
            H_DLT = H_DLT./H_DLT(end, end);
            tempr = [(H_DLT(1, :) * data_Pt(1:3, :)) ./ (H_DLT(3, :) * data_Pt(1:3, :));
                            (H_DLT(2, :) * data_Pt(1:3, :)) ./ (H_DLT(3, :) * data_Pt(1:3, :))] - data_Pt(4:5, :);
            residual = sqrt(sum(tempr.*tempr));
            CS_DLT = sum(residual < epsilon_DLT);
            H_DLT = inv(data_T2)*H_DLT*data_T1;
            H_DLT = H_DLT./H_DLT(end, end); 
        else
            H_DLT = inf(3,3);
            CS_DLT = 0;
        end
                 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%  RESULT =  [method, scene, viewIDL, viewIDR, epsilon, runtime, flagofGlobal, CS, H_GoIA, iter, epsilon_DLT, CS_DLT, H_DLT]  %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%

        temp = temp + 1;
        RESULT(temp).method = currentMethod;
        RESULT(temp).scene = Data(i).scene;
        RESULT(temp).viewIDL = Data(i).Lindex;
        RESULT(temp).viewIDR = Data(i).Rindex;
        RESULT(temp).epsilon = epsilonPixel;
        RESULT(temp).runtime = runtime1 + runtime2;
        RESULT(temp).flagofGlobal = flagofGlobal1 && flagofGlobal2;
        RESULT(temp).CS = CS2;
        RESULT(temp).H_GoIA = H_GoIA;  
        RESULT(temp).iter = iter1 + iter2;
        RESULT(temp).epsilon_DLT = epsilon_DLT;
        RESULT(temp).CS_DLT = CS_DLT;
        RESULT(temp).H_DLT = H_DLT;  
                
        disp([currentMethod, ' ', ' Scene-', num2str(Data(i).scene),  ': Frame ', num2str(Data(i).Lindex), ' and ', num2str(Data(i).Rindex), ' finished Homography estimation!'])
                    
    end
            
   

end
