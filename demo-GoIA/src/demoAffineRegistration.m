function [RESULT] = demoAffineRegistration(Data, epsilonPixel, timeLimit, selectedMethod, methodList)
% demo program for affine registration
% Data: prepocessed data 
% epsilon: inlier threshold
% ETH-Objects Database
% There are 53 groups objects and each of group consits of 5 images. 
% We provide five image pairs for testing.
%% prepare data
    disp(['****************Affine Registration experiments *******************'])
    disp(['----------------------------------------------------------------------------------------'])
    temp = 0;
    dim = 6;
    Repeats = size(Data, 2);
    
    for i = 1: Repeats   

            data_originalPt = Data(i).originaldata;
            data_Pt = Data(i).data;
            data_T1 = Data(i).T1;
            data_T2 = Data(i).T2;
            data_imLrgb = Data(i).imLrgb;
            data_imRrgb = Data(i).imRrgb;
            
          %% excute each method in each scene
            for idxMethod = selectedMethod
                currentMethod = methodList{idxMethod};    
                %solInit
                solInit = rand(dim, 1);                

                %gap  
                gap = 0; %the relative gap between Upper and Lower bound           
                
                data = AffineReg2Quasi(data_Pt);
                
                [runtime, solution, CS, flagofGlobal, iter, outindex] = SloveQuasiProm(currentMethod, data, dim, solInit, epsilonPixel*data_T1(1), timeLimit, gap);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%  RESULT =  [method, scene, Object, viewIDL, viewIDR, epsilon, runtime, flagofGlobal, CS, iter, solution]  %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%
                
                temp = temp + 1;
                RESULT(temp).method = currentMethod;
                RESULT(temp).Object = Data(i).Object;
                RESULT(temp).viewIDL = Data(i).ObjectLindex;
                RESULT(temp).viewIDR = Data(i).ObjectRindex;
                RESULT(temp).epsilon = epsilonPixel;
                RESULT(temp).runtime = runtime;
                RESULT(temp).flagofGlobal = flagofGlobal;
                RESULT(temp).CS = CS;
                RESULT(temp).iter = iter;
                RESULT(temp).solution = solution;  
                
                figure
                imagesc(cat(2, data_imLrgb, data_imRrgb));
                hold on
                plot(data_originalPt(1,:), data_originalPt(2,:), 'r+', 'MarkerSize', 10);
                hold on
                plot(data_originalPt(4,:)+size(data_imLrgb, 2), data_originalPt(5,:), 'b+','MarkerSize', 10);

                flag_in = 1:size(data_originalPt, 2);
                flag_out = outindex;
                flag_in(flag_out)=[];
                xa = data_originalPt(1,:);
                ya = data_originalPt(2,:);
                xb = data_originalPt(4,:)+size(data_imLrgb,2);
                yb = data_originalPt(5,:);

                h = line([xa(flag_in); xb(flag_in)], [ya(flag_in); yb(flag_in)]);
                set(h, 'linewidth', 1, 'color','g');
                ho = line([xa(flag_out); xb(flag_out)], [ya(flag_out); yb(flag_out)]);
                set(ho,'linewidth',1,'color','r');      
                    
            end
            
            disp([currentMethod, ' ', ' Object-', num2str(Data(i).Object),  ': Frame ', num2str(Data(i).ObjectLindex), ' and ', num2str(Data(i).ObjectRindex), ' finished Affine Registration!'])
   
    end

end
