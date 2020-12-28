function [RESULT] = demoTranslationEstimation(Data, epsilon, timeLimit, selectedMethod, methodList)
% demo program for translation estimation
% KITTI Odometry dataset
% Data: prepocessed data 
% epsilon: inlier threshold

%% prepare data
    disp(['****************Translation estimation experiments *******************'])
    disp(['----------------------------------------------------------------------------------------'])
    temp = 0;
    Repeats = size(Data, 2);
    
    for i = 1: Repeats   
            % In one image pair, 
            % originaldata: the location of 2D matching points
            % data: y*Ry' in eq.(20) 
            % t_GT(1)*data(1,:) + t_GT(2)*data(2,:) + t_GT(3)*data(3,:) =
            % 0;
            data = Data{i};
            
          %% excute each method in each scene
            for idxMethod = selectedMethod
                currentMethod = methodList{idxMethod};    
                %solInit
                if strcmp(currentMethod, 'GoIA-UN')
                    solInit = rand(2+1, 1);                
                elseif strcmp(currentMethod, 'GoIA-LR') 
                    solInit = rand(2, 1);
                end       

                %gap  
                gap = 0; %the relative gap between Upper and Lower bound               
                
                [runtime, solution, CS, flagofGlobal, iter, outindex] = SloveLinearProm(currentMethod, data.data, 2, solInit, epsilon, timeLimit, gap);           
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%  RESULT =  [method, viewIDL, viewIDR, epsilon, runtime, flagofGlobal, CS, iter, solution, error_angle]  %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%
                if strcmp(currentMethod, 'GoIA-LR')
                    solution = [solution; -1];
                    solution = solution/norm(solution);         
                end                
                error_angle  = acosd(data.t_GT*solution);
                if error_angle>90
                    error_angle = 180-error_angle;
                end
                
                temp = temp + 1;
                RESULT(temp).method = currentMethod;
                RESULT(temp).viewIDL = data.viewIDL;
                RESULT(temp).viewIDR = data.viewIDR;
                RESULT(temp).epsilon = epsilon;
                RESULT(temp).runtime = runtime;
                RESULT(temp).flagofGlobal = flagofGlobal;
                RESULT(temp).CS = CS;
                RESULT(temp).iter = iter;
                RESULT(temp).solution = solution;  
                RESULT(temp).error_angle = error_angle;  
                
                figure
                imagesc(cat(2, data.imLrgb, data.imRrgb));
                hold on
                plot(data.originaldata(1,:), data.originaldata(2,:), 'r+', 'MarkerSize', 10);
                hold on
                plot(data.originaldata(4,:)+size(data.imLrgb, 2), data.originaldata(5,:), 'b+','MarkerSize', 10);

                flag_in = 1:size(data.originaldata, 2);
                flag_out = outindex;
                flag_in(flag_out)=[];
                xa = data.originaldata(1,:);
                ya = data.originaldata(2,:);
                xb = data.originaldata(4,:)+size(data.imLrgb,2);
                yb = data.originaldata(5,:);

                h = line([xa(flag_in); xb(flag_in)], [ya(flag_in); yb(flag_in)]);
                set(h, 'linewidth', 1, 'color','g');
                ho = line([xa(flag_out); xb(flag_out)], [ya(flag_out); yb(flag_out)]);
                set(ho,'linewidth',1,'color','r');      
                    
            end
            
            disp([currentMethod, ' ', num2str(i), '-scene finished!'])
   
    end

end
