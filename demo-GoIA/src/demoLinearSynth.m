function [RESULT] = demoLinearSynth(N, OutlierRatio, Dim, Repeats, noise_in, noise_out, timeLimit, selectedMethod, methodList)
% demo program for linear synthetic data
% N: number of data points
% Outlier: outlier ratio 
% Dim: dimensionality of linear model (\hat(\theat) in eq.(4))
% Repeats: times of randomly generated data 

    rng(1); 
    seed = randi(10000,1,1000);
%% prepare data
    disp(['****************Linear synthetic experiments *******************'])
    disp(['                Number     dim         Noise       OutlierRatio']);
    disp(['                     ', num2str(N), '           ', num2str(Dim), '            ', num2str(noise_in), '                ', num2str(OutlierRatio*100), '%']);
    disp(['----------------------------------------------------------------------------------------'])
    temp = 0;
    
    for i = 1: Repeats
            [data, para] = PrepareLinearData(N, OutlierRatio, Dim, seed(i), noise_in, noise_out);
   
          %% excute each method in each data
            for idxMethod = selectedMethod
                
                currentMethod = methodList{idxMethod};    
                %solInit
                if strcmp(currentMethod, 'GoIA-UN')
                    solInit = rand(Dim+1, 1);                
                elseif strcmp(currentMethod, 'GoIA-LR') 
                    solInit = rand(Dim, 1);
                end       
                %epsilon
                if strcmp(currentMethod, 'GoIA-UN') 
                    epsilon = 1.5*noise_in;
                else
                    epsilon = 3*noise_in;
                end
                %gap  
                gap = 0; %the relative gap between Upper and Lower bound               
                
                [runtime, solution, CS, flagofGlobal, iter, outindex] = SloveLinearProm(currentMethod, data, Dim, solInit, epsilon, timeLimit, gap);           
                disp([currentMethod, ' ', num2str(i), '-th Repeats finished!'])
                
%                 if strcmp(currentMethod, 'GoIA-UN') 
%                     solution1 = solution/-solution(end-1);
%                     solution2 = [solution1(1:end-2); solution1(end)];
%                     solution = solution2;
%                 end     
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%  RESULT =  [i, method, N, dim, noise, outlierRatio, runtime, flagofGlobal, iter, CS, solution]  %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                temp = temp + 1;
                RESULT(temp).repeats = i;
                RESULT(temp).currentMethod = currentMethod;
                RESULT(temp).N = N;
                RESULT(temp).dim = Dim;
                RESULT(temp).noise = noise_in;
                RESULT(temp).outlier = OutlierRatio;
                RESULT(temp).runtime = runtime;
                RESULT(temp).flagofGlobal = flagofGlobal;
                RESULT(temp).iter = iter;
                RESULT(temp).CS = CS;
                RESULT(temp).solution = solution;   
                
            end
            
    end
    disp('*****************************************************************');
end
