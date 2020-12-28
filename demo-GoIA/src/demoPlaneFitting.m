function [RESULT] = demoPlaneFitting(Data, epsilon, NumofPlane, timeLimit, selectedMethod, methodList)
% demo program for plane fitting
% In one scene, we sequencially estimated two planes.
% data: prepocessed data 
% epsilon: inlier threshold
% NumofPlane: number of planes 

%% prepare data
    disp(['****************Plane fitting experiments *******************'])
    disp(['*********Estimate ', num2str(NumofPlane), ' planes in each scene sequencially.*********']);
    disp(['----------------------------------------------------------------------------------------'])
    temp = 0;
    Repeats = size(Data, 2);
    
    for i = 1: Repeats   
        data = Data{i};
        dataPt = data.originalPtCloud;
        dataImg = data.originalImage;
        dataPtnormalized = data.normalizedPtCloud;
        T = data.T;
        dataPtdownsampled = data.downsamplePtCloud;
            
%% preprocess the orginalPointCloud with your own chose of gridStep
%             removeNaN = isnan(dataPt(1, :));
%             dataPt = [dataPt; ones(1, size(dataPt, 2))];  
%             dataPt(end, removeNaN) = 0;
%             [dataPtnormalized, T] = normalise3dpts(dataPt); %data.normalizedPtCloud, data.T
%             dataPtnormalized(:, removeNaN)=[];
%             ptCloud = pointCloud(dataPtnormalized(1:3, :)');
%             %gridStep
%             gridStep = 0.05;
%             ptCloudD = pcdownsample(ptCloud, 'gridAverage', gridStep);
%             dataPtdownsampled = ptCloudD.Location'; %data.downsamplePtCloud
            
      %% excute each method in each scene
        for idxMethod = selectedMethod
            currentMethod = methodList{idxMethod};    
            %solInit
            if strcmp(currentMethod, 'GoIA-UN')
                solInit = rand(3+1, 1);                
            elseif strcmp(currentMethod, 'GoIA-LR') 
                solInit = rand(3, 1);
            end       

            %gap  
            gap = 0; %the relative gap between Upper and Lower bound               

            datainput = dataPtdownsampled;
            flaginput = zeros(1, size(dataPtdownsampled,2));
            flagdraw = zeros(1, size(dataPtnormalized,2));
            for numP = 1:NumofPlane
                [runtime, solution, CS, flagofGlobal, iter, outindex] = SloveLinearProm(currentMethod, datainput, 3, solInit, epsilon, timeLimit, gap);           
                if strcmp(currentMethod, 'GoIA-LR')
                    flaginput (...
                    abs(dataPtdownsampled(1,:)*solution(1)+dataPtdownsampled(2,:)*solution(2)+solution(3)-dataPtdownsampled(3,:)) < epsilon...
                    ) = numP;

                    flagdraw (...
                    abs(dataPtnormalized(1,:)*solution(1)+dataPtnormalized(2,:)*solution(2)+solution(3)-dataPtnormalized(3,:)) < epsilon...
                    ) = numP;

                elseif strcmp(currentMethod, 'GoIA-UN') 
                    flaginput (...
                    abs(dataPtdownsampled(1,:)*solution(1)+dataPtdownsampled(2,:)*solution(2)+dataPtdownsampled(3,:)*solution(3)+solution(4)) < epsilon...
                    ) = numP;       

                    flagdraw (...
                    abs(dataPtnormalized(1,:)*solution(1)+dataPtnormalized(2,:)*solution(2)+dataPtnormalized(3,:)*solution(3)+solution(4)) < epsilon...
                    ) = numP;                        
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%  RESULT =  [scene, method, numP, epsilon, runtime, flagofGlobal, CS, iter, solution]  %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                temp = temp + 1;
                RESULT(temp).scene = i;
                RESULT(temp).currentMethod = currentMethod;
                RESULT(temp).numP = numP;
                RESULT(temp).epsilon = epsilon;
                RESULT(temp).runtime = runtime;
                RESULT(temp).flagofGlobal = flagofGlobal;
                RESULT(temp).CS = CS;
                RESULT(temp).iter = iter;
                RESULT(temp).solution = solution;  

                dataPtdownsampled(:, flaginput==numP) = Inf;
                dataPtnormalized(:, flagdraw==numP) = Inf;
                datainput = dataPtdownsampled(:, flaginput==0) ;       

            end
            disp([currentMethod, ' ', num2str(i), '-scene finished!'])

            I = dataImg;
            [H, W, D] = size(I);

            M = reshape(flagdraw, H, W);
            color = {'yellow','red','green'};
            II = imoverlay(I, M, 3, color) ;
            figure 
            imshow(II);

        end                 
    end

end
