function [RESULT] = demoTriangulation(Data, epsilon, timeLimit, Frame, selectedMethod, methodList)
% demo program for triangulation 
% Data: the data is the scene of xxx from http://www.maths.lth.se/matematiklth/personal/calle/dataset/dataset.html
% epsilon: inlier threshold

%% prepare data
    disp(['**************** Triangulation experiments *******************'])
    disp(['****************           Prepare data!           *******************'])
    disp(['----------------------------------------------------------------------------------------'])
    
    dim = 3;
    temp = 0;
    P = Data.P; % a series of camera
    NumofCamera = size(P, 2); % number of cameras
    NumofPoints_All = Data.u_uncalib.pointnr; % number of 3D scene points
	ImgArray_r=inf(NumofCamera*2, NumofPoints_All);
    for camera_index=1 : NumofCamera
        ImgArray_r([camera_index*2-1;camera_index*2], Data.u_uncalib.index{camera_index,1})= Data.u_uncalib.points{camera_index, 1}(1:2,:);
    end   
    %remove points that are not visible in two cameras or more
    for i = 1:size(ImgArray_r, 1)
        for j = 1:size(ImgArray_r, 2)
            if ~isfinite(ImgArray_r(i, j))
               ImgArray_r(i, j)=NaN;
            end
        end
    end    
    
    vis = zeros(1, NumofPoints_All); % Each 3D scene point was seen in N views, and vis recorded the number of views (N). 
    for i=1 : NumofCamera
        vis = vis + isfinite(ImgArray_r(i*2, :));
    end

    vis = vis >= Frame;
    ImgArray_o =  ImgArray_r(:, vis);
    NumofPointsAfterProcess = size(ImgArray_o, 2);                   

    %% excute each method in each scene
    
    for idxMethod = selectedMethod
        currentMethod = methodList{idxMethod};
        ScenePt = [];
        disp('Sequencially construct Index3DScenePt : ');
        for Index3DScenePt = 1 : NumofPointsAfterProcess
            fid1 = find( isfinite(ImgArray_o(:, Index3DScenePt)) );  
            fid2 = fid1(2:2:end)/2;
            pt_2d=ImgArray_o(fid1, Index3DScenePt);
            pt_2d=reshape(pt_2d, 2, []);
            data_Pt = [pt_2d; ones(1,size(pt_2d,2))];
            data_P=P(fid2);
            data = Triangulation2Quasi(data_Pt, data_P);
            
             
            %solInit
            solInit = rand(dim, 1);           
            %gap  
            gap = 0; %the relative gap between Upper and Lower bound                   
  
            [runtime, solution, CS, flagofGlobal, iter, outindex] = SloveQuasiProm(currentMethod, data, dim, solInit, epsilon, timeLimit, gap);
            disp([currentMethod,'-','Index3DScenePt: ', num2str(Index3DScenePt),'-All: ', num2str(NumofPointsAfterProcess)])
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%  RESULT =  [method, epsilon, runtime, flagofGlobal, CS, iter, solution]  %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%                

            temp = temp + 1;
            RESULT(temp).method = currentMethod;
            RESULT(temp).epsilon = epsilon;
            RESULT(temp).runtime = runtime;
            RESULT(temp).flagofGlobal = flagofGlobal;
            RESULT(temp).CS = CS;
            RESULT(temp).iter = iter;
            RESULT(temp).solution = solution;  
            ScenePt = [ScenePt, solution];
        end
        
        disp([' Finished! ']);          
        figure
        plot3(ScenePt(1, :), ScenePt(2, :), ScenePt(3, :), 'Marker', '.', 'LineStyle', 'none', 'Color', [1 0 0]);
       
    end

end
                           
%                              temp =  [temp; sum(RESULT(:, 5)), size(RESULT, 1), sum(RESULT(:, 6)), sum(RESULT(:, 4)), 1-sum(RESULT(:, 3))/1111/ size(RESULT, 1),RESULT(1,2), RESULT(1,end-1)];
% 
%                             
%                             if strcmp(currentMethod, 'OurLR-TRian') 
%                                 disp([currentMethod,' is drawing pictures!'])
%                                 Fig_triangulation(RESULT, problemList{i_data});
%                                 saveas(gcf, [problemList{i_data},'hist'], 'fig');
%                                 
%                                 figure
%                                 plot3(RESULT(:, 8),RESULT(:, 9),RESULT(:, 10),'Marker','.','LineStyle','none','Color',[1 0 0])
%                                 hold on
%                                 plot3(U(1,:),U(2,:),U(3,:),'MarkerSize',0.5,'Marker','o','LineStyle','none',...
%     'Color',[0.0784313753247261 0.168627455830574 0.549019634723663])
%                                 saveas(gcf, [problemList{i_data}], 'fig');
% 
%                             end
%                             
%                              RESULT = [];                    
%                             disp([currentMethod,'-',problemList{i_data},'---Finished!'])
%  
% 
% end
