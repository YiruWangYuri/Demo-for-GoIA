function [runtime,sol,CS,outl, flagoflocal,iter] =  SloveHomo(method,XYZ,solInit,th,timeLimit,gapforOur,K)
if nargin < 7
    K = 2;
end
if nargin < 6
    gapforOur = -1;
end
if nargin < 5
    timeLimit = -1;
end

tic;
if strcmp(method, 'OurLR-Homo')
    [sol,CS,iter] = OurLR(XYZ,solInit,th,gapforOur, timeLimit); %%%%%%%%%%3D 4D LR
    inliers = (abs(XYZ(end,:)-(sol'*[XYZ(1:end-1,:);ones(1,size(XYZ,2))]))<th+eps);
    outl = find(inliers == 0);
    CS= sum(inliers);
   % [sol,CS, outl] = OurLRHomo(XYZ,solInit,th,gapforOur, timeLimit); %%%%%%%%%%3D 4D LR
elseif strcmp(method, 'OurUN-Homo')

    [sol,CS,iter] = OurUN_Homo(XYZ,solInit,th,gapforOur, timeLimit); %%%%%%%%%%3D 4D LR
    inliers = (abs(sol'*XYZ)<th+eps); 
    outl = find(inliers == 0);
elseif strcmp(method, 'A*-NAPA') 
        [A,b,c,d] = genMatrixHomographyTT([XYZ(1:2, :); ones(1, size(XYZ, 2))], [XYZ(3:4, :); ones(1, size(XYZ, 2))]);      
        [sol, outl, UniqueNodeNumber, hInit, ubInit, levelMax, NODIBP, runtime] = pseudoConvFit(A,b,c,d,solInit,th,method,timeLimit);
        CS = size(XYZ, 2) - size(outl, 2);
        inliers = true(1, size(XYZ, 2));
        inliers(outl) = 0;
        outl = ~inliers; iter=0;
elseif strcmp(method, 'A*-TOD')
        [A,b,c,d] = genMatrixHomographyTT([XYZ(1:2, :); ones(1, size(XYZ, 2))], [XYZ(3:4, :); ones(1, size(XYZ, 2))]);      
        [sol, outl, UniqueNodeNumber, hInit, ubInit, levelMax, NODIBP, runtime] = pseudoConvFit(A,b,c,d,solInit,th,method,timeLimit);
        CS = size(XYZ, 2) - size(outl, 2);
        inliers = true(1, size(XYZ, 2));
        inliers(outl) = 0;
        outl = ~inliers;iter=0;
elseif strcmp(method, 'A*-NAPA-TOD')
        [A,b,c,d] = genMatrixHomographyTT([XYZ(1:2, :); ones(1, size(XYZ, 2))], [XYZ(3:4, :); ones(1, size(XYZ, 2))]);      
        [sol, outl, UniqueNodeNumber, hInit, ubInit, levelMax, NODIBP, runtime] = pseudoConvFit(A,b,c,d,solInit,th,method,timeLimit);
        CS = size(XYZ, 2) - size(outl, 2);
        inliers = true(1, size(XYZ, 2));
        inliers(outl) = 0;
        outl = ~inliers;iter=0;
elseif strcmp(method, 'A*-NAPA-DIBP') 
        [A,b,c,d] = genMatrixHomographyTT([XYZ(1:2, :); ones(1, size(XYZ, 2))], [XYZ(3:4, :); ones(1, size(XYZ, 2))]);      
        [sol, outl, UniqueNodeNumber, hInit, ubInit, levelMax, NODIBP, runtime] = pseudoConvFit(A,b,c,d,solInit,th,method,timeLimit);
        CS = size(XYZ, 2) - size(outl, 2);
        inliers = true(1, size(XYZ, 2));
        inliers(outl) = 0;
        outl = ~inliers;iter=0;
elseif strcmp(method, 'Li')
    [sol,CS,iter] = LiUN_H(XYZ,solInit,th,gapforOur,timeLimit);  %%%%%3D LR 最大一致集
    inliers = (abs(XYZ(end,:)-(sol'*[XYZ(1:end-1,:)]))<th+eps);
    outl = find(inliers == 0);
    CS= sum(inliers);    
elseif strcmp(method, 'Zheng')
    [sol,CS,iter] = maxconMIDCP(XYZ,solInit,th,gapforOur,timeLimit,K); %%%%3+1D UN
    inliers = (abs(sol'*XYZ)<th+eps); 
    outl = find(inliers == 0);
elseif strcmp(method, 'FLRS') 
%     linearConfig = config;
    linearConfig.QThresh = 1e-12; % numerical accuracy
    [sol, rinls, ransacIters, T1, T2, sampleSet] = ransacfithomography(XYZ(1:2, :), XYZ(3:4, :), th,1, 'F-LORANSAC', [], 0);
    sol = sol./sol(end);sol(end) = [];
    CS = size(rinls, 1);
    outl = true(1, size(XYZ, 2));                     
    outl(rinls) = 0;iter=0;
elseif strcmp(method, 'IBCOFLRS') || strcmp(method, 'IBCO') 
    linearConfig.QThresh = 1e-12; % numerical accuracy
    [A,b,c,d] = genMatrixHomographyTT([XYZ(1:2, :); ones(1, size(XYZ, 2))], [XYZ(3:4, :); ones(1, size(XYZ, 2))]);      
    [sol, runtime] = IBCO_v2(A, b, c, d, solInit, th, linearConfig);
    CS = [];
    outl = [];          iter=0;           
end
runtime = toc;
if runtime>timeLimit
    flagoflocal=1111;
else
    flagoflocal=0000;
end





end