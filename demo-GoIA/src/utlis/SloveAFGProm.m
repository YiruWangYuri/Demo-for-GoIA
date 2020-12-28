function [runtime,sol,CS,outl, flagoflocal] =  SloveAFGProm(method,XYZ,solInit,th,timeLimit,gapforOur,K)
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
if strcmp(method, 'OurLR-AFR')
    [sol,CS, outl] = OurLRAFR(XYZ,solInit,th,gapforOur, timeLimit); %%%%%%%%%%3D 4D LR
elseif strcmp(method, 'A*-NAPA') 
        [A,b,c,d] = genMatrixHomography([XYZ(1,:);XYZ(2,:)], [XYZ(3,:);XYZ(4,:)]);      
        [sol, outl, UniqueNodeNumber, hInit, ubInit, levelMax, NODIBP, runtime] = pseudoConvFit(A,b,c,d,solInit,th,method,timeLimit);
        RESULT_AT = [sol(1) sol(2) sol(3);
                                    sol(4) sol(5) sol(6)];
        temp1 = RESULT_AT(:,1:end-1)*[XYZ(1,:);XYZ(2,:)]+RESULT_AT(:,end)-[XYZ(3,:);XYZ(4,:)];
        temp2 = sqrt(sum(temp1.*temp1));
        CS = sum(temp2 < th+eps);
        outl = find(temp2 >= th + eps);
elseif strcmp(method, 'A*-TOD')
        [A,b,c,d] = genMatrixHomography([XYZ(1,:);XYZ(2,:)], [XYZ(3,:);XYZ(4,:)]);      
        [sol, outl, UniqueNodeNumber, hInit, ubInit, levelMax, NODIBP, runtime] = pseudoConvFit(A,b,c,d,solInit,th,method,timeLimit);
        RESULT_AT = [sol(1) sol(2) sol(3);
                                    sol(4) sol(5) sol(6)];
        temp1 = RESULT_AT(:,1:end-1)*[XYZ(1,:);XYZ(2,:)]+RESULT_AT(:,end)-[XYZ(3,:);XYZ(4,:)];
        temp2 = sqrt(sum(temp1.*temp1));
        CS = sum(temp2 < th+eps);
        outl = find(temp2 >= th + eps);
elseif strcmp(method, 'A*-NAPA-TOD')
        [A,b,c,d] = genMatrixHomography([XYZ(1,:);XYZ(2,:)], [XYZ(3,:);XYZ(4,:)]);      
        [sol, outl, UniqueNodeNumber, hInit, ubInit, levelMax, NODIBP, runtime] = pseudoConvFit(A,b,c,d,solInit,th,method,timeLimit);
        RESULT_AT = [sol(1) sol(2) sol(3);
                                    sol(4) sol(5) sol(6)];
        temp1 = RESULT_AT(:,1:end-1)*[XYZ(1,:);XYZ(2,:)]+RESULT_AT(:,end)-[XYZ(3,:);XYZ(4,:)];
        temp2 = sqrt(sum(temp1.*temp1));
        CS = sum(temp2 < th+eps);
        outl = find(temp2 >= th + eps);
elseif strcmp(method, 'A*-NAPA-DIBP') 
        [A,b,c,d] = genMatrixHomography([XYZ(1,:);XYZ(2,:)], [XYZ(3,:);XYZ(4,:)]);      
        [sol, outl, UniqueNodeNumber, hInit, ubInit, levelMax, NODIBP, runtime] = pseudoConvFit(A,b,c,d,solInit,th,method,timeLimit);
        RESULT_AT = [sol(1) sol(2) sol(3);
                                    sol(4) sol(5) sol(6)];
        temp1 = RESULT_AT(:,1:end-1)*[XYZ(1,:);XYZ(2,:)]+RESULT_AT(:,end)-[XYZ(3,:);XYZ(4,:)];
        temp2 = sqrt(sum(temp1.*temp1));
        CS = sum(temp2 < th+eps);
        outl = find(temp2 >= th + eps);
elseif strcmp(method, 'FLRS') 
%     linearConfig = config;
    linearConfig.QThresh = 1e-12; % numerical accuracy
            
    [sol, CS, outl] = LocalAF(XYZ, th, method, solInit, linearConfig);

elseif strcmp(method, 'IBCOFLRS') || strcmp(method, 'IBCO') 
    linearConfig.QThresh = 1e-12; % numerical accuracy
    [sol, CS, outl] = LocalAF(XYZ, th, method, solInit, linearConfig);

end
runtime = toc;
if runtime>timeLimit
    flagoflocal=1111;
else
    flagoflocal=0000;
end





end