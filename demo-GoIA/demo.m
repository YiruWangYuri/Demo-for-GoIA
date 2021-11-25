%%   
%   Demo program for GoIA:   
%   Y. Wang, Y. Liu, X. Li, C. Wang, M. Wang, Z. Song
%   Practical Globally Optimal Consensus Maximization by Branch-and-Bound based on Interval Arithmetic.   
%   If you encounter any problems or questions please email to yiruwang18@fudan.edu.cn.   
%***********************************************************************
%   The program is free for non-commercial academic use. 
%   Any commercial use is strictly prohibited without the authors' consent. 
%   Please cite the above paper in any academic publications that have made use of this package or part of it.
%***********************************************************************
    clear;
    close all;
    addpath(genpath('data/')); 
    addpath(genpath('src/')); 
%% Problem List --- Method List 

%   Linear problem 
%   Synthetic-N-Outlier-TimeLimit-Dim --- GoIA-LR, GoIA-UN
%   PlaneFitting  --- GoIA-LR, GoIA-UN
%   TranslationEstimation --- GoIA-LR, GoIA-UN
%   HomoLinear --- GoIA-UN

%   Non-linear problem
%   Triangulation --- GoIA
%   AffineRegistration --- GoIA

    problemList = {'Synthetic-N-Outlier-TimeLimit-Dim', 'PlaneFitting', 'TranslationEstimation',...
                             'Triangulation', 'AffineRegistration',...
                             'HomoLinear'};
    methodList = {'GoIA-LR', 'GoIA-UN', 'GoIA'};
                         
    selectedProblem = 1; 
    selectedMethod = 1;

%   You can set the timeLimit to avoid extremely long runtime in each run. 
%   Please notice that the method would not guarantee the global optimality,
%   if it cannot converge within the timeLimit.
%   The output of RESULT.flagofGlobal indicates that if one run guarantees the
%   global optimality.
 
%% excute each problem 

    warning('off', 'all');
    currentProblem = problemList{selectedProblem};
    if strcmp(currentProblem, 'Synthetic-N-Outlier-TimeLimit-Dim')              
        % Synthetic-N-Outlier-TimeLimit-Dim
        
        timeLimit = 10;
        N = 200; % number of data points
        OutlierRatio = 0.2; % outlier ratio
        Dim = 3; % dimensionality of linear model
        Repeats = 50; % times of randomly generated data 
        noise_in = 0.01; % noise varience of inliers
        noise_out = 2; % noise varience of outliers
        RESULT = demoLinearSynth(N, OutlierRatio, Dim, Repeats, noise_in, noise_out, timeLimit, selectedMethod, methodList);
    
    elseif strcmp(currentProblem, 'PlaneFitting')
        % PlaneFitting
        
        timeLimit = 100;
        load('./data/PlaneFittingData.mat');
        NumofPlane = 2; % number of planes 
        disp('In the same task, we recommend that the inlier threshold of GoIA-UN (0.02) is half of that of GoIA-LR (0.04)!');
        epsilon = 0.02; % inlier threshold. We recommend that the inlier threshold of GoIA-UN (0.02) is half of that of GoIA-LR (0.04).
        RESULT = demoPlaneFitting(PlaneFittingData, epsilon, NumofPlane, timeLimit, selectedMethod, methodList);
      
    elseif strcmp(currentProblem, 'TranslationEstimation')
         % TranslationEstimation
         
        timeLimit = 100;
        load('./data/TranslationEstimationData.mat');
        disp('In the same task, we recommend that the inlier threshold of GoIA-UN (0.0025) is half of that of GoIA-LR (0.005)!');
        epsilon = 0.0025; % inlier threshold. We recommend that the inlier threshold of GoIA-UN (0.0025) is half of that of GoIA-LR (0.005).
        RESULT = demoTranslationEstimation(TranslationEstimationData, epsilon, timeLimit, selectedMethod, methodList);

    elseif strcmp(currentProblem, 'Triangulation')
         % Triangulation
         
        timeLimit = 5;
        TriangulationData = load('./data/Triangulation_Buddah_Statue.mat');
        epsilon = 2; 
        Frame = 30; % 3D scene points with more than 30 views were chosen for triangulation.
        RESULT = demoTriangulation(TriangulationData, epsilon, timeLimit, Frame, selectedMethod, methodList);
    
    elseif strcmp(currentProblem, 'AffineRegistration')
         % AffineRegistration
         
        timeLimit = 500; 
        load('./data/AffineRegistrationData.mat');
        epsilon = 4; 
        RESULT = demoAffineRegistration(AffineRegistrationData, epsilon, timeLimit, selectedMethod, methodList);        
    
    elseif strcmp(currentProblem, 'HomoLinear')
         % Our simple parameter decomposition approach converts H
         % estimation from an 8-dimensional model fitting problem to two
         % linear model fitting problems. GoIA-UN is used to estimate
         % decomposed H.
         
        timeLimit = 500; 
        load('./data/HomoLinearData.mat');
        epsilon = 2; 
        selectedMethod = 2;
        RESULT = demoHomoLinear(HomoLinearData, epsilon, timeLimit, selectedMethod, methodList);        
     end

     disp(['Please see the output of RESULT to check the detailed results.']);

 
 






