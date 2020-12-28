function [X_Data, theta, guide_max_thre] = genSyntheticData(Num, noise_in, noise_out, outlierRatio, dim, FlagofBlance, seed)
% Num: number of point, include inliers and outlier
% noise_in: noise varience of inliers
% noise_out: noise varience of outliers
% outlierRatio: number of outliers = outlierRatio*Num
% dim: the dimension of \hat(\theta)  in eq.(4)
    if dim >= 8
        error('We recommend that dim should be less than 8');
    end
    
    rng(seed); 
    t = rand(dim,1); t = [t(1:end-1);-1;t(end)]; t = t/norm(t);
    vector_n = t(1:end-1);    vector_n = vector_n/norm(vector_n);
    temp = null(vector_n');
    
    if dim==3
        vector_nx = temp(:,1);
        vector_ny = temp(:,2);

        x_nxnynz = 2*rand(Num,1)-1;
        y_nxnynz = 2*rand(Num,1)-1;
        z_nxnynz = zeros(Num,1);

        x_vector = x_nxnynz*vector_nx';
        y_vector = y_nxnynz*vector_ny';
        z_vector = z_nxnynz*vector_n';

        X_Data = x_vector + y_vector + z_vector;
        X_Data = [X_Data, ones(Num,1)];
        
    elseif dim==2
        vector_nx = temp(:,1);
        
        x_nxnynz = 2*rand(Num,1)-1;
        z_nxnynz = zeros(Num,1);
 
        x_vector = x_nxnynz*vector_nx';
        z_vector = z_nxnynz*vector_n';

        X_Data = x_vector + z_vector;
        X_Data = [X_Data, ones(Num,1)];       

    elseif dim==4
        vector_nx = temp(:,1);
        vector_ny = temp(:,2);
        vector_np = temp(:,3);

        x_nxnynz = 2*rand(Num,1)-1;
        y_nxnynz = 2*rand(Num,1)-1;
        p_nxnynz = 2*rand(Num,1)-1;
        z_nxnynz = zeros(Num,1);

        x_vector = x_nxnynz*vector_nx';
        y_vector = y_nxnynz*vector_ny';
        p_vector = p_nxnynz*vector_np';
        z_vector = z_nxnynz*vector_n';

        X_Data = x_vector + y_vector + p_vector + z_vector;
        X_Data = [X_Data, ones(Num,1)];

     elseif dim==5
        vector_nx = temp(:,1);
        vector_ny = temp(:,2);
        vector_np = temp(:,3);
        vector_nq = temp(:,4);

        x_nxnynz = 2*rand(Num,1)-1;
        y_nxnynz = 2*rand(Num,1)-1;
        p_nxnynz = 2*rand(Num,1)-1;
        q_nxnynz = 2*rand(Num,1)-1;
        z_nxnynz = zeros(Num,1);

        x_vector = x_nxnynz*vector_nx';
        y_vector = y_nxnynz*vector_ny';
        p_vector = p_nxnynz*vector_np';
        q_vector = q_nxnynz*vector_nq';
        z_vector = z_nxnynz*vector_n';

        X_Data = x_vector + y_vector + p_vector + q_vector + z_vector;
        X_Data = [X_Data, ones(Num,1)];
        
    elseif dim==6
        vector_nx = temp(:,1);
        vector_ny = temp(:,2);
        vector_np = temp(:,3);
        vector_nq = temp(:,4);
        vector_nk = temp(:,5);

        x_nxnynz = 2*rand(Num,1)-1;
        y_nxnynz = 2*rand(Num,1)-1;
        p_nxnynz = 2*rand(Num,1)-1;
        q_nxnynz = 2*rand(Num,1)-1;
        k_nxnynz = 2*rand(Num,1)-1;
        z_nxnynz = zeros(Num,1);

        x_vector = x_nxnynz*vector_nx';
        y_vector = y_nxnynz*vector_ny';
        p_vector = p_nxnynz*vector_np';
        q_vector = q_nxnynz*vector_nq';
        k_vector = k_nxnynz*vector_nk';
        z_vector = z_nxnynz*vector_n';

        X_Data = x_vector + y_vector + p_vector + q_vector + k_vector + z_vector;
        X_Data = [X_Data, ones(Num,1)];       
        
    elseif dim==7        
        vector_nx = temp(:,1);
        vector_ny = temp(:,2);
        vector_np = temp(:,3);
        vector_nq = temp(:,4);
        vector_nk = temp(:,5);
        vector_nl = temp(:,6);

        x_nxnynz = 2*rand(Num,1)-1;
        y_nxnynz = 2*rand(Num,1)-1;
        p_nxnynz = 2*rand(Num,1)-1;
        q_nxnynz = 2*rand(Num,1)-1;
        k_nxnynz = 2*rand(Num,1)-1;
        l_nxnynz = 2*rand(Num,1)-1;
        z_nxnynz = zeros(Num,1);

        x_vector = x_nxnynz*vector_nx';
        y_vector = y_nxnynz*vector_ny';
        p_vector = p_nxnynz*vector_np';
        q_vector = q_nxnynz*vector_nq';
        k_vector = k_nxnynz*vector_nk';
        l_vector = l_nxnynz*vector_nl';
        z_vector = z_nxnynz*vector_n';

        X_Data = x_vector + y_vector + p_vector + q_vector + k_vector + l_vector + z_vector;
        X_Data = [X_Data, ones(Num,1)];    
    
    end
    
    % plane parameter theta 
    
    center = 0.1*rand(1,dim);      
    
    X_Data(:,1:end-1) = X_Data(:,1:end-1)+repmat(center,Num,1);
    theta = zeros(dim+1,1);
    theta(1:end-1) = vector_n;
    theta(end) = -center*vector_n;
    theta = theta/norm(theta);
    
   
    %corrupt inliers by noise_in gaussian noise
    outlier_num =round(Num*outlierRatio);
    inlier_num = Num-outlier_num;
    flag_out = 1:outlier_num;
    flag_in = outlier_num+1:Num;
    
    gNoise_in = -noise_in+2*noise_in*rand(inlier_num, dim);
    X_Data(flag_in,1:end-1) = X_Data(flag_in,1:end-1)+gNoise_in;
    
    gNoise_out = -noise_out+2*noise_out*rand(outlier_num, dim);
    
    if FlagofBlance==1
        gNoise_out = abs(gNoise_out);
    end
    X_Data(flag_out,1:end-1) = X_Data(flag_out,1:end-1)+gNoise_out;
  
    % guide inlier threshold for experiment 
    temp = abs(acos(X_Data*theta./sqrt(sum(X_Data.^2, 2)))-pi/2);
    guide_max_thre = max(temp(flag_in));
end