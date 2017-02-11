function corrected_dwmri = robust_qball_stats_v3(nifti_data,bval,bvec,mask_input)

    % Loading 4d nifti data, bvals and bvecs
    nii = load_untouch_nii(nifti_data);
    %nii = nifti_data;
    bvals = dlmread(bval);
    bvecs = dlmread(bvec);
    mask = load_untouch_nii(mask_input);
    mask = logical(mask.img);
    % Removing the b0 data from the loaded data. Q-ball does not deal with
    % the b0's
    req_4d_dwmri = nii.img(:,:,:,2:97);
    req_bvals = bvals(:,2:97);
    req_bvecs = bvecs(:,2:97);
    % Feeding cleaned data to qball fit
    
    % We can also provide a mask as an input as the last input parameter.
    % For now it has been left out.
    lmax = 6;
    lambda = 0.006;
    
    % L
    P0 = []; Laplac2 = [];
    for L=0:2:lmax
        for m=-L:L
            Pnm = legendre(L, 0); factor1 = Pnm(1);
            P0 = [P0; factor1];
            Laplac2 = [Laplac2; (L^2)*(L + 1)^2];
        end
    end
    L = diag(Laplac2);
    
    % ALLISON -> This is where I am using the construct_SH_basis, this
    % might be at a different place in your laptop, so call it
    % appropriately.
    [basis,~,~] = qball.lib.spherical_harmonics.construct_SH_basis(lmax,req_bvecs',2,'real');
    
    
    counter = 0;
    %output = update(req_4d_dwmri,basis,lambda,L,mask,return_flag,counter);
    counter_array = [];
    dwmri_op = req_4d_dwmri;

    sh_dwi_reg = zeros(size(dwmri_op));
    
    % Zeros matrix for storing the spatial maps
    spatial_map = zeros(78,93,75,7);
    
    for i = 1:size(dwmri_op,1)
        for j = 1:size(dwmri_op,2)
            for k = 1:size(dwmri_op,3)
                if mask(i,j,k)
                    % Reconstruct signal with spherical harmonic coefficients using regularized fit                   
                    C = (basis'*basis + lambda*L)\basis'*squeeze(dwmri_op(i,j,k,:));                  
                    % Reconstruct signal
                    sh_signal = basis*C;
                    % Store it
                    sh_dwi_reg(i,j,k,:) = sh_signal;
                end
            end
        end
    end
    
    residuals = dwmri_op - sh_dwi_reg;
    nonZeroIndices = residuals ~= 0;
    sigma = 1.4826 * (mad(residuals(nonZeroIndices)));
    
    ng_metric = [1:96]';
    ng_metric(:,2) = 0;
    threshold = 3;
    sigma_threshold = 5;
    iteration_count = 0;
    flags = 1;
    while(flags ~= 0)
        display(iteration_count)
        flags = 0;
        for i = 1:size(dwmri_op,1)
            for j = 1:size(dwmri_op,2)
                for k = 1:size(dwmri_op,3)
                    if mask(i,j,k)
                        for ng = 1:size(dwmri_op,4)
                            if (abs(residuals(i,j,k,ng)) > threshold * sigma)
                                diff_check = abs(residuals(i,j,k,ng)) - (threshold * sigma);
                                diff_percent = diff_check/(threshold * sigma);
                                %if (diff_percent > diff_threshold)
                                dwmri_op(i,j,k,ng) = sh_dwi_reg(i,j,k,ng);
                                flags = flags + 1;
                                ng_metric(ng,2) = ng_metric(ng,2) + 1;
                                %end
                            end
                        end
                    end
                end
            end
        end
        
        if (flags >= 1)
            iteration_count = iteration_count + 1; 
            for i = 1:size(dwmri_op,1)
                for j = 1:size(dwmri_op,2)
                    for k = 1:size(dwmri_op,3)
                        if mask(i,j,k)
                        % Reconstruct signal with spherical harmonic coefficients using regularized fit                   
                        C = (basis'*basis + lambda*L)\basis'*squeeze(dwmri_op(i,j,k,:));                  
                        % Reconstruct signal
                        sh_signal = basis*C;
                        % Store it
                        sh_dwi_reg(i,j,k,:) = sh_signal;
                        end
                    end
                end
            end
            residuals = dwmri_op - sh_dwi_reg;
            nonZeroIndices = residuals ~= 0;
            new_sigma = 1.4826 * (mad(residuals(nonZeroIndices)));
            sigma_difference = abs(new_sigma - sigma)/sigma;
            if (sigma_difference > sigma_threshold)
                sigma = new_sigma;
            else
                flags = 0;
            end
            
        
        end                      
    end

%{
clf;
figure(1)
hold on
plot(ng_metric(:,1),ng_metric(:,2));
ng_metric(:,1);
ng_metric(:,2);
xlabel('Gradient direction');
ylabel('Outliers');
text = sprintf('Signal Dropout, No. of Iterations - %d',iteration_count);
title(text);
grid on;

nii.img(:,:,:,2:97) = dwmri_op;
save_untouch_nii(nii,'corrected_dwmri.nii')
%}
corrected_dwmri = dwmri_op;        
end
