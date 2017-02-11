function robust_qball_stats_v2(nifti_data)

    % Loading 4d nifti data, bvals and bvecs
    %nii = load_untouch_nii('dwmri.nii');
    nii = nifti_data;
    bvals = dlmread('dwmri.bval');
    bvecs = dlmread('dwmri.bvec');
    mask = load_untouch_nii('mask.nii');
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
    
    [basis,~,~] = qball.lib.spherical_harmonics.construct_SH_basis(lmax,req_bvecs',2,'real');
    
    
    counter = 0;
    %output = update(req_4d_dwmri,basis,lambda,L,mask,return_flag,counter);
    counter_array = [];
    dwmri_op = req_4d_dwmri;

    sh_dwi_reg = zeros(size(dwmri_op));
    threshold = 50;
    
    % Zeros matrix for storing the spatial maps
    spatial_map = zeros(78,93,75,7);
    
    for i = 1:size(dwmri_op,1)
        for j = 1:size(dwmri_op,2)
            for k = 1:size(dwmri_op,3)
                if mask(i,j,k)
                    flags = 1;
                    counter = 0;
                    while(flags ~=0)
                        % Reconstruct signal with spherical harmonic coefficients using regularized fit                   
                        C = (basis'*basis + lambda*L)\basis'*squeeze(dwmri_op(i,j,k,:));                  
                        % Reconstruct signal
                        sh_signal = basis*C;
                        % Store it
                        sh_dwi_reg(i,j,k,:) = sh_signal;
                        
                        % Calculate residuals
                        residuals = squeeze(dwmri_op(i,j,k,:) - sh_dwi_reg(i,j,k,:));

                        % Take regular standard deviation or median absolute deviation
                        sigma = 1.4826*mad(residuals,1);
                        
                        %outliers_idx = residuals > 3*sigma;
                        
                        %length(find(outliers_idx))
                        
                        flags = 0;
                        dbg_residue = [];
                        directions = [];
                        counter = counter + 1;
                        for ng = 1:size(dwmri_op,4)
                            if (abs(residuals(ng,1)) > threshold*sigma)
                                dwmri_op(i,j,k,ng) = sh_dwi_reg(i,j,k,ng);
                                flags = flags + 1;
                                dbg_residue = [dbg_residue;residuals(ng,1)];
                                directions = [directions;ng];
                            end
                        end
                        %display(counter)
                        %display(dbg_residue)
                        %display(sigma)
                        %display(threshold*sigma)
                        %display(directions)
                        
                        if (counter == 7)
                            display(counter)
                        end
                    end
                    counter_array = [counter_array;counter];
                    spatial_map(i,j,k,counter) = 1;
                end
            end
        end
    end

        
        %{
        % initial flag value/ exit condition as well
        flags = 0;
        
        % Calculate residuals
        residuals = dwmri_op - sh_dwi_reg;

        % Take regular standard deviation or median absolute deviation
        %maddy = std(residuals,0,4);
        maddy = mad(residuals,0,4);

        sig_adjusted = zeros(size(dwmri_op));
        for i = 1:size(dwmri_op,1)
            for j = 1:size(dwmri_op,2)
                for k = 1:size(dwmri_op,3)
                    if mask(i,j,k)
                        for ng = 1:size(dwmri_op,4)
                            if (residuals(i,j,k,ng) > threshold*maddy(i,j,k))
                                sig_adjusted(i,j,k,ng) = threshold*maddy(i,j,k);
                                flags = flags+1;
                            end
                        end
                    end
                end
            end
        end
        
        if (flags ~= 0)
            dwmri_op = sig_adjusted;
        end
        
    end
    %}
end