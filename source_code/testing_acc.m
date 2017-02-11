function testing_acc

    scan1_nii = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_1201/dwmri.nii';
    scan1_bval = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_1201/dwmri.bval';
    scan1_bvec = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_1201/dwmri.bvec';

    scan2_nii = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_501/dwmri.nii';
    scan2_bval = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_501/dwmri.bval';
    scan2_bvec = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_501/dwmri.bvec';

    global_mask = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/mask.nii';
    lmax = 6;
    lambda = 0.006;
        
    % Prepping data for RADISH processing
    [s1_clean,s1_clean_bval,s1_clean_bvec,logic_mask] = cleaning_nii_bval_bvec(scan1_nii,scan1_bval,scan1_bvec,global_mask);
    [s2_clean,s2_clean_bval,s2_clean_bvec,~] = cleaning_nii_bval_bvec(scan2_nii,scan2_bval,scan2_bvec,global_mask);

    % After RADISH
    corrected_scan1 = robust_qball_corrected_volume(s1_clean,s1_clean_bval,s1_clean_bvec,logic_mask);
    corrected_scan2 = robust_qball_corrected_volume(s2_clean,s2_clean_bval,s2_clean_bvec,logic_mask);
    
    [corrected_scan1_sh_coeffs_vol, corrected_scan1_exitcode_vol] = qball.vol_fit(corrected_scan1,s1_clean_bvec,s1_clean_bval,lmax,lambda,logic_mask);
    [corrected_scan2_sh_coeffs_vol, corrected_scan2_exitcode_vol] = qball.vol_fit(corrected_scan2,s2_clean_bvec,s2_clean_bval,lmax,lambda,logic_mask);
    
    % Before RADISH
    [scan1_sh_coeffs_vol, scan1_exitcode_vol] = qball.vol_fit(s1_clean,s1_clean_bvec,s1_clean_bval,lmax,lambda,logic_mask);
    [scan2_sh_coeffs_vol, scan2_exitcode_vol] = qball.vol_fit(s2_clean,s2_clean_bvec,s2_clean_bval,lmax,lambda,logic_mask);
    
    acc = zeros(93,78,75);
    corrected_acc = zeros(93,78,75);
    
    for i = 1:size(s1_clean,1)
        for j = 1:size(s1_clean,2)
            for k = 1:size(s1_clean,3)
                
                acc(i,j,k) = angularCorrCoeff(scan1_sh_coeffs_vol(i,j,k,:),scan2_sh_coeffs_vol(i,j,k,:));
                corrected_acc(i,j,k) = angularCorrCoeff(corrected_scan1_sh_coeffs_vol(i,j,k,:),corrected_scan2_sh_coeffs_vol(i,j,k,:));
                
            end
        end
    end
    
end
