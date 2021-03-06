function landman_test_data_before_after_rodish

    scan1_nii = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_501/dwmri.nii';
    scan1_bval = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_501/dwmri.bval';
    scan1_bvec = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_501/dwmri.bvec';

    scan2_nii = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_1201/dwmri.nii';
    scan2_bval = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_1201/dwmri.bval';
    scan2_bvec = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_1201/dwmri.bvec';
    
    scan3_nii = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_1901/dwmri.nii';
    scan3_bval = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_1901/dwmri.bval';
    scan3_bvec = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_1901/dwmri.bvec';

    scan4_nii = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_2601/dwmri.nii';
    scan4_bval = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_2601/dwmri.bval';
    scan4_bvec = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128374_b3000_2601/dwmri.bvec';
    
    scan5_nii = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128442_b3000_401/dwmri.nii';
    scan5_bval = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128442_b3000_401/dwmri.bval';
    scan5_bvec = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128442_b3000_401/dwmri.bvec';

    scan6_nii = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128442_b3000_1101/dwmri.nii';
    scan6_bval = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128442_b3000_1101/dwmri.bval';
    scan6_bvec = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128442_b3000_1101/dwmri.bvec';
    
    scan7_nii = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128442_b3000_1801/dwmri.nii';
    scan7_bval = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128442_b3000_1801/dwmri.bval';
    scan7_bvec = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128442_b3000_1801/dwmri.bvec';

    scan8_nii = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128457_b3000_401/dwmri.nii';
    scan8_bval = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128457_b3000_401/dwmri.bval';
    scan8_bvec = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128457_b3000_401/dwmri.bvec';
    
    scan9_nii = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128457_b3000_1101/dwmri.nii';
    scan9_bval = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128457_b3000_1101/dwmri.bval';
    scan9_bvec = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128457_b3000_1101/dwmri.bvec';

    scan10_nii = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128457_b3000_1801/dwmri.nii';
    scan10_bval = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128457_b3000_1801/dwmri.bval';
    scan10_bvec = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128457_b3000_1801/dwmri.bvec';
    
    scan11_nii = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128457_b3000_2501/dwmri.nii';
    scan11_bval = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128457_b3000_2501/dwmri.bval';
    scan11_bvec = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/128457_b3000_2501/dwmri.bvec';

    global_mask = '/fs4/masi/nathv/qball_robust_experiment/qball_b3000/mask.nii';
    
    lmax = 6;
    lambda = 0.006;
    
    % Prepping data for RADISH processing
    [s1_clean,s1_clean_bval,s1_clean_bvec,logic_mask] = cleaning_nii_bval_bvec(scan1_nii,scan1_bval,scan1_bvec,global_mask);
    [s2_clean,s2_clean_bval,s2_clean_bvec,~] = cleaning_nii_bval_bvec(scan2_nii,scan2_bval,scan2_bvec,global_mask);
    [s3_clean,s3_clean_bval,s3_clean_bvec,~] = cleaning_nii_bval_bvec(scan3_nii,scan3_bval,scan3_bvec,global_mask);
    [s4_clean,s4_clean_bval,s4_clean_bvec,~] = cleaning_nii_bval_bvec(scan4_nii,scan4_bval,scan4_bvec,global_mask);
    [s5_clean,s5_clean_bval,s5_clean_bvec,~] = cleaning_nii_bval_bvec(scan5_nii,scan5_bval,scan5_bvec,global_mask);
    [s6_clean,s6_clean_bval,s6_clean_bvec,~] = cleaning_nii_bval_bvec(scan6_nii,scan6_bval,scan6_bvec,global_mask);
    [s7_clean,s7_clean_bval,s7_clean_bvec,~] = cleaning_nii_bval_bvec(scan7_nii,scan7_bval,scan7_bvec,global_mask);
    [s8_clean,s8_clean_bval,s8_clean_bvec,~] = cleaning_nii_bval_bvec(scan8_nii,scan8_bval,scan8_bvec,global_mask);
    [s9_clean,s9_clean_bval,s9_clean_bvec,~] = cleaning_nii_bval_bvec(scan9_nii,scan9_bval,scan9_bvec,global_mask);
    [s10_clean,s10_clean_bval,s10_clean_bvec,~] = cleaning_nii_bval_bvec(scan10_nii,scan10_bval,scan10_bvec,global_mask);
    [s11_clean,s11_clean_bval,s11_clean_bvec,~] = cleaning_nii_bval_bvec(scan11_nii,scan11_bval,scan11_bvec,global_mask);
    
    corrected_scan1 = robust_qball_corrected_volume(s1_clean,s1_clean_bval,s1_clean_bvec,logic_mask);
    corrected_scan2 = robust_qball_corrected_volume(s2_clean,s2_clean_bval,s2_clean_bvec,logic_mask);
    corrected_scan3 = robust_qball_corrected_volume(s3_clean,s3_clean_bval,s3_clean_bvec,logic_mask);
    corrected_scan4 = robust_qball_corrected_volume(s4_clean,s4_clean_bval,s4_clean_bvec,logic_mask);
    corrected_scan5 = robust_qball_corrected_volume(s5_clean,s5_clean_bval,s5_clean_bvec,logic_mask);
    corrected_scan6 = robust_qball_corrected_volume(s6_clean,s6_clean_bval,s6_clean_bvec,logic_mask);
    corrected_scan7 = robust_qball_corrected_volume(s7_clean,s7_clean_bval,s7_clean_bvec,logic_mask);
    corrected_scan8 = robust_qball_corrected_volume(s8_clean,s8_clean_bval,s8_clean_bvec,logic_mask);
    corrected_scan9 = robust_qball_corrected_volume(s9_clean,s9_clean_bval,s9_clean_bvec,logic_mask);
    corrected_scan10 = robust_qball_corrected_volume(s10_clean,s10_clean_bval,s10_clean_bvec,logic_mask);
    corrected_scan11 = robust_qball_corrected_volume(s11_clean,s11_clean_bval,s11_clean_bvec,logic_mask);
    
    
    
end
