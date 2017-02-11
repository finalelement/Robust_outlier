function pipeline_acc

    
    %Regular Data
    scan1_nii = '/fs4/masi/nathv/qball_robust_experiment/patient_movement_data/regular_data/splitted_scans/scan1.nii';
    scan1_bval = '/fs4/masi/nathv/qball_robust_experiment/patient_movement_data/regular_data/splitted_scans/scan1.bval';
    scan1_bvec = '/fs4/masi/nathv/qball_robust_experiment/patient_movement_data/regular_data/splitted_scans/scan1.bvec';

    scan2_nii = '/fs4/masi/nathv/qball_robust_experiment/patient_movement_data/regular_data/splitted_scans/scan2.nii';
    scan2_bval = '/fs4/masi/nathv/qball_robust_experiment/patient_movement_data/regular_data/splitted_scans/scan2.bval';
    scan2_bvec = '/fs4/masi/nathv/qball_robust_experiment/patient_movement_data/regular_data/splitted_scans/scan2.bvec';
    
    mask1 = '/fs4/masi/nathv/qball_robust_experiment/patient_movement_data/regular_data/splitted_scans/scan1_strip_mask.nii';
   
    
    % Movement Data
    scan3_nii = '/fs4/masi/nathv/qball_robust_experiment/patient_movement_data/movement_data/splitted_scans/scan1.nii';
    scan3_bval = '/fs4/masi/nathv/qball_robust_experiment/patient_movement_data/movement_data/splitted_scans/scan1.bval';
    scan3_bvec = '/fs4/masi/nathv/qball_robust_experiment/patient_movement_data/movement_data/splitted_scans/scan1.bvec';

    scan4_nii = '/fs4/masi/nathv/qball_robust_experiment/patient_movement_data/movement_data/splitted_scans/scan2.nii';
    scan4_bval = '/fs4/masi/nathv/qball_robust_experiment/patient_movement_data/movement_data/splitted_scans/scan2.bval';
    scan4_bvec = '/fs4/masi/nathv/qball_robust_experiment/patient_movement_data/movement_data/splitted_scans/scan2.bvec';
    
    
    mask2 = '/fs4/masi/nathv/qball_robust_experiment/patient_movement_data/movement_data/splitted_scans/scan1_strip_mask.nii';
    
    [regular_acc,regular_corrected_acc] = acc_before_after_output(scan1_nii,scan1_bval,scan1_bvec,scan2_nii,scan2_bval,scan2_bvec,mask1);
    
    [movement_acc,movement_corrected_acc] = acc_before_after_output(scan3_nii,scan3_bval,scan3_bvec,scan4_nii,scan4_bval,scan4_bvec,mask2);
    
end