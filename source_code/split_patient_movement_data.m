function split_patient_movement_data

    scan_nii = load_untouch_nii('dwmri.nii');
    scan_bval = dlmread('dwmri.bval');
    scan_bvec = dlmread('dwmri.bvec');
    
    % Scan 1 data
    scan1_img = scan_nii.img(:,:,:,1:65);
    scan1_bval = scan_bval(:,1:65);
    scan1_bvec = scan_bvec(:,1:65);
    
    % Scan 2 data
    scan2_img = scan_nii.img(:,:,:,1);
    scan2_bval = scan_bval(:,1);
    scan2_bvec = scan_bvec(:,1);

    scan2_img(:,:,:,2:65) = scan_nii.img(:,:,:,66:129);
    scan2_bval(:,2:65) = scan_bval(:,66:129);
    scan2_bvec(:,2:65) = scan_bvec(:,66:129);
    
    
    % Saving all the data
    scan_nii.hdr.dime.dim = [4    96    96    44   65     1     1     1];
    scan_nii.img = scan1_img;
    dlmwrite('scan1.bval',scan1_bval);
    dlmwrite('scan1.bvec',scan1_bvec);
    save_untouch_nii(scan_nii,'scan1.nii');
    
    scan_nii.img = scan2_img;
    dlmwrite('scan2.bval',scan2_bval);
    dlmwrite('scan2.bvec',scan2_bvec);
    save_untouch_nii(scan_nii,'scan2.nii');
end