function [req_4d_dwmri,req_bvals,req_bvecs,mask] = cleaning_nii_bval_bvec(nifti_data,bval,bvec,mask_input)
    
    nii = load_untouch_nii(nifti_data);
    %nii = nifti_data;
    bvals = dlmread(bval);
    bvecs = dlmread(bvec);
    
    dimensions = size(nii.img);
    
    mask = load_untouch_nii(mask_input);
    mask = logical(mask.img);
    
    req_4d_dwmri = nii.img(:,:,:,2:dimensions(4));
    req_bvals = bvals(:,2:dimensions(4));
    req_bvecs = bvecs(:,2:dimensions(4));

end