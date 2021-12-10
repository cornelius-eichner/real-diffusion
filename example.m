
% Set read/write function for nifti files.
% I personally use the Matlab Toolbox from Jimmy Shem https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

% Add Tools for NIfTI and ANALYZE image to the Matlab Path
PATH_NIITOOLS = '~/Software/niitools/'
addpath(PATH_NIITOOLS)


%%%%%%%%%%%%%%%%%%%%%
% LOAD NII DATA

% Download example dMRI dataset here:
% https://owncloud.gwdg.de/index.php/s/E2K5PjXp6XZCORL
Mag = load_untouch_nii('mag.nii.gz');
Pha = load_untouch_nii('pha.nii.gz');

Img_Pha = double(Pha.img);
Img_Mag = double(Mag.img);


%%%%%%%%%%%%%%%%%%%%%
% PREPRARE DATA 

% Unfortunately, phase data mostly not provided in absolute phase units. 
% For instance, Siemens Systems provide them in an arbitrary 12 bit format (from 0 to 4096). 
% My implementation of Real Diffusion requires a full complex dataset. 
% Therefore, the phase must be scaled to be from 0 to 2*pi

%Process phase data to be within [0 2*pi]
Img_Pha = Img_Pha-min(Img_Pha(:));
Img_Pha = 2*pi/4096 * Img_Pha;

%Make complex dataset
Img_Comp = double(Img_Mag .* exp(1i * Img_Pha));


%%%%%%%%%%%%%%%%%%%%
% CALC REAL DIFFUSION

% First, lets set up the parameters for the background field estimation
lambda = 5; % Regularization
beta = 50;	% Soft threshold

% Please note that those two parameters are a VITAL part of the real diffusion estimation. 
% As always, when it comes to regularized reconstructions, there is no ground truth value. 
% In my experience, the soft threshold (beta) has highest impact on the final reconstruction, 
% as the field is already very smooth in space 

% Run the real diffusion reconstruction
[Img_Comp_corr, BG_field] = real_diffusion(Img_Comp, lambda, beta);


%%%%%%%%%%%%%%%%%%%%
% SAVE RESULTS
Real_Path = 'real.nii.gz';
BG_Field_Path = 'BG_field.nii.gz';

% Just take the header information from the Mag image for the real data and BG phase
Real = Mag; 
BG = Pha; 

Real.hdr.dime.bitpix = 16;
Real.hdr.dime.datatype = 16;
Real.img = real(Img_Comp_corr);
save_untouch_nii(Real, Real_Path);

BG.img = angle(BG_field);
save_untouch_nii(BG, BG_Field_Path);
