Diffusion MRI Phase Data Reconstruction
The here provided steps allow reconstruction of complex-valued diffusion MRI data on Siemens Scanners. Typically, complex phase reconstruction is not selectable as an option for diffusion sequences. Therefore, these phase data need to be collected using retro-reconstruction of previously acquired dMRI data. 


Procedure:
1. Enable “Advanced User Mode” on your Scanner. Your on-site MRI Physicist can support you here. 
2. Start the “Twix Raid Mode” Software by accessing Start > Run > Twix
3. The Twix software displays a list of scans, which have been recently acquired. 
4. Find your dMRI scan of interest in the displayed list. Select your scan, using a single mouse click.
5. In the Twix icon bar, click on “Toggle Reconstruction View”. A new pane, “Retro Recon”, will open below the Twix Scan list.
6. In this new pane, click “Edit”. The Software “XBuilder” will open, displaying all reconstruction information as a structure. 
7. Use the XBuilder find tool (binoculars icon) from the icon bar to search for the reconstruction variable “ucReconstructionMode”. Typically, for magnitude reconstruction, this parameter will be set to “1”.
8. For a reconstruction of both magnitude and phase, change the variable to “8”. Save the file and close XBuilder. 
9. In Twix, you can now click on “Start”. The blue Question Mark Icon will change to a green arrow, indicating the ongoing retro-reconstruction. Once the retro-reconstruction is finished, this will be indicated with a green tick in the “Retro Recon” pane. The retro-reconstructed data are stored in the local patient data folder. They can be distinguished from normal reconstruction by the ‘_RR’ suffix.  

Note:
This document summarizes personal experience with retro-reconstruction on Siemens scanners. There is no guarantee for this strategy producing desired results - especially in future software releases. Please discuss all described steps with your MR Site Scientists before applying them. 
Do not employ the here described retro-reconstruction of phase data for clinical applications
Prevent from using data filters such as Prescan Normalize or Raw Filter. These processing steps are tailored to Magnitude dMRI settings. 
When aiming to retro reconstruct data, it is advisable to not use the vendor provided dMRI reconstructions such as Trace-weighted, ADC, FA, Tensor. Please only save the diffusion-weighted data in Mosaic Format. 
