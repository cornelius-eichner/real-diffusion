function [I_comb_corr, BG_field] = real_diffusion(I_comb, TV_reg, beta)


    %% Initialize output variables 
    I_comb_corr = complex(zeros(size(I_comb)));
    BG_field = complex(zeros(size(I_comb)));

    %% image sizes
    nlin = size(I_comb, 1);
    ncol = size(I_comb, 2);
    nslc = size(I_comb, 3);
    ntr = size(I_comb, 4);

    %% under-sampling mask (full)
    m2d = ones(nlin, ncol);
    msk = ones(nlin, ncol);
    sens = ones(nlin, ncol);

    %% generate point-spread function and gradient sparse matrices
    [Tsp, G_x, G_y] = gen_PSF_sp_mats(sens, m2d, 1e-6);
    G_xpy = (G_x + G_y) / 2;
    G_xmy = (G_x - G_y) / 2;

    %% Iterative TV recon with closed form iterations[
    %TV_reg = 5;  %% TV
    %beta = 5000;  %% soft-thresh for grad

    %%%%%%%%%%%%%%%
    %% The Split-Bregman operator
    %% this is the system matrix
    Cop = Tsp + TV_reg * (G_x' * G_x + G_y' * G_y + G_xpy' * G_xpy + G_xmy' * G_xmy);

    %% form the ilu                
    L = ichol(Cop,struct('type','ict','droptol',1e-02,'michol','off'));        
    LH = L';
   
    for cslc = 1:nslc
        for ctr = 1:ntr

            %% image data
            curr_image = squeeze(I_comb(:,:,cslc,ctr));

            %% define sens map (all ones, no sense)
            x_new = curr_image(:);

            %% zero-fill in FFt
            x0 = x_new(:);

            %% save the error
            x_new = x0;

            %% keep track of the sparsity patterns
            ydx_iter = [];
            ydy_iter = [];


            %% perform the split-bregman
            for t = 1:5


                %% gradient of image
                ydx = G_x * x_new(:);
                ydy = G_y * x_new(:);

                ydxpy = G_xpy * x_new(:);
                ydxmy = G_xmy * x_new(:);

                %% soft-threshold
                ydx = max(abs(ydx) - beta / 2, 0) .* sign(ydx);
                ydy = max(abs(ydy) - beta / 2, 0) .* sign(ydy);

                ydxpy = max(abs(ydxpy) - beta / 2, 0) .* sign(ydxpy);
                ydxmy = max(abs(ydxmy) - beta / 2, 0) .* sign(ydxmy);

                %% save it
                ydx_iter = [ydx_iter, ydx];
                ydy_iter = [ydy_iter, ydy];
                
                %% Solve the matrix problem for the update right-hand side
                b = x0 + TV_reg * (G_x' * ydx + G_y' * ydy + G_xpy' * ydxpy + G_xmy' * ydxmy);

                %% use CG simple no-precon here (could easily use ILU etc)
                x_old = x_new;            
                x_new = pcg(Cop, b, 1e-3, 1000, L, LH, x_old);

            end

            %% subtract off the phase
            I_comb_corr(:,:,cslc,ctr) = I_comb(:,:,cslc,ctr) .* exp(-1i * angle(reshape(x_new, nlin, ncol)));        
            
            %% save the BG field
            BG_field(:,:,cslc,ctr) = exp(1i * angle(reshape(x_new, nlin, ncol)));        

        end  %% Slc
    end  %% Tr
end