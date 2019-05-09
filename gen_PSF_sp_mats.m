
function [Tsp, G_x, G_y] = gen_PSF_sp_mats(S, m2d, sthresh)

nch = size(S,3);
ncol = size(S,1);
nlin = size(S,2);
nvox = nlin * ncol;

im_size = [ncol, nlin];

%% index for mask
indx_m2d = find(m2d(:));
lindx_m2d = length(indx_m2d);


%% define the DFT matrix 
M = ncol * nlin;  %% total size

ts1 = zeros(1, M);

Dcol = fft(eye(ncol));
Dlin = fft(eye(nlin));

%% create the 2D undersampled DFT
II = 1;
for clin = 1:nlin
    for ccol = 1:ncol
        
        t1 = zeros(ncol, nlin);
        
        t1(ccol, clin) = 1;
        t1f = fft2(t1);
        t1fv = t1f(:);
        
        if (II == 1)
            DFM1 = t1fv(indx_m2d);
        end
        
        ts1(II) = DFM1' * t1fv(indx_m2d);
       
        II = II + 1;
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% For the sparse approx
[YY, II] = sort(abs(ts1), 'descend');

% t1s = zeros(length(t1), 1);
% t1s(II(1:sthresh)) = YY(1:sthresh);
% t1s = sparse(t1s);

t1s = sparse(ts1 .* (abs(ts1) > sthresh * max(abs(ts1)))); 
nnz(t1s) / M;

[t1sI] = find(t1s);
nt1s = nnz(t1s);

IndxI = zeros(nt1s * M, 1);
IndxJ = zeros(nt1s * M, 1);
VVals = zeros(nt1s * M, 1);

csp = 0;

for I = 1:M

    cst = csp + 1;
    csp = cst + nt1s - 1;
    
    IndxI(cst:csp) = t1sI + (I-1);
    IndxJ(cst:csp) = I;
    
    VVals(cst:csp) = t1s(t1sI)';
    
end

II = find(IndxI <= M);
IndxI = IndxI(II);
IndxJ = IndxJ(II);
VVals = VVals(II);

sp_mat = sparse(IndxI, IndxJ, VVals, M, M);
spd = diag(sp_mat);
sp_mat = sp_mat + sp_mat' - spdiags(spd, 0, M, M);

% sp_mat = toeplitz(t1s);

%% normalize the sp_mat
sp_mat = sp_mat / M;

%%%%%%%%%%%%%%%%%%%%%%
%% Generate the block tridiagonal G matrix
for diff_typ = [0:1]

    if (diff_typ == 0)
        
        %% fundamental diagonal blocks
        e = ones(ncol, 1);
        A1 = spdiags([-e e], -1:0, ncol, ncol);
        B1 = speye(ncol);
        
        %% total matrix info
        tnz = nnz(A1) * nlin + nnz(B1) * (nlin - 1) * 2;
        Gindx = zeros(tnz, 3);  %% row, col and val
        
        [rindx, cindx, val] = find(A1);
        
        %% loop through to create
        nsp = 0;
        rsp = 0;
        for II = 1:nlin
            
            
            %% diagonal block index
            rst = rsp + 1;
            rsp = rst + ncol - 1;
            
            nst = nsp + 1;
            nsp = nst + length(rindx) - 1;
            
            Gindx(nst:nsp, :) = [rindx + rst - 1, cindx + rst - 1, val];
            
            
        end
        
        %% use index for fast sparse
        G = sparse(Gindx(1:nsp,1), Gindx(1:nsp,2), Gindx(1:nsp,3), nvox, nvox);
        GI  = inv(G);        
        
        GI_x = GI;
        G_x = G;
        
        %%%%%%%%%%%%%%%
        %% form matrix for diagonal inversion
        II = 1;
        Dv = reshape(squeeze(S(:,:,II)), nvox, 1);
        DvM = spdiags(Dv, 0, nvox, nvox);
        Tsp = DvM' * sp_mat * DvM;
        
        for II = 2:nch
                                
            Dv = reshape(squeeze(S(:,:,II)), nvox, 1);
            DvM = spdiags(Dv, 0, nvox, nvox);
            Tsp = Tsp + DvM' * sp_mat * DvM;
            
        end
                
    elseif (diff_typ == 1)
        
        G = spdiags([-ones(nvox,1) ones(nvox,1)], [-nlin,0], nvox, nvox);
        GI = inv(G);
        
        GI_y = GI;
        G_y = G;        
        
    end
    
    
end

