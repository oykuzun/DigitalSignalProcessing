% 1 - Image Sampling

% function to generate a regular and rectangular sampling grid
function [ gridX , gridY ] = generate_rect_grid ( sz , N )
    grid_ = linspace (1 , sz (1) , N ) ;
    [ gridX , gridY ] = meshgrid ( grid_ ) ;
end

function [ im_linear , im_cubic ] = interp_im_rect_grid ( im , gridX , gridY )
    sz = min ( size ( im ) ) ;
    [ imX , imY ] = meshgrid (1: sz ) ; % image is given at these positions
    % now implement the interpolation using linear and cubic interp
    im_linear = interp2 ( imX , imY , im , gridX , gridY , 'linear ') ;
    im_cubic = interp2 ( imX , imY , im , gridX , gridY , 'cubic ') ;
end

function im_rec = reconstruct_from_smaller_image ( im_small , sz , gridX , gridY )
    N = size ( gridX , 1) ;
    [ imX , imY ] = meshgrid (1: sz ) ;
    
    % 1. Apply triscatteredinterp function
    im_interp = TriScatteredInterp ( gridX (:) , gridY (:) , im_small (:) , 'linear') ;
   
    % 2. Then interpolate image to the original image coordinates
    im_rec = im_interp ( imX , imY ) ;
end

function [ hexaX , hexaY ] = generate_hexagonal_grid ()
    sz = 512;
    x_sz = sz -1;
    dx = ( x_sz -1) /(372 -1) ; % horizontal displacement in hexaX
    grid_hx = 1: dx :( x_sz ) ;
    dy = dx * sqrt (3) /2; % vertical displacement in hexaY
    grid_hy = 1: dy : sz (1) ;
    [ hexaX , hexaY ] = meshgrid ( grid_hx , grid_hy ) ;
    % Shift each 2nd row ( instead of separate processing for even /odd):
    hexaX (2:2: end , :) = hexaX (2:2: end , :) + 0 .5 * dx ;
end

function [ im_rec ] = interpolate_reconstruct ( im , hexaX , hexaY )
    sz = min ( size ( im ) ) ;
    [ imX , imY ] = meshgrid (1: sz ) ;
    % Interpoolate to hexagonal grid positions
    im_h = interp2 ( imX , imY , im , hexaX , hexaY , 'linear ') ;
    % Reconstruct the image to its original size
    si_h_rec = TriScatteredInterp ( hexaX (:) , hexaY (:) , im_h (:) , 'linear ') ;
    im_rec = si_h_rec ( imX , imY ) ;
end

function psnr_diff = compute_psnr_diff ( im , im_rec )
    % implement psnr computation according to the formula
    im_diff = im_rec - im ;
    border = 5;
    mse = im_diff ( border +1: end - border , border +1: end - border ) . ^2;
    mse = mean ( mse (:) ) ;
    psnr_diff = 10 * log10 (1/ mse ) ;
end

% 2 -  Zone Plate
function [D , k , Izp ] = calculate_zone_plate ( N )
    t = (0: N -1) ;
    X = repmat (t ,N ,1) ; Y = X ';
    D = sqrt ( X. ^2 + Y. ^2) ;
    dmax = max ( D (:) ) ;
    k = 1 / ( D (end ,end) ^2 - D (end -1 ,end -1) ^2) ;
    cosarg = 2* pi * k /2 * D. ^2;
    Izp = 0.5 + 0.5 * cos ( cosarg ) ;
end

N = 350;
[D , k , Izp ] = calculate_zone_plate ( N ) ;
Izp_21 = Izp (1:2: end ,1:2: end) ;
Izp_41 = Izp (1:4: end ,1:4: end) ;
figure
title ('Zone plate 2:1 subsampling ') ;
imshow ( Izp_21 ) ;
figure
title ('Zone plate 4:1 subsampling ') ;
imshow ( Izp_41 ) ;
imwrite ( Izp_21 , 'res/ Izp_21.png ') ;
imwrite ( Izp_41 , 'res/ Izp_41.png ') ;

% 3 - Image Filters

% Gaussian filters
Hg1 = fspecial ('gaussian ', 3 , 0 .5 ) ;
figure ; freqz2 ( Hg1 ) ; title ('gaussian -3 -0 .5 ') ;
Hg2 = fspecial ('gaussian ', 8 , 0 .5 ) ;
figure ; freqz2 ( Hg2 ) ; title ('gaussian -8 -0 .5 ') ;
Hg3 = fspecial ('gaussian ', 8 , 1 .5 ) ;
figure ; freqz2 ( Hg3 ) ; title ('gaussian -8 -1 .5 ') ;
Hg4 = fspecial ('gaussian ', 8 , 3) ;
figure ; freqz2 ( Hg3 ) ; title ('gaussian -8 -3 ') ;
Hg5 = fspecial ('gaussian ', 64 , 3) ;
figure ; freqz2 ( Hg5 ) ; title ('gaussian -16 -3 ') ;
Hg6 = fspecial ('gaussian ', 16 , 5) ;
figure ; freqz2 ( Hg6 ) ; title ('gaussian -16 -5 ') ;

% Laplacian filter
Hl1 = fspecial ('laplacian ') ;
figure ; freqz2 ( Hl1 ) ; title ('laplacian ') ;

% LoG: Laplacian of Gaussian filter
Hlog1 = fspecial ('log ') ;
figure ; freqz2 ( Hlog1 ) ; title ('log ') ;

% Prewitt
Hp1 = fspecial ('prewitt ') ';
figure ; freqz2 ( Hp1 ) ; title ('prewitt ') ;

function prew_z = prewitt_z ()
    % Define the prewitt filter
    Hp1 = fspecial ('prewitt ') ;
    N =64; % default for: freqz2 (h) uses [Ny Nx] = [64 64]
    % Compute the NxN frequency spectrum with 0 - padding :
    Fp1 = fft2 ( Hp1 , N , N ) ;
    prew_z = abs ( fftshift ( Fp1 ) ) ;
end

function im_filtered = filter_image ( im )
    Hg5 = fspecial ('gaussian ', 16 , 3) ;
    im_filtered = imfilter ( im , Hg5 , 'replicate ') ;
end

% Load original image
im = im2double ( imread ('data / pears.jpg ') ) ;
figure ; imshow ( im ) ; title ('Original image ') ;

% Prewitt
im_p1 = 0 .5 + imfilter ( im , Hp1 ) ;
figure ; imshow ( im_p1 ) ; title ('prewitt horizontal edges ') ;
im_p2 = 0 .5 + imfilter ( im , Hp1 ') ;
figure ; imshow ( im_p2 ) ; title ('prewitt vertical edges ') ;

% Gaussian filters
im_g4 = imfilter ( im , Hg4 ) ;
figure ; imshow ( im_g4 ) ; title ('gaussian -8 -3 ') ;
im_g5 = imfilter ( im , Hg5 ) ;
figure ; imshow ( im_g5 ) ; title ('gaussian -32 -3 ') ;
imf_g5_replicate = imfilter ( im , Hg5 , 'replicate ') ;
figure ; imshow ( imf_g5_replicate ) ; title ('image C') ;

function [ diag_filter_subsampled , diag_nonfilter_subsampled ] = '
    filter_subsample_zone_plate ( Izp )
    F = fspecial ('gaussian ', 16 , 3) ;
    tmp1 = imfilter ( Izp , F , 'symmetric ') ; % filter image
    tmp1 = tmp1 (1:2: end ,1:2: end) ; % subsample
    diag_filter_subsampled = diag ( tmp1 ) ;
    diag_nonfilter_subsampled = diag ( Izp (1:2: end ,1:2: end) ) ;
end