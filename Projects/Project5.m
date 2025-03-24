% 1 - Autoregressive Models

function [ acov , var , cov_h , cov_v ] = autocov ( im )
    neighborhood = 11;
    nh2 = floor ( neighborhood /2) ; % Neighborhood : nh2 elem. forward / backward
    
    % 1. subtract mean
    imm = im - mean ( im (:) ) ;
   
   % 2. exclude borders
    sz = size ( im ) ;
    nh = nh2 * 2 + 1;
    yrange = 1 + nh2 : sz (1) - nh2 ; % Submatrix to process ; exclude borders
    xrange = 1 + nh2 : sz (2) - nh2 ;
    acov = zeros ( nh ) ;
    
    % 3. compute autocovariance matrix
    for t_y = - nh2 : nh2
        for t_x = - nh2 : nh2
        tmp = imm ( yrange , xrange ) . * imm ( yrange + t_y , xrange + t_x ) ;
        acov ( t_y + nh2 +1 , t_x + nh2 +1) = mean ( tmp (:) ) ;
        end
    end
    
    % 4. get variance , covariance in horizontal and vertical directions
    var = acov ( nh2 +1 , nh2 +1) ;
    cov_h = acov ( nh2 +1 , nh2 +1+1) ;
    cov_v = acov ( nh2 +1+1 , nh2 +1) ;
end

function [ ar_gen , rh , rv , std_z ] = generate_ar_image ( cov_h , cov_v , var )
    n = 512; % Size of the image
    b = 100; % " Startup " border
   
    % compute parameters rh ,rv , std_z
    rh = cov_h / var ;
    rv = cov_v / var ;
    std_z = sqrt ( var *(1 - rh ^2) *(1 - rv ^2) ) ;
    
    % allocate and append border
    im = zeros ( b +n , b + n ) ;
    
    % generate AR process
    for i =2: b + n
        for j =2: b + n
        im (i , j ) = rv * im (i -1 , j ) + rh * im (i ,j -1) - rv * rh * im (i -1 ,j -1) + std_z * randn (1) ;
        end
    end
   
    % remove border
    ar_gen = im ( b +1: b +n , b +1: b + n ) ;
end

% 2 - Pyramid Representations

function G = compute_reconstruction_filter ()
    G = 1/4 * [ 1 2 1 ; 2 4 2 ; 1 2 1 ];
end

function laplacian_pyr = generate_laplacian_pyr ( im ,F , G )
    NL = 5; % number of levels
    laplacian_pyr = cell (1 ,5) ; % allocate memory for the Laplacian pyramid
    x = im ;
    
    for level =1:( NL -1)
        % Filter & Downsample
        xf = imfilter (x , F , 'replicate ') ;
        xf = xf (1:2: end , 1:2: end) ;
        % Upsample & filter
        xr = zeros (2* size ( xf ) ) ;
        xr (1:2: end , 1:2: end) = xf ;
        xr = imfilter ( xr , G , 'replicate ') ;
        % Difference image
        laplacian_pyr { level } = x - xr ;
        % For next loop :
        x = xf ;
    end

    laplacian_pyr { level +1} = x ;
end

function [F , G ] = compute_CDF_filter_kernel ( f1d , g1d )
    F = f1d * f1d ;
    G = g1d * g1d;
end

function var_pyramid = calculate_statistics_pyramid ( pyramid_structure )
    NL = 5;
    var_pyramid = zeros ( NL ,1) ;
    b = 4; % border margins
    for i = 1: NL
        cim1 = pyramid_structure { i };
        cim1 = cim1 ( b +1: end -b , b +1: end - b ) ;
        var_pyramid ( i ) = var ( cim1 (:) ) ;
    end
end

function [ im_reconstruct , im_psnr ] = reconstruct_image ( pyramid , im , G )
    NL = numel ( pyramid ) ;
    mode = 'replicate ';
    pyramid {1} = 0* pyramid {1};
    pyramid {2} = 0* pyramid {2};
    b = 4;
    xcur = pyramid { NL };
    
    for level = NL -1: -1:1
        % Upscale & filter previous ( coarser ) level :
        xr = zeros ( size ( pyramid { level }) ) ;
        xr (1:2: end , 1:2: end) = xcur ;
        xr = imfilter ( xr , G , mode ) ;
        % Add current level ( difference image ):
        xcur = xr + pyramid { level };
    end

    im_reconstruct = xcur ;
    mse = ( im - im_reconstruct ) . ^2;
    mse = mse ( b +1: end -b , b +1: end - b ) ; % Ignore border effects
    mse = mean ( mse (:) ) ;
    im_psnr = 10 * log10 (1/ mse ) ;
end

