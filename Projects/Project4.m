% 1 - Discrete Cosine Transform (DCT)

function dct_image = transform_to_dct_domain ( im , dct_matrix , fcnHandleBlkSplitter )
    % DCT on a block : A * block_in * A' where A is the transform matrix
    linear_block_transform = @ ( block_in , A ) A * block_in * A';
    block_sz = size ( dct_matrix ,1) ;
    % Apply block splitter and block processor with the given fcnHandleBlockSplitter
    dct_image = fcnHandleBlkSplitter ( im , block_sz , linear_block_transform , dct_matrix ) ;
end

function im_rec_comp = compress_dct ( dct_image , dct_matrix , fcnHandleBlkSplitter , N )
    block_size = 16;
    block_processor_linear = @ ( block_in , A ) A * block_in * A';
   
    % 1. remove all but first N coefficients
    d = zeros ( block_size ,1) ;
    d (1: sqrt ( N ) ) = 1;
    A_ = diag ( d ) ;
    im_dct_comp = fcnHandleBlkSplitter ( dct_image , block_size , block_processor_linear , A_ ) ;
    
    % 2. back - transform the image using given function handles
    im_rec_comp = fcnHandleBlkSplitter ( im_dct_comp , block_size , '
    block_processor_linear , dct_matrix ') ;
end

function im_dct = apply_3ddct ( im_in , dct_mat )

    % 1. vertical columns
    for i =1: size ( im_in ,3)
        for j =1: size ( im_in ,2)
        vert (: ,j , i ) = dct_mat * im_in (: ,j , i ) ;
        end
    end

    % 2. horizontal rows
    for i =1: size ( im_in ,3)
        for j =1: size ( im_in ,1)
        hor (j ,: , i ) = dct_mat * vert (j ,: , i )';
        end
    end

    % 3. out -of - plane rotations
    for i =1: size ( im_in ,1)
        for j =1: size ( im_in ,2)
        im_dct (i ,j , size ( im_in ,3) : -1:1) = dct_mat * squeeze ( hor (i ,j , end: -1:1) ) ;
        end
    end
end

function [ video_dc , video_hf ] = threed_dct_analysis ( video , dct_mat , fcnHandle3DDCTTransformation )
    sz_b = 2; % block size
    sz_v = size ( video ) ;
    video_out = zeros ( sz_v ) ;
    for i =1: sz_b : sz_v (1) -1
        for j =1: sz_b : sz_v (2) -1
            for k =1: sz_b : sz_v (3) -1
            temp = video ( i : i +1 , j : j +1 , k : k +1) ;
            im_out = fcnHandle3DDCTTransformation ( temp , dct_mat ) ;
            video_out ( i : i +1 , j : j +1 , k : k +1) = im_out ;
            end
        end
    end

    video_dc = video_out (1: sz_b :end , 1: sz_b :end ,2: sz_b :end) ;
    video_hf = video_out (2: sz_b :end , 2: sz_b :end ,1: sz_b :end) ;
end

% Karhunen-Loeve-Transform (KLT)
function im_reshape = reshape_image ( im , block_splitter )
    N =32;
    v = zeros ( N ^2 ,1) ;
    v ([529:544 ,272:32:784]) = 0 .1 ;
    func_handle = @ ( block , add ) reshape ( block (:) + add , N , N ) ;
    im_reshape = block_splitter ( im , N , func_handle , v ) ;
end

function [ blocks , C_im , klt_base ] = compute_klt_basis ( im )
    N =8;

    % 1. Allocate memory for blocks
    sz = size ( im ) ;
    x_arr = 1: N :( sz (2) - N +1) ;
    y_arr = 1: N :( sz (1) - N +1) ;
    blocks = zeros ( length ( x_arr ) * length ( y_arr ) , N * N ) ;

    %2. Regroup the input image into blocks
    for iy = 1: length ( y_arr )
        for ix = 1: length ( x_arr )
        x_start = x_arr ( ix ) ; y_start = y_arr ( iy ) ;
        block_in = im ( y_start : y_start +N -1 , x_start : x_start +N -1) ;
        blocks (( iy -1) * length ( x_arr ) + ix , :) = block_in (:) ;
        end
    end

    %3. Compute autocorrelation over all blocks
    C_im = blocks '* blocks / size ( blocks ,1) ;
   
    %4. Compute eigenvalue decomposition of the autocorrelation matrix and'
    thus get the KLT basis
    [ klt_base_v , ew ] = eig ( C_im ) ;
    [ ew , ew_idx ] = sort ( diag ( ew ) , 'descend ') ;
    klt_base_v = klt_base_v (: , ew_idx ) ;
    klt_base = zeros ( N ^2) ;
    for x = 0: N -1
        for y = 0: N -1
        idx = x * N + y +1;
        klt_base ( y * N +1:( y +1) *N , x * N +1:( x +1) * N ) = ...
        reshape ( klt_base_v (: , idx ) , N , N ) ;
        end
    end
end

function [ im_projected , im_backprojected , klt_psnr ] = project_reconstruct_klt_basis ( im , klt_base_v , block_splitter )
    N = sqrt ( size ( klt_base_v ,1) ) ; % block size

    % 1.a First defined a block processor for orthogonal base for projecting
    function block_out = block_processor_orth_base_proj ( block_in , base )
    block_out = reshape ( base'* reshape ( block_in , N *N ,1) ,N , N ) ;
    end

    % 1.b. Define block processor for othogonal base backprojecting
    function block_out = block_processor_orth_base_backproj ( block_in , base)
    block_out = reshape ( base * reshape ( block_in , N *N ,1) ,N , N ) ;
    end

    % 2. Then using a block splitter project the image
    im_projected = block_splitter ( im, N, @block_processor_orth_base_proj, klt_base_v ) ;

    % 3. Then backproject to the original image
    im_backprojected = block_splitter ( im_projected ,N , @block_processor_orth_base_backproj , klt_base_v ) ;

    % 4. And compute PSNR of the difference
    im_diff = ( im_backprojected - im ) . ^2;
    klt_psnr = 10 * log10 (1/ mean ( im_diff (:) ) ) ;
end

function autocorr = compute_autocorrelation_klt_images ( klt_projected_image , N )
    % 1. Allocate memory for blocks
    sz = size ( klt_projected_image ) ;
    x_arr = 1: N :( sz (2) - N +1) ;
    y_arr = 1: N :( sz (1) - N +1) ;
    blocks = zeros ( length ( x_arr ) * length ( y_arr ) , N * N ) ;
   
    %2. Regroup the input image into blocks
    for iy = 1: length ( y_arr )
        for ix = 1: length ( x_arr )
        x_start = x_arr ( ix ) ; y_start = y_arr ( iy ) ;
        block_in = klt_projected_image ( y_start : y_start +N -1 , x_start : x_start +N -1) ;
        blocks (( iy -1) * length ( x_arr ) + ix , :) = block_in (:) ;
        end
    end

    %3. Compute autocorrelation over all blocks
    autocorr = blocks' * blocks / size ( blocks ,1) ;
end

% Lenna Image
im_lenna = im2double ( imread ('data / lenna_gray.jpg ') ) ;

% KLT basis images of Lenna image
load (' klt_base_v.mat ', 'klt_base_v ') ;
block_splitter = @block_splitter_image ;
[ im_projected_lenna , im_backprojected_lenna , klt_psnr_lenna ] = project_reconstruct_klt_basis ( im_lenna , klt_base_v , block_splitter ) ;

% Calculate Statistics
N = 8;
[ klt_mean_lenna , klt_var_lenna ] = calculate_statistics (im_projected_lenna , N ) ;
figure , surf ( log ( klt_var_lenna ) ) ; title ('Log - variances ') ; view (75.5, 38) ;
xlabel ('X') ; ylabel ('Y') ;

% Mandrill Image
im_mandrill = im2double ( imread ('data / mandrill_gray.png ') ) ;
% KLT basis images of Lenna image
load (' klt_base_v.mat ', 'klt_base_v ') ;
block_splitter = @block_splitter_image ;
[ im_projected_mandrill , im_backprojected_mandrill , klt_psnr_mandrill ] = project_reconstruct_klt_basis ( im_mandrill , klt_base_v , block_splitter ) ;

% Calculate Statistics
N = 8;
[ klt_mean_mandrill , klt_var_mandrill ] = calculate_statistics (im_projected_mandrill , N ) ;
figure , surf ( log ( klt_var_mandrill ) ) ; title ('Log - variances ') ; view (75.5 ,38) ;
xlabel ('X') ; ylabel ('Y') ;

% KLT
im_lenna = im2double ( imread ('data / lenna_gray.jpg ') ) ;
% KLT basis images of Lenna image
load (' klt_base_v.mat ', 'klt_base_v ') ;
block_splitter = @block_splitter_image ;
[ im_projected_lenna , im_backprojected_lenna , klt_psnr_lenna ] = project_reconstruct_klt_basis ( im_lenna , klt_base_v , block_splitter ) ;

% Calculate Statistics
N = 8; % block size
[ klt_mean_lenna , klt_var_lenna ] = calculate_statistics (im_projected_lenna , N ) ;
figure , semilogy ( klt_var_lenna (:) ) ;

% DCT
im_lenna = im2double ( imread ('data / lenna_gray.jpg ') ) ;
M = 8; % block size
dct_matrix = dctmtx ( M ) ;
block_splitter = @block_splitter_image ;
% dct transformed image
dct_image = transform_to_dct_domain ( im_lenna , dct_matrix , block_splitter ) ;

% Calculate Statistics
[ dct_mean_lenna , dct_var_lenna ] = calculate_statistics ( dct_image , M ) ;
dct_var_lenna_ordered = ZigZag8x8 ( dct_var_lenna ) ;
figure , semilogy ( dct_var_lenna_ordered (:) ) ;