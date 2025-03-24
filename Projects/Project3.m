% Discrete Fourier Transform

function [ im_freq_filt_spat , im_freq_filt_freq ] =  filter_frequency_domain ( im , N , sig1 , sig2 )
    % Calculate G in frequency domain
    g = fspecial ('gaussian ', N , sig1 ) - fspecial ('gaussian ', N , sig2 ) ;
    g_pad = padarray (g , ([ size ( im ,1) size ( im ,2) ] - size ( g ) ) , 'post ') ;
    g_pad = circshift ( g_pad , - floor ( size ( g ) /2) ) ;
    G = fft2 ( g_pad ) ;
    % Filter in frequency domain
    imf = fft2 ( im ) ;
    im_freq_filt_freq = imf . * G ;
    im_freq_filt_spat = ifft2 ( im_freq_filt_freq ) ;
end

function conj_point_symm = is_conj_point_symm ( imf )
    % Use fftshift to have the zero frequency component at the center of the spectrum
    ffts = fftshift ( imf ) ;
    center = round ( size ( ffts ) / 2 +1) ; % Point of symmetry

    % Check entire matrix :
    % Copy the first row/col for to the end +1^ th row/col
    ffts (: , end +1) = ffts (: ,1) ;
    ffts ( end +1 , :) = ffts (1 ,:) ;

    % Get upper & lower matrix , check point - symmetry :
    upper = ffts (1: center (1) , 1: end) ;
    lower = ffts ( center (1) :end , 1: end) ;
    is_conj_symm = upper == conj ( rot90 ( lower , 2) ) ;
    conj_point_symm = all ( is_conj_symm (:) ) ;
end

function [x , y ] = find_part_location_image ( im , ims )
    % Pad the ims with zeros to have the same size with the original image
    imt = padarray ( ims , size ( im ) - size ( ims ) , 'post ') ;
    % Use the given formula to calculate the spatial shift
    G1 = fft2 ( imt ) ; G2 = fft2 ( im ) ;
    T = G1 . * conj ( G2 ) . / abs ( G1 . * conj ( G2 ) ) ;
    t = ifft2 ( conj ( T ) ) ;
    [c , i ] = max ( t (:) ) ;
    [y , x ] = ind2sub ( size ( im ) , i ) ;
end

% 2 - Separable Block Transform

function image_out = block_splitter_image ( image_in , block_size , block_func , block_arg )
% 1. Pad the input image based on the block size
    sz = size ( image_in ) ;
    pad = block_size - mod ( sz , block_size ) ;
    pad ( pad == block_size ) = 0;
    image_in = padarray ( image_in , pad , 'replicate ', 'post ') ;
    sz = sz + pad ;

    % 2. Allocate output image variable
    image_out = zeros ( size ( image_in ) ) ;

    % 3. Loop over blocks over x,y with block_size and apply block_func on each block
    for x_start = 1: block_size :( sz (2) - block_size +1)
        for y_start = 1: block_size :( sz (1) - block_size +1)
        block_in = image_in ( y_start : y_start + block_size -1 , x_start :x_start + block_size -1) ;
        block_out = block_func ( block_in , block_arg ) ;
        
        % 4. Save processed block onto output image into appropriate part
        image_out ( y_start : y_start + block_size -1 , x_start : x_start + block_size -1) = block_out ;
        end
    end
   
    % 5. Cut away padded areas
    image_out (end - pad (1) +1: end , :) = [];
    image_out (: , end - pad (2) +1: end) = [];
end

function block_out = block_processor ( block_in )
    block_out = fliplr ( block_in ) -0 .3 ;
    block_out ( block_out <0) = 0;
end

function haar_im = filter_image_haar ( im , block_size , haar_4 , fcnHandleBlockSplitter )
    % Haar transform on a block : A * block_in * A' where A is the transform matrix
    block_processor_linear = @ ( block_in , A ) A * block_in * A ';
    % Apply block splitter and block processor with the given '
    fcnHandleBlockSplitter
    haar_im = fcnHandleBlockSplitter ( im , block_size , block_processor_linear , haar_4 ) ;
end

function im_reorder = reorder_blocks ( im , block_size )
    sz = size ( im ) ; im_reorder = zeros ( sz ) ;
    rows = 1: block_size : sz (1) ; cols = 1: block_size : sz (2) ;
    lr = length ( rows ) ; lc = length ( cols ) ;
    for y = 0: block_size -1
        for x = 0: block_size -1
        block = im ( rows +y , cols + x ) ;
        % Block - based scaling
        % block = block - min( block (:));
        % block = block / max( block (:));
        im_reorder ( y * lr +1:( y +1) * lr , x * lc +1:( x +1) * lc ) = block ;
        end
    end
end

function [ im_mean , im_var ] = calculate_statistics ( haar_im , block_size )
    sz = size ( haar_im ) ;
    % Pad input image
    pad = block_size - mod ( sz , block_size ) ;
    pad ( pad == block_size ) = 0;
    haar_im = padarray ( haar_im , pad , 'replicate ', 'post ') ;
    sz = sz + pad ;
    image_out = zeros ( size ( haar_im ) ) ;
    all_blocks = zeros ( block_size , block_size , sz (1) * sz (2) / block_size'
    ^2) ;
    ibl = 0;

    % Loop over blocks
    for x_start = 1: block_size :( sz (2) - block_size +1)
        for y_start = 1: block_size :( sz (1) - block_size +1)
        ibl = ibl + 1;
        block = haar_im ( y_start : y_start + block_size -1 , x_start : x_start + block_size -1) ;
        all_blocks (1: block_size , 1: block_size , ibl ) = block ;
        end
    end
    % Calc stats
    im_mean = mean ( all_blocks , 3) ;
    im_var = var ( all_blocks , [] , 3) ;
end

function blur_A = return_blur_A_matrix ( n )
    blur_A = 0.5 * diag ( ones (1 , n ) , 0) + 0 .25 * diag ( ones (1 ,n -1) , -1) + 0 .25 * diag ( ones (1 ,n -1) , 1) ;
end

% now apply blur matrix to the image
n = 16;
blur_A = return_blur_A_matrix ( n ) ;
blur_lenna = block_splitter ( im , n , @block_processor_linear , blur_A ) ;
imshow ( blur_lenna ) ;