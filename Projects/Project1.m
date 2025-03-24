% 1 - image manipulation

% we load pears.jpg to variable im and analyze it
im = imread ('data / pears.jpg ') ;
datatype = class ( im ) ; sz = size ( im ) ;
im_min = min ( im (:) ) ; im_max = max ( im (:) ) ;
im_mean = mean ( im (:) ) ; im_stddev = std ( double ( im (:) ) ) ;
imhist ( im ) ;

imd_255 = double ( im ) ;
imd = (1/255) * imd_255 ;
figure (1) ; imshow ( imd ) ; % Works ; range 0 -1 , type double
figure (2) ; imshow ( im ) ; % Works : range 0 -255 , type uint8
figure (3) ; imshow ( imd_255 ) ; % Fails ; range 0 -255 , type double
figure (4) ; colormap ( gray (256) ) ; imagesc ( im /16) ;
figure (5) ; colormap ( jet (256) ) ; imagesc ( imd ) ;
figure (6) ; colormap ( gray (256) ) ; image ( imd_255 ) ;
%Image as 3D surface
figure (7) ; colormap ( gray (256) );
surface ([0 1] , [0 1] , 0 .5 * ones (2) , im , 'FaceColor ','texturemap ', 'CDataMapping ','direct ');
% Image viewer app
imtool ( imd );

% function to draw a 2 pixel wide black boarder
function imd_frame = im_add_frame ( im )
    frame_color = 0;
    imd_frame = im ;
    imd_frame ([1:2 , end -1: end] , :) = frame_color ;
    imd_frame (: , [1:2 , end -1: end ]) = frame_color ;
end

% function to copy every second row and every third column
function im_part = im_cut_part ( im )
    im_part = im (1:2: end , 1:3: end ) ;
end

% we use the Kronecker product to create colorized copies of one grayscale image
function im_kron_rgb = merge_channels_kronecker ( im )
    cR = [ 1 0 1 ; 0 .3 1 0; 0 0 .3 1 ];
    cG = [ 0 1 1 ; 0 .3 0 1; 0 .7 0 .3 0 ];
    cB = [ 0 0 0 ; 1 1 0; 0 .7 1 0 ];
    sz = size ( im ) ;
    im_k_R = zeros (3* sz ) ; im_k_G = zeros (3* sz ) ; im_k_B = zeros (3* sz ) ;
    for r = 1:3
        for c = 1:3
        range_r = (r -1) * sz (1) +1 : r * sz (1) ;
        range_c = (c -1) * sz (2) +1 : c * sz (2) ;
        im_k_R ( range_r , range_c ) = cR (r , c ) * im ;
        im_k_G ( range_r , range_c ) = cG (r , c ) * im ;
        im_k_B ( range_r , range_c ) = cB (r , c ) * im ;
        end
    end
    im_kron_rgb = zeros ([3* sz 3]) ;
    im_kron_rgb (: ,: ,1) = im_k_R ;
    im_kron_rgb (: ,: ,2) = im_k_G ;
    im_kron_rgb (: ,: ,3) = im_k_B ;
end

% 2 - Periodic Sequences
imd = im2double ( imread ('data / texture1.jpg ') ) ;

% corners of the parallelograms
loc = [ 10 10 ]; % loc = [ 1 1 ]; loc = [ 50 20 ];
% we are provided these periodicity matrices
N_A = [ 200 0; 0 200 ];
N_B = [ 200 200; 0 200 ];
N_C = [ 200 200; -200 200 ];

function [ atomic , mask ] = get_atomic ( im_ref , pm , loc_base )
    % function to get atomic element for a given image im_ref
    % Input arguments :
    % loc_base - starting location
    % pm - periodic matrix
    % Output arguments :
    % atomic - atomic element of the image
    % mask - corresponding mask for the atomic element
    % Define polygon :
    poly = zeros (2 ,5) ; poly (: ,1) = loc_base;
    poly (: ,2) = poly (: ,1) + pm (: ,1) ; poly (: ,3) = poly (: ,2) + pm (: ,2);
    poly (: ,4) = poly (: ,3) - pm (: ,1) ; poly (: ,5) = poly (: ,4) - pm (: ,2);
    % Create mask for atomic element ( parallelogram )
    poly = poly + repmat ([0;0 .4 ] - min ( poly ,[] , 2) , 1 , 5);
    sz_mask = round ( max ( poly ,[] ,2) );
    mask = poly2mask ( poly (1 ,:) , poly (2 ,:) , sz_mask (2) , sz_mask (1) );
    % Create atomic element
    atomic = im_ref ( loc_base (2) : loc_base (2) + sz_mask (2) -1 , loc_base (1) : loc_base (1) + sz_mask (1) -1);
    atomic = atomic . * mask;
end

[ atom_A , atom_mask_A ] = get_atomic ( imd , N_A , loc);

% function to generate a periodic image for the given atomic element
function [ im_out ] = create_periodic ( im_atomic , im_mask , pm , loc_base , size_output )
    im_out = zeros ( size_output ) ;
    im_mask_full = zeros ( size_output ) ;
    sz = fliplr ( size ( im_out ) );
    % Calculate required range for nx , ny
    p_ex_ = [ 1 1 ; sz (1) sz (2) ; sz (1) 1 ; 1 sz (2) ]' - repmat ( loc_base , 1 , 4) ;
    p_ex = inv ( pm ) * p_ex_ ;
    nx_range = floor ( min ( p_ex (1 ,:) ) -1) : ceil ( max ( p_ex (1 ,:) ) +1) ;
    ny_range = floor ( min ( p_ex (2 ,:) ) -1) : ceil ( max ( p_ex (2 ,:) ) +1) ;
    % "Copy - Paste " the atomic element multiple times , add to output image
    for nx = nx_range
        for ny = ny_range
        % Calculate pixel positions top - left and bottom - right :
        loc = loc_base + nx * pm (: ,1) + ny * pm (: ,2) ;
        loc2 = loc + flipud (size ( im_atomic )) - 1;
        
        % All positions this repetition addresses :
        all_x = loc (1) : loc2 (1) ;
        all_y = loc (2) : loc2 (2) ;
        atomic_x = 1: size ( im_atomic ,2) ;
        atomic_y = 1: size ( im_atomic ,1) ;
        
        % Remove cols / rows outside matrix :
        valid_x = and ( all_x ≥ 1 , all_x ≤ sz (1) ) ;
        valid_y = and ( all_y ≥ 1 , all_y ≤ sz (2) ) ;
       
        %if isempty ( valid_x ) || isempty ( valid_y ); continue ; end
        if all (¬valid_x ) || all (¬valid_y ) ; continue ; end

        % Copy the atomic element with transparency mask
        at_tran = im_atomic ( atomic_y ( valid_y ) , atomic_x ( valid_x ) ) ;
        im_out ( all_y ( valid_y ) , all_x ( valid_x ) ) = im_out ( all_y ( valid_y ), all_x ( valid_x ) ) + at_tran ;
       
        % And create sum of masks ( should be all 1)
        im_mask_full ( all_y ( valid_y ) , all_x ( valid_x ) ) = im_mask_full ( all_y ( valid_y ) , all_x ( valid_x ) ) + im_mask ( atomic_y ( valid_y ) , atomic_x ( valid_x ) ) ;
        end
    end
    
    assert ( all ( im_mask_full (:) ==1) , 'Your resulting mask is not all ones - you incorrectly sum something ') ;
end

function ncc = get_ncc ( t1 , t2 )
    sz = size ( t1 ) ;
    t1_norm = ( t1 - mean ( t1 (:) ) ) / std ( t1 (:) ) ;
    t2_norm = ( t2 - mean ( t2 (:) ) ) / std ( t2 (:) ) ;
    ncc = 1/( numel ( t1 ) ) * ( t1_norm . * t2_norm ) ;
    ncc = sum ( ncc (:) ) ;
end


function [ ncc_range , ncc_max_loc ] = compute_ncc_find_max ( im , range , fcnHandleNCC )
% Input arguments :
% im - given image (of type double in range 0..1)
% range - vector containing elements from 1 till 300
% fcnHandleNCC - function handle to " get_ncc (i1 ,i2)" -> see C2.3
% Output arguments :
% ncc_range - vector of NCC values for the edge lengths from " range "
% ncc_max - location of maximum in NCC
% Compute ncc in the given range
    sz = size ( im ) ;
    bx = 1; by = 1;
    ncc_range = zeros (1 , max ( range ) ) ;
    
    for si = 1: length ( range )
        s = range ( si ) ;
        if ( by + s > sz (1) || bx + s > sz (2) ) ; continue ; end
        
        % Cut atomic element
        atomic = im ( by : by +s -1 , bx : bx +s -1) ;
        % Create periodic image
        nr = ceil ( max ( sz ) / s ) ;
        im_per = repmat ( atomic , nr , nr ) ;
        im_per = im_per (1: sz (1) , 1: sz (2) ) ;
        % Calc NCC - for full image :
        i1 = im_per ; i2 = im ;
        ncc_range ( s ) = fcnHandleNCC ( i1 , i2 ) ;
    end

    % Find maximum and save the numbers you get
    [¬, ncc_max_loc ] = findpeaks (( ncc_range ) , ' MINPEAKHEIGHT ', 0 .8 ) ;
end