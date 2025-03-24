% Subband Decomposition for DCT

function im_filter_bank = apply_filter_bank ( im , dct8 )
    N =8;
    sz = size ( im ,1) ; szN = sz / N ;
    im_filter_bank = zeros ( sz ) ;

    for r = 0: N -1 % Filter bank level 1 along rows
        kernel1 = dct8 ( r +1 , :) ';
        im_fc = imfilter ( im , kernel1 , 'circular ') ;
        im_fc = im_fc (( N /2) : N :end ,:) ; % Subsample
        for c = 0: N -1 % Filter bank level 2 along columns
            kernel2 = dct8 ( c +1 , :) ;
            im_fr = imfilter ( im_fc , kernel2 , 'circular ') ;
            im_fr = im_fr (: , ( N /2) : N :end) ;
            im_filter_bank (1+ szN * r : szN *( r +1) , 1+ szN * c : szN *( c +1) ) = im_fr ;
        end
    end
end

function pr_rec = verify_condition ( a0 , b0 , c0 , a1 , b1 , c1 )
    syms z
    
    % define here F0 ,F1 in terms of z
    F0 ( z ) = a0 + b0 * z ^ -1+ c0 * z ;
    F1 ( z ) = a1 + b1 * z ^ -1+ c1 * z ;
    F ( z ) =[ F0 ( z ) , F1 ( z ) ; F0 ( - z ) , F1 ( - z ) ];
    G0 ( z ) =(2/ det ( F ( z ) ) ) * F1 ( - z ) ;
    G1 ( z ) =( -2/ det ( F ( z ) ) ) * F0 ( - z ) ;
   
    % Verify the conditions
    cond1 = double ( simplify (( G0 ( z ) * F0 ( z ) + G1 ( z ) * F1 ( z ) ) ) ) ==2;
    cond2 = double ( simplify (( G0 ( z ) * F0 ( - z ) + G1 ( z ) * F1 ( - z ) ) ) ) ==0;
    pr_rec = cond1 && cond2 ;
end


% 2 - Lifting Structure with Motion Compensation

[ video , positions , mask ] = make_video_web () ;
nframes = size ( video , 3) ;
figure ;

for fi = 1: nframes
    imshow ( video (: ,: , fi ) , ' InitialMagnification ', 100) ;
    pause (0 .05 ) ;
end

function [ predict1 , update1 ] = predict_next_frame ( video , positions , mask )
    % 1. put even and odd frames into separate banks
    dE = video (: ,: ,1:2: end) ;
    dO = video (: ,: ,2:2: end) ;
    gray = 0.5 * ones ( size ( mask ) ) ;
    i = 1;
    i2 = 2*( i -1) +1;
    
    % Prediction
    % Cut out the object from frame n and fill the uncovered background with gray
    pred = paste_mask ( dE (: ,: , i ) , positions (: , i2 ) , gray , [1 1] , mask ) ;
     % Paste the object to the new location
    pred = paste_mask ( pred, positions (: , i2 +1), dE (: ,: , i ) , positions (: , i2 ), mask ) ;
   % Negate the prediction
    predict1 = - pred ;
    
    % Update
    % Obtain the difference image
    H = dO (: ,: , i ) + predict1 ;
    % Move the object region in the difference signal back to where it was'
    in the even frame
    upd = paste_mask ( H (: ,: , i ) , positions (: , i2 ) , H (: ,: , i ) , positions (: , i2+1) , mask ) ;
    update1 = 0 .5 * upd ;
end

function [ var_0lp , var_hp , im_lp ] = perform_analysis_multilevel ( video , positions , mask )
    % Lifting analysis function : returns the LP and HP subbands
    function [L , H ] = lift_analysis ( video , positions , mask )
        dE = video (: ,: ,1:2: end) ;
        dO = video (: ,: ,2:2: end) ;
        H = zeros ( size ( dO ) ) ;
        L = zeros ( size ( dE ) ) ;
        % Analysis
        for i = 1: size ( video ,3) /2
            [ L (: ,: , i ) , H (: ,: , i ) ] = predict_next_frame_index ( video , positions , mask , i ) ;
        end 
    end

    % Repeat the analysis hierarchically and compute the variances
    pyr_var_H = cell (1) ;
    pyr_var_L = cell (1) ;
    L = video ;
    
    for j = 1: log2 ( size ( video ,3) )
    [L , H] = lift_analysis (L , positions , mask ) ;
    positions = positions (: , 1:2: end) ;
    s = size ( L ) ;

    if numel ( s ) <3; s (3) =1; end
        H_re = reshape (H , s (1) * s (2) , s (3) ) ;
        L_re = reshape (L , s (1) * s (2) , s (3) ) ;
        pyr_var_H { j } = var ( H_re ) ;
        pyr_var_L { j } = var ( L_re ) ;
    end

    im_lp = L ;
    var_0lp = pyr_var_L {1};
    var_hp = pyr_var_H ;
end

function video = reconstruct_video ( pyr , pyr_vec , mask )
    % Lifting synthesis function : returns the video frame
    function img = lift_synthesis (L , H , vec , mask )
        sz = size ( L ) ;
        if numel ( sz ) <3; sz (3) = 1; end

        img = zeros ( sz (1) , sz (2) , 2* sz (3) ) ;
        gray = 0.5 * ones ( size ( mask ) ) ;
        
        % Synthesis
        for i = 1: size (L ,3)
            i2 = 2*( i -1) +1;
            % Update
            upd = paste_mask ( H (: ,: , i ) , vec (: , i2 ) , H (: ,: , i ) , vec (: , i2 +1), mask ) ;
            L (: ,: , i ) = L (: ,: , i ) - 0 .5 * upd ;
            % Predict
            pred = paste_mask ( L (: ,: , i ) , vec (: , i2 ) , gray ,[1 1] , mask ) ;
            pred = paste_mask ( pred , vec (: , i2 +1) , L (: ,: , i ) , vec (: , i2 ) , mask ) ;
            H (: ,: , i ) = H (: ,: , i ) + pred ;
        end

    img (: ,: ,1:2: end ) = L ;
    img (: ,: ,2:2: end ) = H ;
    end

    % Iterate over the lifting pyramind hierarchically and perform '
    synthesis
    L = pyr {end };
    for li = ( length ( pyr ) -1) : -1:1
        H = pyr { li };
        L = lift_synthesis (L , H , pyr_vec { li } , mask ) ;
    end
    video = L ;
end