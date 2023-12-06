function main()
    %[IM1, IM2] = genTestImages();
    IM1 = phantom('Modified Shepp-Logan', 256);
    IM1_C = cartesian(IM1, 128, 64);
end

% Keep just in case we need for reconstructing image after sampling
function m = magnitude(F)
    mag = sqrt(real(F).^2 + imag(F).^2);
    phase = atan2(imag(F),real(F));
    re = mag .* cos(phase);
    im = mag .* sin(phase);

    m = re + i*im;
end

function [i1, i2] = genTestImages()
    i1 = genPhantom1();
    i2 = genPhantom2();
end

function p1 = genPhantom1() 
    canvasSize = [300, 300];
    mat = zeros(canvasSize,'uint8');
    
    % Calculating the center of the canvas
    centerPoint = [canvasSize(1)/2, canvasSize(2)/2];
    xc = centerPoint(1);
    yc = centerPoint(2);

    % Set radius of the background circle as 40% of the canvas's width
    radius = ceil(0.4 * canvasSize(1));

    %dimensions of rectangle (NOTE: change later for user entry)
    width = 40;
    height = 150;
    
    %coordinates of upper left corner of rectangle
    xr = xc - (width/2);
    yr = yc - (height/2);

    % Draw the background circle on the canvas with a specified color
    IM = insertShape(mat,"filled-circle",[xc yc radius],Color='white');
    IM2 = insertShape(IM,"filled-rectangle",[xr yr width height],Color='white');
    IM3 = rgb2gray(IM2);
    imwrite(IM3, "phantom1.jpg");
    p1 = IM3;
end

function p2 = genPhantom2()
    canvasSize = [300, 300];
    mat = zeros(canvasSize,'uint8');
    
    % Calculating the center of the canvas
    centerPoint = [canvasSize(1)/2, canvasSize(2)/2];
    xc = centerPoint(1);
    yc = centerPoint(2);

    % Set radius of the background circle as 40% of the canvas's width
    radius = ceil(0.4 * canvasSize(1));

    % Draw the background circle on the canvas with a specified color
    IM = insertShape(mat,"filled-circle",[xc yc radius],Color='white');
    
    % x-coordinates and radii for foreground circle series
    innerCircleDistance = ceil(radius/5.4);
    
    c1_x = centerPoint(2)- floor(4.5 * innerCircleDistance);
    c1_r = ceil( min(canvasSize)/55);

    c2_x = centerPoint(2)- 3 * innerCircleDistance;
    c2_r = floor(1.52 * c1_r);

    c3_x = centerPoint(2)- 1 * innerCircleDistance;
    c3_r = floor(1.52 * c2_r);

    c4_x = centerPoint(2) + innerCircleDistance;
    c4_r = floor(1.52 * c3_r);

    c5_x = centerPoint(2) + floor(3.5 * innerCircleDistance);
    c5_r = floor(1.52 * c4_r);

    % Draw series of circles of increasing size on top of background circle
    IM2 = insertShape(IM,"filled-circle",[c1_x yc c1_r],Color='white');
    IM3 = insertShape(IM2,"filled-circle",[c2_x yc c2_r],Color='white');
    IM4 = insertShape(IM3,"filled-circle",[c3_x yc c3_r],Color='white');
    IM5 = insertShape(IM4,"filled-circle",[c4_x yc c4_r],Color='white');
    IM6 = insertShape(IM5,"filled-circle",[c5_x yc c5_r],Color='white');

    IM7 = rgb2gray(IM6);
    imwrite(IM7, 'phantom2.jpg');
    p2 = IM7;
end

function result = cartesian(IM, slines, spoints)
    size_img = length(IM); % size of original phantom
    % ratios of phantom dimensions to sampling grid dimensions
    k = [size_img/slines, size_img/spoints]; 
    % dimensions of new image based on ratios
    new_dim = floor(size_img*k);
    % new padded image
    IM2 = zeros(new_dim(1), new_dim(2), 'uint8');
    IM2(1:size_img, 1:size_img) = IM;
    % fourier transformed image
    FT = fftshift(fft2(IM2));
    F_mag = sqrt(real(FT).^2 + imag(FT).^2);
    F_log = log(1 + FT);
    %imshow(F_log, []);
    % sample k-space
    %smpl = interp2(FT, ,);
    result = IM2;
end