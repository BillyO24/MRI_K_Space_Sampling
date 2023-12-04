function p1 = genPhantom1(N, r, w, h) 
    mat = zeros(250,250,'uint8');
    IM2 = insertShape(mat,"filled-circle",[130 130 90],Opacity=0.5);
    IM2 = insertShape(IM2,"filled-rectangle",[117 60 25 140],Opacity=1);
    IM2 = rgb2gray(IM2);
    imwrite(IM2, 'phantom1.jpg');
    imshow(IM2);
    p1 = IM2; 
end