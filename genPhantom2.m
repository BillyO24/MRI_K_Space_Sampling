function p2 = genPhantom2(N, r1, r2, r3, r4, r5)
    mat = zeros(250,250,'uint8');

    IM2 = insertShape(mat,"filled-circle",[130 130 90],Opacity=0.5);
    IM2 = insertShape(IM2,"filled-circle",[55 130 3],Opacity=1);
    IM2 = insertShape(IM2,"filled-circle",[75 130 6],Opacity=1);
    IM2 = insertShape(IM2,"filled-circle",[100 130 9],Opacity=1);
    IM2 = insertShape(IM2,"filled-circle",[130 130 14],Opacity=1);
    IM2 = insertShape(IM2,"filled-circle",[175 130 20],Opacity=1);

    IM2 = rgb2gray(IM2);
    imwrite(IM2, 'phantom2.jpg');
    imshow(IM2);
    p2 = IM2;
end