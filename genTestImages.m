mat = 255*ones(250,250,'uint8');
IM2 = padarray(mat,[5 5],0,'both');
colorCode = colormap*255;
IM2 = insertShape(IM2,"filled-circle",[130 130 90],ShapeColor="cyan");
IM2 = insertShape(IM2,"circle",[130 130 90],LineWidth=3,ShapeColor="black");
IM2 = insertShape(IM2,"filled-rectangle",[117 60 25 140],LineWidth=3,ShapeColor="white");
IM2 = insertShape(IM2,"rectangle",[117 60 25 140],LineWidth=3,ShapeColor="black");
imshow(IM2);