%% ÕÒsiftÌØÕ÷

[image, descrips, locs] = sift('01lhc_l__Color1200.jpg'); 
showkeys(image, locs);

%% Æ¥ÅäsiftÌØÕ÷
match('01lhc_l__Color1200.jpg','01lhc_m__Color1013.jpg');


%%
i1=imread('01lhc_l__Color1200.jpg');
i2=imread('01lhc_m__Color1013.jpg');
i11=rgb2gray(i1);
i22=rgb2gray(i2);
imwrite(i11,'v1.jpg','quality',80);
imwrite(i22,'v2.jpg','quality',80);
match('v1.jpg','v2.jpg');