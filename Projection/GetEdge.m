function y = GetEdge(I)
%UNTITLED4 此处显示有关此函数的摘要
%   此处显示详细说明
BW1 = edge(I,'sobel');
BW2 = edge(I,'canny');
figure;
imshowpair(BW1,BW2,'montage')
end

