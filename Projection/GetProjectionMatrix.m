function [M, PicturePoints, ModelPoints] = GetProjectionMatrix(start, endd)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
slCharacterEncoding('utf-8');
rawPoints = string(fileread('points.txt'));

rawPoints = rawPoints.split(';').split(',');
[sizeX, sizeY] = rawPoints.size;
rawPoints = rawPoints(2:sizeX, 2:sizeY);
points = double(rawPoints);

points(:,[2,3]) = points(:,[3,2]);
points(:,3) = points(:,3).* (-1.0);

modelPoints = points(:, 1:3);
picturePoints = points(:, 4:5);
HE = ones(sizeX-1, 1);


modelPoints = [modelPoints,HE];
picturePoints = [picturePoints, HE];

PicturePoints = picturePoints(start:endd, :)';
ModelPoints = modelPoints(start: endd, :)';

x=points(start:endd,1);
y=points(start:endd,2);
z=points(start:endd,3);

u=points(start:endd,4);
v=points(start:endd,5);
c = [[x,y,z,ones(endd - start + 1,1),zeros(endd - start + 1,4),-u.*x,-u.*y,-u.*z];...
     [zeros(endd - start + 1,4),x,y,z,ones(endd - start + 1,1),-v.*x,-v.*y,-v.*z]];
b = [u; v];

[M,resnorm,residual,exitflag,output,lambda] = lsqlin(c,b);
M = [M; 1];
M = reshape(M,4,3)';
sum = sqrt(M(3,1) * M(3,1) + M(3,2) * M(3,2) + M(3,3) * M(3,3));
M(3,4) = sqrt(1 - sum)
projectedPicturePoints = M * (modelPoints');
projectedPicturePoints = Normalize(projectedPicturePoints);

sub = projectedPicturePoints - picturePoints';

end

