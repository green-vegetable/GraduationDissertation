modelPoints = ReadObj();
projectionMatrix = GetProjectionMatrix();
number = size(modelPoints, 2);
modelPoints = [modelPoints; ones(1, number)];
picturePoints = projectionMatrix * modelPoints;
% for i = 1:number
%    picturePoints(1,i) = picturePoints(1,i) ./ picturePoints(3,i);
%    picturePoints(2,i) = picturePoints(2,i) ./ picturePoints(3,i);
%    picturePoints(3,i) = picturePoints(3,i) ./ picturePoints(3,i);
% end
X = picturePoints(1, :);
Y = picturePoints(2, :);
scatter(Y, X, 1, 'filled')
axis equal;