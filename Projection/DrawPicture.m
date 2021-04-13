function y = DrawPicture(filteredpicturePoint,originalImage,projectionMatrix,rawModelPoints,rawPicturePoints)
size(filteredpicturePoint)



imshow(originalImage); % 原始图片
hold on;

X = filteredpicturePoint(1, :);
Y = filteredpicturePoint(2, :);
scatter(X, Y, 1, 'filled','b'); % 模型的投影点

rawProjectedPoints = projectionMatrix * rawModelPoints;
rawProjectedPoints = Normalize(rawProjectedPoints);

scatter(rawPicturePoints(1,:), rawPicturePoints(2,:), 60,'filled','g'); % 图片选取的对应点
scatter(rawProjectedPoints(1,:), rawProjectedPoints(2,:), 60,'filled','r'); % 图片选取的对应点的投影点

axis equal;
end

