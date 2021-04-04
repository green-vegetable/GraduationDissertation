function [matrix] = Normalize(matrix)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
number = length(matrix);
for i = 1:number
   matrix(1,i) = matrix(1,i) ./ matrix(3,i);
   matrix(2,i) = matrix(2,i) ./ matrix(3,i);
   matrix(3,i) = matrix(3,i) ./ matrix(3,i);
end


end

