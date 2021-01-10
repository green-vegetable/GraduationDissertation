%% 功能：提取灰度图像的尺度不变特征（SIFT特征�?
% 输入:
% im - 灰度图像，该图像的灰度�?�在0�?1之间（注意：应首先对输入图像的灰度�?�进行归�?化处理）
% octaves - 金字塔的组数：octaves (默认值为4).
% intervals - 该输入参数决定每组金字塔的层�?(默认值为 2).
% object_mask - 确定图像中尺度不变特征点的搜索区域，如果没有特别指出，则算法将搜索整个图�?
% contrast_threshold - 对比度阈�?(默认值为0.03).
% curvature_threshold - 曲率阈�?�（默认值为10.0�?.
% interactive - 函数运行显示标志，将其设定为1,则显示算法运行时间和过程的相关信息；
% 如果将其设定�?2，则仅显示最终运行结�?(default = 1).
%
% 输出:
% pos - Nx2 矩阵，每�?行包括尺度不变特征点的坐�?(x,y)
% scale - Nx3 矩阵，每�?行包括尺度不变特征点的尺度信�?
% (第一列是尺度不变特征点所在的组，
% 第二列是其所在的�?,
% 第三列是尺度不变特征点的sigma).
% orient - Nx1 向量，每个元素是特征点的主方向，其范围在 [-pi,pi)之间.
% desc - N x 128 矩阵，每�?行包含特征点的特征向�?.
 
%% 步骤0：载入图像后归一化，并扩充成 2xN-1 大小
% 用于彻底搞懂SIFT的MATLAB版本代码的一个测试文�?
clear all;clc;close all;
% im = rgb2gray(imread('nest.png'));
im = imread('01lhc_r_0470.jpg');
im = rgb2gray(im);
 
antialias_sigma = 0.5;
signal = im;
 
[X Y] = meshgrid( 1:0.5:size(signal,2), 1:0.5:size(signal,1) );
signal = interp2( signal, X, Y, '*linear' ); % im被扩充成  2 x N - 1  �?
 
subsample = [0.5]; %  降采样率
%% 步骤1：生成高斯和差分高斯(DOG)金字�?
%  这两个金字塔的数据分别存储在名为
%  gauss_pyr{orient,interval}和DOG_pyr{orient,interval}
%  的元胞数组中。高斯金字塔含有s+3层，差分高斯金字塔含有s+2层�??
 
preblur_sigma = sqrt(sqrt(2)^2 - (2*antialias_sigma)^2);
 
g = gaussian_filter( preblur_sigma );
 
gauss_pyr{1,1} = conv2( g, g, signal, 'same' );
clear signal;  %  完成了初始化工作
 
% 第一组第�?层的sigma
initial_sigma = sqrt( (2*antialias_sigma)^2 + preblur_sigma^2 );
 
octaves = 4;                      % �?4层金字塔  
intervals = 2;                    % 每一层有2张图�?    s+3
% 之所�?+3，是因为：生成DOG时会减少�?个图片，在计算极值的时�??
% 第一个和�?后一个图片不能使用，�?以要+3.
% 为了保证差分高斯金字塔DOG在每�?层至少有�?张图片可用来寻找极�?�点
% 必须让高斯金字塔�?+3张图�?
% �?后会有octaves层，每层intervals+3张图�?
 
% 记录每一层和每一个尺度的sigma
absolute_sigma = zeros(octaves,intervals+3);
absolute_sigma(1,1) = initial_sigma * subsample(1);
 
% 记录产生金字塔的滤波器的尺寸和标准差
filter_size = zeros(octaves,intervals+3);
filter_sigma = zeros(octaves,intervals+3);
 
% 生成高斯和差分高斯金字塔
for octave = 1:octaves
   sigma = initial_sigma;
   g = gaussian_filter(sigma);
   filter_size(octave,1) = length(g);
   filter_sigma(octave,1) = sigma;
   DOG_pyr{octave} = ...
       zeros(size(gauss_pyr{octave,1},1),size(gauss_pyr{octave,1},2),intervals+2);
   for interval = 2:(intervals+3)    
      sigma_f = sqrt(2^(2/intervals) - 1)*sigma;
      g = gaussian_filter( sigma_f );
      sigma = (2^(1/intervals))*sigma;
       
      % 记录sigma的�??
      absolute_sigma(octave,interval) = sigma * subsample(octave);
       
      % 记录滤波器的尺寸和标准差
      filter_size(octave,interval) = length(g);
      filter_sigma(octave,interval) = sigma;
       
      gauss_pyr{octave,interval} = conv2(g,g,gauss_pyr{octave,interval-1}, 'same' );   
      DOG_pyr{octave}(:,:,interval-1) = ...
          gauss_pyr{octave,interval} - gauss_pyr{octave,interval-1};                
   end     
   if octave < octaves 
      sz = size(gauss_pyr{octave,intervals+1});
      [X Y] = meshgrid( 1:2:sz(2), 1:2:sz(1) );
      gauss_pyr{octave+1,1} = interp2(gauss_pyr{octave,intervals+1},X,Y,'*nearest');
      absolute_sigma(octave+1,1) = absolute_sigma(octave,intervals+1);
      subsample = [subsample subsample(end)*2];
   end     
end
%% 步骤2：查找差分高斯金字塔中的�?部极值，并�?�过曲率和照度进行检�?
contrast_threshold = 0.02;
curvature_threshold = 10.0;
 
curvature_threshold = ((curvature_threshold + 1)^2)/curvature_threshold;
 
% 二阶微分�?
xx = [ 1 -2  1 ];
yy = xx';
xy = [1 0 -1;
      0 0  0;
      -1 0 1]/4;
 
raw_keypoints = [];
contrast_keypoints = [];
curve_keypoints = [];
 
object_mask = ones(size(im));
 
% 在高斯金字塔中查找局部极�?
loc = cell(size(DOG_pyr));
for octave = 1:octaves
   for interval = 2:(intervals+1)
      keypoint_count = 0;
      contrast_mask = abs(DOG_pyr{octave}(:,:,interval)) >= contrast_threshold;
      loc{octave,interval} = zeros(size(DOG_pyr{octave}(:,:,interval)));
         
      edge = ceil(filter_size(octave,interval)/2);
      
      for y=(1+edge):(size(DOG_pyr{octave}(:,:,interval),1)-edge)        
         for x=(1+edge):(size(DOG_pyr{octave}(:,:,interval),2)-edge)
            if object_mask(round(y*subsample(octave)),round(x*subsample(octave))) == 1   
                interactive = 1;
                if( (interactive >= 2) | (contrast_mask(y,x) == 1) )  % 注意�?           
                  % 通过空间核尺度检测最大�?�和�?小�??
                  tmp = DOG_pyr{octave}((y-1):(y+1),(x-1):(x+1),(interval-1):(interval+1)); 
                  pt_val = tmp(2,2,2);
                  if( (pt_val == min(tmp(:))) | (pt_val == max(tmp(:))) )
                     % 存储对灰度大于对比度阈�?�的点的坐标
                     raw_keypoints = [raw_keypoints; x*subsample(octave) y*subsample(octave)];
                     if abs(DOG_pyr{octave}(y,x,interval)) >= contrast_threshold
                        contrast_keypoints = [contrast_keypoints; raw_keypoints(end,:)];    
                        % 计算�?部极值的Hessian矩阵
                        Dxx = sum(DOG_pyr{octave}(y,x-1:x+1,interval) .* xx);
                        Dyy = sum(DOG_pyr{octave}(y-1:y+1,x,interval) .* yy);
                        Dxy = sum(sum(DOG_pyr{octave}(y-1:y+1,x-1:x+1,interval) .* xy));
                        % 计算Hessian矩阵的直迹和行列�?.
                        Tr_H = Dxx + Dyy;
                        Det_H = Dxx*Dyy - Dxy^2;
                        % 计算主曲�?.
                        curvature_ratio = (Tr_H^2)/Det_H;
                        if ((Det_H >= 0) & (curvature_ratio < curvature_threshold))
                           % 存储主曲率小于阈值的的极值点的坐标（非边缘点�?
                           curve_keypoints = [curve_keypoints; raw_keypoints(end,:)];
                           % 将该点的位置的坐标设�?1，并计算点的数量.
                           loc{octave,interval}(y,x) = 1;
                           keypoint_count = keypoint_count + 1;
                        end
                     end                 
                  end
               end              
            end
         end        
      end
 
         fprintf( 2, '%d keypoints found on interval %d\n', keypoint_count, interval );
 
   end
end
clear raw_keypoints contrast_keypoints curve_keypoints;
%% 步骤3：计算特征点的主方向.
% 在特征点的一个区域内计算其梯度直方图
g = gaussian_filter( 1.5 * absolute_sigma(1,intervals+3) / subsample(1) );
zero_pad = ceil( length(g) / 2 );
 
% 计算高斯金字塔图像的梯度方向和幅�?
mag_thresh = zeros(size(gauss_pyr));
mag_pyr = cell(size(gauss_pyr));
grad_pyr = cell(size(gauss_pyr));
 
for octave = 1:octaves
   for interval = 2:(intervals+1)     
       
      % 计算x,y的差�?
      diff_x = 0.5*(gauss_pyr{octave,interval}(2:(end-1),3:(end))-gauss_pyr{octave,interval}(2:(end-1),1:(end-2)));
      diff_y = 0.5*(gauss_pyr{octave,interval}(3:(end),2:(end-1))-gauss_pyr{octave,interval}(1:(end-2),2:(end-1)));
       
      % 计算梯度幅�??
      mag = zeros(size(gauss_pyr{octave,interval}));     
      mag(2:(end-1),2:(end-1)) = sqrt( diff_x .^ 2 + diff_y .^ 2 );
       
      % 存储高斯金字塔梯度幅�?
      mag_pyr{octave,interval} = zeros(size(mag)+2*zero_pad);
      mag_pyr{octave,interval}((zero_pad+1):(end-zero_pad),(zero_pad+1):(end-zero_pad)) = mag;     
       
      % 计算梯度主方�?
      grad = zeros(size(gauss_pyr{octave,interval}));
      grad(2:(end-1),2:(end-1)) = atan2( diff_y, diff_x );
      grad(find(grad == pi)) = -pi;
       
      % 存储高斯金字塔梯度主方向
      grad_pyr{octave,interval} = zeros(size(grad)+2*zero_pad);
      grad_pyr{octave,interval}((zero_pad+1):(end-zero_pad),(zero_pad+1):(end-zero_pad)) = grad;
   end
end
 
clear mag grad
%% 步骤4：确定特征点的主方向
% 方法：�?�过寻找每个关键点的子区域内梯度直方图的峰�?�（注：每个关键点的主方向可以有不止�?个）
 
% 将灰度直方图分为36等分，每�?10度一�?
num_bins = 36;
hist_step = 2*pi/num_bins;  % 步进�?
hist_orient = [-pi:hist_step:(pi-hist_step)];  % 步进方向
 
% 初始化关键点的位置�?�方向和尺度信息
pos = [];
orient = [];
scale = [];
 
% 给关键点确定主方�?
for octave = 1:octaves
   for interval = 2:(intervals + 1)          
      keypoint_count = 0;
       
      % 构�?�高斯加权掩�?
      g = gaussian_filter( 1.5 * absolute_sigma(octave,interval)/subsample(octave) );
      hf_sz = floor(length(g)/2);
      g = g'*g;     
       
      loc_pad = zeros(size(loc{octave,interval})+2*zero_pad);
      loc_pad((zero_pad+1):(end-zero_pad),(zero_pad+1):(end-zero_pad)) = loc{octave,interval};
       
      [iy ix]=find(loc_pad==1);
      for k = 1:length(iy)
 
         x = ix(k);
         y = iy(k);
         wght = g.*mag_pyr{octave,interval}((y-hf_sz):(y+hf_sz),(x-hf_sz):(x+hf_sz));
         grad_window = grad_pyr{octave,interval}((y-hf_sz):(y+hf_sz),(x-hf_sz):(x+hf_sz));
         orient_hist=zeros(length(hist_orient),1);
         for bin=1:length(hist_orient)    
            diff = mod( grad_window - hist_orient(bin) + pi, 2*pi ) - pi;
            orient_hist(bin) = ...
                orient_hist(bin)+sum(sum(wght.*max(1 - abs(diff)/hist_step,0)));
         end
          
         % 运用非极大抑制法查找主方向直方图的峰�?
         peaks = orient_hist;       
         rot_right = [ peaks(end); peaks(1:end-1) ];
         rot_left = [ peaks(2:end); peaks(1) ];        
         peaks( find(peaks < rot_right) ) = 0;
         peaks( find(peaks < rot_left) ) = 0;
          
         % 提取�?大峰值的值和其索引位�?
         [max_peak_val ipeak] = max(peaks);
          
         % 将大于等于最大峰�?80% 的直方图的也确定为特征点的主方向
         peak_val = max_peak_val;
         while( peak_val > 0.8*max_peak_val )
     
            % �?高峰值最近的三个柱�?��?�过抛物线插值精确得�?
            A = [];
            b = [];
            for j = -1:1
               A = [A; (hist_orient(ipeak)+hist_step*j).^2 (hist_orient(ipeak)+hist_step*j) 1];
                bin = mod( ipeak + j + num_bins - 1, num_bins ) + 1;
               b = [b; orient_hist(bin)];
            end
            c = pinv(A)*b;
            max_orient = -c(2)/(2*c(1));
            while( max_orient < -pi )
               max_orient = max_orient + 2*pi;
            end
            while( max_orient >= pi )
               max_orient = max_orient - 2*pi;
            end           
             
            % 存储关键点的位置、主方向和尺度信�?
            pos = [pos; [(x-zero_pad) (y-zero_pad)]*subsample(octave) ];
            orient = [orient; max_orient];
            scale = [scale; octave interval absolute_sigma(octave,interval)];
            keypoint_count = keypoint_count + 1;  
 
            peaks(ipeak) = 0;
            [peak_val ipeak] = max(peaks);
         end        
      end
         fprintf( 2, '(%d keypoints)\n', keypoint_count );         
   end
end
clear loc loc_pad
%% 步骤5：特征向量生�?
orient_bin_spacing = pi/4;
orient_angles = [-pi:orient_bin_spacing:(pi-orient_bin_spacing)];
 
grid_spacing = 4;
[x_coords y_coords] = meshgrid( [-6:grid_spacing:6] );
feat_grid = [x_coords(:) y_coords(:)]';
[x_coords y_coords] = meshgrid( [-(2*grid_spacing-0.5):(2*grid_spacing-0.5)] );
feat_samples = [x_coords(:) y_coords(:)]';
feat_window = 2*grid_spacing;
 
desc = [];
 
for k = 1:size(pos,1)
   x = pos(k,1)/subsample(scale(k,1));
   y = pos(k,2)/subsample(scale(k,1));  
    
   % 将坐标轴旋转为关键点的方向，以确保旋转不变�??
   M = [cos(orient(k)) -sin(orient(k)); sin(orient(k)) cos(orient(k))];
   feat_rot_grid = M*feat_grid + repmat([x; y],1,size(feat_grid,2));
   feat_rot_samples = M*feat_samples + repmat([x; y],1,size(feat_samples,2));
    
   % 初始化特征向�?.
   feat_desc = zeros(1,128);
    
   for s = 1:size(feat_rot_samples,2)
      x_sample = feat_rot_samples(1,s);
      y_sample = feat_rot_samples(2,s);
       
      % 在采样位置进行梯度插�?
      [X Y] = meshgrid( (x_sample-1):(x_sample+1), (y_sample-1):(y_sample+1) );
      G = interp2( gauss_pyr{scale(k,1),scale(k,2)}, X, Y, '*linear' );  % 耗时太长
      G(find(isnan(G))) = 0;
      diff_x = 0.5*(G(2,3) - G(2,1));
      diff_y = 0.5*(G(3,2) - G(1,2));
      mag_sample = sqrt( diff_x^2 + diff_y^2 );
      grad_sample = atan2( diff_y, diff_x );
      if grad_sample == pi
         grad_sample = -pi;
      end     
       
      % 计算x、y方向上的权重
      x_wght = max(1 - (abs(feat_rot_grid(1,:) - x_sample)/grid_spacing), 0);
      y_wght = max(1 - (abs(feat_rot_grid(2,:) - y_sample)/grid_spacing), 0);
      pos_wght = reshape(repmat(x_wght.*y_wght,8,1),1,128);
       
      diff = mod( grad_sample - orient(k) - orient_angles + pi, 2*pi ) - pi;
      orient_wght = max(1 - abs(diff)/orient_bin_spacing,0);
      orient_wght = repmat(orient_wght,1,16);        
       
      % 计算高斯权重
      g = exp(-((x_sample-x)^2+(y_sample-y)^2)/(2*feat_window^2))/(2*pi*feat_window^2);
       
      feat_desc = feat_desc + pos_wght.*orient_wght*g*mag_sample;
   end
    
   % 将特征向量的长度归一化，则可以进�?步去除光照变化的影响.
   feat_desc = feat_desc / norm(feat_desc);
    
   feat_desc( find(feat_desc > 0.2) ) = 0.2;
   feat_desc = feat_desc / norm(feat_desc);
    
   % 存储特征向量.
   desc = [desc; feat_desc];
   tmp = mod(k,25);
   if ( tmp == 0 )
      fprintf( 2, '.' );
   end
end
 
% 调整采样偏差
sample_offset = -(subsample - 1);
for k = 1:size(pos,1)
   pos(k,:) = pos(k,:) + sample_offset(scale(k,1));
end
 
if size(pos,1) > 0
    scale = scale(:,3);
end
 
hh = display_keypoints( pos, scale, orient);