img_source = double(imread('image/face2.png'));
img_target = double(imread('image/face4.png'));
img_source = imresize(img_source,1.15);
img_target = imresize(img_target,0.5);
% figure,imshow(img_source/255);
% figure,imshow(img_target/255);

%% freedom cropping
% figure;title('target Image'),imshow(img_target/255);
% [mask_target,target_col,target_row]= roipoly(img_target/255);
% line(target_col,target_row);
% figure;title('point a position to paste target image'),imshow(img_source/255);
% [offset_col, offset_row] = ginput(1);

%% fixed cropping 
target_col = [5;231;221;179;56;25;11;3;5];
target_row = [5;6;121;198;200;162;86;53;5];
offset_col =62;
offset_row =2.109999999999999e+02;


%% cropping mask from image in terms of the input 
offset = [abs(target_col(1)-offset_col), abs(target_row(1) -offset_row)];
target_col = int16(target_col); 
target_row = int16(target_row); 
source_col = (target_col)+int16(offset(1));
source_row = (target_row)+int16(offset(2));

result_importing = img_source;
result_mixing = img_source;
result_LocalIllumination = img_source;
for channel = 1:3
imgSingleChannel_target = img_target(:,:,channel);
imgSingleChannel_source = img_source(:,:,channel);

mask_target = roipoly(imgSingleChannel_target/255,target_col,target_row);
mask_source = roipoly(imgSingleChannel_source/255,source_col,source_row);

maskSource_value = imgSingleChannel_source .* mask_source;
maskTarget_value = imgSingleChannel_target .* mask_target;
indexOfBorder_target = cell2mat(bwboundaries(maskTarget_value));
% crop the border from mask to get a smaller mask
small_targetMask = mask_target;
for n =1:size(indexOfBorder_target,1)
    small_targetMask(indexOfBorder_target(n,1), indexOfBorder_target(n,2)) = 0;
end
smallTargetMask_value = imgSingleChannel_target .* small_targetMask;
indexOfSmallTargetMask = find(smallTargetMask_value);
% order all the non-zero pixel in the smaller mask
smallTargetMask_order = zeros(size(smallTargetMask_value));
for n =1:size(indexOfSmallTargetMask);
    smallTargetMask_order(indexOfSmallTargetMask(n)) = n;
end
% create a laplacian operator A based on the smaller mask
A = delsq(smallTargetMask_order);
%% caculate destination function f* defines over (image - cropped_part)
fDestinationBorder_index = zeros(size(indexOfBorder_target));
fDestinationBorder_index(:,1) = int16(indexOfBorder_target(:,1))+ (offset(2));
fDestinationBorder_index(:,2) = int16(indexOfBorder_target(:,2))+ (offset(1));
fDestinationBorder_value = zeros(size(maskTarget_value));
for n =1:size(fDestinationBorder_index,1);
    fDestinationBorder_value(indexOfBorder_target(n,1),indexOfBorder_target(n,2)) =...
        imgSingleChannel_source(fDestinationBorder_index(n,1),fDestinationBorder_index(n,2));
end

f_destination = zeros(size(imgSingleChannel_target));
Vpq_importing  = zeros(size(imgSingleChannel_target));
Vpq_mixing = zeros(size(imgSingleChannel_target));
VDiv = zeros(size(imgSingleChannel_target));
[maskTarget_row, maskTarget_col] = find(mask_target);
[maskSource_row, maskSource_col] = find(maskSource_value);
for n =1:size(maskTarget_row)
    %
    neighbour1_f = fDestinationBorder_value(maskTarget_row(n)-1, maskTarget_col(n));
    neighbour2_f = fDestinationBorder_value(maskTarget_row(n)+1, maskTarget_col(n));
    neighbour3_f = fDestinationBorder_value(maskTarget_row(n), maskTarget_col(n)-1);
    neighbour4_f = fDestinationBorder_value(maskTarget_row(n), maskTarget_col(n)+1);
    f_destination(maskTarget_row(n), maskTarget_col(n)) = neighbour1_f + neighbour2_f + neighbour3_f + neighbour4_f;
    %
    certain_target =  maskTarget_value(maskTarget_row(n), maskTarget_col(n));
    neighbour1_target = certain_target - maskTarget_value(maskTarget_row(n)-1, maskTarget_col(n));
    neighbour2_target = certain_target - maskTarget_value(maskTarget_row(n)+1, maskTarget_col(n));
    neighbour3_target = certain_target - maskTarget_value(maskTarget_row(n), maskTarget_col(n)-1);
    neighbour4_target = certain_target - maskTarget_value(maskTarget_row(n), maskTarget_col(n)+1);
    
    certain_source = maskSource_value(maskSource_row(n),maskSource_col(n));
    neighbour1_source = certain_source - maskSource_value(maskSource_row(n)-1,maskSource_col(n));
    neighbour2_source = certain_source - maskSource_value(maskSource_row(n)+1,maskSource_col(n));
    neighbour3_source = certain_source - maskSource_value(maskSource_row(n),maskSource_col(n)-1);
    neighbour4_source = certain_source - maskSource_value(maskSource_row(n),maskSource_col(n)+1);
    
    if abs(neighbour1_target) < abs(neighbour1_source)
        neighbour1 = neighbour1_source;
    else
        neighbour1 = neighbour1_target;
    end
    if abs(neighbour2_target) < abs(neighbour2_source)
        neighbour2 = neighbour2_source;
    else
        neighbour2 = neighbour2_target;
    end
    if abs(neighbour3_target) < abs(neighbour3_source)
        neighbour3 = neighbour3_source;
    else
        neighbour3 = neighbour3_target;
    end
    if abs(neighbour4_target) < abs(neighbour4_source)
        neighbour4 = neighbour4_source;
    else
        neighbour4 = neighbour4_target;
    end
    
    Vpq_importing(maskTarget_row(n),maskTarget_col(n)) = neighbour1_target + neighbour2_target + neighbour3_target + neighbour4_target;  
    
    Vpq_mixing(maskTarget_row(n),maskTarget_col(n)) = neighbour1 + neighbour2 + neighbour3 + neighbour4;
    
    %
    neighbour1_div = (0.2^0.2)*(abs(neighbour1_source))^(-0.2)* neighbour1_source;
    neighbour2_div = (0.2^0.2)*(abs(neighbour2_source))^(-0.2)* neighbour2_source;
    neighbour3_div = (0.2^0.2)*(abs(neighbour3_source))^(-0.2)* neighbour3_source;
    neighbour4_div = (0.2^0.2)*(abs(neighbour4_source))^(-0.2)* neighbour4_source;
    sum = neighbour1_div + neighbour2_div + neighbour3_div + neighbour4_div;
    VDiv(maskTarget_row(n),maskTarget_col(n)) = sum; 
end
for n =1:size(indexOfBorder_target,1)
    f_destination(indexOfBorder_target(n,1), indexOfBorder_target(n,2)) = 0;
    Vpq_importing(indexOfBorder_target(n,1), indexOfBorder_target(n,2)) = 0;
    Vpq_mixing(indexOfBorder_target(n,1),indexOfBorder_target(n,2)) =0;
    VDiv(indexOfBorder_target(n,1), indexOfBorder_target(n,2)) = 0;
end

fq = f_destination + Vpq_importing;
f_importing = A \ fq(indexOfSmallTargetMask);

fq = f_destination + Vpq_mixing;
f_mixing = A \ fq(indexOfSmallTargetMask);

fq = f_destination + VDiv + Vpq_importing;
f = A \ fq(indexOfSmallTargetMask);

result_importingGradient = imgSingleChannel_source;
result_mixingGradient = imgSingleChannel_source;
result_div = imgSingleChannel_source;

[smallTargetMask_row,smallTargetMask_col] = find(small_targetMask);
for n =1:size(smallTargetMask_row);
    result_importingGradient(int16(smallTargetMask_row(n)+offset(2)),int16(smallTargetMask_col(n)+offset(1))) = f_importing(n);
    result_mixingGradient(int16(smallTargetMask_row(n)+offset(2)),int16(smallTargetMask_col(n)+offset(1))) = f_mixing(n);
    result_div(int16(smallTargetMask_row(n)+offset(2)),int16(smallTargetMask_col(n)+offset(1))) = f(n);
end

result_importing(:,:,channel) = result_importingGradient;
result_mixing(:,:,channel) = result_mixingGradient;
result_LocalIllumination(:,:,channel) = result_div;

end
figure,imshow(result_importing/255);title('Importing Gradients');
figure,imshow(result_mixing/255);title('Mixing Gradients');
figure,imshow(result_LocalIllumination/255);title('Local Illumination Changes');
