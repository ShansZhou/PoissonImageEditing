function result = localIlluminationChanges(resultImg, targetImg,targetCropping_col, targetCropping_row, offsets )
img_source = resultImg;
img_target = targetImg;

%% freedom cropping
% figure;title('target Image'),imshow(img_target/255);
% [mask_target,target_col,target_row]= roipoly(img_target/255);
% figure;title('point a position to paste target image'),imshow(img_source/255);
% [offset_col, offset_row] = ginput(1);

%% fixed cropping 
target_col = targetCropping_col;
target_row = targetCropping_row;
offset_col =offsets(1);
offset_row =offsets(2);


%% cropping mask for image in terms of the input 
offset = [abs(target_col(1)-offset_col), abs(target_row(1) -offset_row)];
target_col = int16(target_col); 
target_row = int16(target_row); 
source_col = (target_col)+int16(offset(1));
source_row = (target_row)+int16(offset(2));

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
[maskTarget_row, maskTarget_col] = find(mask_target);
f_destination = zeros(size(imgSingleChannel_target));
for n =1:size(maskTarget_row)
    neighbour1 = fDestinationBorder_value(maskTarget_row(n)-1, maskTarget_col(n));
    neighbour2 = fDestinationBorder_value(maskTarget_row(n)+1, maskTarget_col(n));
    neighbour3 = fDestinationBorder_value(maskTarget_row(n), maskTarget_col(n)-1);
    neighbour4 = fDestinationBorder_value(maskTarget_row(n), maskTarget_col(n)+1);
    f_destination(maskTarget_row(n), maskTarget_col(n)) = neighbour1 + neighbour2 + neighbour3 + neighbour4;
end
for n =1:size(indexOfBorder_target,1);
    f_destination(indexOfBorder_target(n,1), indexOfBorder_target(n,2)) = 0;
end

vDiv = zeros(size(imgSingleChannel_target));
[maskSource_row, maskSource_col] = find(mask_source);
for n =1:size(maskTarget_row)
%     certain_target =  maskTarget_value(maskTarget_row(n), maskTarget_col(n));
%     neighbour1_target = certain_target - maskTarget_value(maskTarget_row(n)-1, maskTarget_col(n));
%     neighbour2_target = certain_target - maskTarget_value(maskTarget_row(n)+1, maskTarget_col(n));
%     neighbour3_target = certain_target - maskTarget_value(maskTarget_row(n), maskTarget_col(n)-1);
%     neighbour4_target = certain_target - maskTarget_value(maskTarget_row(n), maskTarget_col(n)+1);
    
    certain_source = maskSource_value(maskSource_row(n),maskSource_col(n));
    neighbour1_source = certain_source - maskSource_value(maskSource_row(n)-1,maskSource_col(n));
    neighbour2_source = certain_source - maskSource_value(maskSource_row(n)+1,maskSource_col(n));
    neighbour3_source = certain_source - maskSource_value(maskSource_row(n),maskSource_col(n)-1);
    neighbour4_source = certain_source - maskSource_value(maskSource_row(n),maskSource_col(n)+1);
    %% v = a^b|DelF|.^(-b)DelF <- local illumination changes
    neighbour1_source = (0.2^0.2)*(abs(neighbour1_source))^(-0.2)* neighbour1_source;
    neighbour2_source = (0.2^0.2)*(abs(neighbour2_source))^(-0.2)* neighbour2_source;
    neighbour3_source = (0.2^0.2)*(abs(neighbour3_source))^(-0.2)* neighbour3_source;
    neighbour4_source = (0.2^0.2)*(abs(neighbour4_source))^(-0.2)* neighbour4_source;
    sum = neighbour1_source + neighbour2_source + neighbour3_source + neighbour4_source;
    
    vDiv(maskTarget_row(n),maskTarget_col(n)) = sum;                                                                                   
end
for n =1:size(indexOfBorder_target,1)
    vDiv(indexOfBorder_target(n,1),indexOfBorder_target(n,2)) =0;
end

fq = vDiv + f_destination;
f = A \ fq(indexOfSmallTargetMask);
result_div = imgSingleChannel_source;
[smallTargetMask_row,smallTargetMask_col] = find(small_targetMask);
for n =1:size(smallTargetMask_row);
    result_div(int16(smallTargetMask_row(n)+offset(2)),int16(smallTargetMask_col(n)+offset(1))) = f(n);
end
result_LocalIllumination(:,:,channel) = result_div;
end


result = result_LocalIllumination;
figure,imshow(result/255);title('Local Illumination Changes');

end
