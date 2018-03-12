img = double(rgb2gray(imread('image/apple.jpg')));
disp('Click on the figure to select a region, use "Double Click" to finish');
% x1 = [50,50,150,150];
% y1 = [100,250,250,100];

% for the ROI, find the boundary value of the mask
[mask, x, y]= roipoly(img/255);
line(x,y);
mask_Value = img .* mask;
indexOfBorder = cell2mat(bwboundaries(mask_Value));

% crop the border from mask to get a smaller mask
smallMask = mask;
for n =1:size(indexOfBorder,1)
    smallMask(indexOfBorder(n,1), indexOfBorder(n,2)) = 0;
end
smallMask_value = img .* smallMask;
indexOfSmallMask = find(smallMask_value);

% order all the non-zero pixel in the smaller mask
smallMask_order = zeros(size(smallMask_value));
for n =1:size(indexOfSmallMask);
    smallMask_order(indexOfSmallMask(n)) = n;
end

% create a laplacian operator A based on the smaller mask
A = delsq(smallMask_order);

% caculate destination function f* defines over (image - cropped_part)
fDestinationBorder_index = indexOfBorder;
fDestinationBorder_value = zeros(size(mask_Value));
for n =1:size(indexOfBorder,1);
    fDestinationBorder_value(indexOfBorder(n,1),indexOfBorder(n,2)) = img(indexOfBorder(n,1),indexOfBorder(n,2));
end
[mask_row, mask_col] = find(mask);
f_destination = zeros(size(img));
for n =1:size(mask_row)
    neighbour1 = fDestinationBorder_value(mask_row(n)-1, mask_col(n));
    neighbour2 = fDestinationBorder_value(mask_row(n)+1, mask_col(n));
    neighbour3 = fDestinationBorder_value(mask_row(n), mask_col(n)-1);
    neighbour4 = fDestinationBorder_value(mask_row(n), mask_col(n)+1);
    f_destination(mask_row(n), mask_col(n)) = neighbour1 + neighbour2 + neighbour3 + neighbour4;
end
%assignin('base','old',f_destination);
for n =1:size(indexOfBorder,1);
    f_destination(indexOfBorder(n,1), indexOfBorder(n,2)) = 0;
end

% according to the equation del*F= Fq, F would be the reuslt image in ROI 
fq = f_destination(indexOfSmallMask);
f = A \ fq;
[smallMask_row,smallMask_col] = find(smallMask);
result = img;
for n =1:size(smallMask_row);
    result(smallMask_row(n),smallMask_col(n)) = f(n);
end
figure;imshow((result/255));

