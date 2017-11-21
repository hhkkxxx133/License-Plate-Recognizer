close all
clear;
clc;

%% Load image
im = imread('samples/!!HPIM1232 - tooclose.jpg');
im = imresize(im,[600,800]);
[height, width, channel] = size(im);

%% Preprocess image
img_gray = rgb2gray(im);
img_smooth = medfilt2(img_gray);
%img_smooth = wiener2(img_smooth, [5,5]);
% SE1 = strel('line',10,0);
% img_enhanced = imtophat(~img_smooth,SE1);
%img_contrast = histeq(img_gray);
img_binary = img_smooth > 120;
%img_binary = imbinarize(img_smooth,'adaptive','ForegroundPolarity','dark','Sensitivity',0.3);

%img_binary = imbinarize(img_smooth);
img_binary = bwareaopen(~img_binary,10);

% SE1 = strel('line',10,0);
% edgemap = imtophat(img_binary,SE1);


figure()
subplot(2,2,1)
imshow(im);
title('Original Image');
hold on;

% figure(1)
subplot(2,2,2)
imshow(~img_binary)
title('Processed Binary Image');
hold on;

%% Get image edge
h_mask = fspecial('sobel');
img_vedge = imfilter(img_binary, h_mask');
% figure(2)
subplot(2,2,3)
imshow(~img_vedge)
title('Vertical Edge');

%% First stage - Select rows
avg = zeros([1,height]);
var = zeros([1,height]);
edgepixel = zeros([1,height]);
for row = 1:height
    for col = 1:width
        edgepixel(row) = edgepixel(row) + img_vedge(row,col);
    end
    if edgepixel(row) > 20 % threshold
        avg(row) = mean(img_vedge(row,:));
        var(row) = std(img_vedge(row,:));
    end
end

%% Second stage - Potential license plate regions
band = zeros([1,0]);
maxlength = -10000;
top = zeros([1,0]);
bottom = zeros([1,0]);
for row = 1:height
    if var(row) > 0.15
        band(end+1) = row;
        if length(band) == 1
            top(end+1) = row;%set top=starting row
            bottom(end+1) = row;%set bottom=ending row
            continue
        else if length(band) > 1
                if find(band==(row-1))%check continuous
                    new_len=length(band);
                    if new_len > maxlength
                        maxlength = new_len;
                        bottom(end) = band(end);
                    end
                else
                    band = zeros([1,0]);%start a new band
                    band(end+1) = row;
                    top(end+1) = row;
                    bottom(end+1) = row;
                    maxlength = -10000;
                end
            end
        end
    end
end


%% Third stage - Refine license plate regions
temp_top = top;
temp_bottom = bottom;
top = temp_top((temp_bottom-temp_top)>15 & (temp_bottom-temp_top)<height*0.3);
bottom = temp_bottom((temp_bottom-temp_top)>15 & (temp_bottom-temp_top)<height*0.3);

num_band = length(top);
left = zeros([1,num_band]);
right = zeros([1,num_band]);
for i = 1:num_band
    mean_x = zeros([1,0]);
    var_x = zeros([1,0]);
    for row = top(i):bottom(i)
        mean_x(end+1) = mean(find(img_vedge(row,:)>0));
        var_x(end+1) = std(find(img_vedge(row,:)>0));
    end
    min_x = min(mean_x);
    max_x = max(mean_x);
    max_var = max(var_x);
    
    left(i) = uint32(min_x - max_var);%left=u_xmin-v_xmax
    if left(i)<1
        left(i) = 1;
    end
    right(i) = uint32(max_x + max_var);%right=u_xmax-v_xmax
    if right(i)>width
        right(i) = width;
    end
end

% figure(2)
% subplot(1,3,1)
% imshow(im);
% title('Original Image')
% hold on;

subplot(2,2,4)
%imshow(im);
imshow(~img_vedge);
title('Vertical Edge')

% subplot(1,3,2)
% imshow(im)
hold on;
for i=1:num_band
    rectangle('position', [left(i),top(i),right(i)-left(i),bottom(i)-top(i)],'Edgecolor','r');
end
title('Potential License Plate Region');

%% Fourth stage - Final license plate region(s)
for i = 1:num_band
    hband = bottom(i)-top(i);
    for col = left(i):width
        count = length( find(img_vedge(top(i):bottom(i),col)==1) );  %%find prominent edge
        if count>(hband*0.5)
            break;
        end
    end
    left(i) = col;
    
    for col = right(i):-1:1
        count = length( find(img_vedge(top(i):bottom(i),col)==1) );
        if count>(hband*0.5)
            break;
        end
    end
    right(i) = col;
end

plate_height = bottom - top;
plate_width = right - left;
ratio = plate_width./plate_height;

temp_top = top;
temp_bottom = bottom;
temp_left = left;
temp_right = right;
top = temp_top((temp_right-temp_left)>20 & (temp_right-temp_left)<width*0.6 & ratio>2);
bottom = temp_bottom((temp_right-temp_left)>20 & (temp_right-temp_left)<width*0.6 & ratio>2);
left = temp_left((temp_right-temp_left)>20 & (temp_right-temp_left)<width*0.6 & ratio>2);
right = temp_right((temp_right-temp_left)>20 & (temp_right-temp_left)<width*0.6 & ratio>2);

% figure(4)
% imshow(~img_vedge);
figure()
subplot(1,3,3)
imshow(im)
hold on;
for i=1:length(top)
    rectangle('position', [left(i),top(i),right(i)-left(i),bottom(i)-top(i)],'Edgecolor','r');
end
title('Final License Plate Region');

%% Cut the license plate region
count=0;
for i=1:length(top)
    new_im=im(top(i):bottom(i),left(i):right(i));
    count=count+1;
    filename=strcat(int2str(count),'_plate_found.jpg');
    figure;imshow(new_im);
    imwrite(new_im,filename);
end