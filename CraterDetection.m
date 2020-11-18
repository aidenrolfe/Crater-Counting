clear all;
close all;

%% parameters used for image processing
SEval = 8;
SE2val = 7;

SE = strel('disk',SEval); % structuring element for opening/closing
SE2 = strel('disk',SE2val); % structuring element for morph gradient
range = [500 100000]; % range of crater sizes to find
thresh = 0.29; % threshold for binarized image

%% displaying source image
file = uigetfile('*.*'); % allows user to select image in program folder
moon = imread(file); % read image selected

figure('Position',[100,100,900,600])

[height,width] = size(moon); % dimensions of image

% pixel density of image
pixeldens = 13.10; % in meters/pixel


area = (pixeldens*height)*(pixeldens*width); % area in m^2
area = area/1e6; % m^2 to km^2

subplot(2,4,1)
imshow(moon)
title('Original image')

moon = double(moon);

%% morphological opening and closing on source image
% RGB -> greyscale
moonGrey = sqrt(moon(:,:,1).^2 + moon(:,:,2).^2 + moon(:,:,3).^2);
moonGrey = moonGrey./max(moonGrey(:));

moonOpen = imopen(moonGrey,SE); % morphological opening

subplot(2,4,2)
imshow(moonOpen)
title(['Morphological opening with structuring element ',num2str(SEval)])

moonClose = imclose(moonOpen,SE); % morphological closing

subplot(2,4,3)
imshow(moonClose)
title(['Morphological smoothing with structuring element ',num2str(SEval)])

%% morphological gradient on smoothed image
moonGradient = imdilate(moonClose,SE2) - imerode(moonClose,SE2); % morphological gradient

subplot(2,4,4)
imshow(moonGradient)
title(['Morphological gradient with structuring element ',num2str(SE2val)])

%% binarizing of the gradient image
moonB = imbinarize(moonGradient,thresh); % binarize gradient image

subplot(2,4,5)
imshow(moonB)
title(['Binarized image with threshold ',num2str(thresh)])

%% ellipse fitting on binarized image
% extract area of crater
moonB = imclearborder(moonB);
moonB = bwareafilt(moonB,range);

% calculates properties of ellipse
p = regionprops('table',moonB,{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});
n = numel(p(:,1)); % no. of objects picked up in range
subplot(2,4,6)
imshow(moonB)
title(['Ellipse fitting: ',num2str(n),' craters found in area ',num2str(round(area)),...
    ' km^2 within the range ',num2str(range(1)),' to ',num2str(range(2))])
hold on

theta = linspace(0,2*pi);
col = (p.MajorAxisLength/2).*cos(theta);
row = (p.MinorAxisLength/2).*sin(theta);

for i = 1:n
% calculates the ellipse line
M = makehgtform('translate',[p.Centroid(i,:), 0],'zrotate',deg2rad(-1*p.Orientation(i,:)));
N = M*[col(i,:);row(i,:);zeros(1,numel(row(i,:)));ones(1,numel(row(i,:)))];

% plots ellipses onto binary image
subplot(2,4,6)
plot(N(1,:),N(2,:),'r','LineWidth',2)
hold on

end

%% plotting production function for the moon

T = linspace(0,5,100000);
N = 5.44e-14*(exp(6.93*T)-1) + 8.38e-4*T; % N = cumulative crater count

subplot(2,4,8)
h = plot(T,log10(N));
xlabel('Time (Gyr)')
ylabel('Relative crater count log(N)')
title('Production function for the moon');

%% solve production function for my cumulative crater count N1
N1 = n/area; % craters per km^2

syms X
% interpolate age T1 given the relative cratering rate N1
T1 = vpasolve(N1 == 5.44e-14*(exp(6.93*X)-1) + 8.38e-4*X,X);
T1 = round(double(T1),4);

hold on;
plot(T1,log10(N1),'o');
hold on;
text(T1+0.15,log10(N1),[num2str(T1),' Gyrs']);

%% plotting previous ellipses on the original image
subplot(2,4,7)

moon = imread(file);
imshow(moon)

hold on;

for i = 1:n
    % calculates the ellipse line
    M = makehgtform('translate',[p.Centroid(i,:), 0],'zrotate',deg2rad(-1*p.Orientation(i,:)));
    N = M*[col(i,:);row(i,:);zeros(1,numel(row(i,:)));ones(1,numel(row(i,:)))];
    
    % plots ellipses onto raw image
    subplot(2,4,7)
    plot(N(1,:),N(2,:),'r','LineWidth',2)
    hold on
end

title('Ellipse fitting on original image')

