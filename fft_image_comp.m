clear;
clc;
A = imread('panda.jpg');
m = size(A,1);
n = size(A,2);
B = rgb2gray(A);
figure(1)
subplot(1,2,1)
imagesc(A);

%% FFT
Bt = fft2(A);   % B is grayscale image from above
Bt2 = fft2(B);
Blow2 = log(abs(Bt2)+1); %put fft in log scale
subplot(1,2,2)
imagesc(Blow2)
colorbar
% set(gcf,'Position',[1500 100 size(A,2) size(A,1)]);
sgtitle('Original Image and Fourier Coefficients','FontSize',16)

%% Inv FFT
% Zero out all small coefficients and inverse transform (High Pass Filter)

% Btsort = sort(abs(Bt(:))); %sort by magnitude
Bt2sort = sort(abs(Bt2(:))); %sort by magnitude
counter = 1;
Scale = 2;
for keep=[.99 0.05 0.01 0.002]
    fig = figure(counter+1);
    subplot(1,2,1)
%     thresh = Btsort(floor((1-keep)*length(Btsort)));
    thresh2 = Bt2sort(floor((1-keep)*length(Bt2sort)));
%     ind = abs(Bt)>thresh; %Find small indices
    ind2 = abs(Bt2)>thresh2; %Find small indices
%     Atlow = Bt.*ind*Scale; %Threshold small indices
    Atlow = Bt.*(reshape([ind2,ind2,ind2],size(A)))*Scale;
    Atlow2 = Bt2.*ind2*Scale; %Threshold small indices
    Alow = uint8(ifft2(Atlow)); %Compressed image using inverse fft
    imagesc(Alow); %Plot reconstruction    
    subplot(1,2,2)
    imagesc(log(abs(Atlow2)+1));
    colorbar
%     imagesc(rgb2gray(log(abs(fftshift((Atlow)))+1)));
%     colorbar
    counter = counter + 1;
    sgtitle(['Keeping ',num2str(keep*100),'% of Larger Fourier Coefficients'],'FontSize',16)
end

%% Inv FFT
% Zero out all larger coefficients and inverse transform (Low pass filter)

% Btsort = sort(abs(Bt(:))); %sort by magnitude
Bt2sort = sort(abs(Bt2(:))); %sort by magnitude
counter = 1;
Scale = 4;
for keep=[.9999 0.99 0.95 0.9]
    fig = figure(counter+1);
    subplot(1,2,1)
%     thresh = Btsort(floor((1-keep)*length(Btsort)));
    thresh2 = Bt2sort(floor(keep*length(Bt2sort)));
%     ind = abs(Bt)<thresh; %Find small indices
    ind2 = abs(Bt2)<thresh2; %Find small indices
%     Atlow = Bt.*ind*Scale; %Threshold small indices
    Atlow = Bt.*(reshape([ind2,ind2,ind2],size(A)))*Scale;
    Atlow2 = Bt2.*ind2*Scale; %Threshold small indices
    Alow = uint8(ifft2(Atlow)); %Compressed image using inverse fft
    imagesc(Alow); %Plot reconstruction    
    subplot(1,2,2)
    imagesc(log(abs(Atlow2)+1));
    colorbar
%     imagesc(rgb2gray(log(abs(fftshift((Atlow)))+1)));
%     colorbar
    counter = counter + 1;
    sgtitle(['Including ',num2str(keep*100),'% of Smaller Fourier Coefficients'],'FontSize',16)
end