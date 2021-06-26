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
% Zero out all middle coefficients and inverse transform
counter = 1;
Scale = 1;
for N=[5 10 50 100]
    fig = figure(counter+1);
    subplot(1,2,1)
    ind = ones(m,n)*Scale;
    ind(round(m/N:((N-1)*m)/N),round(n/N:((N-1)*n)/N)) = 0;
    Atlow = Bt.*(reshape([ind,ind,ind],size(A)));
    Atlow2 = Bt2.*ind;
    Alow = uint8(real((ifft2(Atlow)))); %Compressed image using inverse fft
    imagesc(Alow); %Plot reconstruction    
    subplot(1,2,2)
    imagesc(log(abs(Atlow2)+1));
    colorbar
%     imagesc(rgb2gray(log(abs(fftshift((Atlow)))+1)));
%     colorbar
    counter = counter + 1;
    sgtitle(['Including Higher Fourier Coeff.,',' Cutoff(N)=',num2str(N)],'FontSize',16)
end

%% Inv FFT
% Zero out all outer coefficients and inverse transform
counter = 1;
Scale = 4;
for N=[50 100 350 600]
    fig = figure(counter+1);
    subplot(1,2,1)
    ind = zeros(m,n);
    ind(round(m/N:((N-1)*m)/N),round(n/N:((N-1)*n)/N)) = 1*Scale;
    Atlow = Bt.*(reshape([ind,ind,ind],size(A)));
    Atlow2 = Bt2.*ind;
    Alow = uint8(real((ifft2(Atlow)))); %Compressed image using inverse fft
    imagesc(Alow); %Plot reconstruction    
    subplot(1,2,2)
    imagesc(log(abs(Atlow2)+1));
    colorbar
%     imagesc(rgb2gray(log(abs(fftshift((Atlow)))+1)));
%     colorbar
    counter = counter + 1;
    sgtitle(['Including Lower Fourier Coeff.,',' Cutoff(N)=',num2str(N)],'FontSize',16)
end

%% Inv FFT
% Zero out all middle coefficients and inverse transform
counter = 1;
Scale = 1;
for N=[5 10 50 100]
    fig = figure(counter+1);
    subplot(1,2,1)
    ind = ones(m,n)*Scale;
    ind(round(1:((N-1)*m)/N),round(1:end)) = 0;
    Atlow = Bt.*(reshape([ind,ind,ind],size(A)));
    Atlow2 = Bt2.*ind;
    Alow = uint8(real((ifft2(Atlow)))); %Compressed image using inverse fft
    imagesc(Alow); %Plot reconstruction    
    subplot(1,2,2)
    imagesc(log(abs(Atlow2)+1));
    colorbar
%     imagesc(rgb2gray(log(abs(fftshift((Atlow)))+1)));
%     colorbar
    counter = counter + 1;
    sgtitle(['Including Higher Fourier Coeff.,',' Cutoff(N)=',num2str(N)],'FontSize',16)
end
