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

%% Basis Vector
figure(2)
counter = 1;
freq={[2,2];[4,5];[15,10];[20,20]};
for i = 1:4
    subplot(2,2,counter)    
    ind_basis = zeros(m,n);
    pos = cell2mat(freq(i));
    ind_basis(pos(1),pos(2)) = 1;
    basis_t = Bt2.*ind_basis;
    basis = ifft2(basis_t);
    surf(real(basis),'LineStyle','none')
    colorbar
    colormap parula
    counter = counter + 1;
    title_text = sprintf('Freq = (%d,%d)',pos(1),pos(2));
    title(title_text,'FontSize',10)
end
sgtitle('Basis Vectors with different frequencies','FontSize',16)