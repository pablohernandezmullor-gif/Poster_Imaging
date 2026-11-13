% Challenge1_Mi solución
m=256;
n=256;
noise=1e-6;
[B,Ac,Ar,X] = challenge1(m,n,noise);


%Representamos la imagen, que será una imagen RGB
fig=figure(1);
subplot(3,1,1)
imagesc(X);
colormap gray

subplot(3,1,2)
Xnaive=Ac\B/Ar';
imagesc(Xnaive)
title("Deblur con noise: "+noise)
colormap gray

subplot(3,1,3)
noise=1e-7;
[B,Ac,Ar,X] = challenge1(m,n,noise);
Xnaive=Ac\B/Ar';
imagesc(Xnaive)
title("Deblur con noise: "+noise)


carpeta = fileparts(mfilename('fullpath'));
destino=carpeta+"\Output"+"\Naive_Linear_Deblurring.fig";
savefig(fig,destino)
disp("Imagen guardada en  "+ destino);



%% Guardamos las imagenes por separado
fig=figure(2);
imagesc(X); axis off
colormap gray
destino=carpeta+"\Output"+"\Original_image.fig";
savefig(fig,destino)


fig=figure(2);
noise=1e-4;
[B,Ac,Ar,X] = challenge1(m,n,noise);
Xnaive=Ac\B/Ar'; 
imagesc(B); axis off
destino=carpeta+"\Output"+"\Blur_noise_" +noise+ ".fig";
savefig(fig,destino)


fig=figure(2);
imagesc(Xnaive); axis off
destino=carpeta+"\Output"+"\Deblur_Naive_noise_" +noise+ ".fig";
savefig(fig,destino)



%%%% Guardamos la imagen con blur para analizarla con ImageJ
fig2=figure(3);
imagesc(B), colormap gray

B = single(B);                 % tu imagen borrosa (la matriz, no la figura)
destino_2=carpeta+"\Output"+"\Image_Blurred_noise_" + noise + ".tif";
t = Tiff(destino_2,'w');
t.setTag('ImageLength',size(B,1)); t.setTag('ImageWidth',size(B,2));
t.setTag('Photometric',Tiff.Photometric.MinIsBlack);
t.setTag('BitsPerSample',32); t.setTag('SamplesPerPixel',1);
t.setTag('SampleFormat',Tiff.SampleFormat.IEEEFP);
t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
t.write(B); t.close();

disp("Imagen_2 guardada en  "+ destino_2);



%%% Deconvolution Lab nos pide que subamos la PDF. Como 
% estamos trabajando con el libro Hansen-Nagy, que divide la PSF en dos 
% matrices hay que pasarlo a un mapa 2D. Según el capítulo 4.1.3 de Naggy. 
n = size(Ac,1);  m = size(Ar,1);
c = Ac(:, round(n/2));          % 1D column kernel. In the middle (round(n/2)) it has the structure of column [c1,c2,c3.., cn]. 
%the first column would be a column of [c_{n/2},c_{n/2 +1}..., c_{n},c_1,c_2,...,c_{n-1}]
r = Ar(:, round(m/2));          % 1D row kernel
PSF = c * r';                   % 2D PSF
PSF = PSF / sum(PSF(:));  %We normalise the kernel



%We save the PSF as TIFF
t = Tiff('PSF.tif','w');
t.setTag('ImageLength',size(PSF,1)); t.setTag('ImageWidth',size(PSF,2));
t.setTag('Photometric',Tiff.Photometric.MinIsBlack);
t.setTag('BitsPerSample',32); t.setTag('SamplesPerPixel',1);
t.setTag('SampleFormat',Tiff.SampleFormat.IEEEFP);
t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
t.write(single(PSF)); t.close();

%%% Vemos si la PSF es la que efectivamente usa el eejrcico
%%% Diagnóstico de la PSF
kc = round(size(Ac,1)/2);  kr = round(size(Ar,1)/2);
PSFraw = Ac(:,kc) * Ar(:,kr)';                     % PSF SIN normalizar

% 1) Escala: ¿suman 1 las columnas de Ac y Ar?
fprintf('Suma PSF sin normalizar: %.4f\n', sum(PSFraw(:)));

% 2) ¿Son Toeplitz (convolución invariante)?
noToep = @(A) norm(A - toeplitz(A(:,1),A(1,:)),'fro')/norm(A,'fro');
fprintf('Ac no-Toeplitz: %.2e   Ar no-Toeplitz: %.2e\n', noToep(Ac), noToep(Ar));

% 3) Error con la PSF sin normalizar (descarta el problema de escala)
Pc    = circshift(PSFraw, [1-kc, 1-kr]);
Bconv = real(ifft2(fft2(double(X)) .* fft2(Pc)));
Bexact = Ac*double(X)*Ar';
fprintf('Error relativo: %.2e\n', norm(Bconv-Bexact,'fro')/norm(Bexact,'fro'));

% 4) ¿Dónde está el error?
figure; imagesc(abs(Bconv-Bexact)); colorbar; title('|Bconv - Bexact|')
figure; plot(Ac(:,[40 kc 216])); legend('col 40','col centro','col 216')


%%% Este código hace lo mismo que Deconvolution lab excepto que deconvolution lab tiene opciones como padding 
%La siguiente fórmula es para periodic boundary conditions para filtro de mínimos cuadrados restringido
destino=carpeta+"\Output"+"\LSS_Deblurring";
S  = fft2(Pc);                               % Pc ya es PSFraw centrada en (1,1)
[u,v] = ndgrid(0:size(B,1)-1, 0:size(B,2)-1);
L2 = (4 - 2*cos(2*pi*u/size(B,1)) - 2*cos(2*pi*v/size(B,2))).^2;   % |Laplaciano|^2

for lambda = [0 1e-10 1e-8 1e-6]
    Xr = real(ifft2 ( conj(S).*fft2(double(B)) ./ (abs(S).^2 + lambda*L2) ) );

    fig3=figure; imagesc(Xr); colormap gray; axis image
    title("RIF \lambda = "+lambda+",  err = "+norm(Xr-X,'fro')/norm(X,'fro'))
    savefig(fig3, destino +"RIF lambda= "+lambda+".fig")
end



