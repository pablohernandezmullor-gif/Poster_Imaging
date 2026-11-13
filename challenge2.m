%% Challenge 2
load("challenge2.mat")
fig=figure(1);

subplot(2,1,1)
imagesc(B)
title("Blurred image")
colormap gray

subplot(2,1,2)
X_naive=Ac\B/Ar';
imagesc(X_naive)
title("Naive filter deblurring")
colormap gray

carpeta = fileparts(mfilename('fullpath'));
destino=carpeta+"\Output"+"\Naive_Linear_Deblurring_2.fig";
savefig(fig,destino)
disp("Imagen guardada en  "+ destino);