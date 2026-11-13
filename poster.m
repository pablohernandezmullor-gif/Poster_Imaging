%% poster.m
% Deconvolución de una imagen de microscopía de fluorescencia:
% células endoteliales de arteria pulmonar bovina (Datos/FluorescentCells.jpg,
% imagen de ejemplo de ImageJ, dominio público).
%   Azul: núcleos (DAPI) · Verde: microtúbulos (tubulina) · Rojo: actina F
%
% Datos de prueba creados como en Hansen, Nagy & O'Leary, sección 4.6
% ("Creating Realistic Test Data"):
%   1) PSF gaussiana pequeña (psfGauss) rellenada con ceros (padPSF).
%   2) Se emborrona la imagen GRANDE con la FFT y se extrae la parte
%      central. Así la imagen emborronada B es realista: no impone ninguna
%      condición de contorno artificial en sus bordes.
%   3) Se añade un 1% de ruido gaussiano blanco: ||e||/||Ax|| = 0.01.
%
% Reconstrucción (con la FFT, es decir, asumiendo contorno periódico):
%   - Filtro inverso ingenuo (naive):            X = G / H
%   - Filtro inverso regularizado con laplaciano (LSS / RIF):
%                                     X = H* G / (|H|^2 + lambda |L|^2)
% Cada canal de color (R, G, B) se procesa con la misma PSF: fft2 actúa
% sobre cada página de un array m x n x 3.
%
% Solo usa funciones de MATLAB base (no hace falta Image Processing Toolbox).

clear; close all;
rng(0);                                   % ruido reproducible

%% Parámetros
dimPSF   = [31 31];     % tamaño de la PSF pequeña (píxeles)
sigmaPSF = 3;           % anchura de la gaussiana (píxeles)
margen   = 50;          % filas/columnas descartadas en cada borde (>= mitad de la PSF)
nivel    = 0.01;        % nivel de ruido: ||e|| / ||Ax||  (1%)
lambdas  = logspace(-5, 2, 36);           % valores de lambda a explorar
borde    = 30;          % píxeles excluidos al medir el error "interior"

carpeta = fileparts(mfilename('fullpath'));
salida  = fullfile(carpeta, "Output/Poster");
if ~exist(salida, 'dir'); mkdir(salida); end

%% 1. Imagen grande: microscopía de fluorescencia (RGB)
Xbig = double(imread(fullfile(carpeta, "Datos", "FluorescentCells.jpg"))) / 255;

%% 2. PSF gaussiana y emborronado de la imagen grande (Hansen, sec. 4.6)
[P, center] = psfGauss(dimPSF, sigmaPSF);
Pbig = padPSF(P, [size(Xbig,1) size(Xbig,2)]);
Sbig = fft2(circshift(Pbig, 1 - center));
Bbig = real(ifft2(Sbig .* fft2(Xbig)));   % la misma PSF en los 3 canales

% Se extrae la parte central: imagen verdadera X y emborronada B
filas = margen+1 : size(Xbig,1)-margen;
cols  = margen+1 : size(Xbig,2)-margen;
X = Xbig(filas, cols, :);
B = Bbig(filas, cols, :);
[m, n, ~] = size(X);

%% 3. Ruido gaussiano blanco del 1%
E = randn(size(B));
E = E / norm(E(:));
B = B + nivel * norm(B(:)) * E;
G = fft2(B);

%% 4. PSF del tamaño de B para la reconstrucción (contorno periódico)
H = fft2(circshift(padPSF(P, [m n]), 1 - center));

%% 5. Filtro inverso ingenuo:  X = G / H
Xnaive = real(ifft2(G ./ H));

%% 6. Filtro inverso regularizado con laplaciano (LSS), como en solucion_challenge1.m
[u, v] = ndgrid(0:m-1, 0:n-1);
L2 = (4 - 2*cos(2*pi*u/m) - 2*cos(2*pi*v/n)).^2;   % |Laplaciano|^2

% Error relativo en toda la imagen y solo en el interior. Como B no es
% periódica, el filtro (que sí asume contorno periódico) crea artefactos en
% los bordes; el error interior mide la calidad lejos de ellos.
ri = borde+1 : m-borde;   ci = borde+1 : n-borde;
Xi = X(ri, ci, :);
err     = @(Y) norm(Y(:) - X(:)) / norm(X(:));
err_int = @(Y) norm(reshape(Y(ri,ci,:), [], 1) - Xi(:)) / norm(Xi(:));

errores = zeros(size(lambdas));  errores_int = errores;
for k = 1:numel(lambdas)
    Xk = real(ifft2( conj(H) .* G ./ (abs(H).^2 + lambdas(k)*L2) ));
    errores(k)     = err(Xk);
    errores_int(k) = err_int(Xk);
end
[~, kmin] = min(errores_int);             % lambda óptimo según el interior
lambda_opt = lambdas(kmin);
Xlss = real(ifft2( conj(H) .* G ./ (abs(H).^2 + lambda_opt*L2) ));

fprintf('                       error total   error interior\n');
fprintf('Emborronada         :  %.3e     %.3e\n', err(B), err_int(B));
fprintf('Naive               :  %.3e     %.3e\n', err(Xnaive), err_int(Xnaive));
fprintf('LSS (lambda=%.1e) :  %.3e     %.3e\n', lambda_opt, err(Xlss), err_int(Xlss));

%% 7. Figuras
if usejava('jvm')       % las figuras necesitan la JVM de MATLAB
    recorta  = @(Y) min(max(Y, 0), 1);                       % para mostrar en [0,1]
    reescala = @(Y) (Y - min(Y(:))) / (max(Y(:)) - min(Y(:)));

    % Comparación: original, emborronada, naive y LSS
    fig1 = figure('Name', 'Comparación', 'Position', [100 100 950 950]);
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile; image(X); axis image off
    title('Original')

    nexttile; image(recorta(B)); axis image off
    title(sprintf('Emborronada (\\sigma = %g px) + %g%% ruido', sigmaPSF, 100*nivel))

    nexttile; image(reescala(Xnaive)); axis image off
    title(sprintf('Filtro inverso ingenuo, err = %.1e', err(Xnaive)))

    nexttile; image(recorta(Xlss)); axis image off
    title(sprintf('LSS, \\lambda = %.1e, err int. = %.3f', lambda_opt, err_int(Xlss)))

    savefig(fig1, fullfile(salida, "fluor_comparacion.fig"));

    % Cada imagen por separado
    imagenes = {X, recorta(B), reescala(Xnaive), recorta(Xlss)};
    nombres  = ["fluor_original", "fluor_emborronada", "fluor_naive", "fluor_LSS"];
    titulos  = ["Original", ...
                sprintf("Emborronada (\\sigma = %g px) + %g%% ruido", sigmaPSF, 100*nivel), ...
                sprintf("Filtro inverso ingenuo, err = %.1e", err(Xnaive)), ...
                sprintf("LSS, \\lambda = %.1e, err int. = %.3f", lambda_opt, err_int(Xlss))];
    for k = 1:numel(imagenes)
        f = figure('Name', nombres(k));
        image(imagenes{k}); axis image off
        title(titulos(k))
        savefig(f, fullfile(salida, nombres(k) + ".fig"));
    end

    % Error relativo frente a lambda
    fig2 = figure('Name', 'Error frente a lambda');
    loglog(lambdas, errores, 'o-', 'LineWidth', 1.2); hold on
    loglog(lambdas, errores_int, 's-', 'LineWidth', 1.2)
    loglog(lambda_opt, errores_int(kmin), 'r*', 'MarkerSize', 12, 'LineWidth', 1.5)
    yline(err_int(B), '--', 'Emborronada (interior)');
    xlabel('\lambda'); ylabel('Error relativo  ||X_\lambda - X|| / ||X||')
    legend('Imagen completa', sprintf('Interior (sin %d px de borde)', borde), ...
           '\lambda óptimo', 'Location', 'northwest')
    title('Filtro LSS: elección de \lambda'); grid on

    savefig(fig2, fullfile(salida, "fluor_error_lambda.fig"));

    disp("Figuras guardadas en " + salida);
end


%Graficamos la PSF
figure
%Sacamos la imagen de PSF, centrada
centro_img   = floor(size(Pbig)/2) + 1;              % centro de la matriz (257,257)
Pbig_centrada = circshift(Pbig, centro_img - center); % lleva el máximo de la PSF ahí

figure
imagesc(Pbig_centrada); axis image; colorbar
title('PSF gaussiana (\sigma = 3 px)')


%% Funciones locales (obtenidas de psfGauss y padPSF del libro de Hansen)
function [PSF, center] = psfGauss(dim, s)
% PSF gaussiana de tamaño dim = [m n] y desviación típica s (píxeles).
% center es la posición del máximo de la PSF.
m = dim(1); n = dim(2);
x = -fix(n/2) : ceil(n/2)-1;
y = -fix(m/2) : ceil(m/2)-1;
[Xg, Yg] = meshgrid(x, y);
PSF = exp(-(Xg.^2 + Yg.^2) / (2*s^2));
PSF = PSF / sum(PSF(:));                  % normalizada: conserva la intensidad
center = [fix(m/2)+1, fix(n/2)+1];
end


function Pbig = padPSF(P, dim)
% Rellena con ceros la PSF P hasta el tamaño dim, dejándola en la esquina
% superior izquierda (así P y Pbig tienen el mismo centro).
Pbig = zeros(dim);
Pbig(1:size(P,1), 1:size(P,2)) = P;
end
