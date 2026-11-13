# Deconvolución de imágenes biomédicas

Código del póster para la asignatura *Procesado de señales* (Máster en Física Biomédica).
Tema 3: **un algoritmo de deconvolución de imágenes biomédicas**.

El algoritmo estudiado es el **filtro inverso regularizado con laplaciano** (Tikhonov de segundo orden,
también llamado *Constrained Least Squares* o RIF en DeconvolutionLab2), que se compara con el
filtro inverso ingenuo:

$$
\tilde F(\mathbf u) = \frac{H^*(\mathbf u)}{|H(\mathbf u)|^2 + \lambda\,|L(\mathbf u)|^2}\,G(\mathbf u)
$$

donde $H$ es la respuesta en frecuencia de la PSF, $G$ la imagen emborronada y $|L|^2$ el módulo al
cuadrado del laplaciano discreto. Con $\lambda = 0$ se recupera el filtro inverso ingenuo $G/H$.
Es un caso particular del filtro de Wiener con $P_n/P_f = \lambda|L|^2$.

## Archivos

| Archivo | Contenido |
|---|---|
| `poster.m` | Script principal del póster. Emborrona una imagen de microscopía de fluorescencia con una PSF gaussiana, añade ruido y compara el filtro ingenuo con el LSS. |
| `solucion_challenge1.m` | Primer acercamiento con el *Challenge 1* del libro de Hansen: inversión ingenua, extracción de la PSF y filtro LSS. |
| `challenge1.m`, `challenge2.m` | Generadores de los *challenges* del libro de Hansen, Nagy & O'Leary. |
| `Chapter_2.mlx` | Notas del capítulo 2 del libro. |
| `herramientas/otf_to_psf.py` | Convierte la OTF radial medida (formato `.mrc`) del conjunto ZS-DeconvNet en una PSF 2D en TIFF. |
| `Datos/FluorescentCells.jpg` | Imagen de fluorescencia usada en `poster.m`. |

## Cómo ejecutar `poster.m`

1. Abre MATLAB en esta carpeta. Solo se usan funciones de MATLAB base: **no hace falta la Image
   Processing Toolbox** (`psfGauss` y `padPSF` están implementadas como funciones locales).
2. Ejecuta `poster`.
3. Las figuras se guardan con `savefig` en `Outputs/` con el prefijo `fluor_`.

### Qué hace

1. **Datos de prueba realistas** (Hansen et al., sec. 4.6): se emborrona la imagen completa con una PSF
   gaussiana ($31\times31$, $\sigma = 3$ px) usando la FFT y se extrae la parte central. Así la imagen
   emborronada no impone ninguna condición de contorno artificial.
2. **Ruido** gaussiano blanco del 1 %: $\|e\|/\|Ax\| = 0.01$.
3. **Reconstrucción** con el filtro ingenuo y con el LSS para un barrido de $\lambda$, eligiendo el
   $\lambda$ que minimiza el error relativo en el interior de la imagen.

### Resultado principal

| | Error total | Error interior (sin 30 px de borde) |
|---|---|---|
| Emborronada | 0.213 | 0.207 |
| Filtro ingenuo | $2.3\cdot10^{12}$ | $2.2\cdot10^{12}$ |
| LSS ($\lambda = 6.3\cdot10^{-4}$) | 0.373 | **0.127** |

El filtro ingenuo amplifica el ruido sin control. El LSS reduce el error a algo más de la mitad en el
interior, pero aparece *ringing* en los bordes: el filtro asume **condiciones de contorno periódicas**
(FFT) y los datos no lo son. Las mejoras naturales son el contorno reflexivo (DCT), el *padding* o la
apodización, como en las opciones `-pad` y `-apo` de DeconvolutionLab2.

## Datos

| Datos | Fuente | Licencia |
|---|---|---|
| `FluorescentCells.jpg`: células endoteliales de arteria pulmonar bovina (DAPI, tubulina, actina F) | [Wikimedia Commons](https://commons.wikimedia.org/wiki/File:FluorescentCells.jpg), imagen de ejemplo de ImageJ | Dominio público |
| Imágenes TIRF de lisosomas con OTF medida, y microtúbulos con PSF (no incluidas) | [Zenodo 7261163](https://zenodo.org/records/7261163), archivo `2D data.zip` | CC-BY 4.0 |
| Imágenes y datos de los *challenges* (`*.tif`, `*.mat`, no incluidos) | Web del libro de Hansen, Nagy & O'Leary | — |

## Referencias

- P. C. Hansen, J. G. Nagy, D. P. O'Leary. *Deblurring Images: Matrices, Spectra, and Filtering*. SIAM, 2006.
- D. Sage et al. "DeconvolutionLab2: An open-source software for deconvolution microscopy".
  *Methods* 115, 28–41 (2017). [doi:10.1016/j.ymeth.2016.12.015](https://doi.org/10.1016/j.ymeth.2016.12.015)
- C. Qiao et al. "Zero-shot learning enables instant denoising and super-resolution in optical
  fluorescence microscopy". *Nature Communications* 15, 4180 (2024).
  [doi:10.1038/s41467-024-48575-9](https://doi.org/10.1038/s41467-024-48575-9)
