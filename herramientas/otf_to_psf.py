"""Convert the radially averaged TIRF OTF (.mrc) into a 2D PSF on the 512x512 image grid."""
import sys, numpy as np, tifffile
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

base = Path(sys.argv[1])          # .../2D data
outdir = Path(sys.argv[2])        # where to save the PSF
preview = Path(sys.argv[3])       # preview png

with tifffile.TiffFile(base / "Lysosome/train_data/cell01.tif") as t:
    img = t.asarray().astype(float)
    tags = t.pages[0].tags
    xres = tags["XResolution"].value if "XResolution" in tags else None
print("XResolution tag:", xres)

raw = np.fromfile(base / "Lysosome/PSF/TIRF560_cam2_0_z21_OTF2d.mrc", dtype="<f4", offset=1024)
otf_r = (raw[0::2]).astype(float)                 # imaginary part is exactly 0
dk = np.frombuffer(open(base / "Lysosome/PSF/TIRF560_cam2_0_z21_OTF2d.mrc", "rb").read(52)[40:44], "<f4")[0]
N = img.shape[0]
print(f"dk = {dk:.5f} 1/um  ->  pixel size implied for N={N}: {1000/(N*dk):.1f} nm")

# 2D OTF: |k| in units of the radial samples (dk) on an N x N FFT grid
k = np.fft.fftfreq(N) * N                         # 0,1,..,N/2-1,-N/2,..,-1
KX, KY = np.meshgrid(k, k)
kr = np.hypot(KX, KY)
otf2 = np.interp(kr, np.arange(otf_r.size), otf_r, right=0.0)
psf = np.real(np.fft.ifft2(otf2))
psf = np.fft.fftshift(psf)                        # centre at (N/2, N/2) (0-based)
psf[psf < 0] = 0
psf /= psf.sum()
c = N // 2
row = psf[c, :]
fwhm = (row >= row.max() / 2).sum()
print(f"PSF peak at {np.unravel_index(psf.argmax(), psf.shape)}, FWHM ~ {fwhm} px")

outdir.mkdir(parents=True, exist_ok=True)
tifffile.imwrite(outdir / "PSF_Lysosome_TIRF560_512.tif", psf.astype(np.float32))
small = psf[c-32:c+33, c-32:c+33]
tifffile.imwrite(outdir / "PSF_Lysosome_TIRF560_65.tif", (small / small.sum()).astype(np.float32))

fig, ax = plt.subplots(1, 3, figsize=(13, 4.3))
ax[0].imshow(img, cmap="gray"); ax[0].set_title("Lysosome cell01 (raw TIRF)"); ax[0].axis("off")
ax[1].imshow(small, cmap="magma"); ax[1].set_title(f"PSF from measured OTF (65x65 px), FWHM~{fwhm} px"); ax[1].axis("off")
ax[2].semilogy(np.arange(otf_r.size) * dk, np.abs(otf_r) + 1e-6); ax[2].set_xlabel("k (1/um)")
ax[2].set_title("Radial OTF"); ax[2].grid(True)
fig.tight_layout(); fig.savefig(preview, dpi=110)
