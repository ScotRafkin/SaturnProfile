"""Digitization of Smith et al. (1982), Science 215, 504, Fig. 4 (left panel): zonal wind
against planetographic latitude, one point per tracked cloud feature, System III frame.

Method (10 September 2026): page 3 of the journal PDF rendered at 600 dpi; the plot frame
located from full-length black rows and columns; the interior thresholded, opened with a
5x5 structuring element to remove tick marks, and labeled; each connected component
assigned k = round(area / A1) dots with A1 the median area of isolated dots (262 px at this
scale), and k centroids found by k-means on the component's pixels when k > 1; pixel
coordinates mapped linearly with the frame edges at -50 and 500 m/s and at +90 and -90 deg.
The overlay PNG beside this script is the check. Reading error is about one dot radius,
roughly 2.5 m/s and 0.8 deg; merged clusters carry a larger latitude error.
"""
import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy.cluster.vq import kmeans2

PAGE_PNG = "smith_p3-03.png"          # pdftoppm -r 600 -f 3 -l 3 -png <Smith 1982 pdf> smith_p3
CROP = (420, 3720, 2560, 5960)         # Fig. 4 region on the 600 dpi page (5100 x 6600)

im = np.array(Image.open(PAGE_PNG).convert("L").crop(CROP)); b = im < 128
rows = b.sum(axis=1); cols = b.sum(axis=0)
r = np.where(rows > 0.6 * b.shape[1])[0]; c = np.where(cols > 0.6 * b.shape[0])[0]
top, bot, left, right = r.min() + 1, r.max() - 1, c.min() + 1, c.max() - 1
inner = b[top:bot, left:right].copy()
for m in (12,):
    inner[:m, :] = False; inner[-m:, :] = False; inner[:, :m] = False; inner[:, -m:] = False
inner = ndi.binary_opening(inner, structure=np.ones((5, 5)))
lab, n = ndi.label(inner)
areas = np.array(ndi.sum(inner, lab, range(1, n + 1)))
A1 = np.median(areas[(areas > 150) & (areas < 400)])
pts = []
for i in range(1, n + 1):
    a = areas[i - 1]
    if a < 80:
        continue
    rr, cc = np.where(lab == i)
    k = max(1, int(round(a / A1)))
    if k == 1:
        pts.append((rr.mean(), cc.mean()))
    else:
        cent, _ = kmeans2(np.c_[rr, cc].astype(float), k, minit="++", seed=1)
        pts.extend((p[0], p[1]) for p in cent)
pts = np.array(pts)
u = -50 + pts[:, 1] / (right - left) * 550
lat = 90 - pts[:, 0] / (bot - top) * 180
o = np.argsort(-lat)
np.savetxt("smith1982_fig4_points.csv", np.c_[lat, u][o], delimiter=",",
           header="latitude_planetographic_deg,u_ms", comments="", fmt="%.2f")
print(len(pts), "points; A1 =", A1)
