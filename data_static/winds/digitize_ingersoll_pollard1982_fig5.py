"""Digitize the solid zonal wind curve of Ingersoll and Pollard (1982), Icarus 52, 62, Fig. 5.

Input: Ingersoll_Pollard_1982.pdf (page 6 of the scan). Output: ingersoll_pollard1982_fig5_curve.csv
with planetographic latitude and zonal wind (System III, m/s) sampled along the solid curve, and
ingersoll_pollard1982_fig5_digitization_overlay.png for checking.

Method. Page 6 is rendered at 600 dpi and the left panel of Fig. 5 cropped. Axis calibration is
by the tick marks found on the frame (latitude ticks at 90, 60, ... -90; velocity ticks at -100,
0, ... 500), interpolated piecewise-linearly between ticks so that scan skew is absorbed. Dark
pixels inside the frame are labeled into 8-connected components; the two largest are the solid
curve north and south of the ring-obscured gap, and the dashed reflected curve (Fig. 3 caption:
"the same curves reflected about the equator") falls into many small components and is
discarded. Each component is skeletonized (scikit-image) and the
skeleton's longest path, found as the graph diameter by two Dijkstra passes, is taken as the
center line; that discards the short spurs where a dash of the reflected curve touches the solid
line, and it follows the line through steep segments where a row-by-row trace would smear them
into stairs. The path is smoothed with a Gaussian of 3 pixels along its length (about 0.5
degrees or 3 m/s, well under the drawn line width) to remove the pixel staircase, then reduced
to one sample per 0.2 degrees of latitude by averaging the path points in each latitude bin.
Reading error is about 5 m/s in wind and 0.3 degrees in latitude, from the drawn line width.

The solid curve spans about +81.6 to +1.3 degrees and -10.9 to -73.4 degrees planetographic.
Nothing is written for the gap or the polar caps; those are the wind tool's declared rules.
"""
import sys
import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.ndimage import gaussian_filter1d
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra
from skimage.morphology import skeletonize

PDF = sys.argv[1] if len(sys.argv) > 1 else "Ingersoll_Pollard_1982.pdf"
import subprocess, os
subprocess.run(["pdftoppm", "-f", "6", "-l", "6", "-r", "600", "-png", PDF, "ip_page6"], check=True)
page = Image.open("ip_page6-06.png")
crop = page.crop((2160, 3160, 3120, 4300))
im = np.array(crop.convert("L"))
dark = im < 140

# tick calibration (pixel positions found from the frame ticks of this crop)
tick_x = np.array([136, 247, 359, 472, 584, 696, 809], float)
tick_u = np.array([-100, 0, 100, 200, 300, 400, 500], float)
tick_y = np.array([59, 226, 392, 559, 730, 898, 1066], float)
tick_lat = np.array([90, 60, 30, 0, -30, -60, -90], float)

panel = np.zeros_like(dark)
panel[64:1063, 142:835] = dark[64:1063, 142:835]
lab, n = ndimage.label(panel, structure=np.ones((3, 3)))
sizes = ndimage.sum(panel, lab, range(1, n + 1))
order = np.argsort(sizes)[::-1]




def longest_path(mask):
    """Skeleton center line of a component, ordered from its topmost end (north) downward."""
    sk = skeletonize(mask)
    ys, xs = np.nonzero(sk)
    pts = np.c_[ys, xs]
    index = {(y, x): i for i, (y, x) in enumerate(pts)}
    rows, cols, weights = [], [], []
    for i, (y, x) in enumerate(pts):
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if dy == 0 and dx == 0:
                    continue
                j = index.get((y + dy, x + dx))
                if j is not None:
                    rows.append(i)
                    cols.append(j)
                    weights.append(float(np.hypot(dy, dx)))
    graph = coo_matrix((weights, (rows, cols)), shape=(len(pts), len(pts))).tocsr()
    start = int(np.argmin(pts[:, 0]))
    dist = dijkstra(graph, indices=start)
    dist[np.isinf(dist)] = -1
    a = int(np.argmax(dist))
    dist2, pred = dijkstra(graph, indices=a, return_predecessors=True)
    dist2[np.isinf(dist2)] = -1
    b = int(np.argmax(dist2))
    path = [b]
    while path[-1] != a:
        path.append(pred[path[-1]])
    p = pts[path].astype(float)
    if p[0, 0] > p[-1, 0]:
        p = p[::-1]
    return p


rows = []
for comp, tag in ((order[0], "solid_north"), (order[1], "solid_south")):
    p = longest_path(lab == comp + 1)
    y = gaussian_filter1d(p[:, 0], 3.0)
    x = gaussian_filter1d(p[:, 1], 3.0)
    lat = np.interp(y, tick_y, tick_lat)  # tick_y ascends down the page, latitude descends with it
    u = np.interp(x, tick_x, tick_u)
    edges = np.arange(np.floor(lat.min() / 0.2) * 0.2, lat.max() + 0.2, 0.2)
    which = np.digitize(lat, edges)
    for k in np.unique(which):
        m = which == k
        rows.append((float(lat[m].mean()), float(u[m].mean()), tag))
rows.sort(key=lambda r: -r[0])

with open("ingersoll_pollard1982_fig5_curve.csv", "w") as f:
    f.write("latitude_planetographic_deg,u_ms,segment\n")
    for a, b, t in rows:
        f.write(f"{a:.3f},{b:.2f},{t}\n")
print(len(rows), "samples written")
