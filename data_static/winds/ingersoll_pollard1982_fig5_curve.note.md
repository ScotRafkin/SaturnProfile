# Ingersoll and Pollard (1982) Fig. 5: the smoothed Voyager 2 zonal wind profile of Saturn

Source: A. P. Ingersoll and D. Pollard, "Motion in the interiors and atmospheres of Jupiter and
Saturn: scale analysis, anelastic equations, barotropic stability criterion," Icarus 52, 62 to 80
(1982), Fig. 5, left panel. Scan: `Ingersoll_Pollard_1982.pdf`, page 6. Digitized 11 September
2026 by `digitize_ingersoll_pollard1982_fig5.py`; the overlay
`ingersoll_pollard1982_fig5_digitization_overlay.png` shows the samples on the scan and against
the Smith et al. (1982) Fig. 4 points.

## What the figure is

Caption, quoted: "Zonal velocity profile (left) and its curvature or second derivative (right)
for Saturn from Voyager 2 in late August, 1981. The data are from Smith et al. (1982) and are
based on far fewer points than the similar curves for Jupiter (Figs. 3 and 4). The two solid
curves on the right have the same meaning as in Fig. 4. Effects of the metallic core become
important at about 65 degrees latitude. Data are missing at southern equatorial latitudes
because that region was obscured by rings and ring shadow at the time of the Voyager
encounters."

The solid curve is the authors' smooth profile through the Smith et al. (1982) cloud-tracking
data; how it was drawn is not stated. The dashed curve is defined in the Fig. 3 caption for
Jupiter, which Fig. 5 follows: "The dashed curves show the same curves reflected about the
equator to test for symmetry." Latitude convention: the Fig. 3 caption says "the latitudes are
planetographic"; Fig. 5 does not restate it, and the digitized curve overlays the planetographic
Smith points without a latitude offset, which is taken as confirmation. Wind is with respect to
System III (Smith et al. 1982).

Lindal et al. (1985) cite this paper alongside Smith et al. (1982) as the source of the zonal
winds used in their geodetic calculations. How they combined the two is not stated.

## What was digitized

The solid curve only: 708 samples at 0.2 degree spacing, from +81.1 to +1.3 degrees planetographic
(north of the ring gap) and from −10.9 to −72.8 degrees (south of it). The dashed reflected
curve was not digitized: by definition it is the solid curve mirrored, which a tool can
construct exactly from these samples. Nothing is recorded in the ring gap (−10.9 to +1.3
degrees) or poleward of the curve's ends; filling those is the wind tool's declared rule, not
this file's.

Method: 600 dpi render, tick-mark calibration interpolated piecewise-linearly between ticks,
8-connected component labeling (the solid curve forms the two largest components, the dashes
fall into small ones and are discarded), skeletonization of each component with the skeleton's
longest path (graph diameter) taken as the center line, which drops the spurs where a dash
touches the solid line and follows steep segments without smearing them; a 3 pixel Gaussian
along the path removes the pixel staircase; one sample per 0.2 degrees by averaging path points
in latitude bins. A first version traced row by row and produced stair artifacts on the steep
flanks near 10 to 13 degrees north; that method was replaced. Reading error is about 5 m/s in
wind and 0.3 degrees in latitude, from the drawn line width.

## Checks

Against the 323 Smith et al. (1982) Fig. 4 points digitized separately: 314 points lie within
the curve's latitude range and outside the gap; their residuals about the curve have RMS 21.3
m/s and mean −4.3 m/s. The pooled within-bin scatter of the same points about their own 2
degree bin means is 21.4 m/s, so the curve is, to the precision of the digitizations, the mean
curve of the Smith data. Values: 490.5 m/s peak at +7.4 degrees; 1.9 m/s at +36.3 degrees;
75 m/s at +30.8 degrees; about 33 m/s at the northern end (+81.1); about 8 m/s at the southern
end (−72.8).

## Caveat on the southern equatorial gap

The ten Smith points on the southern flank, −7.5 to −10.9 degrees, read 306 to 385 m/s, while
the northern curve reflected across the equator gives 420 to 480 m/s at those latitudes.
Measured point by point against the assembled curve of the wind tool (REPORT_01_step7), the
nine points inside the gap lie 95 to 160 m/s below it. The reflection therefore overstates
the southern flank by of order 100 m/s, or, read the other way,
the southern jet is a few degrees narrower than the northern, or the flank points are affected
by ring shadow. The peak itself, −7 to 0 degrees, has no points at all. The reflection is the
declared fill because it is the authors' own construction and the only estimate of the peak;
the flank discrepancy is recorded here so that the sensitivity can be sized.
