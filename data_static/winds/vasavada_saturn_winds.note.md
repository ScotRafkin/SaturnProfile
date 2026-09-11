# `vasavada_saturn_winds.txt`: the Sanchez-Lavega et al. (2000) zonal wind profile

Received 11 September 2026 from J. Moses, who obtained it from A. Vasavada, who obtained it
from A. Sanchez-Lavega (per J. Moses); placed in `data_static/winds/` by the author. The file header reads "Saturn Zonal Winds from
Sanchez-Lavega et al. (2000)", with columns planetocentric latitude, planetographic latitude,
u (m/s), u_rms (m/s). 262 rows at 0.5 degrees planetocentric from +78.6 to −67.0, plus two rows
at ±90.0 marked in the file itself "fake data added here" (u = 0.0, rms 8.7); those two rows
are not data and are to be dropped on ingest, the polar rule of the wind tool applying instead.

**Epoch not yet confirmed.** Sanchez-Lavega, Rojas and Sada (2000, Icarus 147, 405) reanalyzed
the Voyager 1 and 2 images and compared them with 1990s Hubble data; whether this table is the
Voyager profile alone or a blend is to be read from the paper before the file is used. It was
described in the accompanying email as a Cassini profile; the header and the provenance chain
(Sanchez-Lavega to Vasavada to Moses) say it is the Sanchez-Lavega et al. (2000) table. The meaning
of u_rms (scatter of tracked features per latitude bin, or a formal error) is also to be read
from the paper. Frame: System III assumed, to be confirmed.

Checks made on receipt: the two latitude columns imply a flattening of 0.0983 ± 0.0012 (from
tan φ_c = (1 − f)² tan φ_g), consistent with the 100 mbar oblateness 0.09822, so the file's own
conversion is consistent with the reduction's. Against the digitized Ingersoll and Pollard
(1982) curve over 253 overlapping rows: mean difference −5.8 m/s, RMS 21.4 m/s, the same scatter
the Smith (1982) points show about that curve. Peak 467.1 m/s at +5.6° planetographic (curve:
490.5 at +7.4°). Coverage gaps in planetographic latitude: 80.7 to 90; −2.1 to −10.7 (the ring
gap, narrower than in Ingersoll and Pollard); −14.8 to −17.4; −70.9 to −90. Rows exist at −1.0
to −2.1° (about 436 m/s) and on the southern flank at −10.7 to −11.8° (377 to 349 m/s), which
agree with the Smith flank points and sit about 90 m/s below the reflected northern curve used
as the Step 7 gap fill. See `sanchezlavega2000_vs_ip1982.png` in the notes.

Intended use: a second declared reduction wind source (`wind_source = "sanchezlavega2000"`)
through the wind tool's curve path, with u_rms as its uncertainty, for the Step 7 sensitivity
and the SPEC_02 Monte Carlo; not a change to the accepted Step 7 build. Median u_rms 10.6 m/s.
