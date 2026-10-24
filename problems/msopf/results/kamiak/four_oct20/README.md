# Four-substation regions, Kamiak, 2026-10-20

Stages `four_ieee123a` (ieee123_5poi sources 1-4), `four_ieee123b` (sources 1,3,4,5) and
`four_ieee2522` (deck source plus substations at buses 947, 996, 1120, 0.311 ohm), T = 24, PV
reactive power dispatched, 26 directions (a cube's faces, edges, corners) in the three angle
differences, idle and cyclic, by the secant search (jobs 31328463-76, commit c0587aa1); plus the
affine import model at 30 random idle points per system. Merged into
`../../multiperiod/multi_windows.csv` and `multi_linearity.csv`; faces checked by
`root_level/check_four.jl` -> `../../multiperiod/four_faces_check.txt`.

- Affine import model: largest hourly error 0.36% / 0.71% (ieee123), 3.6% (ieee2522).
- ieee123: four faces, one per substation, planar to 0-0.8% and within 0-1.1 deg of the predicted
  normals in 13 of 16 face fits; faces at 1.01-1.06x the predicted distance, except substation 4's,
  pushed out 1.19-1.47x by PV reactive dispatch (as in the three-substation case). Cyclic faces
  2.0-2.3x farther out than idle ones, same normals: nested tetrahedra.
- ieee2522: the 947 and 1120 faces planar (1.4-4.6%) and at 1.04-1.11x; the 996 face not planar
  (12-15%) and 0.92-0.97x, where the voltage limit binds (10 of the 52 directions fell back from the
  secant to the bracket).
- Time: ieee123 jobs 15-25 min; ieee2522 cyclic 1.5-2.5 h per 4-5 directions.
