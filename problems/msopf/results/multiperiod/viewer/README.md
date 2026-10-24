# Four-substation region viewer

`four_regions.html` is a self-contained page (data included): the four-substation feasible angle
regions, predicted from the angle-to-import sensitivities, against the measured boundary points,
for ieee123 (two source sets), the 2522-bus feeder and the 2705-bus feeder (from the IEEE 9500). Rotate, view each face edge-on (is it
flat?) or face-on (where is the prediction off?), and read each face's residuals.

Online: <https://realife-brahmin.github.io/MultiPeriodDistOPF/four_regions.html> (GitHub Pages,
served from `docs/` on master; after a rebuild, copy the page to
`docs/four_regions.html`). Offline: download the file and double-click it (it needs an internet connection only for the
plotting library), or view it online through
<https://raw.githack.com/Realife-Brahmin/MultiPeriodDistOPF/msopf-oct09-further-investigations/envs/multi_poi/results/multiperiod/viewer/four_regions.html>.

Rebuilt from `viewer_head.html` + `../four_viewer_data_fib78.json` + `viewer_tail.html`; the data
come from `envs/multi_poi/root_level/four_viewer_data.jl fib78`.
