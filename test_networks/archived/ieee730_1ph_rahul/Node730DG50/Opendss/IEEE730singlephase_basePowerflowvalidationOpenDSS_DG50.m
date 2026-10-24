clear all;
clc;

[DSSStartOK, DSSObj, DSSText] = DSSStartup;

if DSSStartOK
    

DSSObj.AllowForms = false;    


DSSText.Command = 'Set DataPath = C:\Users\rahul\Dropbox\730nodesystem\lossmin\OpenDSS730node';

DSSText.Command = 'compile singlephase730.dss';      

%% defining PV in opendss 
DSSText.Command = 'New PVSystem.41 phases=1 bus1=41 kV=0.24 kVA=0.997 irradiance=1 Pmpp=0.831 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.44 phases=1 bus1=44 kV=0.24 kVA=1.051 irradiance=1 Pmpp=0.876 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.46 phases=1 bus1=46 kV=0.24 kVA=0.8304 irradiance=1 Pmpp=0.692 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.48 phases=1 bus1=48 kV=0.24 kVA=1.814 irradiance=1 Pmpp=1.512 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.50 phases=1 bus1=50 kV=0.24 kVA=0.9864 irradiance=1 Pmpp=0.822 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.52 phases=1 bus1=52 kV=0.24 kVA=1.061 irradiance=1 Pmpp=0.884 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.59 phases=1 bus1=59 kV=0.24 kVA=1.246 irradiance=1 Pmpp=1.038 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.61 phases=1 bus1=61 kV=0.24 kVA=1.812 irradiance=1 Pmpp=1.51 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.65 phases=1 bus1=65 kV=0.24 kVA=3.1032 irradiance=1 Pmpp=2.586 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.63 phases=1 bus1=63 kV=0.24 kVA=2.273 irradiance=1 Pmpp=1.894 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.68 phases=1 bus1=68 kV=0.24 kVA=2.340 irradiance=1 Pmpp=1.95 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.74 phases=1 bus1=74 kV=0.24 kVA=2.441 irradiance=1 Pmpp=2.034 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.77 phases=1 bus1=77 kV=0.24 kVA=1.589 irradiance=1 Pmpp=1.324 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.79 phases=1 bus1=79 kV=0.24 kVA=2.1696 irradiance=1 Pmpp=1.808 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.81 phases=1 bus1=81 kV=0.24 kVA=2.418 irradiance=1 Pmpp=2.015 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.87 phases=1 bus1=87 kV=0.24 kVA=1.664 irradiance=1 Pmpp=1.387 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.90 phases=1 bus1=90 kV=0.24 kVA=1.607 irradiance=1 Pmpp=1.339 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.94 phases=1 bus1=94 kV=0.24 kVA=3.527 irradiance=1 Pmpp=2.939 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.92 phases=1 bus1=92 kV=0.24 kVA=1.4472 irradiance=1 Pmpp=1.206 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.100 phases=1 bus1=100 kV=0.24 kVA=1.18 irradiance=1 Pmpp=0.983 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.103 phases=1 bus1=103 kV=0.24 kVA=1.795 irradiance=1 Pmpp=1.496 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.105 phases=1 bus1=105 kV=0.24 kVA=2.593 irradiance=1 Pmpp=2.161 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.107 phases=1 bus1=107 kV=0.24 kVA=0.830 irradiance=1 Pmpp=0.692 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.115 phases=1 bus1=115 kV=0.24 kVA=1.873 irradiance=1 Pmpp=1.561 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.118 phases=1 bus1=118 kV=0.24 kVA=2.848 irradiance=1 Pmpp=2.373 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.120 phases=1 bus1=120 kV=0.24 kVA=0.9072 irradiance=1 Pmpp=0.756 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.122 phases=1 bus1=122 kV=0.24 kVA=0.828 irradiance=1 Pmpp=0.690 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.124 phases=1 bus1=124 kV=0.24 kVA=0.956  irradiance=1 Pmpp=0.797 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.126 phases=1 bus1=126 kV=0.24 kVA=1.351 irradiance=1 Pmpp=1.126 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.134 phases=1 bus1=134 kV=0.24 kVA=1.86 irradiance=1 Pmpp=1.55 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.136 phases=1 bus1=136 kV=0.24 kVA=1.0764 irradiance=1 Pmpp=0.897 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.138 phases=1 bus1=138 kV=0.24 kVA=1.074 irradiance=1 Pmpp=0.895 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.140 phases=1 bus1=140 kV=0.24 kVA=0.7752 irradiance=1 Pmpp=0.646 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.142 phases=1 bus1=142 kV=0.24 kVA=1.014 irradiance=1 Pmpp=0.845 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.145 phases=1 bus1=145 kV=0.24 kVA=0.872 irradiance=1 Pmpp=0.727 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.153 phases=1 bus1=153 kV=0.24 kVA=0.808 irradiance=1 Pmpp=0.673 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.155 phases=1 bus1=155 kV=0.24 kVA=1.68 irradiance=1 Pmpp=1.4 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.157 phases=1 bus1=157 kV=0.24 kVA=2.581 irradiance=1 Pmpp=2.151 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.160 phases=1 bus1=160 kV=0.24 kVA=0.778 irradiance=1 Pmpp=0.648 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.162 phases=1 bus1=162 kV=0.24 kVA=0.888 irradiance=1 Pmpp=0.74 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.164 phases=1 bus1=164 kV=0.24 kVA=0.895 irradiance=1 Pmpp=0.746 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.172 phases=1 bus1=172 kV=0.24 kVA=1.523 irradiance=1 Pmpp=1.269 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.174 phases=1 bus1=174 kV=0.24 kVA=1.7892 irradiance=1 Pmpp=1.491 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.176 phases=1 bus1=176 kV=0.24 kVA=0.994 irradiance=1 Pmpp=0.828 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.178 phases=1 bus1=178 kV=0.24 kVA=0.7644 irradiance=1 Pmpp=0.637 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.180 phases=1 bus1=180 kV=0.24 kVA=1.374 irradiance=1 Pmpp=1.145 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.183 phases=1 bus1=183 kV=0.24 kVA=1.432 irradiance=1 Pmpp=1.193 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.189 phases=1 bus1=189 kV=0.24 kVA=0.932 irradiance=1 Pmpp=0.777 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.192 phases=1 bus1=192 kV=0.24 kVA=1.466 irradiance=1 Pmpp=1.222 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.194 phases=1 bus1=194 kV=0.24 kVA=2.6688 irradiance=1 Pmpp=2.224 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.196 phases=1 bus1=196 kV=0.24 kVA=2.474 irradiance=1 Pmpp=2.062 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.202 phases=1 bus1=202 kV=0.24 kVA=2.628 irradiance=1 Pmpp=2.19 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.204 phases=1 bus1=204 kV=0.24 kVA=2.243 irradiance=1 Pmpp=1.869 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.206 phases=1 bus1=206 kV=0.24 kVA=2.588 irradiance=1 Pmpp=2.157 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.209 phases=1 bus1=209 kV=0.24 kVA=3.491 irradiance=1 Pmpp=2.909 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.214 phases=1 bus1=214 kV=0.24 kVA=1.578 irradiance=1 Pmpp=1.315 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.216 phases=1 bus1=216 kV=0.24 kVA=2.586 irradiance=1 Pmpp=2.155 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.219 phases=1 bus1=219 kV=0.24 kVA=4.961 irradiance=1 Pmpp=4.134 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.225 phases=1 bus1=225 kV=0.24 kVA=3.719 irradiance=1 Pmpp=3.099 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.228 phases=1 bus1=228 kV=0.24 kVA=3.541 irradiance=1 Pmpp=2.951 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.230 phases=1 bus1=230 kV=0.24 kVA=1.6908 irradiance=1 Pmpp=1.409 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.232 phases=1 bus1=232 kV=0.24 kVA=3.637 irradiance=1 Pmpp=3.031 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.238 phases=1 bus1=238 kV=0.24 kVA=3.197 irradiance=1 Pmpp=2.664 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.241 phases=1 bus1=241 kV=0.24 kVA=2.495 irradiance=1 Pmpp=2.079 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.243 phases=1 bus1=243 kV=0.24 kVA=2.6124 irradiance=1 Pmpp=2.177 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.245 phases=1 bus1=245 kV=0.24 kVA=1.783 irradiance=1 Pmpp=1.486 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.253 phases=1 bus1=253 kV=0.24 kVA=1.693 irradiance=1 Pmpp=1.411 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.255 phases=1 bus1=255 kV=0.24 kVA=1.5468 irradiance=1 Pmpp=1.289 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.257 phases=1 bus1=257 kV=0.24 kVA=2.376 irradiance=1 Pmpp=1.98 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.259 phases=1 bus1=259 kV=0.24 kVA=2.7804 irradiance=1 Pmpp=2.317 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.261 phases=1 bus1=261 kV=0.24 kVA=2.762 irradiance=1 Pmpp=2.302 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.264 phases=1 bus1=264 kV=0.24 kVA=2.447 irradiance=1 Pmpp=2.039 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.270 phases=1 bus1=270 kV=0.24 kVA=3.936 irradiance=1 Pmpp=3.28 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.273 phases=1 bus1=273 kV=0.24 kVA=2.507 irradiance=1 Pmpp=2.089 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.275 phases=1 bus1=275 kV=0.24 kVA=2.412 irradiance=1 Pmpp=2.01 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.277 phases=1 bus1=277 kV=0.24 kVA=3.672 irradiance=1 Pmpp=3.06 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.283 phases=1 bus1=283 kV=0.24 kVA=2.783 irradiance=1 Pmpp=2.319 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.286 phases=1 bus1=286 kV=0.24 kVA=2.178 irradiance=1 Pmpp=1.815 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.288 phases=1 bus1=288 kV=0.24 kVA=2.9028 irradiance=1 Pmpp=2.419 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.290 phases=1 bus1=290 kV=0.24 kVA=2.196 irradiance=1 Pmpp=1.83 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.296 phases=1 bus1=296 kV=0.24 kVA=1.903 irradiance=1 Pmpp=1.586 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.298 phases=1 bus1=298 kV=0.24 kVA=2.8692 irradiance=1 Pmpp=2.391 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.300 phases=1 bus1=300 kV=0.24 kVA=2.627 irradiance=1 Pmpp=2.189 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.303 phases=1 bus1=303 kV=0.24 kVA=2.677 irradiance=1 Pmpp=2.231 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.310 phases=1 bus1=310 kV=0.24 kVA=2.104 irradiance=1 Pmpp=1.753 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.312 phases=1 bus1=312 kV=0.24 kVA=3.3516 irradiance=1 Pmpp=2.793 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.314 phases=1 bus1=314 kV=0.24 kVA=4.469 irradiance=1 Pmpp=3.724 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.316 phases=1 bus1=316 kV=0.24 kVA=1.6548 irradiance=1 Pmpp=1.379 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.319 phases=1 bus1=319 kV=0.24 kVA=1.667 irradiance=1 Pmpp=1.389 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.326 phases=1 bus1=326 kV=0.24 kVA=2.364 irradiance=1 Pmpp=1.97 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.328 phases=1 bus1=328 kV=0.24 kVA=1.5228 irradiance=1 Pmpp=1.269 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.330 phases=1 bus1=330 kV=0.24 kVA=3.145 irradiance=1 Pmpp=2.621 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.332 phases=1 bus1=332 kV=0.24 kVA=2.586 irradiance=1 Pmpp=2.155 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.335 phases=1 bus1=335 kV=0.24 kVA=1.813 irradiance=1 Pmpp=1.511 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.342 phases=1 bus1=342 kV=0.24 kVA=3.071 irradiance=1 Pmpp=2.559 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.344 phases=1 bus1=344 kV=0.24 kVA=2.0904 irradiance=1 Pmpp=1.742 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.347 phases=1 bus1=347 kV=0.24 kVA=3.731 irradiance=1 Pmpp=3.109 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.349 phases=1 bus1=349 kV=0.24 kVA=1.8756 irradiance=1 Pmpp=1.563 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.351 phases=1 bus1=351 kV=0.24 kVA=2.486 irradiance=1 Pmpp=2.072 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.358 phases=1 bus1=358 kV=0.24 kVA=3.028 irradiance=1 Pmpp=2.523 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.360 phases=1 bus1=360 kV=0.24 kVA=1.8756 irradiance=1 Pmpp=1.563 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.363 phases=1 bus1=363 kV=0.24 kVA=2.186 irradiance=1 Pmpp=1.822 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.365 phases=1 bus1=365 kV=0.24 kVA=2.2512 irradiance=1 Pmpp=1.876 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.367 phases=1 bus1=367 kV=0.24 kVA=2.318 irradiance=1 Pmpp=1.932 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.374 phases=1 bus1=374 kV=0.24 kVA=3.378 irradiance=1 Pmpp=2.815 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.376 phases=1 bus1=376 kV=0.24 kVA=10.1988 irradiance=1 Pmpp=8.499 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.379 phases=1 bus1=379 kV=0.24 kVA=8.419 irradiance=1 Pmpp=7.016 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.381 phases=1 bus1=381 kV=0.24 kVA=2.178 irradiance=1 Pmpp=1.815 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.383 phases=1 bus1=383 kV=0.24 kVA=2.276 irradiance=1 Pmpp=1.897 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.390 phases=1 bus1=390 kV=0.24 kVA=3.166 irradiance=1 Pmpp=2.638 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.392 phases=1 bus1=392 kV=0.24 kVA=2.022 irradiance=1 Pmpp=1.685 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.394 phases=1 bus1=394 kV=0.24 kVA=5.06 irradiance=1 Pmpp=4.217 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.396 phases=1 bus1=396 kV=0.24 kVA=4.8972 irradiance=1 Pmpp=4.081 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.399 phases=1 bus1=399 kV=0.24 kVA=4.07 irradiance=1 Pmpp=3.392 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.407 phases=1 bus1=407 kV=0.24 kVA=2.605 irradiance=1 Pmpp=2.171 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.409 phases=1 bus1=409 kV=0.24 kVA=1.6488 irradiance=1 Pmpp=1.374 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.411 phases=1 bus1=411 kV=0.24 kVA=1.330 irradiance=1 Pmpp=1.108 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.413 phases=1 bus1=413 kV=0.24 kVA=1.2288 irradiance=1 Pmpp=1.024 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.415 phases=1 bus1=415 kV=0.24 kVA=0.888 irradiance=1 Pmpp=0.74 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.418 phases=1 bus1=418 kV=0.24 kVA=1.465 irradiance=1 Pmpp=1.221 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.426 phases=1 bus1=426 kV=0.24 kVA=1.172 irradiance=1 Pmpp=0.977 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.430 phases=1 bus1=430 kV=0.24 kVA=1.2324 irradiance=1 Pmpp=1.027 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.430 phases=1 bus1=430 kV=0.24 kVA=1.566 irradiance=1 Pmpp=1.305 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.432 phases=1 bus1=432 kV=0.24 kVA=2.148 irradiance=1 Pmpp=1.79 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.434 phases=1 bus1=434 kV=0.24 kVA=1.028 irradiance=1 Pmpp=0.857 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.437 phases=1 bus1=437 kV=0.24 kVA=1.136 irradiance=1 Pmpp=0.947 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.445 phases=1 bus1=445 kV=0.24 kVA= 1.7388 irradiance=1 Pmpp=1.499 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.447 phases=1 bus1=447 kV=0.24 kVA=1.9572 irradiance=1 Pmpp=1.631 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.449 phases=1 bus1=449 kV=0.24 kVA= 1.7388 irradiance=1 Pmpp=1.449 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.451 phases=1 bus1=451 kV=0.24 kVA=1.7832 irradiance=1 Pmpp=1.486 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.453 phases=1 bus1=453 kV=0.24 kVA=1.122 irradiance=1 Pmpp=0.935 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.456 phases=1 bus1=456 kV=0.24 kVA=1.706 irradiance=1 Pmpp=1.422 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.464 phases=1 bus1=464 kV=0.24 kVA=1.615 irradiance=1 Pmpp=1.346 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.466 phases=1 bus1=466 kV=0.24 kVA=1.824 irradiance=1 Pmpp=1.52 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.468 phases=1 bus1=468 kV=0.24 kVA=1.435 irradiance=1 Pmpp=1.196 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.470 phases=1 bus1=470 kV=0.24 kVA=3.0576 irradiance=1 Pmpp=2.548 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.472 phases=1 bus1=472 kV=0.24 kVA=1.73 irradiance=1 Pmpp=1.442 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.475 phases=1 bus1=475 kV=0.24 kVA=1.481 irradiance=1 Pmpp=1.234 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.483 phases=1 bus1=483 kV=0.24 kVA=1.487 irradiance=1 Pmpp=1.239 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.486 phases=1 bus1=486 kV=0.24 kVA=2.414 irradiance=1 Pmpp=2.012 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.488 phases=1 bus1=488 kV=0.24 kVA=1.4472 irradiance=1 Pmpp=1.206 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.490 phases=1 bus1=490 kV=0.24 kVA=1.702 irradiance=1 Pmpp=1.418 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.492 phases=1 bus1=492 kV=0.24 kVA=1.4964 irradiance=1 Pmpp=1.247 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.494 phases=1 bus1=494 kV=0.24 kVA=1.429 irradiance=1 Pmpp=1.191 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.502 phases=1 bus1=502 kV=0.24 kVA=1.466 irradiance=1 Pmpp=1.222 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.504 phases=1 bus1=504 kV=0.24 kVA=1.6932 irradiance=1 Pmpp=1.411 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.506 phases=1 bus1=506 kV=0.24 kVA=1.465 irradiance=1 Pmpp=1.221 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.509 phases=1 bus1=509 kV=0.24 kVA=1.988 irradiance=1 Pmpp=1.657 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.511 phases=1 bus1=511 kV=0.24 kVA=1.5348 irradiance=1 Pmpp=1.279 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.513 phases=1 bus1=513 kV=0.24 kVA=1.632 irradiance=1 Pmpp=1.36 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.521 phases=1 bus1=521 kV=0.24 kVA=1.582 irradiance=1 Pmpp=1.318 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.523 phases=1 bus1=523 kV=0.24 kVA=2.041 irradiance=1 Pmpp=1.701 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.525 phases=1 bus1=525 kV=0.24 kVA=1.517 irradiance=1 Pmpp=1.264 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.527 phases=1 bus1=527 kV=0.24 kVA=1.517 irradiance=1 Pmpp=1.264 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.529 phases=1 bus1=529 kV=0.24 kVA=1.776 irradiance=1 Pmpp=1.48 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.532 phases=1 bus1=532 kV=0.24 kVA=1.742 irradiance=1 Pmpp=1.452 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.540 phases=1 bus1=540 kV=0.24 kVA=1.817 irradiance=1 Pmpp=1.514 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.542 phases=1 bus1=542 kV=0.24 kVA=3.6408 irradiance=1 Pmpp=3.034 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.544 phases=1 bus1=544 kV=0.24 kVA=1.835 irradiance=1 Pmpp=1.529 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.546 phases=1 bus1=546 kV=0.24 kVA=1.687 irradiance=1 Pmpp=1.406 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.548 phases=1 bus1=548 kV=0.24 kVA=1.754 irradiance=1 Pmpp=1.462 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.551 phases=1 bus1=551 kV=0.24 kVA=2.044 irradiance=1 Pmpp=1.703 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.557 phases=1 bus1=557 kV=0.24 kVA=2.011 irradiance=1 Pmpp=1.676 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.559 phases=1 bus1=559 kV=0.24 kVA=1.782 irradiance=1 Pmpp=1.485 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.561 phases=1 bus1=561 kV=0.24 kVA=1.945 irradiance=1 Pmpp=1.621 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.564 phases=1 bus1=564 kV=0.24 kVA=2.04 irradiance=1 Pmpp=1.7 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.570 phases=1 bus1=570 kV=0.24 kVA=1.597 irradiance=1 Pmpp=1.331 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.573 phases=1 bus1=573 kV=0.24 kVA=1.646 irradiance=1 Pmpp=1.372 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.575 phases=1 bus1=575 kV=0.24 kVA=2.1552 irradiance=1 Pmpp=1.796 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.577 phases=1 bus1=577 kV=0.24 kVA=4.856 irradiance=1 Pmpp=4.047 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.583 phases=1 bus1=583 kV=0.24 kVA=1.801 irradiance=1 Pmpp=1.501 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.586 phases=1 bus1=586 kV=0.24 kVA=3.119 irradiance=1 Pmpp=2.599 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.588 phases=1 bus1=588 kV=0.24 kVA=2.1552 irradiance=1 Pmpp=1.796 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.590 phases=1 bus1=590 kV=0.24 kVA=2.659 irradiance=1 Pmpp=2.216 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.598 phases=1 bus1=598 kV=0.24 kVA=1.848 irradiance=1 Pmpp=1.54 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.600 phases=1 bus1=600 kV=0.24 kVA=2.3076 irradiance=1 Pmpp=1.923 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.602 phases=1 bus1=602 kV=0.24 kVA=3.064 irradiance=1 Pmpp=2.553 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.604 phases=1 bus1=604 kV=0.24 kVA=3.8124 irradiance=1 Pmpp=3.177 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.606 phases=1 bus1=606 kV=0.24 kVA=1.75 irradiance=1 Pmpp=1.458 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.609 phases=1 bus1=609 kV=0.24 kVA=1.055 irradiance=1 Pmpp=0.879 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.617 phases=1 bus1=617 kV=0.24 kVA=0.881 irradiance=1 Pmpp=0.734 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.619 phases=1 bus1=619 kV=0.24 kVA=1.0512 irradiance=1 Pmpp=0.876 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.621 phases=1 bus1=621 kV=0.24 kVA=0.830 irradiance=1 Pmpp=0.692 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.623 phases=1 bus1=623 kV=0.24 kVA=1.814 irradiance=1 Pmpp=1.512 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.625 phases=1 bus1=625 kV=0.24 kVA=0.986 irradiance=1 Pmpp=0.822 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.628 phases=1 bus1=628 kV=0.24 kVA=1.246 irradiance=1 Pmpp=1.038 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.636 phases=1 bus1=636 kV=0.24 kVA=2.848 irradiance=1 Pmpp=2.373 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.638 phases=1 bus1=638 kV=0.24 kVA=1.9404 irradiance=1 Pmpp=1.617 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.641 phases=1 bus1=641 kV=0.24 kVA=3.103 irradiance=1 Pmpp=2.586 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.643 phases=1 bus1=643 kV=0.24 kVA=1.53 irradiance=1 Pmpp=1.275 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.645 phases=1 bus1=645 kV=0.24 kVA=2.441 irradiance=1 Pmpp=2.034 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.647 phases=1 bus1=647 kV=0.24 kVA=2.015 irradiance=1 Pmpp=1.679 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.654 phases=1 bus1=654 kV=0.24 kVA=1.589 irradiance=1 Pmpp=1.324 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.656 phases=1 bus1=656 kV=0.24 kVA=2.169 irradiance=1 Pmpp=1.808 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.658 phases=1 bus1=658 kV=0.24 kVA=2.418 irradiance=1 Pmpp=2.015 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.660 phases=1 bus1=660 kV=0.24 kVA=1.6356 irradiance=1 Pmpp=1.363 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.663 phases=1 bus1=663 kV=0.24 kVA=1.651 irradiance=1 Pmpp=1.376 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.670 phases=1 bus1=670 kV=0.24 kVA=1.447 irradiance=1 Pmpp=1.206 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.672 phases=1 bus1=672 kV=0.24 kVA=3.5268 irradiance=1 Pmpp=2.939 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.674 phases=1 bus1=674 kV=0.24 kVA=3.098 irradiance=1 Pmpp=2.582 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.676 phases=1 bus1=676 kV=0.24 kVA=1.7952 irradiance=1 Pmpp=1.496 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.679 phases=1 bus1=679 kV=0.24 kVA=3.428 irradiance=1 Pmpp=2.857 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.687 phases=1 bus1=687 kV=0.24 kVA=0.830 irradiance=1 Pmpp=0.692 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.690 phases=1 bus1=690 kV=0.24 kVA=0.954 irradiance=1 Pmpp=0.795 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.692 phases=1 bus1=692 kV=0.24 kVA=2.4792 irradiance=1 Pmpp=2.066 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.694 phases=1 bus1=694 kV=0.24 kVA=0.840 irradiance=1 Pmpp=0.7 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.696 phases=1 bus1=696 kV=0.24 kVA=1.1796 irradiance=1 Pmpp=0.983 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.698 phases=1 bus1=698 kV=0.24 kVA=1.614 irradiance=1 Pmpp=1.345 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.706 phases=1 bus1=706 kV=0.24 kVA=1.351 irradiance=1 Pmpp=1.126 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.709 phases=1 bus1=709 kV=0.24 kVA=1.076 irradiance=1 Pmpp=0.897 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.713 phases=1 bus1=713 kV=0.24 kVA=1.074 irradiance=1 Pmpp=0.895 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.713 phases=1 bus1=713 kV=0.24 kVA=0.775 irradiance=1 Pmpp=0.646 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.715 phases=1 bus1=715 kV=0.24 kVA=1.014 irradiance=1 Pmpp=0.845 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.717 phases=1 bus1=717 kV=0.24 kVA=1.142 irradiance=1 Pmpp=0.952 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.723 phases=1 bus1=723 kV=0.24 kVA=0.872 irradiance=1 Pmpp=0.727 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.726 phases=1 bus1=726 kV=0.24 kVA=1.680 irradiance=1 Pmpp=1.4 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.728 phases=1 bus1=728 kV=0.24 kVA=2.5812 irradiance=1 Pmpp=2.151 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.730 phases=1 bus1=730 kV=0.24 kVA=0.868 irradiance=1 Pmpp=0.723 kVAr = 0 %cutin=0.1 %cutout=0.1';


% 
%% defining PV with optimzation in opendss 

% DSSText.Command = 'New PVSystem.41 phases=1 bus1=41 kV=0.24 kVA=0.997 irradiance=1 Pmpp=0.831 kVAr = 0.5512 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.59 phases=1 bus1=59 kV=0.24 kVA=1.246 irradiance=1 Pmpp=1.038 kVAr = 0.2161 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.74 phases=1 bus1=74 kV=0.24 kVA=2.441 irradiance=1 Pmpp=2.034 kVAr = 1.3492 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.87 phases=1 bus1=87 kV=0.24 kVA=1.664 irradiance=1 Pmpp=1.387 kVAr = 0.4161 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.100 phases=1 bus1=100 kV=0.24 kVA=1.18 irradiance=1 Pmpp=0.983 kVAr = 0.1618 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.115 phases=1 bus1=115 kV=0.24 kVA=1.873 irradiance=1 Pmpp=1.561 kVAr = 1.0355 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.134 phases=1 bus1=134 kV=0.24 kVA=1.86 irradiance=1 Pmpp=1.55 kVAr = 1.0282 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.153 phases=1 bus1=153 kV=0.24 kVA=0.808 irradiance=1 Pmpp=0.673 kVAr = 0.4464 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.172 phases=1 bus1=172 kV=0.24 kVA=1.523 irradiance=1 Pmpp=1.269 kVAr = 0.8418 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.189 phases=1 bus1=189 kV=0.24 kVA=0.932 irradiance=1 Pmpp=0.777 kVAr = 0.5154 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.202 phases=1 bus1=202 kV=0.24 kVA=2.628 irradiance=1 Pmpp=2.19 kVAr = 1.0675 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.214 phases=1 bus1=214 kV=0.24 kVA=1.578 irradiance=1 Pmpp=1.315 kVAr = 0.164 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.225 phases=1 bus1=225 kV=0.24 kVA=3.719 irradiance=1 Pmpp=3.099 kVAr = 1.0063 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.238 phases=1 bus1=238 kV=0.24 kVA=3.197 irradiance=1 Pmpp=2.664 kVAr = 0.8878 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.253 phases=1 bus1=253 kV=0.24 kVA=1.693 irradiance=1 Pmpp=1.411 kVAr = 0.936 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.270 phases=1 bus1=270 kV=0.24 kVA=3.936 irradiance=1 Pmpp=3.28 kVAr = 2.1757 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.283 phases=1 bus1=283 kV=0.24 kVA=2.783 irradiance=1 Pmpp=2.319 kVAr = 1..5383 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.296 phases=1 bus1=296 kV=0.24 kVA=1.903 irradiance=1 Pmpp=1.586 kVAr = 1.051 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.310 phases=1 bus1=310 kV=0.24 kVA=2.104 irradiance=1 Pmpp=1.753 kVAr = 1.1628 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.326 phases=1 bus1=326 kV=0.24 kVA=2.364 irradiance=1 Pmpp=1.97 kVAr = 1.2865 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.342 phases=1 bus1=342 kV=0.24 kVA=3.071 irradiance=1 Pmpp=2.559 kVAr = 0.8526 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.358 phases=1 bus1=358 kV=0.24 kVA=3.028 irradiance=1 Pmpp=2.523 kVAr = 1.6736 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.374 phases=1 bus1=374 kV=0.24 kVA=3.378 irradiance=1 Pmpp=2.815 kVAr = 1.8673 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.390 phases=1 bus1=390 kV=0.24 kVA=3.166 irradiance=1 Pmpp=2.638 kVAr = 1.7499 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.407 phases=1 bus1=407 kV=0.24 kVA=2.605 irradiance=1 Pmpp=2.171 kVAr = 1.4401 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.426 phases=1 bus1=426 kV=0.24 kVA=1.172 irradiance=1 Pmpp=0.977 kVAr = 0.6481 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.445 phases=1 bus1=445 kV=0.24 kVA=1.799 irradiance=1 Pmpp=1.499 kVAr = 0.9943 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.464 phases=1 bus1=464 kV=0.24 kVA=1.615 irradiance=1 Pmpp=1.346 kVAr = 0.8503 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.483 phases=1 bus1=483 kV=0.24 kVA=1.487 irradiance=1 Pmpp=1.239 kVAr = 0.754 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.502 phases=1 bus1=502 kV=0.24 kVA=1.466 irradiance=1 Pmpp=1.222 kVAr = 0.8106 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.521 phases=1 bus1=521 kV=0.24 kVA=1.582 irradiance=1 Pmpp=1.318 kVAr = 0.8743 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.540 phases=1 bus1=540 kV=0.24 kVA=1.817 irradiance=1 Pmpp=1.514 kVAr = 0.0267 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.557 phases=1 bus1=557 kV=0.24 kVA=2.011 irradiance=1 Pmpp=1.676 kVAr = 1.117 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.570 phases=1 bus1=570 kV=0.24 kVA=1.597 irradiance=1 Pmpp=1.331 kVAr = 0.8735 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.583 phases=1 bus1=583 kV=0.24 kVA=1.801 irradiance=1 Pmpp=1.501 kVAr = 0.4513 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.598 phases=1 bus1=598 kV=0.24 kVA=1.848 irradiance=1 Pmpp=1.54 kVAr = 1.0002 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.617 phases=1 bus1=617 kV=0.24 kVA=0.881 irradiance=1 Pmpp=0.734 kVAr = 0.4733 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.636 phases=1 bus1=636 kV=0.24 kVA=2.848 irradiance=1 Pmpp=2.373 kVAr = 1.5741 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.654 phases=1 bus1=654 kV=0.24 kVA=1.589 irradiance=1 Pmpp=1.324 kVAr = 0.8175 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.670 phases=1 bus1=670 kV=0.24 kVA=1.447 irradiance=1 Pmpp=1.206 kVAr = 0.7946 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.687 phases=1 bus1=687 kV=0.24 kVA=0.830 irradiance=1 Pmpp=0.692 kVAr = 0.459 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.706 phases=1 bus1=706 kV=0.24 kVA=1.351 irradiance=1 Pmpp=1.126 kVAr = 0.7469 %cutin=0.1 %cutout=0.1';
% DSSText.Command = 'New PVSystem.723 phases=1 bus1=723 kV=0.24 kVA=0.872 irradiance=1 Pmpp=0.727 kVAr = 0.4822 %cutin=0.1 %cutout=0.1';


    
DSSCircuit      = DSSObj.ActiveCircuit;
DSSSolution     = DSSCircuit.Solution;
DSSSolution.Solve;                          %% solving the circuit without PV
Buses=DSSCircuit.Allbusnames;               %% knowing all the buses in the system
DSSCircuit.SetActiveElement('Vsource.SOURCE');                  %% power generated by source
Power_subsopend = -1*(DSSCircuit.ActivecktElement.Powers);
Vallnodes = DSSCircuit.AllBusVmagPU;
MyLosses = DSSCircuit.Losses;   

Buses1=str2double(Buses);



          a = 'DSS start properly';
else
    
    a = 'DSS did not start properly';
    disp(a)
end

