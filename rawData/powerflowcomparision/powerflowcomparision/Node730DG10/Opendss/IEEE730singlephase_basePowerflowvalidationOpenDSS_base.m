clear all;
clc;

[DSSStartOK, DSSObj, DSSText] = DSSStartup;

if DSSStartOK
    

   
DSSObj.AllowForms = false;    


DSSText.Command = 'Set DataPath = C:\Users\rjha\Dropbox\UNCC-MTU-WSU\From WSU\730nodesystem\powerflowcomparision\Basecase\OpenDSSmodelbase';

DSSText.Command = 'compile singlephase730.dss';          
    
%% defining PV in opendss 
DSSText.Command = 'New PVSystem.41 phases=1 bus1=41 kV=0.24 kVA=0.997 irradiance=1 Pmpp=0.831 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.59 phases=1 bus1=59 kV=0.24 kVA=1.246 irradiance=1 Pmpp=1.038 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.74 phases=1 bus1=74 kV=0.24 kVA=2.441 irradiance=1 Pmpp=2.034 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.87 phases=1 bus1=87 kV=0.24 kVA=1.664 irradiance=1 Pmpp=1.387 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.100 phases=1 bus1=100 kV=0.24 kVA=1.18 irradiance=1 Pmpp=0.983 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.115 phases=1 bus1=115 kV=0.24 kVA=1.873 irradiance=1 Pmpp=1.561 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.134 phases=1 bus1=134 kV=0.24 kVA=1.86 irradiance=1 Pmpp=1.55 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.153 phases=1 bus1=153 kV=0.24 kVA=0.808 irradiance=1 Pmpp=0.673 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.172 phases=1 bus1=172 kV=0.24 kVA=1.523 irradiance=1 Pmpp=1.269 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.189 phases=1 bus1=189 kV=0.24 kVA=0.932 irradiance=1 Pmpp=0.777 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.202 phases=1 bus1=202 kV=0.24 kVA=2.628 irradiance=1 Pmpp=2.19 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.214 phases=1 bus1=214 kV=0.24 kVA=1.578 irradiance=1 Pmpp=1.315 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.225 phases=1 bus1=225 kV=0.24 kVA=3.719 irradiance=1 Pmpp=3.099 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.238 phases=1 bus1=238 kV=0.24 kVA=3.197 irradiance=1 Pmpp=2.664 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.253 phases=1 bus1=253 kV=0.24 kVA=1.693 irradiance=1 Pmpp=1.411 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.270 phases=1 bus1=270 kV=0.24 kVA=3.936 irradiance=1 Pmpp=3.28 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.283 phases=1 bus1=283 kV=0.24 kVA=2.783 irradiance=1 Pmpp=2.319 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.296 phases=1 bus1=296 kV=0.24 kVA=1.903 irradiance=1 Pmpp=1.586 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.310 phases=1 bus1=310 kV=0.24 kVA=2.104 irradiance=1 Pmpp=1.753 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.326 phases=1 bus1=326 kV=0.24 kVA=2.364 irradiance=1 Pmpp=1.97 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.342 phases=1 bus1=342 kV=0.24 kVA=3.071 irradiance=1 Pmpp=2.559 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.358 phases=1 bus1=358 kV=0.24 kVA=3.028 irradiance=1 Pmpp=2.523 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.374 phases=1 bus1=374 kV=0.24 kVA=3.378 irradiance=1 Pmpp=2.815 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.390 phases=1 bus1=390 kV=0.24 kVA=3.166 irradiance=1 Pmpp=2.638 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.407 phases=1 bus1=407 kV=0.24 kVA=2.605 irradiance=1 Pmpp=2.171 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.426 phases=1 bus1=426 kV=0.24 kVA=1.172 irradiance=1 Pmpp=0.977 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.445 phases=1 bus1=445 kV=0.24 kVA=1.799 irradiance=1 Pmpp=1.499 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.464 phases=1 bus1=464 kV=0.24 kVA=1.615 irradiance=1 Pmpp=1.346 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.483 phases=1 bus1=483 kV=0.24 kVA=1.487 irradiance=1 Pmpp=1.239 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.502 phases=1 bus1=502 kV=0.24 kVA=1.466 irradiance=1 Pmpp=1.222 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.521 phases=1 bus1=521 kV=0.24 kVA=1.582 irradiance=1 Pmpp=1.318 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.540 phases=1 bus1=540 kV=0.24 kVA=1.817 irradiance=1 Pmpp=1.514 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.557 phases=1 bus1=557 kV=0.24 kVA=2.011 irradiance=1 Pmpp=1.676 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.570 phases=1 bus1=570 kV=0.24 kVA=1.597 irradiance=1 Pmpp=1.331 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.583 phases=1 bus1=583 kV=0.24 kVA=1.801 irradiance=1 Pmpp=1.501 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.598 phases=1 bus1=598 kV=0.24 kVA=1.848 irradiance=1 Pmpp=1.54 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.617 phases=1 bus1=617 kV=0.24 kVA=0.881 irradiance=1 Pmpp=0.734 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.636 phases=1 bus1=636 kV=0.24 kVA=2.848 irradiance=1 Pmpp=2.373 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.654 phases=1 bus1=654 kV=0.24 kVA=1.589 irradiance=1 Pmpp=1.324 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.670 phases=1 bus1=670 kV=0.24 kVA=1.447 irradiance=1 Pmpp=1.206 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.687 phases=1 bus1=687 kV=0.24 kVA=0.830 irradiance=1 Pmpp=0.692 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.706 phases=1 bus1=706 kV=0.24 kVA=1.351 irradiance=1 Pmpp=1.126 kVAr = 0 %cutin=0.1 %cutout=0.1';
DSSText.Command = 'New PVSystem.723 phases=1 bus1=723 kV=0.24 kVA=0.872 irradiance=1 Pmpp=0.727 kVAr = 0 %cutin=0.1 %cutout=0.1';

%%

DSSCircuit      = DSSObj.ActiveCircuit;
DSSSolution     = DSSCircuit.Solution;
DSSSolution.Solve;                          %% solving the circuit without PV
Buses=DSSCircuit.Allbusnames;               %% knowing all the buses in the system
DSSCircuit.SetActiveElement('Vsource.SOURCE');                  %% power generated by source
Power_subsopend = -1*(DSSCircuit.ActivecktElement.Powers);
Vallnodes = DSSCircuit.AllBusVmagPU;
    

Buses1=str2double(Buses);

%% Voltage as per ascending node number

% bus_voltage = [busname Vallnodes'];         %% busname and corresponding volatge
% Vasced = sort(busname);                     %% sorting the bus name in ascending order
% VoltOpenDSS = [];
% for i = 1:size(Vasced,1)
% indx = find(i ==bus_voltage(:,1));
% Vo = bus_voltage(indx,2);
% VoltOpenDSS = [VoltOpenDSS Vo];
% end


          a = 'DSS start properly';
else
    
    a = 'DSS did not start properly';
    disp(a)
end

