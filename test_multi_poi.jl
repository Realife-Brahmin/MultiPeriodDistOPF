# multi_poi.jl
include("./src/setupMultiPeriodDistOPF.jl") 

data = Dict()
V_1_pu = 1.07
delta_1_deg = 0.0
V_2_pu = 1.06
delta_2_deg = -10.0
alpha_share = 0.5
r1_ohm = 15.00
x1_ohm = 10.00
r2_ohm = 15.00
x2_ohm = 10.00
P_L_kW = 500
Q_L_kW = 500*0.75
kVA_B = 1000
kV_B = 20.00

data = Dict(
    :V_1_pu => V_1_pu,
    :delta_1_deg => delta_1_deg,
    :V_2_pu => V_2_pu,
    :delta_2_deg => delta_2_deg,
    :alpha_share => alpha_share,
    :r1_ohm => r1_ohm,
    :x1_ohm => x1_ohm,
    :r2_ohm => r2_ohm,
    :x2_ohm => x2_ohm,
    :P_L_kW => P_L_kW,
    :Q_L_kW => Q_L_kW,
    :kVA_B => kVA_B,
    :kV_B => kV_B
)