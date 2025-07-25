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
r2_ohm = 10.00
x2_ohm = 15.00
P_L_kW = 500
Q_L_kW = 500*0.75
kVA_B = 1000
kV_B = 20.00
C_1_dollar_per_kWh = 0.50
C_2_dollar_per_kWh = 0.25

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
    :kV_B => kV_B,
    :C_1_dollar_per_kWh => C_1_dollar_per_kWh,
    :C_2_dollar_per_kWh => C_2_dollar_per_kWh
)

function process_data!(data)
    # Base values
    S_base = data[:kVA_B]
    V_base = data[:kV_B]

    # Per-unit conversions
    P_L_pu = data[:P_L_kW] / S_base
    Q_L_pu = data[:Q_L_kW] / S_base
    r1_pu = data[:r1_ohm] / ((V_base^2) / S_base)
    x1_pu = data[:x1_ohm] / ((V_base^2) / S_base)
    r2_pu = data[:r2_ohm] / ((V_base^2) / S_base)
    x2_pu = data[:x2_ohm] / ((V_base^2) / S_base)

    # Cost per unit energy (pu)
    C_1_dollar_pu = data[:C_1_dollar_per_kWh] / S_base
    C_2_dollar_pu = data[:C_2_dollar_per_kWh] / S_base

    # Angles in radians
    delta_1_rad = data[:delta_1_deg] * pi / 180
    delta_2_rad = data[:delta_2_deg] * pi / 180

    # Voltage limits
    Vminpu = 0.75
    Vmaxpu = 1.15

    # Update data dict with all required fields
    data[:P_L_pu] = P_L_pu
    data[:Q_L_pu] = Q_L_pu
    data[:r1_pu] = r1_pu
    data[:x1_pu] = x1_pu
    data[:r2_pu] = r2_pu
    data[:x2_pu] = x2_pu
    data[:C_1_dollar_pu] = C_1_dollar_pu
    data[:C_2_dollar_pu] = C_2_dollar_pu
    data[:delta_1_rad] = delta_1_rad
    data[:delta_2_rad] = delta_2_rad
    data[:Vminpu] = Vminpu
    data[:Vmaxpu] = Vmaxpu

end

process_data!(data);