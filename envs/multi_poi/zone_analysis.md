# Zone Assignment for IEEE 123-Bus 5-POI System

## Substation-Zone Mapping

Each distribution bus is assigned to a zone based on which substation feeds it (determined by tracing parent connections):

| Zone | Substation | Load Pattern | Characteristics |
|------|------------|--------------|-----------------|
| 1 | 1s (SW) | Smooth Sinusoid | Residential - peak midday, Range: [0.75, 1.00] |
| 2 | 2s (NW) | Smooth Sinusoid | Commercial - peak afternoon, Range: [0.75, 1.00] |
| 3 | 3s (N)  | Multi-level Square | Industrial with 3 shifts (morning/afternoon/evening) |
| 4 | 4s (NE) | Two-level Square | Office building - on/off (business hours 8-18) |
| 5 | 5s (SE) | Ramping | Charging station - gradual ramp up then down |

## Resource Distribution by Zone

### PV Systems (17 total)
```
Zone 1 (1s): Buses [3, 9, 18, 26, 35]         = 5 PV units
Zone 2 (2s): Buses [41, 48, 53]                = 3 PV units  
Zone 3 (3s): Buses [61, 67, 73, 79]            = 4 PV units
Zone 4 (4s): Buses [87, 94, 101, 108, 116]     = 5 PV units
Zone 5 (5s): Buses []                          = 0 PV units
```

### Battery Systems (26 total)
```
Zone 1 (1s): Buses [3, 7, 11, 18, 22, 30, 34, 37]         = 8 batteries
Zone 2 (2s): Buses [41, 47, 50, 53]                       = 4 batteries
Zone 3 (3s): Buses [58, 62, 67, 71, 75, 79]               = 6 batteries
Zone 4 (4s): Buses [84, 87, 92, 97, 101, 106, 111, 116]   = 8 batteries
Zone 5 (5s): Buses []                                      = 0 batteries
```

## Explanation of Power Dispatch Fluctuations

### Why Subs 1 has huge fluctuations:
- **Zone 1 (Subs 1s)** has **5 PV units** and **8 batteries**
- Smooth sinusoidal load pattern with peak midday
- PV generation peaks during day (hours 7-18) → large power EXPORT
- Night time (hours 1-6, 19-24) → large power IMPORT for load
- Battery dispatch actively balances this → extreme swings

### Why Subs 2 and 3 have moderate fluctuations:
- **Zone 2 (Subs 2s)**: 3 PV, 4 batteries - moderate resources
- **Zone 3 (Subs 3s)**: 4 PV, 6 batteries - moderate resources  
- Load patterns (sinusoid + industrial shifts) combined with PV create moderate swings

### Why Subs 4 and 5 have minimal fluctuations:
- **Zone 4 (Subs 4s)**: Despite having 5 PV and 8 batteries, the two-level square load pattern (office building) is relatively flat during hours 8-18 when PV is active → PV directly serves local load → minimal net power exchange
- **Zone 5 (Subs 5s)**: **NO PV and NO batteries** → purely load-driven with ramping pattern → very stable power import

## Key Insight
The fluctuations in substation power dispatch are primarily driven by:
1. **PV generation patterns** (day/night cycle)
2. **Battery charging/discharging** to optimize costs
3. **Mismatch between local load and generation** in each zone

Zones with high renewable penetration but mismatched load profiles (Zone 1) show the largest swings.
Zones with no DERs (Zone 5) or well-matched load-generation (Zone 4) show stable dispatch.
