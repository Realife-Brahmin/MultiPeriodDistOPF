# evaluate_voltage_limits.py

def evaluate_voltage_limits(load_data, pv_data, battery_data):
    """
    Evaluates voltage limits for a set of components by calculating the most restrictive
    voltage limits based on load, PV, and battery data.
    """

    # Create Compset as the union of the sets from load, PV, and battery data
    Compset = set(load_data.get("NLset", set())) | set(pv_data.get("Dset", set())) | set(battery_data.get("Bset", set()))

    # Initialize dictionaries to store voltage limits
    Vminpu_Comp = {}
    Vmaxpu_Comp = {}

    # Iterate over each bus in Compset to calculate voltage limits
    for j in Compset:
        # Get Vminpu and Vmaxpu for load, PV, and battery components if they exist
        Vminpu_L = load_data.get("Vminpu_L", {}).get(j, float("-inf"))
        Vmaxpu_L = load_data.get("Vmaxpu_L", {}).get(j, float("inf"))

        Vminpu_D = pv_data.get("Vminpu_D", {}).get(j, float("-inf"))
        Vmaxpu_D = pv_data.get("Vmaxpu_D", {}).get(j, float("inf"))

        Vminpu_B = battery_data.get("Vminpu_B", {}).get(j, float("-inf"))
        Vmaxpu_B = battery_data.get("Vmaxpu_B", {}).get(j, float("inf"))

        # Calculate the most restrictive voltage limits for this bus
        Vminpu_Comp[j] = max(Vminpu_L, Vminpu_D, Vminpu_B)  # Choose the highest lower limit
        Vmaxpu_Comp[j] = min(Vmaxpu_L, Vmaxpu_D, Vmaxpu_B)  # Choose the lowest upper limit

    # Create the output dictionary
    component_data = {
        "Compset": Compset,
        "Vminpu_Comp": Vminpu_Comp,
        "Vmaxpu_Comp": Vmaxpu_Comp
    }

    return component_data
