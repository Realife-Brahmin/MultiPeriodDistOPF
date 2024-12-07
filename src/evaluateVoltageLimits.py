def evaluate_voltage_limits(load_data, pv_data, battery_data):
    """
    Returns:
        dict: A dictionary with the component set and evaluated voltage limits.
    """
    # Create Compset as the union of the sets from the three components
    Compset = set(load_data["NLset"]).union(pv_data["Dset"], battery_data["Bset"])

    # Initialize dictionaries to store Vminpu_Comp and Vmaxpu_Comp
    Vminpu_Comp = {}
    Vmaxpu_Comp = {}

    # Iterate over each bus in Compset to calculate voltage limits
    for j in Compset:
        # Get Vminpu and Vmaxpu for load, PV, and battery components if they exist
        Vminpu_L = load_data["Vminpu_L"].get(j, float("-inf"))  # Default to -inf if not in set
        Vmaxpu_L = load_data["Vmaxpu_L"].get(j, float("inf"))   # Default to inf if not in set

        Vminpu_D = pv_data["Vminpu_D"].get(j, float("-inf"))
        Vmaxpu_D = pv_data["Vmaxpu_D"].get(j, float("inf"))

        Vminpu_B = battery_data["Vminpu_B"].get(j, float("-inf"))
        Vmaxpu_B = battery_data["Vmaxpu_B"].get(j, float("inf"))

        # Calculate the most restrictive voltage limits for this bus
        Vminpu_Comp[j] = max(Vminpu_L, Vminpu_D, Vminpu_B)  # Choose the highest lower limit
        Vmaxpu_Comp[j] = min(Vmaxpu_L, Vmaxpu_D, Vmaxpu_B)  # Choose the lowest upper limit

    # Create the output dictionary
    component_data = {
        "Compset": Compset,
        "Vminpu_Comp": Vminpu_Comp,
        "Vmaxpu_Comp": Vmaxpu_Comp,
    }

    return component_data
