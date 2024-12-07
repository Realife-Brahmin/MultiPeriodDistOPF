import os
import re
from collections import defaultdict

def myprintln(verbose, message):
    if verbose:
        print(message)

def parse_branch_data(systemName,
    kVA_B=1000,
    kV_B=2.4018,
    Z_B=5.768643240000001,
    verbose=False):

    # Now let's take a small detour to ensure that we know exactly what base impedance to compute pu values of impedances:
    MVA_B = kVA_B / 1000

    # Rule: If Z_B is not provided and kVA_B is specified but not kV_B, throw an error
    if Z_B is None and kVA_B is not None and kV_B is None:
        raise ValueError("Error: You must specify both kV_B and kVA_B to calculate Z_B, or provide Z_B directly.")

    # Rule: If Z_B is not specified, and kV_B is provided, compute Z_B using the default or given kVA_B
    if Z_B is None and kV_B is not None:
        Z_B = (kV_B ** 2) / MVA_B
        myprintln(verbose, f"Computed Z_B = (kV_B^2) / MVA_B = {Z_B} using kV_B = {kV_B} and kVA_B = {kVA_B}")

    # Rule: If Z_B is provided and kV_B/kVA_B are not specified, just use the provided Z_B
    if Z_B is not None:
        myprintln(verbose, f"Using user-specified Z_B = {Z_B}")

    # Rule: If neither Z_B nor kV_B/kVA_B are provided, use the default value of Z_B
    if Z_B is None and kV_B is None and kVA_B is None:
        Z_B = 5.768643240000001  # Default value
        myprintln(verbose, f"Using default Z_B = {Z_B}")

    # Now Z_B is guaranteed to be set to a valid value at this point
    myprintln(verbose, f"Final Z_B = {Z_B}")

    # Todo: Ensure that substation bus being equal to 1 is not taken for granted, have some kwarg or something to ensure that even bus 153 can be the substation bus

    wd = os.path.dirname(os.path.abspath(__file__))
    # Construct the file path using wd
    filename = os.path.join(wd, "..", "..", "rawData", systemName, "BranchData.dss")

    # Initialize data structures
    Nset = set()
    Lset = set()
    rdict = {}
    xdict = {}

    rdict_pu = {}
    xdict_pu = {}

    parent = defaultdict(lambda: None)
    children = defaultdict(list)

    N1set = set()
    Nm1set = set()
    Nc1set = set()
    Nnc1set = set()
    L1set = set()
    Lm1set = set()

    # Regular expression to match key=value pairs with optional spaces
    kv_pattern = re.compile(r"(\w+)\s*=\s*([\S]+)")

    # Open and read the file
    with open(filename, "r") as file:
        for line in file:
            # Remove comments and strip whitespace
            line = line.split("!")[0].strip()

            # Skip empty lines
            if not line:
                continue

            # Parse lines starting with "New Line."
            if line.startswith("New Line."):
                branch_info = {}
                for match in kv_pattern.finditer(line):
                    key, value = match.groups()
                    branch_info[key.strip()] = value.strip()

                # Extract from_bus and to_bus
                if "Bus1" in branch_info and "Bus2" in branch_info:
                    bus1_str = branch_info["Bus1"]
                    bus2_str = branch_info["Bus2"]

                    bus1_parts = bus1_str.split(".")
                    bus2_parts = bus2_str.split(".")

                    from_bus = int(bus1_parts[0])
                    to_bus = int(bus2_parts[0])

                    Nset.update([from_bus, to_bus])
                    Lset.add((from_bus, to_bus))

                    parent[to_bus] = from_bus
                    children[from_bus].append(to_bus)

                    # Ensure to_bus is a key in the children dictionary
                    if to_bus not in children:
                        children[to_bus] = []

                    if from_bus == 1:
                        N1set.add(from_bus)
                        Nc1set.add(to_bus)
                        L1set.add((from_bus, to_bus))
                    elif to_bus == 1:
                        N1set.add(from_bus)
                        Nc1set.add(from_bus)
                        L1set.add((from_bus, to_bus))
                    else:
                        Nm1set.update([from_bus, to_bus])
                        Lm1set.add((from_bus, to_bus))

                    Nm1set = Nset - N1set
                    Nnc1set = Nm1set - Nc1set
                else:
                    raise ValueError("Bus1 or Bus2 not specified for a line in BranchData.dss")

                rdict[(from_bus, to_bus)] = float(branch_info.get("r1", 0.0))
                xdict[(from_bus, to_bus)] = float(branch_info.get("x1", 0.0))

                rdict_pu[(from_bus, to_bus)] = rdict[(from_bus, to_bus)] / Z_B
                xdict_pu[(from_bus, to_bus)] = xdict[(from_bus, to_bus)] / Z_B

    N = len(Nset)
    m = len(Lset)

    N1 = len(N1set)
    Nm1 = len(Nm1set)
    Nc1 = len(Nc1set)
    Nnc1 = len(Nnc1set)
    m1 = len(L1set)
    mm1 = len(Lm1set)

    Nset = sorted(Nset)
    Lset = sorted(Lset)
    N1set = sorted(N1set)
    Nm1set = sorted(Nm1set)
    Nc1set = sorted(Nc1set)
    Nnc1set = sorted(Nnc1set)
    L1set = sorted(L1set)
    Lm1set = sorted(Lm1set)

    branchData = {
        "Nset": Nset,
        "Lset": Lset,
        "rdict": rdict,
        "xdict": xdict,
        "rdict_pu": rdict_pu,
        "xdict_pu": xdict_pu,
        "parent": parent,
        "children": children,
        "N": N,
        "m": m,
        "N1set": N1set,
        "Nm1set": Nm1set,
        "Nc1set": Nc1set,
        "Nnc1set": Nnc1set,
        "L1set": L1set,
        "Lm1set": Lm1set,
        "N1": N1,
        "Nm1": Nm1,
        "Nc1": Nc1,
        "Nnc1": Nnc1,
        "m1": m1,
        "mm1": mm1,
    }

    return branchData
