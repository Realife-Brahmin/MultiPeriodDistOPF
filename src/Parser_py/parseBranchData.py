import os
import re
from collections import defaultdict

# Assuming `myprintln` and `generateBinaryLoadShape` are in a helper module
from src.helperFunctions import myprintln

def parse_branch_data(system_name, kVA_B=1000, kV_B=2.4018, Z_B=5.768643240000001, verbose=False):
    """
    Parses branch data from a `.dss` file for the specified system.
    """
    
    # Calculate MVA_B
    MVA_B = kVA_B / 1000

    # Rule: If Z_B is not provided and kVA_B is specified but not kV_B, raise an error
    if Z_B is None and kVA_B is not None and kV_B is None:
        raise ValueError("Error: You must specify both kV_B and kVA_B to calculate Z_B, or provide Z_B directly.")

    # Rule: If Z_B is not specified, compute it using kV_B
    if Z_B is None and kV_B is not None:
        Z_B = (kV_B ** 2) / MVA_B
        myprintln(verbose, f"Computed Z_B = (kV_B^2) / MVA_B = {Z_B} using kV_B = {kV_B} and kVA_B = {kVA_B}")

    # Rule: Use provided Z_B if specified
    if Z_B is not None:
        myprintln(verbose, f"Using user-specified Z_B = {Z_B}")

    # Rule: Use default value of Z_B if none provided
    if Z_B is None and kV_B is None and kVA_B is None:
        Z_B = 5.768643240000001
        myprintln(verbose, f"Using default Z_B = {Z_B}")

    # Final Z_B confirmation
    myprintln(verbose, f"Final Z_B = {Z_B}")

    # Construct file path
    wd = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(wd, "..", "..", "rawData", system_name, "BranchData.dss")

    # Initialize data structures
    Nset = set()
    Lset = set()
    rdict, xdict = {}, {}
    rdict_pu, xdict_pu = {}, {}
    parent = defaultdict(lambda: None)
    children = defaultdict(list)

    # Initialize additional sets and parameters
    N1set = set()
    Nm1set = set()
    Nc1set = set()
    Nnc1set = set()
    L1set = set()
    Lm1set = set()

    # Regular expression to match key=value pairs
    kv_pattern = re.compile(r"(\w+)\s*=\s*([\S]+)")

    # Open and parse the file
    with open(filename, "r") as file:
        for line in file:
            line = line.split("!")[0].strip()  # Remove comments and whitespace
            if not line:
                continue

            if line.startswith("New Line."):
                branch_info = {}
                for match in kv_pattern.finditer(line):
                    key, value = match.groups()
                    branch_info[key] = value

                # Extract buses
                if "Bus1" in branch_info and "Bus2" in branch_info:
                    from_bus = int(branch_info["Bus1"].split(".")[0])
                    to_bus = int(branch_info["Bus2"].split(".")[0])

                    # Update sets and dictionaries
                    Nset.update([from_bus, to_bus])
                    Lset.add((from_bus, to_bus))
                    parent[to_bus] = from_bus
                    children[from_bus].append(to_bus)
                    parent.setdefault(from_bus, None)
                    children.setdefault(to_bus, [])

                    # Update specific sets based on substation node (assumed to be 1)
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

                    # Nm1set: all nodes except the substation
                    Nm1set = Nset - N1set

                    # Nnc1set: nodes not connected to the substation
                    Nnc1set = Nm1set - Nc1set
                else:
                    raise ValueError("Bus1 or Bus2 not specified for a line in BranchData.dss")

                # Extract resistance and reactance values
                rdict[(from_bus, to_bus)] = float(branch_info.get("r1", 0.0))
                xdict[(from_bus, to_bus)] = float(branch_info.get("x1", 0.0))

                # Calculate per-unit values
                rdict_pu[(from_bus, to_bus)] = rdict[(from_bus, to_bus)] / Z_B
                xdict_pu[(from_bus, to_bus)] = xdict[(from_bus, to_bus)] / Z_B

    # Summary statistics
    N = len(Nset)
    m = len(Lset)
    N1 = len(N1set)
    Nm1 = len(Nm1set)
    Nc1 = len(Nc1set)
    Nnc1 = len(Nnc1set)
    m1 = len(L1set)
    mm1 = len(Lm1set)

    # Sort lists for deterministic order
    sorted_data = {
        "Nset": sorted(Nset),
        "Lset": sorted(Lset),
        "rdict": rdict,
        "xdict": xdict,
        "rdict_pu": rdict_pu,
        "xdict_pu": xdict_pu,
        "parent": dict(parent),
        "children": dict(children),
        "N": N,
        "m": m,
        "N1set": sorted(N1set),
        "Nm1set": sorted(Nm1set),
        "Nc1set": sorted(Nc1set),
        "Nnc1set": sorted(Nnc1set),
        "L1set": sorted(L1set),
        "Lm1set": sorted(Lm1set),
        "N1": N1,
        "Nm1": Nm1,
        "Nc1": Nc1,
        "Nnc1": Nnc1,
        "m1": m1,
        "mm1": mm1
    }

    return sorted_data
