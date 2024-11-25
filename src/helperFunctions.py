# # helper_functions.py

# import os

# def myprintln(verbose, msg):
#     """
#     Prints the message only if verbose is True.
#     """
#     if verbose:
#         print(msg)


# def generateBinaryLoadShape(T, filename_load_shape="LoadShapePSubsCostDefault.dss",
#                                hi=None, lo=None, peak_hours_fraction=0.3, verbose=False):
#     """
#     Generates a binary load shape for a given period `T`, based on a file or default peak/off-peak values.
#     """
#     # Get the directory path of this script
#     wd = os.path.dirname(os.path.abspath(__file__))
#     filepath = os.path.join(wd, "..", "rawData", filename_load_shape)

#     myprintln(verbose, f"Attempting to read LoadShape from: {filepath}")

#     # Initialize an empty list for prices
#     prices = []

#     # Read the LoadShape file and extract the prices
#     with open(filepath, "r") as file:
#         myprintln(verbose, "Reading file content:")
#         for line in file:
#             line = line.strip()
#             myprintln(verbose, f"Read line: {line}")

#             if not line or line.startswith("!"):
#                 myprintln(verbose, "Skipping comment or empty line")
#                 continue

#             if "mult=" in line:
#                 price_str = line.split("mult=[")[1].split("]")[0]
#                 prices = [float(p) for p in price_str.split()]
#                 myprintln(verbose, f"Parsed prices: {prices}")
#             else:
#                 myprintln(verbose, "Skipping line, does not contain mult=")

#     if not prices:
#         raise ValueError(f"No data found in {filename_load_shape}. The prices array is empty.")

#     # Determine hi and lo based on user input or defaults to the max/min of prices
#     hi = max(prices) if hi is None else hi
#     lo = min(prices) if lo is None else lo
#     myprintln(verbose, f"Computed hi: {hi}, lo: {lo}")

#     # Adjust prices to match period T through subsampling or supersampling
#     if T < 24:
#         start_idx = (24 - T) // 2
#         lambda_vals0 = prices[start_idx:start_idx + T]
#         myprintln(verbose, f"Subsampling: lambdaVals0 = {lambda_vals0}")
#     elif T > 24:
#         lambda_vals0 = []
#         for i in range(T):
#             scaled_idx = (i * (24 - 1)) / (T - 1)
#             lower_idx = int(scaled_idx)
#             upper_idx = min(lower_idx + 1, 23)
#             frac = scaled_idx - lower_idx
#             interpolated_value = (1 - frac) * prices[lower_idx] + frac * prices[upper_idx]
#             lambda_vals0.append(interpolated_value)
#         myprintln(verbose, f"Supersampled: lambdaVals0 = {lambda_vals0}")
#     else:
#         lambda_vals0 = prices
#         myprintln(verbose, f"T == 24, using original prices: lambdaVals0 = {lambda_vals0}")

#     # Calculate the number of peak hours
#     num_peak_hours = max(1, int(peak_hours_fraction * T))
#     myprintln(verbose, f"Number of peak hours: {num_peak_hours}")

#     # Initialize the cost array to the low value
#     cost_array = [lo] * T
#     myprintln(verbose, f"Initialized costArray with low values: {cost_array}")

#     # Assign the high value to the peak hours
#     sorted_indices = sorted(range(T), key=lambda i: lambda_vals0[i], reverse=True)
#     for idx in sorted_indices[:num_peak_hours]:
#         cost_array[idx] = hi
#     myprintln(verbose, f"Updated costArray with peak values: {cost_array}")

#     # Create the output dictionary
#     cost_data = {
#         "LoadShapeCost": cost_array,
#         "peakCost": hi,
#         "offPeakCost": lo,
#         "peakHoursFraction": peak_hours_fraction
#     }

#     return cost_data


# def generateLoadShape(T, filename_load_shape=None):
#     """
#     Generates a load shape based on a specified file or default values.
#     """
#     wd = os.path.dirname(os.path.abspath(__file__))

#     # Default load shape file if none is provided
#     if filename_load_shape is None:
#         filename_load_shape = "LoadShapePVDefault.dss"
    
#     loadshape_filepath = os.path.join(wd, "..", "rawData", filename_load_shape)

#     # Read the load shape file
#     default_load_shape = []
#     with open(loadshape_filepath, "r") as file:
#         for line in file:
#             line = line.split("!")[0].strip()
#             if not line:
#                 continue
#             if "mult=[" in line:
#                 shape_str = line.split("mult=[")[1].split("]")[0]
#                 default_load_shape = [float(val) for val in shape_str.split()]

#     # Adjust load shape values to match period T through subsampling or supersampling
#     if T < 24:
#         start_idx = (24 - T) // 2
#         load_shape = default_load_shape[start_idx:start_idx + T]
#     elif T > 24:
#         load_shape = []
#         for i in range(T):
#             scaled_idx = (i * (24 - 1)) / (T - 1)
#             lower_idx = int(scaled_idx)
#             upper_idx = min(lower_idx + 1, 23)
#             frac = scaled_idx - lower_idx
#             interpolated_value = (1 - frac) * default_load_shape[lower_idx] + frac * default_load_shape[upper_idx]
#             load_shape.append(interpolated_value)
#     else:
#         load_shape = default_load_shape

#     return load_shape

## new helper functions

import os

def myprintln(verbose, msg):
    if verbose:
        print(msg)

def generate_binary_load_shape(T, filename_load_shape="LoadShapePSubsCostDefault.dss",
                               hi=None, lo=None, peak_hours_fraction=0.3,
                               verbose=False):
    # Get the working directory of this script
    wd = os.path.dirname(os.path.abspath(__file__))

    # Construct the full file path using the filename
    filepath = os.path.join(wd, "..", "rawData", filename_load_shape)

    # Print the file path to confirm it's the correct file
    myprintln(verbose, f"Attempting to read LoadShape from: {filepath}")

    # Initialize a list for prices
    prices = []

    # Read the LoadShapePSubsCostDefault.dss file
    myprintln(verbose, "Reading file content:")
    with open(filepath, "r") as file:
        for line in file:
            line = line.strip()
            # Print each line read from the file
            myprintln(verbose, f"Read line: {line}")
            # Skip empty lines and comments
            if not line or line.startswith("!"):
                myprintln(verbose, "Skipping comment or empty line")
                continue
            # Check if the line contains the 'mult=' string
            if "mult=" in line:
                # Extract the values inside the brackets
                price_str = line.split("mult=")[-1].strip("[]")
                # Parse the prices from the string
                prices = [float(p) for p in price_str.split()]
                myprintln(verbose, f"Parsed prices: {prices}")
            else:
                myprintln(verbose, "Skipping line, does not contain mult=")

    # Check if prices is empty
    if not prices:
        raise ValueError(f"No data found in {filename_load_shape}. The prices array is empty.")

    # Compute hi and lo based on the user input or default to the max/min of prices
    hi = max(prices) if hi is None else hi
    lo = min(prices) if lo is None else lo

    myprintln(verbose, f"Computed hi: {hi}, lo: {lo}")

    # Subsampling or supersampling the 24-step price array to get lambdaVals0
    if T < 24:
        # Subsample the middle T points
        start_idx = (24 - T) // 2
        lambda_vals0 = prices[start_idx:(start_idx + T)]
        myprintln(verbose, f"Subsampling: lambdaVals0 = {lambda_vals0}")
    elif T > 24:
        # Supersample using interpolation to create T values from the 24-step array
        lambda_vals0 = []
        myprintln(verbose, f"Supersampling to {T} values:")
        for i in range(T):
            scaled_idx = i * (24 - 1) / (T - 1)
            lower_idx = int(scaled_idx)
            upper_idx = min(lower_idx + 1, 23)
            frac = scaled_idx - lower_idx
            interpolated_value = (1 - frac) * prices[lower_idx] + frac * prices[upper_idx]
            lambda_vals0.append(interpolated_value)
        myprintln(verbose, f"Supersampled: lambdaVals0 = {lambda_vals0}")
    else:
        # If T == 24, use the original prices
        lambda_vals0 = prices
        myprintln(verbose, f"T == 24, using original prices: lambdaVals0 = {lambda_vals0}")

    # Now that we have lambdaVals0, compute the number of peak hours
    num_peak_hours = max(1, int(peak_hours_fraction * T))

    myprintln(verbose, f"Number of peak hours: {num_peak_hours}")

    # Initialize the cost array to the low value (lo)
    cost_array = [lo] * T
    myprintln(verbose, f"Initialized costArray with low values: {cost_array}")

    # Sort lambdaVals0 in descending order and get the sorted indices
    sorted_indices = sorted(range(T), key=lambda x: lambda_vals0[x], reverse=True)
    myprintln(verbose, f"Sorted indices for peak hours: {sorted_indices}")

    # Assign the high value (hi) to the peak hours
    for idx in sorted_indices[:num_peak_hours]:
        cost_array[idx] = hi
    myprintln(verbose, f"Updated costArray with peak values: {cost_array}")

    # Create the output dictionary
    cost_data = {
        "LoadShapeCost": cost_array,
        "peakCost": hi,
        "offPeakCost": lo,
        "peakHoursFraction": peak_hours_fraction
    }

    return cost_data

def generate_load_shape(T, filename_load_shape=None):
    wd = os.path.dirname(os.path.abspath(__file__))

    default_load_shape = []
    if filename_load_shape is None:
        # Default to LoadShapePVDefault.dss if no filename_load_shape is provided
        filename_load_shape = "LoadShapePVDefault.dss"

    # Construct the full file path using the provided filename_load_shape
    loadshape_filepath = os.path.join(wd, "..", "rawData", filename_load_shape)

    # Read the LoadShape file
    with open(loadshape_filepath, "r") as file:
        for line in file:
            line = line.split("!")[0].strip()  # Remove comments
            if not line:
                continue
            # Parse the load shape values
            if "mult=[" in line:
                shape_str = line.split("mult=")[-1].strip("[]")
                default_load_shape = [float(val) for val in shape_str.split()]

    # Generate load shape values based on T
    if T < 24:
        # Take the middle T values
        start_idx = (24 - T) // 2
        load_shape = default_load_shape[start_idx:(start_idx + T)]
    elif T > 24:
        # Supersample the values
        load_shape = []
        for i in range(T):
            scaled_idx = i * (24 - 1) / (T - 1)
            lower_idx = int(scaled_idx)
            upper_idx = min(lower_idx + 1, 23)
            frac = scaled_idx - lower_idx
            interpolated_value = (1 - frac) * default_load_shape[lower_idx] + frac * default_load_shape[upper_idx]
            load_shape.append(interpolated_value)
    else:
        # T == 24, use the default load shape
        load_shape = default_load_shape

    load_shape_load = [float(x) for x in load_shape]  # Ensure all values are floats
    return load_shape_load
