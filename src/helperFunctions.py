import os
import re
import math
import numpy as np

def myprintln(verbose: bool, msg: str):
    if verbose:
        print(msg)

def generateBinaryLoadShape(T: int, filenameLoadShape="LoadShapePSubsCostDefault.dss", hi=None, lo=None, peakHoursFraction=0.3, verbose=False):
    # Get the working directory of this script
    wd = os.path.dirname(__file__)

    # Construct the full file path using the filename
    filepath = os.path.join(wd, "..", "rawData", filenameLoadShape)

    # Print the file path to confirm it's the correct file
    myprintln(verbose, f"Attempting to read LoadShape from: {filepath}")

    # Initialize a list for prices
    prices = []

    # Read the LoadShapePSubsCostDefault.dss file
    myprintln(verbose, "Reading file content:")
    with open(filepath, "r") as file:
        for line in file:
            line = line.strip()
            myprintln(verbose, f"Read line: {line}")
            # Skip empty lines and comments
            if not line or line.startswith("!"):
                myprintln(verbose, "Skipping comment or empty line")
                continue
            # Check if the line contains the 'mult=' string
            if "mult=" in line:
                # Extract the values inside the brackets
                match = re.search(r"mult=\[(.*)\]", line)
                if match:
                    price_str = match.group(1)
                    prices = [float(p) for p in price_str.split()]
                    myprintln(verbose, f"Parsed prices: {prices}")

    # Check if prices is empty
    if not prices:
        raise ValueError(f"No data found in {filenameLoadShape}. The prices array is empty.")

    # Compute hi and lo based on the user input or default to the max/min of prices
    hi = max(prices) if hi is None else hi
    lo = min(prices) if lo is None else lo

    myprintln(verbose, f"Computed hi: {hi}, lo: {lo}")

    # Subsampling or supersampling the 24-step price array to get lambdaVals0
    if T < 24:
        start_idx = int(math.floor((24 - T) / 2))
        lambdaVals0 = prices[start_idx:(start_idx + T)]
        myprintln(verbose, f"Subsampling: lambdaVals0 = {lambdaVals0}")
    elif T > 24:
        lambdaVals0 = []
        myprintln(verbose, f"Supersampling to {T} values:")
        for i in range(1, T + 1):
            scaled_idx = (i - 1) * (24 - 1) / (T - 1)
            lower_idx = int(math.floor(scaled_idx))
            upper_idx = min(lower_idx + 1, 23)
            frac = scaled_idx - lower_idx
            interpolated_value = (1 - frac) * prices[lower_idx] + frac * prices[upper_idx]
            lambdaVals0.append(interpolated_value)
        myprintln(verbose, f"Supersampled: lambdaVals0 = {lambdaVals0}")
    else:
        lambdaVals0 = prices
        myprintln(verbose, f"T == 24, using original prices: lambdaVals0 = {lambdaVals0}")

    num_peak_hours = max(1, int(math.floor(peakHoursFraction * T)))

    myprintln(verbose, f"Number of peak hours: {num_peak_hours}")

    costArray = [lo] * T
    myprintln(verbose, f"Initialized costArray with low values: {costArray}")

    sortedIndices = np.argsort(lambdaVals0)[::-1]
    myprintln(verbose, f"Sorted indices for peak hours: {sortedIndices}")

    for idx in sortedIndices[:num_peak_hours]:
        costArray[idx] = hi
    myprintln(verbose, f"Updated costArray with peak values: {costArray}")

    costData = {
        "LoadShapeCost": costArray,
        "peakCost": hi,
        "offPeakCost": lo,
        "peakHoursFraction": peakHoursFraction,
    }

    return costData

def generateLoadShape(T: int, filenameLoadShape=None):
    wd = os.path.dirname(__file__)

    if filenameLoadShape is None:
        filenameLoadShape = "LoadShapePVDefault.dss"

    loadshape_filepath = os.path.join(wd, "..", "rawData", filenameLoadShape)

    defaultLoadShape = []
    with open(loadshape_filepath, "r") as file:
        for line in file:
            line = line.split("!")[0].strip()
            if not line:
                continue
            if "mult=[" in line:
                shape_str = re.search(r"mult=\[(.*)\]", line).group(1)
                defaultLoadShape = [float(val) for val in shape_str.split()]

    LoadShape = []
    if T < 24:
        start_idx = int(math.floor((24 - T) / 2))
        LoadShape = defaultLoadShape[start_idx:(start_idx + T)]
    elif T > 24:
        for i in range(1, T + 1):
            scaled_idx = (i - 1) * (24 - 1) / (T - 1)
            lower_idx = int(math.floor(scaled_idx))
            upper_idx = min(lower_idx + 1, 23)
            frac = scaled_idx - lower_idx
            interpolated_value = (1 - frac) * defaultLoadShape[lower_idx] + frac * defaultLoadShape[upper_idx]
            LoadShape.append(interpolated_value)
    else:
        LoadShape = defaultLoadShape

    return [float(x) for x in LoadShape]

def trim_number_for_printing(number):
    if number < 1:
        # Round to 4 decimal places
        trimmed_number = round(number, 4)
        trimmed_str = str(trimmed_number)
    elif number < 1e5:
        # Round to 5 significant digits
        trimmed_str = f"{number:.5g}"
    else:
        # For large numbers, also 5 significant digits (scientific notation if needed)
        trimmed_str = f"{number:.5g}"

    formatted_number = trimmed_str.replace('.', '_')
    return formatted_number


