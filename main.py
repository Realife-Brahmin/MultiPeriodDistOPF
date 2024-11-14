# solve_mpopf.py
import sys
import os
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src'))

# Import necessary modules
from src.MultiPeriodDistOPF import parse_all_data

# Define system parameters
system_name = "ads10_1ph"  # Switch between systems if needed
# system_name = "ieee123_1ph"
T0 = 24
factor = 1
T = int(T0 * factor)

# Parse all data
data = parse_all_data(system_name, T)
print(5)
