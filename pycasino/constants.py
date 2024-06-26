"""
A set of useful constants
"""

# Avogadro number 
N_a = 6.02214076 * 1E023

elements = [
    {"symbol": "H", "atomic_number": 1, "density": 0.00008988, "atomic_mass": 1.008},
    {"symbol": "He", "atomic_number": 2, "density": 0.0001785, "atomic_mass": 4.002602},
    {"symbol": "Li", "atomic_number": 3, "density": 0.534, "atomic_mass": 6.94},
    {"symbol": "Be", "atomic_number": 4, "density": 1.85, "atomic_mass": 9.0122},
    {"symbol": "B", "atomic_number": 5, "density": 2.34, "atomic_mass": 10.81},
    {"symbol": "C", "atomic_number": 6, "density": 2.267, "atomic_mass": 12.011},
    {"symbol": "N", "atomic_number": 7, "density": 0.0012506, "atomic_mass": 14.007},
    {"symbol": "O", "atomic_number": 8, "density": 0.001429, "atomic_mass": 15.999},
    {"symbol": "F", "atomic_number": 9, "density": 0.001696, "atomic_mass": 18.998},
    {"symbol": "Ne", "atomic_number": 10, "density": 0.0008999, "atomic_mass": 20.180},
    {"symbol": "Na", "atomic_number": 11, "density": 0.971, "atomic_mass": 22.989769},
    {"symbol": "Mg", "atomic_number": 12, "density": 1.738, "atomic_mass": 24.305},
    {"symbol": "Al", "atomic_number": 13, "density": 2.70, "atomic_mass": 26.981538},
    {"symbol": "Si", "atomic_number": 14, "density": 2.33, "atomic_mass": 28.085},
    {"symbol": "P", "atomic_number": 15, "density": 1.82, "atomic_mass": 30.973762},
    {"symbol": "S", "atomic_number": 16, "density": 2.067, "atomic_mass": 32.06},
    {"symbol": "Cl", "atomic_number": 17, "density": 0.003214, "atomic_mass": 35.45},
    {"symbol": "Ar", "atomic_number": 18, "density": 0.0017837, "atomic_mass": 39.948},
    {"symbol": "K", "atomic_number": 19, "density": 0.862, "atomic_mass": 39.0983},
    {"symbol": "Ca", "atomic_number": 20, "density": 1.54, "atomic_mass": 40.078},
    {"symbol": "Sc", "atomic_number": 21, "density": 2.985, "atomic_mass": 44.955908},
    {"symbol": "Ti", "atomic_number": 22, "density": 4.506, "atomic_mass": 47.867},
    {"symbol": "V", "atomic_number": 23, "density": 6.11, "atomic_mass": 50.9415},
    {"symbol": "Cr", "atomic_number": 24, "density": 7.15, "atomic_mass": 51.9961},
    {"symbol": "Mn", "atomic_number": 25, "density": 7.44, "atomic_mass": 54.938044},
    {"symbol": "Fe", "atomic_number": 26, "density": 7.874, "atomic_mass": 55.845},
    {"symbol": "Co", "atomic_number": 27, "density": 8.86, "atomic_mass": 58.933194},
    {"symbol": "Ni", "atomic_number": 28, "density": 8.912, "atomic_mass": 58.6934},
    {"symbol": "Cu", "atomic_number": 29, "density": 8.96, "atomic_mass": 63.546},
    {"symbol": "Zn", "atomic_number": 30, "density": 7.134, "atomic_mass": 65.38},
    {"symbol": "Ga", "atomic_number": 31, "density": 5.91, "atomic_mass": 69.723},
    {"symbol": "Ge", "atomic_number": 32, "density": 5.323, "atomic_mass": 72.63},
    {"symbol": "As", "atomic_number": 33, "density": 5.727, "atomic_mass": 74.921595},
    {"symbol": "Se", "atomic_number": 34, "density": 4.81, "atomic_mass": 78.971},
    {"symbol": "Br", "atomic_number": 35, "density": 3.12, "atomic_mass": 79.904},
    {"symbol": "Kr", "atomic_number": 36, "density": 0.003733, "atomic_mass": 83.798},
    {"symbol": "Rb", "atomic_number": 37, "density": 1.532, "atomic_mass": 85.4678},
    {"symbol": "Sr", "atomic_number": 38, "density": 2.64, "atomic_mass": 87.62},
    {"symbol": "Y", "atomic_number": 39, "density": 4.469, "atomic_mass": 88.90584},
    {"symbol": "Zr", "atomic_number": 40, "density": 6.506, "atomic_mass": 91.224},
    {"symbol": "Nb", "atomic_number": 41, "density": 8.57, "atomic_mass": 92.90637},
    {"symbol": "Mo", "atomic_number": 42, "density": 10.22, "atomic_mass": 95.95},
    {"symbol": "Tc", "atomic_number": 43, "density": 11.5, "atomic_mass": 98},
    {"symbol": "Ru", "atomic_number": 44, "density": 12.37, "atomic_mass": 101.07},
    {"symbol": "Rh", "atomic_number": 45, "density": 12.41, "atomic_mass": 102.90550},
    {"symbol": "Pd", "atomic_number": 46, "density": 12.02, "atomic_mass": 106.42},
    {"symbol": "Ag", "atomic_number": 47, "density": 10.501, "atomic_mass": 107.8682},
    {"symbol": "Cd", "atomic_number": 48, "density": 8.69, "atomic_mass": 112.414},
    {"symbol": "In", "atomic_number": 49, "density": 7.31, "atomic_mass": 114.818},
    {"symbol": "Sn", "atomic_number": 50, "density": 7.287, "atomic_mass": 118.710},
    {"symbol": "Sb", "atomic_number": 51, "density": 6.685, "atomic_mass": 121.760},
    {"symbol": "Te", "atomic_number": 52, "density": 6.232, "atomic_mass": 127.60},
    {"symbol": "I", "atomic_number": 53, "density": 4.93, "atomic_mass": 126.90447},
    {"symbol": "Xe", "atomic_number": 54, "density": 0.005887, "atomic_mass": 131.293},
    {"symbol": "Cs", "atomic_number": 55, "density": 1.873, "atomic_mass": 132.90545196},
    {"symbol": "Ba", "atomic_number": 56, "density": 3.594, "atomic_mass": 137.327},
    {"symbol": "La", "atomic_number": 57, "density": 6.145, "atomic_mass": 138.90547},
    {"symbol": "Ce", "atomic_number": 58, "density": 6.77, "atomic_mass": 140.116},
    {"symbol": "Pr", "atomic_number": 59, "density": 6.77, "atomic_mass": 140.90766},
    {"symbol": "Nd", "atomic_number": 60, "density": 7.01, "atomic_mass": 144.242},
    {"symbol": "Pm", "atomic_number": 61, "density": 7.26, "atomic_mass": 145},
    {"symbol": "Sm", "atomic_number": 62, "density": 7.52, "atomic_mass": 150.36},
    {"symbol": "Eu", "atomic_number": 63, "density": 5.243, "atomic_mass": 151.964},
    {"symbol": "Gd", "atomic_number": 64, "density": 7.895, "atomic_mass": 157.25},
    {"symbol": "Tb", "atomic_number": 65, "density": 8.229, "atomic_mass": 158.92535},
    {"symbol": "Dy", "atomic_number": 66, "density": 8.55, "atomic_mass": 162.500},
    {"symbol": "Ho", "atomic_number": 67, "density": 8.795, "atomic_mass": 164.93033},
    {"symbol": "Er", "atomic_number": 68, "density": 9.066, "atomic_mass": 167.259},
    {"symbol": "Tm", "atomic_number": 69, "density": 9.321, "atomic_mass": 168.93422},
    {"symbol": "Yb", "atomic_number": 70, "density": 6.965, "atomic_mass": 173.045},
    {"symbol": "Lu", "atomic_number": 71, "density": 9.84, "atomic_mass": 174.9668},
    {"symbol": "Hf", "atomic_number": 72, "density": 13.31, "atomic_mass": 178.49},
    {"symbol": "Ta", "atomic_number": 73, "density": 16.654, "atomic_mass": 180.94788},
    {"symbol": "W", "atomic_number": 74, "density": 19.25, "atomic_mass": 183.84},
    {"symbol": "Re", "atomic_number": 75, "density": 21.02, "atomic_mass": 186.207},
    {"symbol": "Os", "atomic_number": 76, "density": 22.59, "atomic_mass": 190.23},
    {"symbol": "Ir", "atomic_number": 77, "density": 22.56, "atomic_mass": 192.217},
    {"symbol": "Pt", "atomic_number": 78, "density": 21.45, "atomic_mass": 195.084},
    {"symbol": "Au", "atomic_number": 79, "density": 19.32, "atomic_mass": 196.966569},
    {"symbol": "Hg", "atomic_number": 80, "density": 13.534, "atomic_mass": 200.59},
    {"symbol": "Tl", "atomic_number": 81, "density": 11.85, "atomic_mass": 204.38},
    {"symbol": "Pb", "atomic_number": 82, "density": 11.34, "atomic_mass": 207.2},
    {"symbol": "Bi", "atomic_number": 83, "density": 9.78, "atomic_mass": 208.98040},
    {"symbol": "Po", "atomic_number": 84, "density": 9.32, "atomic_mass": 209},
    {"symbol": "At", "atomic_number": 85, "density": 7, "atomic_mass": 210},
    {"symbol": "Rn", "atomic_number": 86, "density": 0.00973, "atomic_mass": 222},
    {"symbol": "Fr", "atomic_number": 87, "density": 1.87, "atomic_mass": 223},
    {"symbol": "Ra", "atomic_number": 88, "density": 5.5, "atomic_mass": 226},
    {"symbol": "Ac", "atomic_number": 89, "density": 10.07, "atomic_mass": 227},
    {"symbol": "Th", "atomic_number": 90, "density": 11.72, "atomic_mass": 232.0377},
    {"symbol": "Pa", "atomic_number": 91, "density": 15.37, "atomic_mass": 231.03588},
    {"symbol": "U", "atomic_number": 92, "density": 18.95, "atomic_mass": 238.02891}
]