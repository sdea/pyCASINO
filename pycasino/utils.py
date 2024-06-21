""" 
A set of utilities functions
"""

import numpy as np 
from matplotlib import pyplot as plt

from pycasino import constants

def get_element_properties(symbol):
    """
    This function returns the properties of an element by its symbol.

    Parameters
    ----------
    symbol : str
        The symbol of the element.

    Returns
    -------
    dict
        A dictionary containing the atomic number, atomic mass, and density of the element.
    """
    for element in constants.elements:
        if element["symbol"] == symbol:
            return {
                "atomic_number": element["atomic_number"],
                "atomic_mass": element["atomic_mass"],
                "density": element["density"]
            }
        else:
            print(f'Element {symbol} not found! Returning None...')
            return None
        
def plot_simulation_results(x_traj_list, y_traj_list, is_bse_list, plot_interval=10, figsize=(10,10)):
    
    """
    A convenience function to plot the results of a simulation

    Parameters
    ----------
    x_traj_list : list
        The list of x-coordinates for the electrons trajectory
     y_traj_list : list
        The list of y-coordinates for the electrons trajectory
    is_bse_list : list 
        The list of Bools indicating the the simulated electron is BSE
    plot_interval : int, optional
        The plot interval

    """

    plt.figure(figsize = figsize)
    num_electrons = len(x_traj_list)
    for n_el in np.arange(1, num_electrons, plot_interval):

        # Invert the coordinates to have the electron beam on top
        if is_bse_list[n_el] is True:
            plt.plot(np.array(y_traj_list[n_el]), -np.array(x_traj_list[n_el]), 'r-', linewidth = 2.0)
        else:
            plt.plot(np.array(y_traj_list[n_el]), -np.array(x_traj_list[n_el]), 'b-', linewidth = 1.5)

    plt.grid(True)
    plt.show()
