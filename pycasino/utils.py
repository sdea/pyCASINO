""" 
A set of utilities functions
"""

import numpy as np 
from matplotlib import pyplot as plt

import constants

class ElectronTrajectory:

    _can_instantiate = False
    
    def __init__(self, x_list, y_list, theta_list, E_list, is_BSE, is_trasmitted):
        self.x_list = x_list
        self.y_list = y_list
        self.theta_list = theta_list 
        self.E_list = E_list
        self.is_BSE = is_BSE
        self.is_trasmitted = is_trasmitted

    # Check that the users cannot init this class
    def __new__(cls, *args, **kwargs):
        if not cls._can_instantiate:
            raise TypeError("This class cannot be instantiated directly")
        return super().__new__(cls)



def generate_electron_position(center, radius=30, std_dev=15):
    """
    A function to generate the initial y position of an electron within an electron beam, just before interacting with the sample.
    We assume that the beam has an intensity profile that follows a gaussian distribution. 

    Parameters:
        center (float): The initial central position for the Gaussian distribution.
        radius (float): The radius of the beam in nm (default is 30 nm).
        
    Returns:
        float: The position of the electron.
    """
    
    radius = radius * 1E-07
    std_dev = radius / 2  # Standard deviation of the Gaussian distribution
    while True:
        position = np.random.normal(center, std_dev)  # Sample around the initial position
        if -radius <= position <= radius:  # Ensure the position is within the range
            return position

class SimulationParameters:
    """""
    A class to define and store the parameters for the simulation. 
    The default values correspond to a 30 keV beam interacting with a bulk silicon sample.  

    Attributes:
    -----------
    E_beam : float
        The energy of the electron beam in kiloelectronvolts (keV). Default is 30 keV.
    Z : int
        The atomic number of the target material. Default is 14 (Silicon).
    A : float
        The atomic weight of the target material in grams per mole (g/mol). Default is 28.0855 g/mol.
    rho : float
        The density of the target material in grams per cubic centimeter (g/cm³). Default is 2.33 g/cm³.
    x : float
        The x-coordinate of the initial beam position in simulation units. Default is 0.
    y : float
        The y-coordinate of the initial beam position in simulation units. Default is 0.
    theta : float
        The angle of incidence of the electron beam in degrees. Default is 0.
    N_electrons : int
        The number of electrons to simulate in the simulation. Default is 10,000.
    max_events : int
        The maximum number of scattering events allowed per electron in the simulation. Default is 1,000.

    Usage:
    ------
    The `SimulationParameters` class acts as a container for simulation parameters and is meant to be used
    as a configuration object. Simply access its attributes to retrieve or modify the values as needed.
    """""

    def __init__(self):
        
        self.simulation_type = 'bulk'
        self.E_beam = 30     # keV 
        self.Z = 14          # Silicon 
        self.A = 28.0855     # Atomic weight (g/mol) 
        self.rho = 2.33      # Density (g/cm3)
        self.thickness = -1  # Thickness of the sample in nm

        self.radius = 30    # The radius of the electron beam (nm)
        self.x_ini = 0
        self.y_ini = 0
        self.theta = 0

        self.N_electrons = 10000 # The number of electrons in simulation
        self.max_events = 1000   # The maximum number of scattering events in the simulation


    def print_params(self):
        """
        Prints all the simulation parameters in a nicely formatted way.
        """
        params = {
            "Simulation type": self.simulation_type,
            "E_beam (keV)": self.E_beam,
            "Z (Atomic number)": self.Z,
            "A (Atomic weight, g/mol)": self.A,
            "rho (Density, g/cm³)": self.rho,
            "x (Initial x-coordinate)": self.x_ini,
            "y (Initial y-coordinate)": self.y_ini,
            "theta (Incidence angle, degrees)": self.theta,
            "N_electrons (Number of electrons)": self.N_electrons,
            "max_events (Maximum scattering events)": self.max_events,
        }

        print("Simulation Parameters:")
        print("-" * 30)
        for key, value in params.items():
            print(f"{key:<35}: {value}")


class SimulationResults:
    
    # Simulation results should not be created manually by the user
    _can_instantiate = False
    traj_list = []

    # Check that the users cannot init this class
    def __new__(cls, *args, **kwargs):
        if not cls._can_instantiate:
            raise TypeError("This class cannot be instantiated directly. It is only meant as output of a simulation.")
        return super().__new__(cls)

    def __init__(self, traj_list):
        
        # This is enough
        self.traj_list = traj_list

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
        A dictionary containing the atomic number, atomic mass, and mass density of the element.
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
        
def plot_simulation_results(sim_results, plot_interval=10, figsize=(10, 10),
                             beam_color='b', bse_color='r', te_color='g', grid=True):
    """
    A convenience function to plot the results of a simulation

    Parameters
    ----------
    x_traj_list : list
        The list of x-coordinates for the electrons trajectory
    y_traj_list : list
        The list of y-coordinates for the electrons trajectory
    is_bse_list : list 
        The list of Bools indicating if the simulated electron is BSE
    plot_interval : int, optional
        The plot interval (default is 10)
    figsize : tuple, optional
        The figure size (default is (10, 10))
    beam_color : str, optional
        The color for non-BSE trajectories (default is 'b')
    bse_color : str, optional
        The color for BSE trajectories (default is 'r')
    grid : bool, optional
        Whether to display the grid (default is True)
    """

    plt.figure(figsize=figsize)
    num_electrons = len(sim_results.traj_list)

    for n_el in range(0, num_electrons, plot_interval):
        # Invert the coordinates to have the electron beam on top
        x_coords = -np.array(sim_results.traj_list[n_el].x_list) * 1E07
        y_coords = np.array(sim_results.traj_list[n_el].y_list) * 1E07
        if sim_results.traj_list[n_el].is_BSE:
            color = bse_color
        elif sim_results.traj_list[n_el].is_trasmitted:
            color = te_color
        else:
            color = beam_color
        
        linewidth = 2.0 if sim_results.traj_list[n_el].is_BSE else 1.5
        #color = bse_color if sim_results.traj_list[n_el].is_BSE else beam_color
        #linewidth = 2.0 if sim_results.traj_list[n_el].is_BSE else 1.5

        plt.plot(y_coords, x_coords, color=color, linewidth=linewidth)

    if grid:
        plt.grid(True)

    plt.xlabel('Y-axis', fontsize=18)
    plt.ylabel('X-axis', fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.title('Electron Trajectories')
    plt.show()

def find_y_extremes(trajectories):
    """
    Finds the global minimum and maximum y value considering all electron trajectories resulting from a simulation.

    Parameters
    ----------
    nested_list : list of list
        A list of lists containing numerical values.

    Returns
    -------
    list
        A list containing two elements:
        - The global minimum value (float or int).
        - The global maximum value (float or int).

    Raises
    ------
    ValueError
        If the input list is empty or contains no numerical values.
    """
    list_max = []
    list_min = []
    for traj in trajectories:
        y_coors = np.array(traj.y_list)
        if len(y_coors) > 0:
            list_max.append(np.max(y_coors))
            list_min.append(np.min(y_coors))
        else:
            list_max.append(0)
            list_min.append(0)
                

    
    global_min = min(list_min)
    global_max = max(list_max)
    
    return [global_min, global_max]
