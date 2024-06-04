"""
A set of functions needed to compute the electron trajectories  
"""

import numpy as np
import constants

def compute_alpha(E, Z):

    """
    Compute the screening factor, taking into account the electron cloud surrounding each atom.
    
    Parameters:
    E (float): The energy of the electron.
    Z (int): The atomic number of the element.
    
    Returns:
    float: The screening factor.
    """

    alpha = 3.4E-03 * Z ** (0.67) / E
    return alpha

def compute_sigma(E, Z, alpha):

    """
    Compute the Rutherford scattering cross-section.
    
    Parameters:
    E (float): The energy of the electron.
    Z (int): The atomic number of the element.
    alpha (float): The screening factor.
    
    Returns:
    float: The Rutherford scattering cross-section in cm^2.
    """

    numeric_factor = 5.21E-21 * (Z**2)/(E**2) 
    alpha_factor = (4 * np.pi) / (alpha * (1 + alpha)) 
    energy_factor = ((E + 511) / (E + 1024))**2

    sigma = numeric_factor * alpha_factor * energy_factor
    return sigma
 
def compute_step(lambda_mean):

    """
    Compute the random path length taken by the electron.
    
    Parameters:
    lambda_mean (float): The mean free path length.
    
    Returns:
    float: The random path length.
    """

    random_step = -lambda_mean * np.log(np.random.uniform(0, 1))
    return random_step  

def compute_scatt_angle(alpha):

    """
    Compute the cosine of the solid angle.
    
    Parameters:
    alpha (float): The screening factor.
    
    Returns:
    float: The scattering angle in radians.
    """

    rnd_number = np.random.uniform(1, 0)
    cos_theta = 1 - ((2 * alpha * rnd_number))/(1 + alpha - rnd_number)
    theta = np.arccos(cos_theta)
    return theta

def compute_lambda(A, rho, sigma):

    """
    Compute the mean free path length.
    
    Parameters:
    A (float): The atomic mass number.
    rho (float): The density of the material.
    sigma (float): The scattering cross-section.
    
    Returns:
    float: The mean free path length.
    """

    lambda_mean = A / (constants.N_a * rho * sigma)
    return lambda_mean 
    
def compute_J(Z):

    """
    Compute the mean ionization potential.
    
    Parameters:
    Z (int): The atomic number of the element.
    
    Returns:
    float: The mean ionization potential in keV.
    """

    J = (9.76 * Z + 58.5 / Z**0.19)*1E-03
    return J


def compute_energy_loss(Z, E, rho, A, step):

    """
    Compute the energy loss per unit distance (keV/cm).
    
    Parameters:
    Z (int): The atomic number of the element.
    E (float): The energy of the electron.
    rho (float): The density of the material.
    A (float): The atomic mass number.
    step (float): The step length.
    
    Returns:
    float: The energy loss in keV/cm.
    """

    J = compute_J(Z)
    dEdS = -78500 * ((rho * Z) / (A * E)) * np.log((1.166 * E) / J + 1)
    E_loss = step * dEdS

    return E_loss

def compute_electron_trajectory(E, Z, A, rho, x_ini = 0, y_ini = 0, max_steps = 1000, energy_cutoff = 0.2):

    """
    This function simulates the trajectory of a single electron, starting a x_ini and y_ini, with initial energy E.

    Parameters:
    E (float): Initial energy of the electron.
    Z (int): Atomic number of the element.
    A (float): Atomic mass number of the element.
    rho (float): Density of the material.
    x_ini (float, optional): Initial x-coordinate of the electron. Default is 0.
    y_ini (float, optional): Initial y-coordinate of the electron. Default is 0.
    max_steps (int, optional): Maximum number of steps to simulate. Default is 1000.
    energy_cutoff (float, optional): Energy cutoff to stop the simulation. Default is 0.2 keV.

    Returns:
    tuple: 
        - x_list (list of float): List of x-coordinates of the electron trajectory.
        - y_list (list of float): List of y-coordinates of the electron trajectory.
        - is_backscattered (bool): True if the electron is backscattered, False otherwise.
    """

    x_list = []
    y_list = []
    is_backscattered = False
    x = x_ini
    y = y_ini

    # Include tilt in the simulations
    # It should be as easy as including a theta ini in the parameters
    theta = 0

    for i in range(0, max_steps):

        alpha = compute_alpha(E, Z)
        sigma = compute_sigma(E, Z, alpha)
        lambda_mean = compute_lambda(A, rho, sigma)
        step = compute_step(lambda_mean)
        theta_new = compute_scatt_angle(alpha)

        # Get the random on the sign
        rand_sign = 2 * np.random.randint(0, 2) - 1
        theta_new = theta_new * rand_sign

        theta = theta + theta_new

        # Update the new trajectory
        x = x + step * np.cos(theta);
        y = y + step * np.sin(theta);
        
        # We check if the electron leaves the material
        if (x <= 0):
            is_backscattered = True
            break
        else:
            is_backscattered = False
        
        if (E < 0.2):
            break
        else:
            x_list.append(x)
            y_list.append(y)
    
    return x_list, y_list, is_backscattered


def get_element_properties(symbol):
    """
    Get the properties of an element by its symbol.

    Parameters:
    symbol (str): The symbol of the element.

    Returns:
    dict: A dictionary containing the atomic number, atomic mass, and density of the element.
    """
    for element in constants.elements:
        if element["symbol"] == symbol:
            return {
                "atomic_number": element["atomic_number"],
                "atomic_mass": element["atomic_mass"],
                "density": element["density"],
            }
        else:
            print(f'Element {symbol} not found! Returning None...')
            return None












