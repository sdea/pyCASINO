"""
A set of functions needed to compute the electron trajectories
"""
import sys
sys.path.append('/Users/sdea/Coding/pyCASINO')

import numpy as np
from pycasino import constants

def compute_alpha(E, Z):
    """
    Compute the screening factor, taking into account the electron cloud surrounding each atom.
    
    Parameters
    ----------
    E : float
        The energy of the electron.
    Z : int
        The atomic number of the element.
    
    Returns
    -------
    float
        The screening factor.
    """
    alpha = 3.4E-03 * Z ** (0.67) / E
    return alpha

def compute_sigma(E, Z, alpha):
    """
    Compute the Rutherford scattering cross-section.
    
    Parameters
    ----------
    E : float
        The energy of the electron.
    Z : int
        The atomic number of the element.
    alpha : float
        The screening factor.
    
    Returns
    -------
    float
        The Rutherford scattering cross-section in cm^2.
    """
    numeric_factor = 5.21E-21 * (Z**2)/(E**2)
    alpha_factor = (4 * np.pi) / (alpha * (1 + alpha))
    energy_factor = ((E + 511) / (E + 1024))**2

    sigma = numeric_factor * alpha_factor * energy_factor
    return sigma

def compute_step(lambda_mean):
    """
    Compute the random path length taken by the electron.
    
    Parameters
    ----------
    lambda_mean : float
        The mean free path length.
    
    Returns
    -------
    float
        The random path length.
    """
    random_step = -lambda_mean * np.log(np.random.uniform(0, 1))
    return random_step

def compute_scatt_angle(alpha):
    """
    Compute the cosine of the solid angle.
    
    Parameters
    ----------
    alpha : float
        The screening factor.
    
    Returns
    -------
    float
        The scattering angle in radians.
    """
    rnd_number = np.random.uniform(1, 0)
    cos_theta = 1 - ((2 * alpha * rnd_number))/(1 + alpha - rnd_number)
    theta = np.arccos(cos_theta)
    return theta

def compute_lambda(A, rho, sigma):
    """
    Compute the mean free path length.
    
    Parameters
    ----------
    A : float
        The atomic mass number.
    rho : float
        The density of the material.
    sigma : float
        The scattering cross-section.
    
    Returns
    -------
    float
        The mean free path length in cm.
    """
    lambda_path = A / (rho * constants.N_a * sigma * 1E-21)
    return lambda_path

def compute_single_trajectory(E, Z, A, rho, x_ini, y_ini, max_steps=1000):
    """
    Simulate the trajectory of an electron through a material.
    
    Parameters
    ----------
    E : float
        The initial energy of the electron.
    Z : int
        The atomic number of the element.
    A : float
        The atomic mass number.
    rho : float
        The density of the material.
    x_ini : float
        The initial x-coordinate of the electron.
    y_ini : float
        The initial y-coordinate of the electron.
    max_steps : int, optional
        The maximum number of steps to simulate (default is 1000).
    
    Returns
    -------
    tuple
        A tuple containing:
        - x_list (list of float): List of x-coordinates of the electron trajectory.
        - y_list (list of float): List of y-coordinates of the electron trajectory.
        - is_backscattered (bool): True if the electron is backscattered, False otherwise.
    """
    x_list = []
    y_list = []
    is_backscattered = False
    x = x_ini
    y = y_ini

    theta = 0

    for i in range(max_steps):
        alpha = compute_alpha(E, Z)
        sigma = compute_sigma(E, Z, alpha)
        lambda_mean = compute_lambda(A, rho, sigma)
        step = compute_step(lambda_mean)
        theta_new = compute_scatt_angle(alpha)

        rand_sign = 2 * np.random.randint(0, 2) - 1
        theta_new = theta_new * rand_sign

        theta = theta + theta_new

        x = x + step * np.cos(theta)
        y = y + step * np.sin(theta)
        
        if x <= 0:
            is_backscattered = True
            break
        else:
            is_backscattered = False
        
        if E < 0.2:
            break
        else:
            x_list.append(x)
            y_list.append(y)
    
    return x_list, y_list, is_backscattered

def simulate_bulk_interaction(E, Z, A, rho, radius, center, num_electrons = 5000):

    """
    Simulate the trajectory of an electron through a material.
    
    Parameters
    ----------
    E : float
        The initial energy of the electron.
    Z : int
        The atomic number of the element.
    A : float
        The atomic mass number.
    rho : float
        The density of the material.
    radius : float
        The radius of the electron beam
    center : float
        The y-coordinate for the center of the electron beam
    num_electrons : int, optional
        The number of electrons to simulate (default is 5000).
    
    Returns
    -------
    
    """
    x_list_final = [] 
    y_list_final = []
    bse_list = []
    for n_el in range(0, num_electrons):
    
        # Reset parameters  
        x = 0

        # For the radius we need a gauss distribution
        y = center
        x_list = []
        y_list = []
        is_backscattered = False
        
        x_list, y_list, is_backscattered = compute_single_trajectory(E, Z, A, rho, x, y) 
            
        # The main loop
        bse_list.append(is_backscattered)
        x_list_final.append(x_list)
        y_list_final.append(y_list)
    
    return bse_list, x_list_final, y_list_final