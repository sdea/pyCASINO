"""
A set of functions needed to compute the electron trajectories
"""
import sys
sys.path.append('/Users/sdea/Coding/pyCASINO')

import numpy as np
from pycasino import constants

def compute_alpha(E, Z):
    """
    This function computes the screening factor, which accounts for the sorrounding electron cloud of the nucleus. The expression used is empirically estimated (Bishop 1976).
    
    Parameters
    ----------
    E : float
        The energy of the electron.
    Z : int
        The atomic number of the target.
    
    Returns
    -------
    float
        The screening factor.
    """
    alpha = 3.4E-03 * Z ** (0.67) / E
    return alpha

def compute_sigma(E, Z, alpha):
    """
    This functions computes the relativistic corrected, screened Rutherford elastic cross-section. E is the energy of the electron in keV, Z is the atomic number of the target and alpha is a screening factor. The screening factor accounts for the fact that the incident electron does not see all the charge of the nucleus because of the sorrounding electron cloud.
    
    Parameters
    ----------
    E : float
        The energy of the electron.
    Z : int
        The atomic number of the target.
    alpha : float
        The screening factor.
    
    Returns
    -------
    float
        The Rutherford elastic scattering cross-section in cm^2.
    """
    numeric_factor = 5.21E-21 * (Z**2)/(E**2)
    alpha_factor = (4 * np.pi) / (alpha * (1 + alpha))
    energy_factor = ((E + 511) / (E + 1024))**2

    sigma = numeric_factor * alpha_factor * energy_factor
    return sigma

def compute_step(lambda_mean):
    """
    This function computes the random path length traveled by the electron between two scattering events.
    
    Parameters
    ----------
    lambda_mean : float
        The mean free path length.
    
    Returns
    -------
    float
        The random path length expressed in cm.
    """
    random_step = -lambda_mean * np.log(np.random.uniform(0, 1))
    return random_step

def compute_scatt_angle(alpha):
    """
    This function computes the scattering angle related to a scattering event.
    
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
    This function computes the mean free path of the electron in the target material, given the Rutherford cross-section sigma.
    
    Parameters
    ----------
    A : float
        The atomic mass number of the target (g/mol).
    rho : float
        The density of the target (g/cm^3).
    sigma : float
        The scattering cross-section (cm^2).
    
    Returns
    -------
    float
        The mean free path length in cm.
    """
    # lambda_path = A / (rho * constants.N_a * sigma * 1E-21)
    lambda_path = A / (rho * constants.N_a * sigma)
    return lambda_path

def compute_J(Z):
    """
    This function computes the mean ionization potential for a given element with atomic number Z.
    
    Parameters
    ----------

    Z : int 
        Atomic number of the target material.
    
    Returns
    -------
    float
        The mean ionization potential in keV.
    """

    J = (9.76 * Z + 58.5 / Z**0.19) * 1E-03
    return J

def compute_energy_loss(E, Z, A, rho, step):
    """
    This function computes the energy loss experienced by the electron as it moves through the target.
    
    Parameters
    ----------
    E : float 
        Energy of the electron in keV.
    Z : int
        Atomic number of the target material.
    A : float 
        Atomic mass of the target in g/mol.
    rho : float 
        Density of the target in g/cm^3.
    step : float 
        The random step length between scattering events in cm.
    
    Returns
    ------- 
    float
        Energy loss of the electron over the given step size. The energy loss is measured in MeV.
    """

    J = compute_J(Z)
    dEdS = -78500 * ((rho * Z) / (A * E)) * np.log((1.166 * E) / J + 1)
    E_loss = step * dEdS
    return E_loss

def compute_single_trajectory(E, Z, A, rho, x_ini, y_ini, max_steps=1000):
    """
    This function computes the trajectory of a single electron inside a uniform, bulk material. 
    
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
            - x_list (list of float): List of x-coordinates of the single electron trajectory.
            - y_list (list of float): List of y-coordinates of the single electron trajectory.
            - is_backscattered (bool): True if the electron is backscattered, False otherwise.
    """
    x_list = []
    y_list = []
    is_backscattered = False
    x = x_ini
    y = y_ini

    # This is temporary, to include tilt in a future version
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
        
        E_loss = compute_energy_loss(E, Z, A, rho, step)
        E = E + E_loss
        
        if E < 0.2:
            break
        else:
            x_list.append(x)
            y_list.append(y)

        
    
    return x_list, y_list, is_backscattered

def simulate_bulk_interaction(E_beam, Z, A, rho, radius, center, num_electrons = 5000):

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
     A tuple containing:
            A list of bools - for each electron simulated, True if backscattered, False otherwise
            A list of trajectories - each element of the list contains the x-coordinates for the trajectory of a simulated electron  
            A list of trajectories - each element of the list contains the y-coordinates for the trajectory of a simulated electron
    """
    x_list_final = [] 
    y_list_final = []
    bse_list = []
    for n_el in range(0, num_electrons):
    
        # Reset parameters  
        x = 0
        E = E_beam

        # For the radius we need a gauss distribution
        y = center
        x_list = []
        y_list = []
        is_backscattered = False
        
        x_list, y_list, is_backscattered = compute_single_trajectory(E, Z, A, rho, x, y) 
            
        bse_list.append(is_backscattered)
        x_list_final.append(x_list)
        y_list_final.append(y_list)
    
    return bse_list, x_list_final, y_list_final