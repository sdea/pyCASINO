"""
A set of functions needed to compute the electron trajectories
"""

import numpy as np
from tqdm import tqdm
import constants, utils

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

def compute_single_trajectory_thin(E, Z, A, rho, x_ini, y_ini, theta_ini, thickness, max_steps=1000):
    """
    This function computes the trajectory of a single electron across a thin-film sample. These types of sample are used in TEM analysis 
    
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
            - is_trasmitted (bool): True if the electron is backscattered, False otherwise.
    """
    x_list = []
    y_list = []
    theta_list = []
    E_list =  []
    is_backscattered = False
    is_trasmitted = False
    x = x_ini
    y = y_ini

    theta = np.deg2rad(theta_ini)

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
        
        # Check if the electron is backscattered
        if x <= 0:
            is_backscattered = True
            break
        else:
            if thickness >= 0:
                thickness_cm = thickness * 1E-07
                if x >= thickness_cm:
                    is_trasmitted = True
                    break
                else:
                    is_trasmitted = False
            
        E_loss = compute_energy_loss(E, Z, A, rho, step)
        E = E + E_loss
        
        if E < 0.2:
            break
        else:
            x_list.append(x)
            y_list.append(y)
            theta_list.append(theta)
            E_list.append(E)

    utils.ElectronTrajectory._can_instantiate = True
    trajectory = utils.ElectronTrajectory(x_list, y_list, theta_list, E_list, is_backscattered, is_trasmitted)    
    
    return trajectory

def compute_single_trajectory(E, Z, A, rho, x_ini, y_ini, theta_ini, max_steps=1000):
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
            - is_trasmitted (bool): True if the electron is backscattered, False otherwise.
    """
    x_list = []
    y_list = []
    theta_list = []
    E_list =  []
    is_backscattered = False
    is_trasmitted = False
    x = x_ini
    y = y_ini

    theta = np.deg2rad(theta_ini)

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
        
        # Check if the electron is backscattered
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
            theta_list.append(theta)
            E_list.append(E)

    utils.ElectronTrajectory._can_instantiate = True
    trajectory = utils.ElectronTrajectory(x_list, y_list, theta_list, E_list, is_backscattered, is_trasmitted)    
    
    return trajectory

def compute_single_trajectory_interface(E, Z, A, rho, x_ini, y_ini, theta_ini, max_steps=1000):
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
            - is_trasmitted (bool): True if the electron is backscattered, False otherwise.
    """
    
    # E, Z, A, rho are in this case a list of two elements
    x_list = []
    y_list = []
    theta_list = []
    E_list =  []
    is_backscattered = False
    is_trasmitted = False
    x = x_ini
    y = y_ini

    theta = np.deg2rad(theta_ini)

    for i in range(max_steps):
        
        if y <= 0:
            alpha = compute_alpha(E, Z[0])
            sigma = compute_sigma(E, Z[0], alpha)
            lambda_mean = compute_lambda(A[0], rho[0], sigma)
        else:
            alpha = compute_alpha(E, Z[1])
            sigma = compute_sigma(E, Z[1], alpha)
            lambda_mean = compute_lambda(A[1], rho[1], sigma)
        
        step = compute_step(lambda_mean)
        theta_new = compute_scatt_angle(alpha)
        rand_sign = 2 * np.random.randint(0, 2) - 1
        theta_new = theta_new * rand_sign

        theta = theta + theta_new

        x = x + step * np.cos(theta)
        y = y + step * np.sin(theta)
        
        # Check if the electron is backscattered
        if x <= 0:
            is_backscattered = True
            break
        else:
            is_backscattered = False
        
        if y<=0:
            E_loss = compute_energy_loss(E, Z[0], A[0], rho[0], step)
        else:
            E_loss = compute_energy_loss(E, Z[1], A[1], rho[1], step)
        
        E = E + E_loss
        
        if E < 0.2:
            break
        else:
            x_list.append(x)
            y_list.append(y)
            theta_list.append(theta)
            E_list.append(E)

    utils.ElectronTrajectory._can_instantiate = True
    trajectory = utils.ElectronTrajectory(x_list, y_list, theta_list, E_list, is_backscattered, is_trasmitted)    
    
    return trajectory

def get_bulk_results(sim_param):
    """
    A function to compute all electron trajectories in the material.
    
    Parameters
    ----------
    sim_param : SimulationParameters
        The SimulationParameters class storing the parameters for the simulation
    
    Returns
    -------
     A simulation result struct *to write it better*        
    """
    traj_list = []
    
    for n_el in tqdm(range(0, sim_param.N_electrons), desc = 'N Electrons'):
    
        # Reset parameters  
        x = 0 
        E = sim_param.E_beam

        # For the radius we need a gaussian distribution
        y = utils.generate_electron_position(sim_param.y_ini, sim_param.radius)
        
        x_list = []
        y_list = []
        is_backscattered = False
        
        if sim_param.simulation_type == 'bulk':
            trajectory = compute_single_trajectory(E, sim_param.Z, sim_param.A, sim_param.rho, x, y, sim_param.theta)
        elif sim_param.simulation_type == 'interface':
            trajectory = compute_single_trajectory_interface(E, sim_param.Z, sim_param.A, sim_param.rho, x, y, sim_param.theta) 
        elif sim_param.simulation_type == 'thin-film':
            trajectory = compute_single_trajectory_thin(E, sim_param.Z, sim_param.A, sim_param.rho, x, y, sim_param.theta, sim_param.thickness)         
        
        traj_list.append(trajectory)


    # Instantiate the simulation results class 
    utils.SimulationResults._can_instantiate = True
    sim_results = utils.SimulationResults(traj_list)

    return sim_results