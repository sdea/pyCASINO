"""
A set of functions needed to compute the electron trajectories  
"""

import numpy as np
import constants

# Screening factor, taking into account the electron clound sorrounding each atom
def compute_alpha(E, Z):

    alpha = 3.4E-03 * Z**(0.67)/E
    return alpha

# Rutherford scattering cross-section (cm2)
def compute_sigma(E, Z, alpha):

    numeric_factor = 5.21E-21 * (Z**2)/(E**2) 
    alpha_factor = (4 * np.pi) / (alpha * (1 + alpha)) 
    energy_factor = ((E + 511)/(E + 1024))**2

    sigma = numeric_factor * alpha_factor * energy_factor
    return sigma

# The random path lenght takend by the electron 
def compute_step(lambda_mean):

    random_step = -lambda_mean * np.log(np.random.uniform(0, 1))
    return random_step  

# Compute the cosine of the solid angle (radomic)
def compute_scatt_angle(alpha):

    rnd_number = np.random.uniform(1, 0)
    cos_theta = 1 - ((2 * alpha * rnd_number))/(1 + alpha - rnd_number)
    return cos_theta

def compute_lambda(A, rho, sigma):
    lambda_mean = A / (constants.N_a * rho * sigma)
    return lambda_mean 



