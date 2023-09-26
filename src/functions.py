"""
A set of functions needed to compute the electron trajectories  
"""

import numpy as np
import constants
from datamodels import Pos

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

# Given the angle, compute the next position (basde on current P)
def compute_next_P(P_current: Pos, theta: np.float32) -> Pos:
    
    P_new = Pos(0, 0)
    # We include a random (left or right scattering)
    choices = np.array([1, -1]) 
    R = np.random.choice(choices)

    mod_coor = (P_current.rx**2 + P_current.ry**2)**(1/2)
    P_new.rx = (P_current.ry * np.sin(theta) * R) / mod_coor + P_current.rx * np.cos(theta)
    P_new.ry = (-P_current.rx * np.sin(theta) * R) / mod_coor + P_current.ry * np.cos(theta)

    return P_new



