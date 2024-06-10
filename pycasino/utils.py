""" 
A set of utilities functions
"""

from pycasino import constants

def get_element_properties(symbol):
    """
    Get the properties of an element by its symbol.

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
                "density": element["density"],
            }
        else:
            print(f'Element {symbol} not found! Returning None...')
            return None
        

