import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import itertools
from scipy.linalg import solve

def ellipitcal_lift_distribution(b_calc, CL_calc):
    
    """ 
        Creates a perfect elliptical curve for comparison

    Returns:
        floats:
                CL_el: Is the array that contains the elliptical lift coefficients.
                y_el:  Is the array that contains the y_coord of the elliptical curve. 
    """
    
    if b_calc is None:
        raise ValueError("You need to calculate 'b' using lifting_line_theory first.")

    # Determine CL at points A and B
    CL_A = CL_calc[N]                                                                     
    y_el = np.linspace(0, b_calc / 2, N+1)
    x = (2 * y_el) / b_calc
    CL_el = CL_A * np.sqrt(1 - x**2)
    return y_el, CL_el
