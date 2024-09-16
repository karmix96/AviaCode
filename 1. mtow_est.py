import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import itertools
from scipy.linalg import solve

def mtow_estimation():
    
    """ 
    Calculates the estimated mtom for payload.
    Equations used were derived from file "Database.xlsx".

    # Returns:
    #     int: mtom: Maximum Take Off Mass in kg
    # """
    
    # Get user input for payload
    payload = float(input("Enter payload: "))
    
    # Calculate estimated mtom for a given payload
    if payload <= 10:
        mtom = 3.6522 * math.exp(0.274 * payload)
    elif payload <=100:
        mtom = 40.514 * math.exp(0.0332 * payload)
    elif payload > 100:
        mtom = 714.87 * math.exp(0.0013 * payload)
        
    mtow = mtom*g
    # print (f"mtow is: ", mtow)
    return mtom, mtow