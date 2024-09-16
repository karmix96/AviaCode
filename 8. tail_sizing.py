import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import itertools
from scipy.linalg import solve


def horizontal_tail_calculator(CLw_calc, S, AR, C_mean):
    """ 
        Inside this function, we define the parameters of the horizontal tail.
        
    Returns:
        float:  h_tail, ho_tail, Sh, Vh, CLh, a_h, l_opt_h, L
    """
    Vh = 0.8 # [0.4 - 1.2]
    ho_tail = 0.25 # ho denotes the non-dimensional wing/fuselage aerodynamic center position in relation to MAC [0.2 - 0.25]
    h_tail = 0.2 # The parameter h denotes the non-dimensional aircraft cg position in relation to MAC [0.1 - 0.3]
    Cmowf = -0.06
    CLh = (Cmowf + CLw_calc*(h_tail - ho_tail ))/Vh # desired tail lift coefficient
    a_h = -1 # Horizontal tail setting angle
    Sh = Vh*C_mean*S
    Kc = 1.2 # [1 - 1.4] correction factor
    Df = 0.3  # max fuselage diameter
    l_opt_h = Kc*np.sqrt((4*C_mean*S*Vh)/(np.pi*Df))   # horizontal tail moment arm, the distance between the tail aerodynamic center and the aircraft center of gravity
    L = 0.6/l_opt_h   # l/L = 0.6 The ratio between tail arm and Fuselage Length
    ARh = (2/3)*AR    # [3 - 5]
    a_set_eng = 3 # [2 - 4], the typical engine setting angle is about 2–4 deg
    dh = 0.2 # r. The aircraft non-dimensional center of gravity limit (h) is the difference between the most forward and the most aft position of the aircraft cg. [0.1 - 0.3].
    
    return  h_tail, ho_tail, Sh, Vh, CLh, a_h, l_opt_h, L

def vertical_tail_calculator(l_opt_h):
    """ 
        Inside this function, we define the parameters of the vertical tail.
        
    Returns:
        float:  Vv, Sv, i_v
    """
    l_opt_v = l_opt_h # vertical tail moment arm
    Kf = 0.7 # [0.65 - 0.78]
    Cnv = 0.3 # [0.2 - 0.4]
    Vv = 0.075 # [0.05 - 0.1]
    Sv = 0.125 * S_final # [0.1 - 0.15]
    i_v = 1.5 # [1 - 2] vertical tail incidence to prevent roll from propeller revolution
    ARv = 1.5 # [1 - 2]
    
    return Vv, Sv, i_v