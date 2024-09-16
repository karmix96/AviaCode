import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import itertools
from scipy.linalg import solve

def wing_aerodynamics_calculator():
    
    """ 
        Application of the Lifting Line Theory to calculate the Wing Lift Coefficient.

    Returns:
        float: y_calc, CL_calc, CLw_calc, b_calc
    """
    b_calc = math.sqrt(AR * S)  # wing spaN(m)
    MAC = S / b_calc  # MeaNAerodynamic Chord (m)
    Croot = (1.5 * (1 + taper) * MAC) / (1 + taper + taper ** 2)  # root chord (m)

    theta = np.linspace((math.pi / (2 * N)), (math.pi / 2), N, endpoint=True)
    alpha = np.linspace(i_w + alpha_twist, i_w, N)
    z = (b_calc / 2) * np.cos(theta)
    c = Croot * (1 - (1 - taper) * np.cos(theta))  # Mean Aerodynamics
    mu = c * a_2d / (4 * b_calc)

    LHS = mu * (np.array(alpha) - alpha_0) / 57.3  # .reshape((N-1),1)# Left Hand Side

    RHS = []
    for i in range(1, 2 * N + 1, 2):
        RHS_iter = np.sin(i * theta) * (1 + (mu * i) / (np.sin(list(theta))))  # .reshape(1,N)
        # print(RHS_iter,"RHS_iter shape")
        RHS.append(RHS_iter)

    test = np.asarray(RHS)
    x = np.transpose(test)
    inv_RHS = np.linalg.inv(x)

    ans = np.matmul(inv_RHS, LHS)

    mynum = np.divide((4 * b_calc), c)

    CL = (np.sin((1) * theta)) * ans[0] * mynum
    for i in range(1,N-1):
        CL = CL + (np.sin((2*i+1) * theta)) * ans[i] * mynum
        
    CL_calc = np.append(0, CL)

    y_calc = np.concatenate(([b_calc / 2], z))

    # The Available Wing Lift Coefficient for Given Geometrical Parameters
    CLw_calc = (math.pi * AR * ans[0])
    Lw = 0.5*rho*(Vc**2)*S*CLw_calc
    # print(f"The Available Wing Lift Coef. is: ", CLw_calc,
        # f"The Available Wing Lift is: ", Lw )
    return  y_calc, CL_calc, CLw_calc, b_calc
