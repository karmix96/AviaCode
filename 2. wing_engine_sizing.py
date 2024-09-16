import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import itertools
from scipy.linalg import solve

def wing_and_engine_sizing(mtom, mtow, S, P_design, AR):
    
    """ 
        Inside this function, we define the permitable limits base on the combination of our demands.
        We are also able to find the optimal design point, if the motor is NOT already known but can be chosen.
        
    Returns:
        float: opt_WL 
        float: wl_design 
        float: opt_PL
        float: pl_design
    """
    
    # General
    data = pd.DataFrame()
    k = 1/(pi*e*AR)                                                             

    # Wing Loading Based on Stall Speed (Vs)
    WL0 = 0.5*rho0*(Vs**2)*CLmax 

    # Power Loading Based on Max Speed (Vmax)
    sigma = rho/rho0
    a = 0.5*rho0*Cd0
    n = (2*k)/(rho*sigma)
    Vmax = 1.2*Vc
    WL = np.arange(1, 150, 0.01)
    data['WL'] = WL
    data['PL1'] = (hp_max_speed)/((a*(Vmax**3)/WL)+((n*WL)/Vmax)) 

    # Power Loading Based on Take Off Runway (Sto)
    Cdolg = 0.009
    Cdohld_to = 0.005
    Cdoto = Cd0 + Cdolg + Cdohld_to
    Clc = 0.3
    DCl_flap_to = 0
    Clto = Clc + DCl_flap_to
    Cdto = Cdoto + k*(Clto**2)
    mi = 0.07
    Cdg = Cdto - mi*Clto
    Vr = 1.1*Vs                                                 # Vr is the aircraft speed at rotation, which is about 1.1Vs to 1.2Vs
    Clr = (2*mtom*g)/(rho*S*(Vr**2))
    Vto = 1.2*Vs
    vector = np.vectorize(math.exp)
    x = vector(0.6*rho*g*Cdg*Sto*(1/WL))
    y = vector(0.6*rho*g*Cdg*Sto*(1/WL))
    data['PL2'] = ((1-x)*hp_take_off)/((mi-(mi+(Cdg/Clr))*y)*Vto)

    # Power Loading Based on Rate Off Climb (ROC)
    vector = np.vectorize(math.sqrt)
    x = vector(2*WL/(rho0*math.sqrt((3*Cd0)/k)))
    ROC = 1.67 # 0.762                                                 # The rate of climb run is required to be greater than 150 fpm (or 2.5ft/s) at sea level or 0.762 m/sec
    data['PL3']=1/((ROC/hp_climb)+x*(1.155/(LoverDmax*hp_climb)))

    # Power Loading Based on Ceiling
    h = 110 
    sigma_ac = (1-6.873*(10**(-6))*h)**4.26                     # ac=absolute ceiling
    z = math.sqrt((3*Cd0)/k)
    vector = np.vectorize(math.sqrt)
    f = vector(2*WL/(rho*z))
    data['PL4'] = sigma_ac/(f*(1.155/(LoverDmax*hp_ceiling)))

    # find intersection
    common_area_start = np.max(np.minimum.reduce([data['PL1'], data['PL2'], data['PL3'], data['PL4']]))
    opt_PL = common_area_start

    # Retrieve the corresponding opt_WL value at the opt_index
    row, col = np.where(data.values == common_area_start)
    opt_WL = data.loc[row[0], 'WL']

    common_area1 = np.minimum(data['PL1'], data['PL2'])  
    common_area2 = np.minimum(data['PL1'], data['PL3'])  
    common_area3 = np.minimum(data['PL1'], data['PL4'])  
    common_area4 = np.minimum(data['PL2'], data['PL3'])
    common_area5 = np.minimum(data['PL2'], data['PL4'])
    common_area6 = np.minimum(data['PL3'], data['PL4'])

    # Find the minimum common area among the three pairs
    common_area = np.minimum.reduce([common_area1, common_area2, common_area3, common_area4, common_area5, common_area6])

    # Find the Vertical line
    WL0 = 0.5 * rho0 * (Vs ** 2) * CLmax  

    # Design values of Wing Loading and Power Loading
    PL0 = mtow/P_design
    wl_design = WL0 
    pl_design = PL0

    # # Real Values of Wing Loading and Power Loading for our given data.
    # PLreal = mtow/Pmax
    # WLreal = WL0

    max_value_2 = data['PL4'].max()
    
    ratio1 = opt_PL/opt_WL
    y_coord1 = ratio1*wl_design
    last_PL2_value = data['PL2'].iloc[-1]                               # it is the end of the graph for Take Off runway
    last_PL2_value = 1.4*last_PL2_value                                 # 1.4 is a correction because it is not a straight line but a paravola
    # print(last_PL2_value)
    ratio2 = (last_PL2_value-opt_PL)/(np.max(WL)-opt_WL)
    y_coord2 = opt_PL + ratio2*wl_design
    valid_1 = 0
    # valid_2 = 0
    valid = "zero"
    
    # # Create the plot
    # plt.xlabel("W/S (N/m2)")
    # plt.ylabel("W/P (N/W)")
    # plt.plot(WL, data['PL1'], label="Vmax")
    # plt.plot(WL, data['PL2'], label="Take Off run")
    # plt.plot(WL, data['PL3'], label="Rate of Climb")
    # plt.plot(WL, data['PL4'], label="Ceiling")
    # plt.scatter(opt_WL, opt_PL, color='red', label='Optimal Point')
    # plt.fill_between(WL, common_area, color='blue', alpha=0.5, label='Permissible Area')
    # plt.vlines(x = WL0, ymin = 0, ymax = max_value_2, colors='magenta', label='Vstall')
    # plt.hlines(y = PL0, xmin = 0, xmax = WL0, colors='yellow', label='Available Power of given motor for MTOW')    # Only if the engine is known
    
    # For the design point
    if 0 < wl_design <= np.max(WL) and wl_design <= WL0 and 0 < pl_design < opt_PL:
        if (wl_design <= opt_WL):
            if pl_design<=y_coord1: 
                plt.scatter(wl_design, pl_design, color='yellow', label='Design Point (Inside Permissible Area)')
                valid_1 = 1
            else:
                plt.scatter(wl_design, pl_design, color='orange', label='Design Point (Outside Permissible Area)')
        else:
            if pl_design <= y_coord2:
                plt.scatter(wl_design, pl_design, color='yellow', label='Design Point (Inside Permissible Area)')
                valid_1 = 1
            else:
                plt.scatter(wl_design, pl_design, color='orange', label='Design Point (Outside Permissible Area)')                  
    else:
        plt.scatter(wl_design, pl_design, color='orange', label='Design Point (Outside Permissible Area)')
        
    # # For the real point    
    # if 0 < WLreal <= np.max(WL) and WLreal <= WL0 and 0 < PLreal < opt_PL:    
    #     if (WLreal <= opt_WL): 
    #         if PLreal<=y_coord1: 
    #             plt.scatter(WLreal, PLreal, color='magenta', label='Real Design Point (Inside Permissible Area)')
    #             valid_2 = 1
    #         else:
    #             plt.scatter(WLreal, PLreal, color='blue', label='Real Design Point (Outside Permissible Area)')
    #     else:
    #         if PLreal <= y_coord2:
    #             plt.scatter(WLreal, PLreal, color='magenta', label='Real Design Point (Inside Permissible Area)')
    #             valid_2 = 1
    #         else:
    #             plt.scatter(WLreal, PLreal, color='blue', label='Real Design Point (Outside Permissible Area)')                  
    # else:
    #     plt.scatter(WLreal, PLreal, color='blue', label='Real Design Point (Outside Permissible Area)')
    
    # leg = plt.legend(loc='upper center')
    # plt.show()

    if valid_1 == 1: # +valid_2 == 2:
        valid = "Inside Permissable Area"
    else:
        valid = "Outside Permissable Area"

    return (opt_WL, wl_design, opt_PL, pl_design, valid)
