import numpy as np
import pandas as pd
from functools import reduce
import math
import matplotlib.pyplot as plt
import itertools
from scipy.linalg import solve
from scipy.interpolate import interp1d
from ConstantsOOP import Const
import aerosandbox as asb
import aerosandbox.numpy as np
import aerosandbox.tools.pretty_plots as p
from typing import Union
from shutil import which
import random
from deap import base, creator, tools, algorithms

class WingEngineSizing: 

    def __init__(self, opt_WL, wl_design, opt_PL, pl_design, flag_1):
        
        self.opt_WL = opt_WL
        self.wl_design = wl_design
        self.opt_PL = opt_PL
        self.pl_design = pl_design
        self.flag_1 = flag_1

    def wing_and_engine_sizing(mtom, mtow, wing_area_ga, P_design, aspect_ratio_ga):

        """ 
        Inside this function, we define the permitable limits base on the combination of our demands.
        We are also able to find the optimal design point, if the motor is NOT already known but can be chosen.
        
        Returns:
            float:  opt_WL, wl_design, opt_PL, pl_design
            string: valid
        """

        # General
        data = pd.DataFrame()
        k = 1/(Const.pi*Const.e*aspect_ratio_ga)                                                             

        # Wing Loading Based on Stall Speed (Vs)
        WL0 = 0.5*Const.rho0*(Const.Vs**2)*Const.CLmax 

        # Power Loading Based on Max Speed (Vmax)
        sigma = Const.rho/Const.rho0
        a = 0.5*Const.rho0*Const.Cd0
        n = (2*k)/(Const.rho*sigma)
        Vmax = 1.2*Const.Vc
        WL = np.arange(1, 150, 0.01)
        data['WL'] = WL
        data['PL1'] = (Const.hp_max_speed)/((a*(Vmax**3)/WL)+((n*WL)/Vmax)) 

        # Power Loading Based on Take Off Runway (Sto)
        Cdolg = 0.009
        Cdohld_to = 0.005
        Cdoto = Const.Cd0 + Cdolg + Cdohld_to
        Clc = 0.3
        DCl_flap_to = 0
        Clto = Clc + DCl_flap_to
        Cdto = Cdoto + k*(Clto**2)
        mi = 0.07
        Cdg = Cdto - mi*Clto
        Vr = 1.1*Const.Vs                                                 # Vr is the aircraft speed at rotation, which is about 1.1Vs to 1.2Vs
        Clr = (2*mtom*Const.g)/(Const.rho*wing_area_ga*(Vr**2))
        Vto = 1.2*Const.Vs
        vector = np.vectorize(math.exp)
        x = vector(0.6*Const.rho*Const.g*Cdg*Const.Sto*(1/WL))
        y = vector(0.6*Const.rho*Const.g*Cdg*Const.Sto*(1/WL))
        data['PL2'] = ((1-x)*Const.hp_take_off)/((mi-(mi+(Cdg/Clr))*y)*Vto)

        # Power Loading Based on Rate Off Climb (ROC)
        vector = np.vectorize(math.sqrt)
        x = vector(2*WL/(Const.rho0*math.sqrt((3*Const.Cd0)/k)))
        ROC = 1.67 # 0.762                                                 # The rate of climb run is required to be greater than 150 fpm (or 2.5ft/s) at sea level or 0.762 m/sec
        data['PL3']=1/((ROC/Const.hp_climb)+x*(1.155/(Const.LoverDmax*Const.hp_climb)))

        # Power Loading Based on Ceiling
        h = 110 
        sigma_ac = (1-6.873*(10**(-6))*h)**4.26                     # ac=absolute ceiling
        z = math.sqrt((3*Const.Cd0)/k)
        vector = np.vectorize(math.sqrt)
        f = vector(2*WL/(Const.rho*z))
        data['PL4'] = sigma_ac/(f*(1.155/(Const.LoverDmax*Const.hp_ceiling)))

        # find intersection
        common_area_start = np.max(reduce(np.minimum, [data['PL1'], data['PL2'], data['PL3'], data['PL4']]))
        opt_PL = common_area_start

        # Retrieve the corresponding opt_WL value at the opt_index
        opt_index = np.argwhere(data.values == common_area_start)
        if opt_index.size > 0:
            row = opt_index[0]
        else:
            row = None
        opt_WL = data.loc[row[0], 'WL']     # Select the first occurance

        common_area1 = np.minimum(data['PL1'], data['PL2'])  
        common_area2 = np.minimum(data['PL1'], data['PL3'])  
        common_area3 = np.minimum(data['PL1'], data['PL4'])  
        common_area4 = np.minimum(data['PL2'], data['PL3'])
        common_area5 = np.minimum(data['PL2'], data['PL4'])
        common_area6 = np.minimum(data['PL3'], data['PL4'])

        # Find the minimum common area among the three pairs
        common_area = reduce(np.minimum, [common_area1, common_area2, common_area3, common_area4, common_area5, common_area6])

        # Find the Vertical line
        WL0 = 0.5 * Const.rho0 * (Const.Vs ** 2) * Const.CLmax  

        # Design values of Wing Loading and Power Loading
        PL0 = mtow/P_design
        wl_design = WL0 
        pl_design = PL0

        max_value_2 = data['PL4'].max()

        ratio1 = opt_PL/opt_WL
        y_coord1 = ratio1*wl_design
        last_PL2_value = data['PL2'].iloc[-1]                               # it is the end of the graph for Take Off runway
        last_PL2_value = 1.4*last_PL2_value                                 # 1.4 is a correction because it is not a straight line but a paravola

        # Calculate ratio2 and y_coord2, handling division by zero
        if np.max(WL) != opt_WL:
            ratio2 = (last_PL2_value - opt_PL) / (np.max(WL) - opt_WL)
        else:
            ratio2 = 0  # or any other appropriate value
        y_coord2 = opt_PL + ratio2*wl_design
        valid_1 = 0
        # valid_2 = 0
        flag_1 = 0                 # Point outside permissable area. flag_1 = 1: point inside permissable area.

        # Create the plot
        plt.xlabel("W/S (N/m2)")
        plt.ylabel("W/P (N/W)")
        plt.plot(WL, data['PL1'], label="Vmax")
        plt.plot(WL, data['PL2'], label="Take Off run")
        plt.plot(WL, data['PL3'], label="Rate of Climb")
        plt.plot(WL, data['PL4'], label="Ceiling")
        plt.scatter(opt_WL, opt_PL, color='red', label='Optimal Point')
        plt.fill_between(WL, common_area, color='blue', alpha=0.5, label='Permissible Area')
        plt.vlines(x = WL0, ymin = 0, ymax = max_value_2, colors='magenta', label='Vstall')
        plt.hlines(y = PL0, xmin = 0, xmax = WL0, colors='yellow', label='Available Power of given motor for MTOW')    # Only if the engine is known

        # For the design point
        if 0 < wl_design <= np.max(WL) and wl_design <= WL0 and 0 < pl_design < opt_PL:
            if (wl_design <= opt_WL):
                if pl_design<=y_coord1: 
                    plt.scatter(wl_design, pl_design, color='yellow', label='Design Point (Inside Permissible Area)')
                    flag_1 = 1
                else:
                    plt.scatter(wl_design, pl_design, color='orange', label='Design Point (Outside Permissible Area)')
            else:
                if pl_design <= y_coord2:
                    plt.scatter(wl_design, pl_design, color='yellow', label='Design Point (Inside Permissible Area)')
                    flag_1 = 1
                else:
                    plt.scatter(wl_design, pl_design, color='orange', label='Design Point (Outside Permissible Area)')                  
        else:
            plt.scatter(wl_design, pl_design, color='orange', label='Design Point (Outside Permissible Area)')

        leg = plt.legend(loc='upper center')
        # plt.show()

        return opt_WL, wl_design, opt_PL, pl_design, flag_1

class Wing:

    def __init__(self, Cl_i, Cl_max_gross, Cl_net_max, CL_cr, CL_to_req, CL_to_fl_req, optimal_airfoil, first_wing_incidence_angle, a_stall_optimal_airfoil, y_calc, CL_calc, CLw_calc, b_calc, y_el, CL_el, CL_TO_calc, i_set):
        self.Cl_i = Cl_i
        self.Cl_max_gross = Cl_max_gross
        self.Cl_net_max = Cl_net_max
        self.CL_cr = CL_cr
        self.CL_to_req = CL_to_req
        self.CL_to_fl_req = CL_to_fl_req
        self.optimal_airfoil = optimal_airfoil
        self.first_wing_incidence_angle = first_wing_incidence_angle
        self.a_stall_optimal_airfoil = a_stall_optimal_airfoil
        self.y_calc = y_calc
        self.CL_calc = CL_calc
        self.CLw_calc = CLw_calc
        self.b_calc = b_calc
        self.y_el = y_el
        self.CL_el = CL_el
        self.CL_TO_calc = CL_TO_calc
        self.i_set = i_set

    # instance method
    def airfoil_calc(mtow, wing_area_ga):
        
        """ 
        Calculates basic aerodynamic lift coef. for the aircraft,
        the wing and the airfoil. Those are used in def: "airfoil-selection"
        to choose the best airfoil for our needs.

        Returns:
            float: Cl_i, Cl_max_gross, Cl_net_max, CL_cr, CL_to_req, CL_to_fl_req
        """
        CL_cr = 2*mtow/(Const.rho*(Const.Vc**2)*wing_area_ga)                # CL aircraft required
        CLw_cr = CL_cr/0.95
        Cl_i = CLw_cr/0.9
        CL_max = 2*mtow/(Const.rho0*(Const.Vs**2)*wing_area_ga)
        CL_w_max = CL_max/0.95
        Cl_max_gross = CL_w_max/0.9                   # Max Cl for the airfoil (flap-up) 
        Cl_net_max = Cl_max_gross - Const.D_CL_HLD          # Max Cl for the airfoil (flap down)
        CL_to_req = 0.85*2*mtow/(Const.rho0*(Const.Vc**2)*wing_area_ga)      # CL_aircraft_take_off_required
        V_takeoff = 1.2*Const.Vs
        CL_to_fl_req = 2*mtow/(Const.rho0*(V_takeoff**2)*wing_area_ga) # CL_aircraft_take_off_required_with_flap
        
        return  Cl_i, Cl_max_gross, Cl_net_max, CL_cr, CLw_cr, CL_to_req, CL_to_fl_req
    
    # instance method
    def airfoil_selection(Cl_i, Cl_max_gross, Cl_net_max):

        """ 
        Obtain airfoil data from the file:
        'D:\AviaCode\Airfoil Selection\Paper\Table S2. Main findings.xlsx'.
        Based on the airfoil database and the aerodynamic calculations of 
        the def: airfoil_calc, we choose airfoil that fits better our needs.
        We examine two cases:
        1. HLD not employed, 
        2. HLD employed.

        - column 9:  Airfoil max lift coef, flap up, Re = 100.000 take off
        - column 15: Airfoil ideal lift coef, flap up, Re = 200.000 cruise    
        - column 33: Airfoil max lift coef, flap down, Re = 100.000 landing, take off
        
        For more info refer to 'Fast Airfoil Selection' paper.
    
        Returns:
            string: optimal airfoil
            float:  first_wing_incidence_angle, a_stall_optimal_airfoil
        """

        # file_path = "/home/mike/Desktop/Aviacode/Table S2. Main findings.xlsx"
        file_path = r'D:\AviaCode\Airfoil Selection\Paper\Table S2. Main findings.xlsx'

        # Read the Excel file skipping header rows
        airfoil_data = pd.read_excel(file_path, header=None, skiprows=6)

        ######################################## Flap - Up ################################################

        # Specify the columns of interest
        columns_ind1 = [9, 15]                                                                      
        selected_columns1 = airfoil_data[columns_ind1]
        
        # # Plot the data
        # plt.scatter(selected_columns1[9], selected_columns1[15])

        # # Add labels and title
        # plt.xlabel('Airfoil gross max lift coefficient (Clmax) Re=100.000, zero flap deflection')
        # plt.ylabel('Ideal lift coefficient (Cli) Re=200')
        # plt.title('Fast Airfoil Selection - Flap Up')
        # plt.plot(Cl_max_gross, Cl_i, color='red', marker='x')
        # plt.show()

        ######################################## Flap - Down ################################################

        # Specify the columns of interest
        columns_ind2 = [33, 15]                                                                      
        selected_columns2 = airfoil_data[columns_ind2]    

        # # Plot the data
        # plt.scatter(selected_columns2[33], selected_columns2[15])

        # # Add labels and title
        # plt.xlabel('Airfoil net maximum lift coefficient (Clmax) Re=100, 30 deg flap deflection')
        # plt.ylabel('Ideal lift coefficient (Cli) Re=200')
        # plt.title('Fast Airfoil Selection - Flap Down')
        # plt.plot(Cl_net_max, Cl_i, color='green', marker='+')
        # plt.show()

        # Sum of distances for both cases. The min sum corresponds to the optimal airfoil
        airfoil_data['Score'] = (np.sqrt((selected_columns1[9] - Cl_max_gross)**2 + (selected_columns1[15] - Cl_i)**2) + 
                                np.sqrt((selected_columns2[33] - Cl_net_max)**2 + (selected_columns2[15] - Cl_i)**2))
        min_distance_index = np.argmin(airfoil_data['Score'])
        optimal_airfoil= airfoil_data.iloc[min_distance_index, 1]
        # print(f"Optimal Airfoil is: ", optimal_airfoil)

        # Locate the angle where cl_ideal occurs. That is at the same time the wing_incidence_ angle first aproximation for def: wing_variables_exploration
        first_wing_incidence_angle = airfoil_data.iloc[min_distance_index, 16]

        # find stall angle for chosen airfoil
        a_stall_optimal_airfoil = airfoil_data.iloc[min_distance_index, 8]

        return optimal_airfoil, first_wing_incidence_angle, a_stall_optimal_airfoil

    # instance method
    def wing_aerodynamics_calculator(wing_area_ga, aspect_ratio_ga, taper_ratio_ga, wing_incidence_ga, alpha_twist_ga, CLw_cr):

        """ 
        Application of the Lifting Line Theory to calculate the Wing Lift Coefficient.

        Returns:
            float: y_calc, CL_calc, CLw_calc, b_calc
        """
        
        b_calc = math.sqrt(aspect_ratio_ga * wing_area_ga)                      # wing spaN(m) (includes both wings)
        # MAC = Const.S / b_calc                                                # MeaNAerodynamic Chord (m)
        MAC = b_calc/aspect_ratio_ga
        Croot = (1.5 * (1 + taper_ratio_ga) * MAC) / (1 + taper_ratio_ga + taper_ratio_ga ** 2)  # root chord (m)

        theta = np.linspace((math.pi / (2 * Const.N)), (math.pi / 2), Const.N)
        alpha = np.linspace(wing_incidence_ga + alpha_twist_ga, wing_incidence_ga, Const.N)
        z = (b_calc / 2) * np.cos(theta)
        c = Croot * (1 - (1 - taper_ratio_ga) * np.cos(theta))  # Mean Aerodynamics
        mu = c * Const.a_2d / (4 * b_calc)

        LHS = mu * (np.array(alpha) - Const.alpha_0) / 57.3  # .reshape((N-1),1)# Left Hand Side

        RHS = []
        for i in range(1, 2 * Const.N + 1, 2):
            RHS_iter = np.sin(i * theta) * (1 + (mu * i) / (np.sin(list(theta))))  # .reshape(1,N)
            # print(RHS_iter,"RHS_iter shape")
            RHS.append(RHS_iter)

        test = np.asarray(RHS)
        x = np.transpose(test)
        inv_RHS = np.linalg.inv(x)

        ans = np.matmul(inv_RHS, LHS)

        mynum = np.divide((4 * b_calc), c)

        CL = (np.sin((1) * theta)) * ans[0] * mynum
        for i in range(1,Const.N-1):
            CL = CL + (np.sin((2*i+1) * theta)) * ans[i] * mynum
            
        CL_calc = np.append(0, CL)

        y_calc = np.concatenate(([b_calc / 2], z))

        # The Available Wing Lift Coefficient for Given Geometrical Parameters
        CLw_calc = (math.pi * aspect_ratio_ga * ans[0])
        Lw = 0.5*Const.rho*(Const.Vc**2)*wing_area_ga*CLw_calc
        # print(f"The Available Wing Lift Coef. is: ", CLw_calc,
        # f"The Available Wing Lift is: ", Lw )

        flag_2 = 0

        if CLw_calc >= CLw_cr:
            flag_2 = 1

        return  y_calc, CL_calc, CLw_calc, b_calc, Lw, MAC, Croot, flag_2

    def elliptical_lift_distribution(b_calc, CL_calc, y_calc):
        CL_el_max = CL_calc[Const.N]
        CL_el = CL_el_max * np.sqrt(1 - (y_calc**2 / ((b_calc/2)**2)))
        return CL_el

    def calculate_spacing(y_calc, CL_calc, CL_el):
        y_actual, CL_actual = y_calc, CL_calc
        y_elliptical, CL_elliptical = y_calc, CL_el
        
        # Interpolate the actual lift distribution onto the spanwise locations of the elliptical lift distribution
        interp_CL_actual = interp1d(y_actual, CL_actual, kind='linear', fill_value="extrapolate")
        CL_actual_interp = interp_CL_actual(y_elliptical)
        
        # Calculate the absolute distance between the curves at each spanwise location
        distances = np.abs(CL_actual_interp - CL_elliptical)
        
        # Calculate the total absolute distance
        total_distance = np.sum(distances)
        
        return total_distance

    # instance method
    def llt_flap(a_stall_optimal_airfoil, Croot, CL_to_fl_req, wing_area_ga, aspect_ratio_ga, taper_ratio_ga, wing_incidence_ga, alpha_twist_ga):

        """
        It is a function that calculates the  wing-generated take-off lift coefficient, while flaps are employed.
        All values are calculated based on deflection of 30 degrees of flaps.
        It receives the wing design parameters as well as the required take-off lift coefficient of the aircraft,
        that was calculated in "airfoil_calc" function. As i_set we define as: (i_set = a_stall_optimal_airfoil-2).

        Returns:
            float: CL_TO_calc, i_set
        """

        
        i_w_set_take_off = a_stall_optimal_airfoil-2    # This is 2 the wing angle of attack during take-off. This angle is assumed to be as high as possible. Based on the airfoil stall angle.
        N_flap = 10
        flap_df = pd.DataFrame(columns=['Aspect Ratio', 'Taper Ratio', 'Alpha Twist', 'Wing Set Angle', 'CL_flap_calc'])

        b = np.sqrt(aspect_ratio_ga * wing_area_ga)  # wing span
        # MAC = Const.S / b  # Mean Aerodynamic Chord
        MAC = b/aspect_ratio_ga
        Croot = (1.5 * (1 + taper_ratio_ga) * MAC) / (1 + taper_ratio_ga + taper_ratio_ga ** 2)  # root chord

        theta = np.linspace(np.pi / (2 * N_flap), np.pi / 2, N_flap)
        alpha = np.linspace(i_w_set_take_off + alpha_twist_ga, i_w_set_take_off, N_flap)  
        # alpha = i_w + alpha_twist - np.arange(0, N_flap) * alpha_twist / (N_flap - 1)

        # Segment’s angle of attack
        alpha_zero = np.zeros(N_flap)
        for i in range(N_flap):
            if i / N_flap > (1 - Const.bf_b):
                alpha_zero[i] = Const.a_0_fd  # flap down zero lift AOA
            else:
                alpha_zero[i] = Const.alpha_0     # flap up zero lift AOA

        z = (b / 2) * np.cos(theta)
        c = Croot * (1 - (1 - taper_ratio_ga) * np.cos(theta))  # MAC at each segment
        mu = c * Const.a_2d / (4 * b)
        LHS = mu * (alpha - Const.alpha_0) / 57.3  # Left Hand Side

        # Solving N equations to find coefficients A(i):
        B = np.zeros((N_flap, N_flap))
        for i in range(N_flap):
            for j in range(N_flap):
                B[i, j] = np.sin((2 * j + 1) * theta[i]) * (1 + (mu[i] * (2 * j + 1)) / np.sin(theta[i]))

        A = solve(B, LHS)

        sum1 = np.zeros(N_flap)
        sum2 = np.zeros(N_flap)

        for i in range(N_flap):
            for j in range(N_flap):
                sum1[i] += (2 * j + 1) * A[j] * np.sin((2 * j + 1) * theta[i])
                sum2[i] += A[j] * np.sin((2 * j + 1) * theta[i])

        CL_TO_calc = np.pi * aspect_ratio_ga * A[0]
        # print(f"Produced Take Off Lift Coef with Flaps is: ", CL_TO_calc)

        flag_3 = 0

        if CL_TO_calc >= CL_to_fl_req:
            flag_3 = 1        

        fuselage_pitching_angle = i_w_set_take_off - wing_incidence_ga
        b_flap = b*Const.bf_b
        c_flap = Croot*Const.flap_chord_ratio

        return CL_TO_calc, i_w_set_take_off, fuselage_pitching_angle, b_flap, c_flap, flag_3

class Tail:
    
    def __init__(self, h_tail, ho_tail, Sh, Vh, CLh, a_h, l_opt_h, L):
        self.h_tail = h_tail
        self.ho_tail = ho_tail
        self.Sh = Sh
        self.Vh = Vh
        self.CLh = CLh
        self.a_h = a_h
        self.v = l_opt_h
        self.L = L
    
    def horizontal_tail_calculator(CLw_calc, wing_area_ga, aspect_ratio_ga, MAC):
        
        """ 
        Inside this function, we define the parameters of the horizontal tail.
        
        Returns:
            float:  h_tail, ho_tail, Sh, Vh, CLh, a_h, l_opt_h, L
        """
        
        tail_efficiency = 0.85
        taper_h = 0.7
        Vh = 0.8 # [0.4 - 1.2]
        ho_tail = 0.2 # ho denotes the non-dimensional wing/fuselage aerodynamic center position in relation to MAC [0.2 - 0.25]
        h_tail = 0.25 # The parameter h denotes the non-dimensional aircraft cg position in relation to MAC [0.1 - 0.3]
        wing_pos = 1.8*MAC
        Cmowf = -0.06
        CLh = (Cmowf + CLw_calc*(h_tail - ho_tail ))/Vh # desired tail lift coefficient
        a_h = 0 # Horizontal tail setting angle
        Kc = 1.2 # [1 - 1.4] correction factor
        Df = 0.35  # max fuselage diameter
        l_opt_h = Kc*np.sqrt((4*MAC*wing_area_ga*Vh)/(np.pi*Df))   # horizontal tail moment arm, the distance between the tail aerodynamic center and the aircraft center of gravity
        Sh = Vh*MAC*wing_area_ga/l_opt_h
        L_fus = 0.6/l_opt_h   # l/L = 0.6 The ratio between tail arm and Fuselage Length
        ARh = (2/3)*aspect_ratio_ga    # [3 - 5]
        a_set_eng = 3 # [2 - 4], the typical engine setting angle is about 2–4 deg
        
        return  tail_efficiency, h_tail, ho_tail, Sh, Vh, CLh, a_h, l_opt_h, L_fus, ARh, a_set_eng, taper_h, wing_pos
    
    def tail_aerodynamics_calculator(Sh, a_h, ARh, taper_h, alpha_twist_ga):
        
        """ 
        Application of the Lifting Line Theory to calculate the Wing Lift Coefficient.

        Returns:
            float: CL_calc_tail
        """
        
        
        b_calc_tail = math.sqrt(ARh * Sh)  
        # MAC_tail = Sh / b_calc_tail
        MAC_tail = b_calc_tail/ARh
        Croot_tail = (1.5 * (1 + taper_h) * MAC_tail) / (1 + taper_h + taper_h ** 2)  # root chord (m)

        theta_tail = np.linspace((math.pi / (2 * Const.N)), (math.pi / 2), Const.N)
        alpha_tail = np.linspace(a_h + alpha_twist_ga, a_h, Const.N)
        z_tail = (b_calc_tail / 2) * np.cos(theta_tail)
        c_tail = Croot_tail * (1 - (1 - taper_h) * np.cos(theta_tail))  # Mean Aerodynamics
        mu = c_tail * Const.a_2d / (4 * b_calc_tail)

        LHS = mu * (np.array(alpha_tail) - Const.alpha_0) / 57.3  # .reshape((N-1),1)# Left Hand Side

        RHS = []
        for i in range(1, 2 * Const.N + 1, 2):
            RHS_iter = np.sin(i * theta_tail) * (1 + (mu * i) / (np.sin(list(theta_tail))))  # .reshape(1,N)
            # print(RHS_iter,"RHS_iter shape")
            RHS.append(RHS_iter)

        test = np.asarray(RHS)
        x = np.transpose(test)
        inv_RHS = np.linalg.inv(x)

        ans = np.matmul(inv_RHS, LHS)

        mynum = np.divide((4 * b_calc_tail), c_tail)

        CL = (np.sin((1) * theta_tail)) * ans[0] * mynum
        for i in range(1,Const.N-1):
            CL = CL + (np.sin((2*i+1) * theta_tail)) * ans[i] * mynum
            
        CL_calc_tail = np.append(0, CL)

        y_calc_tail = np.concatenate(([b_calc_tail / 2], z_tail))

        # The Available Wing Lift Coefficient for Given Geometrical Parameters
        CLw_calc_tail = (math.pi * ARh * ans[0])
        L_tail = 0.5*Const.rho*(Const.Vc**2)*Sh*CLw_calc_tail
        # print(f"The Available Wing Lift Coef. is: ", CLw_calc,
        # f"The Available Wing Lift is: ", Lw )
        
        return  CL_calc_tail, b_calc_tail, CLw_calc_tail, L_tail, MAC_tail, Croot_tail

    def vertical_tail_calculator(l_opt_h, wing_area_ga, ARh):
        """ 
            Inside this function, we define the parameters of the vertical tail.
            
        Returns:
            float:  Vv, Sv, i_v
        """

        l_opt_v = l_opt_h # vertical tail moment arm
        Kf = 0.7 # [0.65 - 0.78]
        Cnv = 0.3 # [0.2 - 0.4]
        Vv = 0.075 # [0.05 - 0.1]
        Sv = 0.15 * wing_area_ga # [0.1 - 0.15]
        i_v = 1.5 # [1 - 2] vertical tail incidence to prevent roll from propeller revolution
        ARv = 0.6*ARh #3 # [1 - 2]
        bv_calc = np.sqrt(ARv*Sv)
        # C_mean_v = Sv/bv_calc
        C_mean_v = bv_calc/ARv
        Croot_v = (1.5 * (1 + Const.taper_v) * C_mean_v) / (1 + Const.taper_v + Const.taper_v ** 2)
        
        return l_opt_v, Vv, Sv, i_v, ARv, bv_calc, Croot_v

class Airplane:   

    def __init__(self):
        pass

    @staticmethod    
    def mtow_estimation(payload_ga):
        """ 
        Calculates the estimated mtom for payload based on statistical data.
        Equations used were derived from file "Database.xlsx".

        Returns:
            float: mtom, mtow
        """
        
        # Calculate estimated mtom for a given payload
        if payload_ga <= 10:
            mtom = 3.6522 * math.exp(0.274 * payload_ga)
        elif payload_ga <=100:
            mtom = 40.514 * math.exp(0.0332 * payload_ga)
        elif payload_ga > 100:
            mtom = 714.87 * math.exp(0.0013 * payload_ga)
            
        mtow = mtom*Const.g
        
        return mtom, mtow

    def create_airplane(self, individual): # payload_ga, wing_area_ga, aspect_ratio_ga, taper_ratio_ga, wing_incidence_ga, alpha_twist_ga, alpha_sweep_rad_ga, dihedral_rad_ga):
        """
        It combines the derived data to obtain the full dataset of the airplane configuration_
        
        Returns:
            list: overall_airplane_charact
        """
        payload_ga, wing_area_ga, aspect_ratio_ga, taper_ratio_ga, wing_incidence_ga, alpha_twist_ga, alpha_sweep_rad_ga, dihedral_rad_ga = individual
        flag = 0
        overall_airplane_charact = (0,) * 53    # Initialize a tuple with zeros
        
        mtom, mtow = Airplane.mtow_estimation(payload_ga)
        opt_WL, wl_design, opt_PL, pl_design, flag_1 = WingEngineSizing.wing_and_engine_sizing(mtom, mtow, wing_area_ga, Const.P_design, aspect_ratio_ga)
        
        # Check that the combination of power loading and the wing loading is inside permissable area
        if flag_1 != 0:
            flag += 0.1
            # print(f"flag = ", flag)
        else:
            return overall_airplane_charact, flag
            
        Cl_i, Cl_max_gross, Cl_net_max, CL_cr, CLw_cr, CL_to_req, CL_to_fl_req  = Wing.airfoil_calc(mtow, wing_area_ga)
        optimal_airfoil, first_wing_incidence_angle, a_stall_optimal_airfoil = Wing.airfoil_selection(Cl_i, Cl_max_gross, Cl_net_max)
        y_calc, CL_calc, CLw_calc, b_calc, Lw, MAC, Croot, flag_2 = Wing.wing_aerodynamics_calculator(wing_area_ga, aspect_ratio_ga, taper_ratio_ga, wing_incidence_ga, alpha_twist_ga, CLw_cr)
        
        # Check that the produced lift coef. is greater or equal than the required
        if flag_2 != 0:
            flag += 0.1
            # print(f"flag = ", flag)
        else:
            return overall_airplane_charact, flag
        
        CL_el = Wing.elliptical_lift_distribution(b_calc, CL_calc, y_calc)
        total_distance = Wing.calculate_spacing(y_calc, CL_calc, CL_el)
        CL_TO_calc, i_w_set_take_off, fuselage_pitching_angle, b_flap, c_flap, flag_3 = Wing.llt_flap(a_stall_optimal_airfoil, Croot, CL_to_fl_req, wing_area_ga, aspect_ratio_ga, taper_ratio_ga, wing_incidence_ga, alpha_twist_ga)
        
        # Check that the produced lift coef. ata take off with flaps is greater or equal than the required
        if flag_3 != 0:
            flag += 0.1
            # print(f"flag = ", flag)
        else:
            return individual, flag
        
        tail_efficiency, h_tail, ho_tail, Sh, Vh, CLh, a_h, l_opt_h, L_fus, ARh, a_set_eng, taper_h, wing_pos= Tail.horizontal_tail_calculator(CLw_calc, wing_area_ga, aspect_ratio_ga, MAC)
        CL_calc_tail, b_calc_tail, CLw_calc_tail, L_tail, MAC_tail, Croot_tail = Tail.tail_aerodynamics_calculator(Sh, a_h, ARh, taper_h, alpha_twist_ga)
        l_opt_v, Vv, Sv, i_v, ARv, bv_calc, Croot_v = Tail.vertical_tail_calculator(l_opt_h, wing_area_ga, ARh)
        
        # Combine characteristics
        overall_airplane_charact = {
                                        "mtom": mtom,
                                        "mtow": mtow,
                                        "opt_WL": opt_WL,
                                        "wl_design": wl_design,
                                        "opt_PL": opt_PL,
                                        "pl_design": pl_design,
                                        "Cl_i": Cl_i,
                                        "Cl_max_gross": Cl_max_gross,
                                        "Cl_net_max": Cl_net_max,
                                        "CL_cr": CL_cr,
                                        "CLw_cr": CLw_cr,
                                        "CL_to_req": CL_to_req,
                                        "CL_to_fl_req": CL_to_fl_req,
                                        "optimal_airfoil": optimal_airfoil,
                                        "first_wing_incidence_angle": first_wing_incidence_angle,
                                        "a_stall_optimal_airfoil": a_stall_optimal_airfoil,
                                        "y_calc": y_calc,
                                        "CL_calc": CL_calc,
                                        "CLw_calc": CLw_calc,
                                        "b_calc": b_calc,
                                        "Lw": Lw,
                                        "MAC": MAC,
                                        "Croot": Croot,
                                        "total_distance": total_distance,
                                        "CL_TO_calc": CL_TO_calc,
                                        "i_w_set_take_off": i_w_set_take_off,
                                        "fuselage_pitching_angle": fuselage_pitching_angle,
                                        "b_flap": b_flap,
                                        "c_flap": c_flap,
                                        "tail_efficiency": tail_efficiency,
                                        "h_tail": h_tail,
                                        "ho_tail": ho_tail,
                                        "Sh": Sh,
                                        "Vh": Vh,
                                        "CLh": CLh,
                                        "a_h": a_h,
                                        "l_opt_h": l_opt_h,
                                        "L_fus": L_fus,
                                        "ARh": ARh,
                                        "a_set_eng": a_set_eng,
                                        "taper_h": taper_h,
                                        "CL_calc_tail": CL_calc_tail,
                                        "b_calc_tail": b_calc_tail,
                                        "CLw_calc_tail": CLw_calc_tail,
                                        "L_tail": L_tail,
                                        "MAC_tail": MAC_tail,
                                        "Croot_tail": Croot_tail,
                                        "l_opt_v": l_opt_v,
                                        "Vv": Vv,
                                        "Sv": Sv,
                                        "i_v": i_v,
                                        "ARv": ARv,
                                        "bv_calc": bv_calc,
                                        "Croot_v": Croot_v
                                    }

        return (overall_airplane_charact, flag)

class GeneticAlgorithm:
    def __init__(self, population_size, num_generations, mutation_prob, crossover_prob):
        self.population_size = population_size
        self.num_generations = num_generations
        self.mutation_prob = mutation_prob
        self.crossover_prob = crossover_prob

        # Set up DEAP framework
        self.toolbox = base.Toolbox()
        self.setup_deap()

    def setup_deap(self):
        # Define the problem as a maximization (positive fitness values)
        creator.create("FitnessMax", base.Fitness, weights=(1.0,))
        creator.create("Individual", list, fitness=creator.FitnessMax)

        # Register the necessary functions to the toolbox
        self.toolbox.register("individual", tools.initIterate, creator.Individual, self.create_individual)
        self.toolbox.register("population", tools.initRepeat, list, self.toolbox.individual)
        self.toolbox.register("evaluate", self.evaluate_individual)
        self.toolbox.register("mate", tools.cxBlend, alpha=0.5)
        self.toolbox.register("mutate", tools.mutGaussian, mu=0, sigma=1, indpb=0.2)
        self.toolbox.register("select", tools.selTournament, tournsize=3)

    def create_individual(self):
        payload_ga = np.random.uniform(0.1, 0.5)
        wing_area_ga = np.random.uniform(0.3, 0.5)
        aspect_ratio_ga = np.random.uniform(6.0, 12.0)
        taper_ratio_ga = np.random.uniform(0.3, 1.0)
        wing_incidence_ga = np.random.uniform(-2.0, 2.0)
        alpha_twist_ga = np.random.uniform(-2.0, 0.0)
        alpha_sweep_rad_ga = np.random.uniform(0.0, 0.2)
        dihedral_rad_ga = np.random.uniform(0.0, 0.2)

        individual = payload_ga, wing_area_ga, aspect_ratio_ga, taper_ratio_ga, wing_incidence_ga, alpha_twist_ga, alpha_sweep_rad_ga, dihedral_rad_ga
        
        return individual

    def evaluate_individual(self, individual):
        overall_airplane_charact, flag = self.calculate_individual_aerodynamics(individual)

        score = self.score_function(individual, overall_airplane_charact, flag)
        return (score,)

    def calculate_individual_aerodynamics(self, individual):
        
        # Extract parameters from the individual list
        payload_ga, wing_area_ga, aspect_ratio_ga, taper_ratio_ga, wing_incidence_ga, alpha_twist_ga, alpha_sweep_rad_ga, dihedral_rad_ga = individual

        my_plane = Airplane()  

        overall_airplane_charact, flag = my_plane.create_airplane(individual)
        
        if flag >= 0.3:
            
            # Convert the dictionary to a pandas DataFrame
            df = pd.DataFrame(list(overall_airplane_charact.items()), columns=["Variables", "Values"])

            # Export the DataFrame to an Excel file
            output_file = "airplane_characteristics.xlsx"
            df.to_excel(output_file, index=False)
            
            b_calc = overall_airplane_charact["b_calc"]
            Croot = overall_airplane_charact["Croot"]
            h_tail = overall_airplane_charact["h_tail"]
            MAC = overall_airplane_charact["MAC"]
            b_calc_tail = overall_airplane_charact["b_calc_tail"]
            Sh = overall_airplane_charact["Sh"]
            Croot_tail = overall_airplane_charact["Croot_tail"]
            bv_calc = overall_airplane_charact["bv_calc"]
            Sv = overall_airplane_charact["Sv"]
            Croot_v = overall_airplane_charact["Croot_v"]
            optimal_airfoil = overall_airplane_charact["optimal_airfoil"]
            a_h = overall_airplane_charact["a_h"]
            l_opt_h = overall_airplane_charact["l_opt_h"]
            l_opt_v = overall_airplane_charact["l_opt_v"]

            # Calculations for wing
            d_1 = np.sin(alpha_sweep_rad_ga)*b_calc
            beff = np.sqrt((b_calc**2)-(d_1**2))
            AReff = (beff**2)/wing_area_ga
            C_tip = wing_area_ga*Croot
            d_2 = d_1 - C_tip/2
            dist_1 = Croot/2 + d_2
            dist_2 = np.tan(dihedral_rad_ga)*beff/2
            dist_3 = h_tail*MAC                                        # Aerodynamic Center of the aircraft
            dist_4 = np.sin(wing_incidence_ga)*Croot

            # Calculations for horizontal_tail
            d_1_tail = np.sin(alpha_sweep_rad_ga)*b_calc_tail
            beff_tail = np.sqrt((b_calc_tail**2)-(d_1_tail**2))
            AReff_tail = ((beff_tail/2)**2)/Sh
            C_tip_tail = taper_ratio_ga*Croot_tail
            d_2_tail = d_1_tail - C_tip_tail/2
            dist_1_tail = Croot_tail/2 + d_2_tail
            dist_2_tail = np.tan(dihedral_rad_ga)*beff_tail/2
            
            # Calculations for vertical_tail
            d_1_tail_v = np.sin(alpha_sweep_rad_ga)*bv_calc
            beff_tail_v = np.sqrt((bv_calc**2)-(d_1_tail_v**2))
            AReff_tail_v = ((beff_tail_v/2)**2)/Sv
            C_tip_tail_v = Const.taper_v*Croot_v
            d_2_tail_v = d_1_tail_v - C_tip_tail_v/2
            dist_1_tail_v = Croot_v/2 + d_2_tail_v
            dist_2_tail_v = np.tan(Const.dihedral_v_rad)*beff_tail_v/2
            
            ## The coordinates of the Leading Edge points: A:[0,0,0]  ->  [dist_1, beff, dist_2]
            ## The coordinates of the Trailing Edge points: A:[Croot,0,0]  ->  [dist_1 + Ctip, beff, dist_2]

            # Run Aerosandbox
            wing_airfoil = asb.Airfoil(optimal_airfoil)
            tail_airfoil = asb.Airfoil("naca0010")
                
            airplane = asb.Airplane(
                name="Mike's Airplane",                        
                xyz_ref=[0.25*MAC, 0, 0],  # Reference for moments
                s_ref=wing_area_ga,
                c_ref=MAC,
                b_ref=beff/2,
                wings=[
                    asb.Wing(
                        name="Main Wing",
                        symmetric=True,  # Should this wing be mirrored across the XZ plane?
                        xsecs=[  # The wing's cross ("X") sections
                            asb.WingXSec(  # Root
                                xyz_le=[0-h_tail*Croot, 0, 0.05],  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                                chord=Croot,
                                twist=wing_incidence_ga,  # degrees
                                airfoil=wing_airfoil,  # Airfoils are blended between a given XSec and the next one.
                                # control_surface_is_symmetric=True,
                                # control_surface_deflection=0,  # degrees
                            ),
                            # asb.WingXSec(  # Mid
                            #     xyz_le=[0-h_tail*Croot, 0, 0],
                            #     chord=Croot,
                            #     twist=Const.i_w,
                            #     airfoil=wing_airfoil,
                            #     # control_surface_is_symmetric=True,  # Aileron
                            #     # control_surface_deflection=0,
                            # ),
                            asb.WingXSec(  # Tip
                                xyz_le= [dist_1-h_tail*Croot, beff/2, dist_2+0.05],
                                chord=C_tip,
                                twist=wing_incidence_ga + alpha_twist_ga,
                                airfoil=wing_airfoil,
                            ),
                        ]
                    ),
                    asb.Wing(
                        name="Horizontal Stabilizer",
                        symmetric=True,
                        xsecs=[
                            asb.WingXSec(  # root
                                xyz_le=[0-h_tail*Croot_tail, 0, 0],
                                chord=Croot_tail,
                                twist=a_h,
                                airfoil=tail_airfoil,
                                control_surface_is_symmetric=True,  # Elevator
                                control_surface_deflection=0,
                            ),
                            asb.WingXSec(  # tip
                                xyz_le=[dist_1_tail-h_tail*Croot_tail, beff_tail/2, dist_2_tail],
                                chord=C_tip_tail,
                                twist=alpha_twist_ga + a_h,
                                airfoil=tail_airfoil,
                            )
                        ]
                    ).translate([l_opt_h, 0, 0.01]),
                    asb.Wing(
                        name="Vertical Stabilizer",
                        symmetric=False,
                        xsecs=[
                            asb.WingXSec( # Root
                                xyz_le=[0-Croot_tail/4, 0, 0],
                                chord=Croot_v,
                                twist=0,
                                airfoil=tail_airfoil,
                                control_surface_is_symmetric=True,  # Rudder
                                control_surface_deflection=0,
                            ),
                            asb.WingXSec( # Tip
                                xyz_le=[dist_1_tail_v, 0, beff_tail_v/2],
                                chord=C_tip_tail_v,
                                twist=0,
                                airfoil=tail_airfoil
                            )
                        ]
                    ).translate([l_opt_v, 0, 0.02])
                ],
                    fuselages=[
                        asb.Fuselage(
                            name="Fuselage",
                            xsecs=[
                                asb.FuselageXSec(
                                    xyz_c=[xi * (0.5 + l_opt_h+3/4*Croot_tail) - 0.5, 0, 0],
                                    radius=0.35*asb.Airfoil("naca0024").local_thickness(x_over_c=xi)
                                )
                                for xi in np.cosspace(0, 1, 20)
                        ]
                    )
                ]
            )
            
            alpha = np.linspace(-10, 15, 100)
            velocity = np.linspace(0,30,100)

            aero = asb.AeroBuildup(
                airplane=airplane,
                op_point=asb.OperatingPoint(
                    velocity=Const.Vc,
                    alpha=alpha,
                    beta=0
                ),
            ).run_with_stability_derivatives()

            # Check lengths of all lists
            lengths = [len(v) for v in aero.values()]
            max_length = max(lengths)

            # Pad shorter lists/tuples with None
            for key, value in aero.items():
                if len(value) < max_length:
                    if isinstance(value, tuple):
                        aero[key] = value + (None,) * (max_length - len(value))  # Pad tuples
                    elif isinstance(value, list):
                        aero[key] = value + [None] * (max_length - len(value))  # Pad lists

            # Convert the dictionary to a DataFrame
            aero_df = pd.DataFrame(aero)

            # Add AoA column in the DataFrame
            aero_df['AoA'] = alpha
            
            # Calculate CL/CD ratio for each AoA
            aero_df['CL/CD'] = aero_df['CL'] / aero_df['CD']

            # Find the row with the maximum CL/CD ratio
            max_cl_cd_row = aero_df.loc[aero_df['CL/CD'].idxmax()]

            # Extract the maximum CL/CD ratio
            max_cl_cd_ratio = max_cl_cd_row['CL/CD']

            # Extract the AoA corresponding to the maximum CL/CD ratio
            max_cl_cd_aoa = max_cl_cd_row['AoA']
            
            # Find the index of the row with the maximum CL/CD ratio
            max_cl_cd_index = aero_df['CL/CD'].idxmax()
            
            # Find velocity at max CL/CD
            max_cl_cd_vel = velocity[max_cl_cd_index]
            
            # Export the DataFrame to Excel
            output_path = 'aero.xlsx'
            aero_df.to_excel(output_path, index=False)
            
            overall_airplane_charact["CL/CD max"] = max_cl_cd_ratio
            overall_airplane_charact["AoA CL/CD max"] = max_cl_cd_aoa
            overall_airplane_charact["Vel CL/CD max"] = max_cl_cd_vel

        return overall_airplane_charact, flag

    def score_function(self, individual, overall_airplane_charact, flag):
        if flag >= 0.3:
            score = (np.abs(individual[0]) / Const.payload_ref) + (overall_airplane_charact["CL/CD max"] / Const.L_over_d_ref) + Const.wing_area_ref / individual[1]  # Adjust as per your constants
        else:
            score = flag 
            if flag > 0:
                if (individual[0] >= 0.1) & (individual[0] <= 2):
                    if (individual[1] >= 0.3) & (individual[1] <= 1):
                        if (individual[2] >= 6) & (individual[2] <= 14):
                            if (individual[3] >= 0.3) & (individual[3] <= 1):
                                if (individual[6] >= 0) & (individual[6] <= 0.5):
                                    if (individual[7] >= 0) & (individual[7] <= 0.5):
                                        score = score + np.abs(individual[0]/Const.payload_ref) + Const.wing_area_ref / individual[1]
            
        print ("Flag = ", flag, " and score = ", score)
        return score

    def run(self):
        population = self.toolbox.population(n=self.population_size)

        # Define statistics object to track desired statistics
        stats = tools.Statistics(lambda ind: ind.fitness.values)
        stats.register("avg", np.mean)
        stats.register("std", np.std)
        stats.register("min", np.min)
        stats.register("max", np.max)

        # Create a Logbook to record the statistics
        logbook = tools.Logbook()
        logbook.header = ["gen", "nevals"] + stats.fields

        # List to keep track of the best 100 individuals
        best_individuals = []

        # Evaluate the initial population
        print("Initial Population:")
        for i, ind in enumerate(population):
            ind.fitness.values = self.evaluate_individual(ind)
            print(f"Individual {i}: {ind}, Fitness: {ind.fitness.values[0]}")

        # Record the statistics for the initial population
        record = stats.compile(population)
        logbook.record(gen=0, nevals=len(population), **record)
        print(logbook.stream)

        for gen in range(1, self.num_generations + 1):
            print(f"\nGeneration {gen}:")

            # Perform the evolutionary algorithm (selection, crossover, mutation)
            result_population = algorithms.varAnd(population, self.toolbox, cxpb=self.crossover_prob, mutpb=self.mutation_prob)

            # Check every variable of every individual in the new population
            for individual in result_population:
                for var_idx, var_value in enumerate(individual):
                    if var_value < Const.lb[var_idx] or var_value > Const.ub[var_idx]:
                        # If a variable is out of bounds, set the individual's fitness to a very low value
                        individual.fitness.values = (0.0,)
                        break  # Move on to the next individual

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in result_population if not ind.fitness.valid]
            for i, ind in enumerate(invalid_ind):
                ind.fitness.values = self.evaluate_individual(ind)
                print(f"Evaluating Individual {i} in Generation {gen}: {[round(val, 2) for val in ind]}, Fitness: {ind.fitness.values[0]:.2f}")

            # Update the best individuals list
            best_individuals.extend(result_population)
            best_individuals = sorted(best_individuals, key=lambda ind: ind.fitness.values[0], reverse=True)[:100]

            # Select the next generation population
            population[:] = self.toolbox.select(result_population, len(population))

            # Extract the best individual of this generation
            best_individual = tools.selBest(population, k=1)[0]
            print("Best Individual of Generation:", gen, " is: ", best_individual)
            print("Best Score:", best_individual.fitness.values[0])

            # Record the statistics for the current population
            record = stats.compile(population)
            logbook.record(gen=gen, nevals=len(invalid_ind), **record)
            print(logbook.stream)

        # Export the logbook to a CSV file
        df_logbook = pd.DataFrame(logbook)
        df_logbook.to_csv("evolution_statistics.csv", index=False)

        # Export the best 100 individuals to a CSV file
        df_best_individuals = pd.DataFrame(
            [{"Individual": ind, "Fitness": ind.fitness.values[0]} for ind in best_individuals]
        )
        df_best_individuals.to_csv("best_100_individuals.csv", index=False)

        # Print the best individual found
        print("Best Individual:", best_individual)
        
        # After all generations, return the best individual found
        return best_individual

class AirplaneAnalysis:
    
    def _init_(self):
        pass

    def airplane_analysis(self, best_individual):
        
        # Get Characteristics from best (generic) individual
        generic_payload = best_individual[0]
        generic_wing_area = best_individual[1]
        generic_aspect_ratio = best_individual[2]
        generic_taper_ratio = best_individual[3]
        generic_wing_incidence = best_individual[4]
        generic_alpha_twist = best_individual[5]
        generic_sweep_rad = best_individual[6]
        generic_dihedral_rad = best_individual[7]
        
        # Run analysis for generic individual
        mtom, mtow = Airplane.mtow_estimation(generic_payload)
        opt_WL, wl_design, opt_PL, pl_design, valid = WingEngineSizing.wing_and_engine_sizing(mtom, mtow, generic_wing_area, Const.P_design, generic_aspect_ratio)
        Cl_i, Cl_max_gross, Cl_net_max, CL_cr, CLw_cr, CL_to_req, CL_to_fl_req = Wing.airfoil_calc(mtow, generic_wing_area)
        optimal_airfoil, first_wing_incidence_angle, a_stall_optimal_airfoil = Wing.airfoil_selection(Cl_i, Cl_max_gross, Cl_net_max)
        y_calc, CL_calc, CLw_calc, b_calc, Lw, MAC, Croot, flag_2 = Wing.wing_aerodynamics_calculator(generic_wing_area, generic_aspect_ratio, generic_taper_ratio, generic_wing_incidence, generic_alpha_twist, CLw_cr)
        CL_el = Wing.elliptical_lift_distribution(b_calc, CL_calc, y_calc)
        total_distance = Wing.calculate_spacing(y_calc, CL_calc, CL_el)
        CL_TO_calc, i_w_set_take_off, fuselage_pitching_angle, b_flap, c_flap, flag_3 = Wing.llt_flap(a_stall_optimal_airfoil, Croot, CL_to_fl_req, generic_wing_area, generic_aspect_ratio, generic_taper_ratio,  generic_wing_incidence, generic_alpha_twist)
        tail_efficiency, h_tail, ho_tail, Sh, Vh, CLh, a_h, l_opt_h, L_fus, ARh, a_set_eng, taper_h, wing_pos = Tail.horizontal_tail_calculator(CLw_calc, generic_wing_area, generic_aspect_ratio, MAC)
        CL_calc_tail, b_calc_tail, CLw_calc_tail, L_tail, MAC_tail, Croot_tail = Tail.tail_aerodynamics_calculator(Sh, a_h, ARh, taper_h, generic_alpha_twist)
        l_opt_v, Vv, Sv, i_v, ARv, bv_calc, Croot_v = Tail.vertical_tail_calculator(l_opt_h, generic_wing_area, ARh)

        # Calculations for generic wing
        d_1 = np.sin(generic_sweep_rad)*b_calc
        beff = np.sqrt((b_calc**2)-(d_1**2))
        AReff = (beff**2)/generic_wing_area
        C_tip =generic_taper_ratio*Croot
        d_2 = d_1 - C_tip/2
        dist_1 = Croot/2 + d_2
        dist_2 = np.tan(generic_dihedral_rad)*beff/2
        dist_3 = h_tail*MAC                                        
        dist_4 = np.sin(generic_wing_incidence)*Croot

        # Calculations for generic horizontal_tail
        d_1_tail = np.sin(generic_sweep_rad)*b_calc_tail
        beff_tail = np.sqrt((b_calc_tail**2)-(d_1_tail**2))
        AReff_tail = ((beff_tail/2)**2)/Sh
        C_tip_tail =generic_taper_ratio*Croot_tail
        d_2_tail = d_1_tail - C_tip_tail/2
        dist_1_tail = Croot_tail/2 + d_2_tail
        dist_2_tail = np.tan(generic_dihedral_rad)*beff_tail/2

        # Calculations for generic vertical_tail
        d_1_tail_v = np.sin(Const.sweep_v_rad)*bv_calc
        beff_tail_v = np.sqrt((bv_calc**2)-(d_1_tail_v**2))
        AReff_tail_v = ((beff_tail_v/2)**2)/Sv
        C_tip_tail_v =Const.taper_v*Croot_v
        d_2_tail_v = d_1_tail_v - C_tip_tail_v/2
        dist_1_tail_v = Croot_v/2 + d_2_tail_v
        dist_2_tail_v = np.tan(Const.dihedral_v_rad)*beff_tail_v/2

        # Run Aerosandbox
        wing_airfoil = asb.Airfoil(optimal_airfoil)
        tail_airfoil = asb.Airfoil("naca0010")
            
        airplane = asb.Airplane(
            name="Mike's Airplane",                        
            xyz_ref=[h_tail*MAC, 0, 0],  # Reference for moments
            s_ref=generic_wing_area,
            c_ref=MAC,
            b_ref=beff/2,
            wings=[
                asb.Wing(
                    name="Main Wing",
                    symmetric=True,  # Should this wing be mirrored across the XZ plane?
                    xsecs=[  # The wing's cross ("X") sections
                        asb.WingXSec(  # Root
                            xyz_le=[0-h_tail*Croot, 0, 0.05],  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                            chord=Croot,
                            twist=generic_wing_incidence,  # degrees
                            airfoil=wing_airfoil,  # Airfoils are blended between a given XSec and the next one.
                            # control_surface_is_symmetric=True,
                            # control_surface_deflection=0,  # degrees
                        ),
                        # asb.WingXSec(  # Mid
                        #     xyz_le=[0-h_tail*Croot, 0, 0],
                        #     chord=Croot,
                        #     twist=Const.i_w,
                        #     airfoil=wing_airfoil,
                        #     # control_surface_is_symmetric=True,  # Aileron
                        #     # control_surface_deflection=0,
                        # ),
                        asb.WingXSec(  # Tip
                            xyz_le= [dist_1-h_tail*Croot, beff/2, dist_2+0.05],
                            chord=C_tip,
                            twist=generic_wing_incidence + generic_alpha_twist,
                            airfoil=wing_airfoil,
                        ),
                    ]
                ),
                asb.Wing(
                    name="Horizontal Stabilizer",
                    symmetric=True,
                    xsecs=[
                        asb.WingXSec(  # root
                            xyz_le=[0-h_tail*Croot_tail, 0, 0],
                            chord=Croot_tail,
                            twist=a_h,
                            airfoil=tail_airfoil,
                            control_surface_is_symmetric=True,  # Elevator
                            control_surface_deflection=0,
                        ),
                        asb.WingXSec(  # tip
                            xyz_le=[dist_1_tail-h_tail*Croot_tail, beff_tail/2, dist_2_tail],
                            chord=C_tip_tail,
                            twist=generic_alpha_twist + a_h,
                            airfoil=tail_airfoil,
                        )
                    ]
                ).translate([l_opt_h, 0, 0.01]),
                asb.Wing(
                    name="Vertical Stabilizer",
                    symmetric=False,
                    xsecs=[
                        asb.WingXSec( # Root
                            xyz_le=[0-Croot_tail/4, 0, 0],
                            chord=Croot_v,
                            twist=0,
                            airfoil=tail_airfoil,
                            control_surface_is_symmetric=True,  # Rudder
                            control_surface_deflection=0,
                        ),
                        asb.WingXSec( # Tip
                            xyz_le=[dist_2_tail_v, 0, beff_tail_v/2],
                            chord=C_tip_tail_v,
                            twist=0,
                            airfoil=tail_airfoil
                        )
                    ]
                ).translate([l_opt_v, 0, 0.02])
            ],
                fuselages=[
                    asb.Fuselage(
                        name="Fuselage",
                        xsecs=[
                            asb.FuselageXSec(
                                xyz_c=[xi * (2*MAC + l_opt_h+1/4*Croot_tail) - wing_pos, 0, 0],
                                radius=0.35*asb.Airfoil("naca0024").local_thickness(x_over_c=xi)
                            )
                            for xi in np.cosspace(0, 1, 50)
                    ]
                )
            ]
        )

        # airplane.draw_three_view()
        airplane.draw()

        alpha = np.linspace(-10, 15, 100)
        velocity = np.linspace(0,30,100)

        aero = asb.AeroBuildup(
            airplane=airplane,
            op_point=asb.OperatingPoint(
                velocity=Const.Vc,
                alpha=alpha,
                beta=0
            ),
        ).run_with_stability_derivatives()

        # Check lengths of all lists
        lengths = [len(v) for v in aero.values()]
        max_length = max(lengths)

        # Pad shorter lists/tuples with None
        for key, value in aero.items():
            if len(value) < max_length:
                if isinstance(value, tuple):
                    aero[key] = value + (None,) * (max_length - len(value))  # Pad tuples
                elif isinstance(value, list):
                    aero[key] = value + [None] * (max_length - len(value))  # Pad lists

        # Convert the dictionary to a DataFrame
        aero_df = pd.DataFrame(aero)

        # Add AoA column in the DataFrame
        aero_df['AoA'] = alpha

        # Calculate CL/CD ratio for each AoA
        aero_df['CL/CD'] = aero_df['CL'] / aero_df['CD']

        # Find the row with the maximum CL/CD ratio
        max_cl_cd_row = aero_df.loc[aero_df['CL/CD'].idxmax()]

        # Extract the maximum CL/CD ratio
        max_cl_cd_ratio = max_cl_cd_row['CL/CD']

        # Extract the AoA corresponding to the maximum CL/CD ratio
        max_cl_cd_aoa = max_cl_cd_row['AoA']

        # Find the index of the row with the maximum CL/CD ratio
        max_cl_cd_index = aero_df['CL/CD'].idxmax()

        # Find velocity at max CL/CD
        max_cl_cd_vel = velocity[max_cl_cd_index]

        # Export the DataFrame to Excel
        output_path = 'aero.xlsx'
        aero_df.to_excel(output_path, index=False)

        # print(aero)   
        return aero_df

if __name__ == "__main__":    
# def generic_airplane_creation(self):

    # Create an instance of GeneticAlgorithm with desired parameters
    ga = GeneticAlgorithm(population_size=5, num_generations=5, mutation_prob=0.2, crossover_prob=0.5)

    # Run the genetic algorithm
    best_individual = ga.run()
    print(best_individual)
    
    # Run analysis for the generic airplane
    best_ga_airplane = AirplaneAnalysis()
    best_ga_characteristics = best_ga_airplane.airplane_analysis(best_individual)



