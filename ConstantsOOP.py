import numpy as np

class Const:
    pi = 3.14
    g = 9.81                                                                    # gravity:  m/sec2
    rho0 = 1.225                                                                # air density at sea level: kg/m3  

    ##############################  AERODYNAMIC REQUIREMENTS   #################################
    CLmax = 1.2
    Vs = 10                                                                     # Vstall (m/sec)
    Cd0 = 0.035                                                                 # 
    e = 0.7                                                                     #
    rho = 1.213                                                                 # kg/m3 Air density at 100m
    hp_max_speed = 0.7                                                          # hp_max_speed: prop efficiency at max speed is [0.7 - 0.85]
    hp_take_off = 0.5                                                           # ηp_take_off: prop efficiency for a fixed-pitch propeller at take_off accelerating motion is: 0.5
    hp_climb = 0.6                                                              # ηp_climb: prop efficiency for a fixed-pitch propeller at climbing is [0.5 - 0.7]
    hp_ceiling = 0.75  
    hp_cruise = 0.8                                                             # hp_cruise: prop efficiency for a fixed-pitch propeller at cruise is [0.7 - 0.85]
    Vc = 20                                                                     # m/sec
    Sto = 50                                                                    # (1.65*Wto)/(rho*g*S*Cdg)*np.log((PL2-μ)/(PL2-μ-(Cdg/Clr)))
    LoverDmax = 4 
    Pmax = 400
    P_design = Pmax * 0.75                                                      # P_design = P at cruise is calculated at [0.75 - 0.8] of the max Power
    D_CL_HLD = 0.2                                                              # HLD contribution
    a_2d: float = 6.3                                                          # Must be revised. It is different for every airfoil!!!!!!!!!!!!!!!!!!!!!!!!! BOTH WING AND TAIL     
    alpha_0: float = -1.5                                                      # Flap up zero lift AoA. Must be revised. It is different for every airfoil!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    #######################################  FLAP PARAMETERS   ######################################
    flap_chord_ratio = 0.2
    df = 30                                                                    # assuming flap deflection during take off = 30 degrees
    a_0_fd =  -1.15*flap_chord_ratio*df                                        # Flap down zero lift AoA
    bf_b = 0.2

#########################################     WING PARAMETERS    ####################################
    # AR = 12
    # S = 0.45
    N = 10
    # taper = 1
    taper_v = 0.7
    # alpha_twist = -3
    # S_final = S
    # i_w = -2    # WING SETTING ANGLE CRUISE
    # payload = 2
    # dihedral = 2 # deg
    # dihedral_rad = dihedral*np.pi/180
    # dihedral_h = 2
    # dihedral_h_rad = dihedral_h*np.pi/180
    dihedral_v = 0
    dihedral_v_rad = dihedral_v*np.pi/180
    # sweep = 0 # deg
    sweep_v = 0
    # sweep_rad = sweep*np.pi/180
    sweep_v_rad = sweep_v*np.pi/180
    tail_airfoil = "naca0010"

#########################################    genetic parameters  ########################################
    payload_ref = 10
    L_over_d_ref = 30
    wing_area_ref = 0.1
    total_distance_ref = 0.1
    lb = [0.1, 0.3, 6.0, 0.3, -2.0, -2.0, 0.0, 0.0]
    ub = [2.0, 0.8, 12.0, 1.0, 2.0, 0.0, 0.5, 0.5]
    
    names = ("mtom", "mtow",                                                                                                                                # weight_charact
            "opt_WL", "wl_design", "opt_PL", "pl_design",                                                                                                   # wing_engine sizing
            "Cl_i", "Cl_max_gross", "Cl_net_max", "CL_cr", "CLw_cr", "CL_to_req", "CL_to_fl_req", "optimal_airfoil", "first_wing_incidence_angle",          # wing_charact
            "a_stall_optimal_airfoil", "y_calc", "CL_calc", "CLw_calc", "b_calc", "Lw", "MAC", "Croot", "total_distance",                                   # wing_charact
            "CL_TO_calc", "i_w_set_take_off", "fuselage_pitching_angle", "b_flap", "c_flap",                                                                # tail_charact
            "tail_efficiency", "h_tail", "ho_tail", "Sh", "Vh", "CLh", "a_h", "l_opt_h", "L_fus", "ARh", "a_set_eng", "taper_h",
            "CL_calc_tail", "b_calc_tail", "CLw_calc_tail", "L_tail", "MAC_tail", "Croot_tail",
            "l_opt_v", "Vv", "Sv", "i_v", "ARv", "bv_calc", "Croot_v")
