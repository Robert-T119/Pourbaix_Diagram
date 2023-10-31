import numpy as np

# define useful values
R = 8.31434  # J/mol.K
F = 96484.56  # C/mol

pO2 = 1  # partial pressures, both in atm
pH2 = 1  # for now =1

#at 298.15k
delta_G_cu1p = 50300
delta_G_cu2p = 65700
delta_G_cu2o = -148100
delta_G_cuo = -134000
delta_G_cuoh2 = -359500
delta_G_hcuo2 = -258900
delta_G_cuo2 = -183900
delta_G_cus = -53200
delta_G_cu2s = -87600
delta_G_OH = -157293
delta_G_w = -237178


def R1(cu2p, T): #Cu(2+)+ 2e- = Cu
    C1 = -np.log(10)*R*T #-5707.7
    C2 = -C1/F #0.0592
    delta_G_R1 = -delta_G_cu2p
    E_theta_1 = -delta_G_R1/(2*F)
    V1 = E_theta_1 + (C2/2)*np.log10(cu2p)
    return V1


def R2(pH, T, cu2p): #2cu2+ +H2o+ +2e = Cu2O +2H+
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R2 = delta_G_cu2o - 2*delta_G_cu2p -delta_G_w
    E_theta_2 = -delta_G_R2/(2*F)
    V2 = []
    if type(pH) == list:  # for list(vector), below is only for single value
        for pHval in pH:
            V2.append(E_theta_2 +(C2*2*pHval)/2 +(C2*np.log10(cu2p)))
    else:
        V2 = E_theta_2 +(C2*2*pH)/2 +(C2*np.log10(cu2p))
    return V2


def R3(pH, T): #2CuO +2H+ +2e- =Cu2O +H2O
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R3 = delta_G_cu2o +delta_G_w -2*delta_G_cuo
    E_theta_3 = -delta_G_R3/(2*F)
    V3 = []
    for pHval in pH:
        V3.append(E_theta_3 -(C2*2*pHval)/2)
    return V3


def R4(pH, T, cuo2):  #2CuO2- + 6H+ +2e- =Cu2O +3H2O
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R4 = delta_G_cu2o +3*delta_G_w -2*delta_G_cuo2
    E_theta_4 = -delta_G_R4/(2*F)
    V4 = []
    for pHval in pH:
        V4.append(E_theta_4 -(C2*6*pHval)/2 +(C2*np.log10(cuo2)))
    return V4

def B1(pH, T): #Cu2O +2H+ +2e- =2Cu +H2O
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_T1 = delta_G_w - delta_G_cu2o
    E_theta_T1 = -delta_G_T1/(2*F)
    if type(pH) == list:
        vt1 = []
        for pHval in pH:
            vt1.append(E_theta_T1 - C2*pHval)
    else:
        vt1 = E_theta_T1 - C2*pH
    return vt1

def W1(pH, T): #O2+ 4H+ +4e- =2H20
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltahw1 = -571660
    delta_G_W1_0 = 2*delta_G_w
    delta_G_W1 = (T*delta_G_W1_0/298)+T*deltahw1*((1/T)-(1/298))
    E_theta_W1 = -delta_G_W1/(4*F)
    VW1 = []
    for pHval in pH:
        VW1.append(E_theta_W1 - (C2*4*pHval/4) + (C2*1*np.log10(pO2)/4))
    return VW1


def W2(pH, T):   #2H+ 2e- = H2
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    VW2 = []
    for pHval in pH:
        VW2.append(-(C2*2*pHval/2) - (C2*1/2))
    return VW2


