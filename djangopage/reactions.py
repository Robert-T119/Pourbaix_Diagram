import numpy as np

# define useful values
R = 8.31434  # J/mol.K
F = 96484.56  # C/mol

pO2 = 1  # partial pressures, both in atm
pH2 = 1  # for now =1

#at 298k
delta_G_nio3 = -541800
delta_G_nip2 = -46400
delta_G_nio2 = -453100
delta_G_w = -237178
delta_G_nin6p2 = -250000
delta_G_nin4p2 = -192000
delta_G_o2a = 16300
delta_G_n = -26600
delta_G_ng = -16500
delta_G_nhp = -79370


def R1(x, T): #NI(2+)+ 2e- = Ni
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R1 = -delta_G_nip2
    E_theta_1 = -delta_G_R1/(2*F)
    V1 = E_theta_1 + (C2/2)*np.log10(x)
    return V1

def R2(pH, T): #Ni(OH)2 +2H+ +2e- =Ni + 2H2O
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R2 = 2*delta_G_w - delta_G_nio2
    E_theta_2 = -delta_G_R2/(2*F)
    V2 = []
    if type(pH) == list:  # for list(vector), below is only for single value
        for pHval in pH:
            V2.append(E_theta_2 - (C2*2*pHval)/2)
    else:
        V2 = E_theta_2 - (C2*2*pH)/2
    return V2


def R3(pH, nhp, nin4p2, T): #Ni(NH3)4(2+) + 4H+ +2e- = Ni + 4NH4+
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R3 = 4*delta_G_nhp - delta_G_nin4p2
    E_theta_3 = - delta_G_R3/(2*F)
    V3 = []
    if type(pH) == list:  # for list(vector), below is only for single value
        for pHval in pH:
            V3.append(E_theta_3 - (C2*4*pHval/2) - (C2*4*np.log10(nhp)/2) + (C2*np.log10(nin4p2)/2))
    else:
        V3 = (E_theta_3 - (C2*4*pH/2) - (C2*4*np.log10(nhp)/2) + (C2*np.log10(nin4p2)/2))
    return V3


def R4(pH, nhp, nin6p2, T): #Ni(NH3)+6H+ +2e-=Ni +6NH4+
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R4 = 6*delta_G_nhp - delta_G_nin6p2
    E_theta_4 = -delta_G_R4/(2*F)
    V4 = []
    if type(pH) == list:  # for list(vector), below is only for single value
        for pHval in pH:
            V4.append(E_theta_4 - (C2*6*pHval/2) - (C2*6*np.log10(nhp)/2) + (C2*np.log10(nin6p2)/2))
    else:
        V4 = (E_theta_4 - (C2*6*pH/2) - (C2*6*np.log10(nhp)/2) + (C2*np.log10(nin6p2)/2))
    return V4

def R5(n, nin6p2, T): # Ni(NH3)+6H+ +2e- =Ni + 6NH3
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R5 = 6*delta_G_n - delta_G_nin6p2
    E_theta_6 = -delta_G_R5/(F*2)
    V5 = E_theta_6 - (C2*6*np.log10(n)/2) + (C2*np.log10(nin6p2)/2)
    return V5

def T1(nip2, pH, T): #Ni(OH)3 +3H+ +e- =Ni(2+) +3H2O
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_T1 = 3*delta_G_w + delta_G_nip2 - delta_G_nio3
    E_theta_T1 = -delta_G_T1/(F)
    if type(pH) == list:
        vt1 = []
        for pHval in pH:
            vt1.append(E_theta_T1 - C2*np.log10(nip2) - C2*3*pHval)
    else:
        vt1 = E_theta_T1 - C2*np.log10(nip2) - C2*3*pH
    return vt1


def T2(pH, T): #Ni(OH)3 +H+ +e- =Ni(OH)2 +H2O
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_T2 = delta_G_nio2 + delta_G_w - delta_G_nio3
    E_theta_T2 = -delta_G_T2/F
    if type(pH) == list:
        vt2 = []
        for pHval in pH:
            vt2.append(E_theta_T2 - (C2*1*pHval)/1)
    else:
        vt2 = E_theta_T2 - (C2*1*pH)/1
    return vt2


def T3(nhp, nin4p2, pH, T): #Ni(OH)3 +4NH4+ +e- =Ni(Nh3)4(2+) +3H2O +H+
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_T3 = 3*delta_G_w + delta_G_nin4p2 - 4*delta_G_nhp - delta_G_nio3
    E_theta_T3 = - delta_G_T3/F
    if type(pH) == list:
        vt3 = []
        for pHval in pH:
            vt3.append(E_theta_T3 + C2*pHval + C2*4*np.log10(nhp) - C2*np.log10(nin4p2))
    else:
        vt3 = E_theta_T3 + C2*pH + C2*4*np.log10(nhp) - C2*np.log10(nin4p2)
    return vt3


def T4(pH, nhp, nin6p2, T): #Ni(OH)3 +6NH4+ +e- =Ni(Nh3)6(2+) +3H2O +3H+
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_T4 = 3*delta_G_w + delta_G_nin6p2 - 6*delta_G_nhp - delta_G_nio3
    E_theta_T4 = -delta_G_T4/F
    if type(pH) == list:
        vt4 = []
        for pHval in pH:
            vt4.append(E_theta_T4 + C2*3*pHval + C2*6*np.log10(nhp) - C2*np.log10(nin6p2))
    else:
        vt4 = E_theta_T4 + C2*3*pH + C2*6*np.log10(nhp) - C2*np.log10(nin6p2)
    return vt4


def T5(n, pH, nin6p2, T): #Ni(OH)3 +6NH3 +3H +e- =Ni(Nh3)6(2+) +3H2O
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_T5 = 3*delta_G_w + delta_G_nin6p2 - 6*delta_G_n - delta_G_nio3
    E_theta_T5 = -delta_G_T5/F
    if type(pH) == list:
        vt5 = []
        for pHval in pH:
            vt5.append(E_theta_T5 + 6*C2*np.log10(n) - 3*C2*pHval - C2*np.log10(nin6p2))
    else:
        vt5 = E_theta_T5 + 6*C2*np.log10(n) - 3*C2*pH - C2*np.log10(nin6p2)
    return vt5

def W1(pH, T): #O2+ 4H+ +4e- =2H20
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_W1 = 2*delta_G_w
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
