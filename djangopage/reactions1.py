import numpy as np

# define useful values
R = 8.31434  # J/mol.K
F = 96484.56  # C/mol

pO2 = 1  # partial pressures, both in atm
pH2 = 1  # for now =1

#at 298k
delta_G_nip2 = -46400
delta_G_nio2 = -453100
delta_G_nio3 = -541800
delta_G_H3cit = -1243400
delta_G_H2cit = -1226300
delta_G_Hcit = -1199200
delta_G_cit3 = -1162700
delta_G_Nicit = -1239906.03
delta_G_NiHcit = -1264425.91
delta_G_NiH2cit = -1282683.44
delta_G_OH = -157293
delta_G_w = -237178
delta_G_o2a = 16300
delta_G_o2 = 16300

def R1(nip2, T): #NI(2+)+ 2e- = Ni
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R1 = -delta_G_nip2
    E_theta_1 = -delta_G_R1/(2*F)
    V1 = E_theta_1 + (C2/2)*np.log10(nip2)
    return V1


def R2(pH, T, NiH2cit, H3cit): #NiH2cit +H+ +2e = Ni +H3cit
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R2 = delta_G_H3cit - delta_G_NiH2cit
    E_theta_2 = -delta_G_R2/(2*F)
    V2 = []
    if type(pH) == list:  # for list(vector), below is only for single value
        for pHval in pH:
            V2.append(E_theta_2 -(C2*pHval)/2 +(C2*np.log10(NiH2cit))/2 -(C2*np.log10(H3cit))/2)
    else:
        V2 = E_theta_2 -(C2*pH)/2 +(C2*np.log10(NiH2cit))/2 -(C2*np.log10(H3cit))/2
    return V2


def R3(pH, T, NiHcit, H3cit): #NiHcit +2H+ +2e = Ni +H3cit
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R3 = delta_G_H3cit - delta_G_NiHcit
    E_theta_3 = -delta_G_R3/(2*F)
    V3 = []
    if type(pH) == list:  # for list(vector), below is only for single value
        for pHval in pH:
            V3.append(E_theta_3 -(C2*2*pHval)/2 +(C2*np.log10(NiHcit))/2 -(C2*np.log10(H3cit))/2)
    else:
        V3 = E_theta_3 -(C2*2*pH)/2 +(C2*np.log10(NiHcit))/2 -(C2*np.log10(H3cit))/2
    return V3


def R4(pH, T, NiHcit, H2cit):  #NiHcit +H+ +2e = Ni +H2cit-
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R4 = delta_G_H2cit - delta_G_NiHcit
    E_theta_4 = -delta_G_R4/(2*F)
    V4 = []
    if type(pH) == list:  # for list(vector), below is only for single value
        for pHval in pH:
            V4.append(E_theta_4 -(C2*pHval)/2 +(C2*np.log10(NiHcit))/2 -(C2*np.log10(H2cit))/2)
    else:
        V4 = (E_theta_4 -(C2*pH)/2 +(C2*np.log10(NiHcit))/2 -(C2*np.log10(H2cit))/2)
    return V4


def R5(pH, T, Nicit, H2cit):  #Nicit +2H+ +2e = Ni +H2cit-
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R5 = delta_G_H2cit - delta_G_Nicit
    E_theta_5 = -delta_G_R5/(2*F)
    V5 = []
    if type(pH) == list:  # for list(vector), below is only for single value
        for pHval in pH:
            V5.append(E_theta_5 -(C2*2*pHval)/2 +(C2*np.log10(Nicit))/2 -(C2*np.log10(H2cit))/2)
    else:
        V5 = (E_theta_5 -(C2*2*pH)/2 +(C2*np.log10(Nicit))/2 -(C2*np.log10(H2cit))/2)
    return V5


def R6(pH, T, Nicit, Hcit):  #Nicit +H+ +2e = Ni +Hcit-
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R6 = delta_G_Hcit - delta_G_Nicit
    E_theta_6 = -delta_G_R6/(2*F)
    V6 = []
    if type(pH) == list:  # for list(vector), below is only for single value
        for pHval in pH:
            V6.append(E_theta_6 -(C2*pHval)/2 +(C2*np.log10(Nicit))/2 -(C2*np.log10(Hcit))/2)
    else:
        V6 = (E_theta_6 -(C2*pH)/2 +(C2*np.log10(Nicit))/2 -(C2*np.log10(Hcit))/2)
    return V6


def R7(T, Nicit, cit3):  #Nicit +2e = Ni +cit3-
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R7 = delta_G_cit3 - delta_G_Nicit
    E_theta_7 = -delta_G_R7/(2*F)
    V7 = E_theta_7 +(C2*np.log10(Nicit))/2 -(C2*np.log10(cit3))/2
    return V7


def R8(pH,T):  #Nio2 +2H+ +2e = Ni+ 2H2O
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_R8 = 2*delta_G_w - delta_G_nio2
    E_theta_8 = -delta_G_R8/(2*F)
    V8=[]
    if type(pH) == list:  # for list(vector), below is only for single value
        for pHval in pH:
            V8.append(E_theta_8 -(C2*2*pH) / 2)
    else:
        V8 = E_theta_8 -(C2*2*pH) / 2
    return V8


def T1(nip2, pH, T): #Ni(OH)3 +3H+ +e- =Ni(2+) +3H2O
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    delta_G_T1 = 3*delta_G_w + delta_G_nip2 - delta_G_nio3
    E_theta_T1 = -delta_G_T1/F
    if type(pH) == list:
        vt1 = []
        for pHval in pH:
            vt1.append(E_theta_T1 - C2*np.log10(nip2) - C2*3*pHval)
    else:
        vt1 = E_theta_T1 - C2*np.log10(nip2) - C2*3*pH
    return vt1


def T2(pH,NiH2cit,H3cit,T): #Ni(OH)3 +H3cit +2H +e- =NiH2cit+ +3H2O
    C1 = -np.log(10) * R * T
    C2 = -C1 / F
    delta_G_T2 = delta_G_NiH2cit + 3*delta_G_w - delta_G_H3cit - delta_G_nio3
    E_theta_T2 = -delta_G_T2/F
    if type(pH) == list:
        vt2 = []
        for pHval in pH:
            vt2.append(E_theta_T2 +C2*np.log10(H3cit)-C2*2*pHval-C2*np.log10(NiH2cit))
    else:
        vt2 = E_theta_T2 +C2*np.log10(H3cit)-C2*2*pH-C2*np.log10(NiH2cit)
    return vt2


def T3(pH,NiHcit,H3cit,T): #Ni(OH)3 +H3cit +H +e- =NiHcit+ +3H2O
    C1 = -np.log(10)*R*T
    C2 = -C1 / F
    delta_G_T3 = delta_G_NiHcit + 3*delta_G_w - delta_G_H3cit - delta_G_nio3
    E_theta_T3 = -delta_G_T3/F
    if type(pH) == list:
        vt3 = []
        for pHval in pH:
            vt3.append(E_theta_T3 +C2*np.log10(H3cit)-C2*pHval-C2*np.log10(NiHcit))
    else:
        vt3 = E_theta_T3 +C2*np.log10(H3cit)-C2*pH-C2*np.log10(NiHcit)
    return vt3


def T4(pH,NiHcit,H2cit,T): #Ni(OH)3 +H2cit +2H +e- =NiHcit +3H2O
    C1 = -np.log(10)*R*T
    C2 = -C1 / F
    delta_G_T4 = delta_G_NiHcit + 3*delta_G_w - delta_G_H2cit - delta_G_nio3
    E_theta_T4 = -delta_G_T4/F
    if type(pH) == list:
        vt4 = []
        for pHval in pH:
            vt4.append(E_theta_T4 +C2*np.log10(H2cit)-C2*2*pHval-C2*np.log10(NiHcit))
    else:
        vt4 = E_theta_T4 +C2*np.log10(H2cit)-C2*2*pH-C2*np.log10(NiHcit)
    return vt4


def T5(pH,Nicit,H2cit,T): #Ni(OH)3 +H2cit +H +e- =Nicit +3H2O
    C1 = -np.log(10)*R*T
    C2 = -C1 / F
    delta_G_T5 = delta_G_Nicit + 3*delta_G_w - delta_G_H2cit - delta_G_nio3
    E_theta_T5 = -delta_G_T5/F
    if type(pH) == list:
        vt5 = []
        for pHval in pH:
            vt5.append(E_theta_T5 +C2*np.log10(H2cit)-C2*pHval-C2*np.log10(Nicit))
    else:
        vt5 = E_theta_T5 +C2*np.log10(H2cit)-C2*pH-C2*np.log10(Nicit)
    return vt5


def T6(pH,Nicit,Hcit,T): #Ni(OH)3 +Hcit +2H +e- =Nicit +3H2O
    C1 = -np.log(10)*R*T
    C2 = -C1 / F
    delta_G_T6 = delta_G_Nicit + 3*delta_G_w - delta_G_Hcit - delta_G_nio3
    E_theta_T6 = -delta_G_T6/F
    if type(pH) == list:
        vt6 = []
        for pHval in pH:
            vt6.append(E_theta_T6 +C2*np.log10(Hcit)-C2*2*pHval-C2*np.log10(Nicit))
    else:
        vt6 = E_theta_T6 +C2*np.log10(Hcit)-C2*2*pH-C2*np.log10(Nicit)
    return vt6

def T7(pH,Nicit,cit3,T): #Ni(OH)3 +cit3 +3H +e- =Nicit +3H2O
    C1 = -np.log(10)*R*T
    C2 = -C1 / F
    delta_G_T7 = delta_G_Nicit + 3*delta_G_w - delta_G_cit3- delta_G_nio3
    E_theta_T7 = -delta_G_T7/F
    if type(pH) == list:
        vt7 = []
        for pHval in pH:
            vt7.append(E_theta_T7 +C2*np.log10(cit3)-C2*3*pHval-C2*np.log10(Nicit))
    else:
        vt7 = E_theta_T7 +C2*np.log10(cit3)-C2*3*pH-C2*np.log10(Nicit)
    return vt7

def T8(pH,T): #Ni(OH)3 +H +e- =Ni02 +H2O
    C1 = -np.log(10)*R*T
    C2 = -C1 / F
    delta_G_T8 = delta_G_nio2 + delta_G_w -delta_G_nio3
    E_theta_T8 = -delta_G_T8/F
    if type(pH) == list:
        vt8 = []
        for pHval in pH:
            vt8.append(E_theta_T8 -C2*pHval)
    else:
        vt8 = E_theta_T8 -C2*pH
    return vt8

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
