import numpy as np

from ase import Atoms
from scipy.constants import Boltzmann, Planck


def kcal_per_mol_to_J(value):
    return value*(4184/(6.02214076*(10**23)))

def J_to_kcal_per_mol(value):
    return value*((6.02214076*(10**23))*0.000239006)

def J_to_eV(value):
    return value*(6.241509*(10**18))

def find_intersection_point(angle,y):
    slope = np.tan(np.deg2rad(angle))
    x_inter = y / slope
    inter_point = np.array([x_inter,y])
    return inter_point

def find_angle_rad_2d(v1, v2):
    norm_v1 = np.linalg.norm(v1) 
    norm_v2 = np.linalg.norm(v2) 
    angle_rad = np.arccos(np.dot(v1, v2) / (norm_v1*norm_v2))   
    return angle_rad

def find_angle_deg_2d(v1, v2):
    norm_v1 = np.linalg.norm(v1) 
    norm_v2 = np.linalg.norm(v2) 
    angle_deg = np.arccos(np.dot(v1, v2) / (norm_v1*norm_v2)) * (180/np.pi) 
    return angle_deg

def evaluate_equal(coord1, coord2):
    return round(coord1[0], 2) == round(coord2[0], 2) and round(coord1[1], 2) == round(coord2[1], 2) and round(coord1[2], 2) == round(coord2[2], 2)

def moment_of_inertia(molecule):
    """
    molecule : ase.Atoms object
    """
    mol = molecule.copy()
    center_of_mass = mol.get_center_of_mass()
    mol_pos = mol.get_positions().copy()
    vectors = [mol_pos[i]-center_of_mass for i in range(len(mol_pos))]
    molecule = Atoms(mol.symbols, positions = vectors)
    pos = molecule.get_positions().copy()
    amu_to_kg = 1.66054*(10**(-27))
    ang_to_m = (10**(-10))
    Ixx = sum([molecule[i].mass * (pos[i][1]**2 + pos[i][2]**2) for i in range(len(molecule))])
    Iyy = sum([molecule[i].mass * (pos[i][0]**2 + pos[i][2]**2) for i in range(len(molecule))])
    Izz = sum([molecule[i].mass * (pos[i][0]**2 + pos[i][1]**2) for i in range(len(molecule))])
    Ixy = -sum([molecule[i].mass * (pos[i][0] * pos[i][1]) for i in range(len(molecule))])
    Ixz = -sum([molecule[i].mass * (pos[i][0] * pos[i][2]) for i in range(len(molecule))])
    Iyz = -sum([molecule[i].mass * (pos[i][1] * pos[i][2]) for i in range(len(molecule))])

    I = np.array([[Ixx, Ixy, Ixz],
                  [Ixy, Iyy, Iyz],
                  [Ixz, Iyz, Izz]])
    adj = 10
    if round(Ixy, adj) != 0 or round(Ixz, adj) != 0 or round(Iyz, adj) != 0:
        eigenvalues, eigenvectors = np.linalg.eig(I)
        D = np.diag(eigenvalues)
        P = eigenvectors
        I_diag = np.dot(np.dot(np.linalg.inv(P), I), P)
        I = I_diag

    print(f'\n=== Moment of Inertia tensor (amu*angstrom) ===\n{I}')
    I_list = [I[0][0], I[1][1], I[2][2]]
    I_C = max(I_list)
    I_C_index = I_list.index(I_C)
    I_C = I_C * amu_to_kg * (ang_to_m**(2))
    del I_list[I_C_index]
    I_B = max(I_list)
    I_B_index = I_list.index(I_B)
    I_B = I_B * amu_to_kg * (ang_to_m**(2))
    del I_list[I_B_index]
    I_A = I_list[0] * amu_to_kg * (ang_to_m**(2))
    print(f'\n==> I_A : {I_A}\n==> I_B : {I_B}\n==> I_C : {I_C}   (kg*m)')
    print('===============================================\n')
    return I_A, I_B, I_C

def calculate_molecular_mass(molecule):
    """
    molecule : ase.Atoms object
    """
    return sum([i.mass*1.66054*(10**(-27)) for i in molecule])

def zero_point_enery(v_list):
    return (1/2)*(sum([Planck*i for i in v_list]))

def electronic_contributions(w0):
    """
    Calculate contributions from electronic motion.
    This calculation is for cases 
    where the excited state energy levels are much higher than kB*T
    """
    kB = Boltzmann
    q_e = w0
    S_e = kB * np.log(q_e)
    H_e = 0

    print('\n###### Electronic Contribution ###############\n')
    print(f"q_e = {q_e:5.5f}")
    print('----------------------------------------------')
    print(f"S_e = {S_e:5.5f} J/K")
    print(f"    = {J_to_kcal_per_mol(S_e):5.5f} kcal/K*mol")
    print(f"    = {J_to_eV(S_e):5.5f} eV/K")
    print('----------------------------------------------')
    print(f"H_e = {H_e:5.5f} J")
    print(f"    = {J_to_kcal_per_mol(H_e):5.5f} kcal/mol")
    print(f"    = {J_to_eV(H_e):5.5f} eV")

    return q_e, S_e, H_e

def vibrational_contribution(v_list, T):
    def q_vK(vK, T):
        """
        Calculate partition function for each vibrational mode
        vK is K_th vibrational frequency
        T is temperature(K)
        """
        return 1/(1-(np.exp(-theta_vK(vK)/T)))

    def theta_vK(vK):
        """
        Calculate theta(vK) for each vibrational mode.
        theta(vK) is characteristic vibrational temperature.
        vK is K_th vibrational frequency.
        """
        return (Planck*vK)/kB

    kB = Boltzmann
    first_term_list = [np.log(q_vK(vK, T)) for vK in v_list]
    last_term_list = [(theta_vK(vK)/T)*(1/(np.exp(theta_vK(vK)/T)-1)) for vK in v_list]
    enthalpy_list = [theta_vK(vK)*(1/(np.exp(theta_vK(vK)/T)-1)) for vK in v_list]
    q_v_list = [q_vK(vK, T) for vK in v_list]

    q_v = np.prod(q_v_list)
    S_v = kB * (sum(first_term_list) + sum(last_term_list))
    H_v = kB * sum(enthalpy_list)

    print('\n###### Vibrational Contribution ##############\n')
    print(f"q_v = {q_v:5.5f}")
    print('----------------------------------------------')
    print(f"S_v = {S_v:5.5f}     J/K")
    print(f"    = {J_to_kcal_per_mol(S_v):5.5e} kcal/K*mol")
    print(f"    = {J_to_eV(S_v):5.5f}     eV/K")
    print('----------------------------------------------')
    print(f"H_v = {H_v:5.5f}     J")
    print(f"    = {J_to_kcal_per_mol(H_v):5.5f}     kcal/mol")
    print(f"    = {J_to_eV(H_v):5.5f}     eV")
    
    return q_v, S_v, H_v

def rotational_contirubution(molecule, sigma, T, linear=False):
    kB = Boltzmann
    def theta_r(I):
        """ 
        Calculate theta(r) for moment of inertia.
        theta(r) is characteristic rotational temperature.
        I is moment of inertia.
        """
        return (Planck**2)/(8*(np.pi**2)*I*kB)
    
    print('\n###### Rotational Contribution ###############\n')    
    I_A, I_B, I_C = moment_of_inertia(molecule)
    #print(round(I_A, 50), round(I_B, 50))
    if linear or round(I_A, 50) == 0:
        print(f'  ** This molecule is considered "linear" **\n')
        I = I_B
        q_r = T/(sigma*theta_r(I))
        S_r = kB * ((np.log(q_r)) + 1)
        H_r = kB * T
    else:
        print(f'  ** This molecule is considered "non-linear" **')
        print(f'  ** If not, please set it as "linear = True" **\n')
        theta_A = theta_r(I_A)
        theta_B = theta_r(I_B)
        theta_C = theta_r(I_C)
        q_r = (np.sqrt(np.pi)/sigma)*(np.sqrt((T**3)/(theta_A*theta_B*theta_C)))
        S_r = kB * ((np.log(q_r)) + 3/2)
        H_r = (3/2) * kB * T

    print(f"q_r = {q_r:5.5f}")
    print('----------------------------------------------')
    print(f"S_r = {S_r:5.5f}     J/K")
    print(f"    = {J_to_kcal_per_mol(S_r):5.5e} kcal/K*mol")
    print(f"    = {J_to_eV(S_r):5.5f}     eV/K")
    print('----------------------------------------------')
    print(f"H_r = {H_r:5.5f}     J")
    print(f"    = {J_to_kcal_per_mol(H_r):5.5f}     kcal/mol")
    print(f"    = {J_to_eV(H_r):5.5f}     eV")

    return q_r, S_r, H_r

def translational_contribution(molecule, T, P=None, M=None):
    if M is None:
        M = calculate_molecular_mass(molecule)
    if P is None:
        P = 1*101325 # 1 atm to Pa

    kB = Boltzmann
    q_t = (((2*np.pi*M*kB*T)/(Planck**2))**(3/2)) * (kB*T/P)
    S_t = kB * (np.log(q_t) + (5/2))
    H_t = (5/2)*kB*T

    print('\n###### Translational Contribution ############\n')
    print(f"q_t = {q_t:5.5e}")
    print('----------------------------------------------')
    print(f"S_t = {S_t:5.5f}     J/K")
    print(f"    = {J_to_kcal_per_mol(S_t):5.5e} kcal/K*mol")
    print(f"    = {J_to_eV(S_t):5.5f}     eV/K")
    print('----------------------------------------------')
    print(f"H_t = {H_t:5.5f}     J")
    print(f"    = {J_to_kcal_per_mol(H_t):5.5f}     kcal/mol")
    print(f"    = {J_to_eV(H_t):5.5f}     eV")

    return q_t, S_t, H_t

def total_contribution(v_list, w0, T,
                       molecule, sigma,
                       P=None, M=None,
                       linear=False):
    """
    - molecule : ase.Atoms object
    - sigma    : Rotational symmetry number
    - w0       : Spin multiplicity
    - v_list   : Vibrational frequencies list 
    """
    zpe = zero_point_enery(v_list)
    q_e, S_e, H_e = electronic_contributions(w0)
    q_v, S_v, H_v = vibrational_contribution(v_list, T)
    q_r, S_r, H_r = rotational_contirubution(molecule, sigma, T, linear)
    q_t, S_t, H_t = translational_contribution(molecule, T, P, M)

    S = S_e + S_v + S_r + S_t   #+ Boltzman
    Cp = H_e + H_v + H_r + H_t  #+ Boltzmann*T
    Cp_TS_ZPE = Cp - T*S + zpe

    print('\n###### Total Contribution ####################\n')
    print(f"ZPE        = {zpe:9.5f}     J")
    print(f"           = {J_to_kcal_per_mol(zpe):9.5f}     kcal/mol")
    print(f"           = {J_to_eV(zpe):9.5f}     eV")
    print('----------------------------------------------')
    print(f"-TS        = {-T*S:9.5f}     J")
    print(f"           = {-J_to_kcal_per_mol(T*S):9.5f}     kcal/mol")
    print(f"           = {-J_to_eV(T*S):9.5f}     eV")
    print('----------------------------------------------')
    print(f"Cp         = {Cp:9.5f}     J")
    print(f"           = {J_to_kcal_per_mol(Cp):9.5f}     kcal/mol")
    print(f"           = {J_to_eV(Cp):9.5f}     eV")
    print('----------------------------------------------')
    print(f"Cp-TS+ZPE  = {Cp_TS_ZPE:9.5f}     J")
    print(f"           = {J_to_kcal_per_mol(Cp_TS_ZPE):9.5f}     kcal/mol")
    print(f"           = {J_to_eV(Cp_TS_ZPE):9.5f}     eV\n")

    return zpe, S, Cp, Cp_TS_ZPE
                           

