# coding: utf-8
'''
This function read_struct_get_ewald(df,df_name) reads structure file where df=train/test and output Ewald energies.
'''

class Atom:

    # Equality epsilon
    ee = 1e-6

    def __init__(self,
                 x=0.0,
                 y=0.0,
                 z=0.0,
                 t="",
                 c=0):

        self.x = x
        self.y = y
        self.z = z
        self.t = t
        self.c = c

def read_geometry_file(path_to_file):

    with open(path_to_file) as f:
        lines = f.readlines()

    vec_x = lines[3].split()
    vec_y = lines[4].split()
    vec_z = lines[5].split()

    vec_x = [float(vec_x[i]) for i in range(1, len(vec_x))]
    vec_y = [float(vec_y[i]) for i in range(1, len(vec_y))]
    vec_z = [float(vec_z[i]) for i in range(1, len(vec_z))]

    ga_charge = periodic_table.Element['Ga'].Z # 31.0
    al_charge = periodic_table.Element['Al'].Z # 13.0
    in_charge = periodic_table.Element['In'].Z # 49.0
    o_charge = periodic_table.Element['O'].Z # 8.0

    vectors = [vec_x, vec_y, vec_z]
    uc_atoms = []
    for i in range(6, len(lines)):
        sl = lines[i].split()
        x = float(sl[1])
        y = float(sl[2])
        z = float(sl[3])
        t = sl[4]

        if sl[4] == "Ga":
            c = ga_charge
        elif sl[4] == "Al":
            c = al_charge
        elif sl[4] == "In":
            c = in_charge
        elif sl[4] == "O":
            c = o_charge

        a = Atom(x, y, z, t, c)
        uc_atoms.append(a)

    return vectors, uc_atoms

def convert_uc_atoms_to_input_for_pymatgen(uc_atoms):

    n = len(uc_atoms)
    atom_coords = []
    atom_labels = []
    charge_list = []
    for i in range(n):
        x = uc_atoms[i].x
        y = uc_atoms[i].y
        z = uc_atoms[i].z
        t = uc_atoms[i].t
        c = uc_atoms[i].c

        vec = [x, y, z]

        atom_coords.append(vec)
        atom_labels.append(t)
        charge_list.append(c)
    site_properties = {"charge": charge_list}

    return atom_coords, atom_labels, site_properties

def read_struct_get_ewald(df,df_name):
    r_ene = [] 
    k_ene = []
    p_ene = []
    eta = []
    for i in range(len(df)):
        filename = str(str(df_name)+'/'+str(i+1)+'/geometry.xyz')
        vectors, uc_atoms = read_geometry_file(filename)
        atom_coords, atom_labels, site_properties = convert_uc_atoms_to_input_for_pymatgen(uc_atoms)
        lattice = Lattice.from_parameters(a=df.lattice_vector_1_ang[i], b=df.lattice_vector_2_ang[i], c=df.lattice_vector_3_ang[i], 
                                          alpha=df.lattice_angle_alpha_degree[i], beta=df.lattice_angle_beta_degree[i], gamma=df.lattice_angle_gamma_degree[i])
        structure = Structure(lattice, atom_labels, atom_coords, site_properties=site_properties)
        ewald_sum = ewald.EwaldSummation(structure)
        r_ene.append(ewald_sum.real_space_energy)
        k_ene.append(ewald_sum.reciprocal_space_energy)
        p_ene.append(ewald_sum.point_energy)
        eta.append(ewald_sum.eta)
    df['r_ene'] = pd.DataFrame(r_ene)
    df['k_ene'] = pd.DataFrame(k_ene)
    df['p_ene'] = pd.DataFrame(p_ene)
    df['eta'] = pd.DataFrame(eta)
    return df['r_ene'], df['k_ene'], df['p_ene'], df['eta']
