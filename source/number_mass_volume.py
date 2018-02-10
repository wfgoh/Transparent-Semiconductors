'''
These functions:
1. get_vol(a,b,c,alpha,beta,gamma) convert the lattice parameters to volume.
2. get_n_atom(df,df_name) where df=train/test output the total number of each atom for each train/test set.
3. get_mass(df) where df=train/test output the total mass of each atoms for each train/test set.
'''

def get_vol(a,b,c,alpha,beta,gamma):
    
    vol = a*b*c*np.sqrt(1+2*np.cos(np.deg2rad(alpha))*np.cos(np.deg2rad(beta))*np.cos(np.deg2rad(gamma))-np.cos(np.deg2rad(alpha))**2-np.cos(np.deg2rad(beta))**2-np.cos(np.deg2rad(gamma))**2)
    
    return vol
    
def get_n_atom(df,df_name):
    
    N_Ga = []
    N_Al = []
    N_In = []
    N_O = []
    for i in df.id:
        n_Ga = 0
        n_Al = 0
        n_In = 0
        n_O = 0
        with open(str(df_name) + '/' + str(i) + '/geometry.xyz','r') as finp: 
            for line in finp:
                if line.split()[0] == 'atom':
                    if line.split()[4] == 'Ga':
                        n_Ga += 1
                    elif line.split()[4] == 'Al':
                        n_Al += 1
                    elif line.split()[4] == 'In':
                        n_In += 1
                    elif line.split()[4] == 'O':
                        n_O += 1
        N_Ga.append(n_Ga)
        N_Al.append(n_Al)
        N_In.append(n_In)
        N_O.append(n_O)
    df['n_ga'] = pd.DataFrame(N_Ga)
    df['n_al'] = pd.DataFrame(N_Al)
    df['n_in'] = pd.DataFrame(N_In)
    df['n_o'] = pd.DataFrame(N_O)
    
    return df['n_al'], df['n_ga'], df['n_in'], df['n_o']

def get_mass(df):
    
    mass_ga = periodic_table.Element['Ga'].atomic_mass
    mass_al = periodic_table.Element['Al'].atomic_mass
    mass_in = periodic_table.Element['In'].atomic_mass
    mass_o = periodic_table.Element['O'].atomic_mass
    tot_mass_ga = df['n_ga']*mass_ga
    tot_mass_al = df['n_al']*mass_al
    tot_mass_in = df['n_in']*mass_in
    tot_mass_o = df['n_o']*mass_o
    
    return tot_mass_al, tot_mass_ga, tot_mass_in, tot_mass_o

