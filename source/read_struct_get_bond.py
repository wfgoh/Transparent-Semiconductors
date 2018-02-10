'''
These functions 
1. read_struct_get_bond(df,df_name) reads structure file where df=train/test and output min and mean bond length and inverse and inverse_square of bond length for metal oxides.
2. read_struct_get_mean(df,df_name) reads structure file where df=train/test and output mean bond length for all metal oxides.
'''

def calculate_rij(struct,rcut=3.0):
    rgao = np.empty(0)
    ralo = np.empty(0)
    rino = np.empty(0)
    invrgao = np.empty(0)
    invralo = np.empty(0)
    invrino = np.empty(0)
    inv2rgao = np.empty(0)
    inv2ralo = np.empty(0)
    inv2rino = np.empty(0)
    if not struct.loc[struct['atom'] == 'O'].empty:
        for _,O in struct.loc[struct['atom'] == 'O'].iterrows():
            if not struct.loc[struct['atom'] == 'Ga'].empty:
                for _,Ga in struct.loc[struct['atom'] == 'Ga'].iterrows():
                    rij = np.linalg.norm([Ga.x-O.x, Ga.y-O.y, Ga.z-O.z])
                    if rij < rcut:
                        rgao = np.append(rgao,rij)
                        invrgao = np.append(invrgao,1/rij)
                        inv2rgao = np.append(inv2rgao,1/np.square(rij))
            if not struct.loc[struct['atom'] == 'Al'].empty:
                for _,Al in struct.loc[struct['atom'] == 'Al'].iterrows():
                    rij = np.linalg.norm([Al.x-O.x, Al.y-O.y, Al.z-O.z])
                    if rij < rcut:
                        ralo = np.append(ralo,rij)
                        invralo = np.append(invralo,1/rij)
                        inv2ralo = np.append(inv2ralo,1/np.square(rij))

            if not struct.loc[struct['atom'] == 'In'].empty:
                for _,In in struct.loc[struct['atom'] == 'In'].iterrows():
                    rij = np.linalg.norm([In.x-O.x, In.y-O.y, In.z-O.z])
                    if rij < rcut:
                        rino = np.append(rino,rij)
                        invrino = np.append(invrino,1/rij)
                        inv2rino = np.append(inv2rino,1/np.square(rij))
        if not struct.loc[struct['atom'] == 'Ga'].empty:
            rgao_ave = rgao.mean()
            rgao_min = rgao.min()
            invrgao_sum = invrgao.sum()
            inv2rgao_sum = inv2rgao.sum()
        else: 
            rgao_ave = np.inf
            rgao_min = np.inf
            invrgao_sum = 0
            inv2rgao_sum = 0
        if not struct.loc[struct['atom'] == 'Al'].empty:
            ralo_ave = ralo.mean()
            ralo_min = ralo.min()
            invralo_sum = invralo.sum()
            inv2ralo_sum = inv2ralo.sum()
        else: 
            ralo_ave = np.inf
            ralo_min = np.inf
            invralo_sum = 0
            inv2ralo_sum = 0

        if not struct.loc[struct['atom'] == 'In'].empty:
            rino_ave = rino.mean()
            rino_min = rino.min()
            invrino_sum = invrino.sum()
            inv2rino_sum = inv2rino.sum()
        else: 
            rino_ave = np.inf
            rino_min = np.inf
            invrino_sum = 0
            inv2rino_sum = 0
    return rgao_ave, rgao_min, ralo_ave, ralo_min, rino_ave, rino_min, invrgao_sum, inv2rgao_sum, invralo_sum, inv2ralo_sum, invrino_sum, inv2rino_sum

def read_struct_get_bond(df,df_name):
    rgaoave = []
    rgaomin = []
    raloave = []
    ralomin = []
    rinoave = []
    rinomin = []
    invrgaosum = []
    inv2rgaosum = []
    invralosum = []
    inv2ralosum = []
    invrinosum = []
    inv2rinosum = []
    for i in range(len(df)):
        filename = str(str(df_name)+'/'+str(i+1)+'/geometry.xyz')
        struct = pd.read_csv(filename, sep=' ',skiprows=6, names = ['-','x','y','z','atom'])
        rgao_ave, rgao_min, ralo_ave, ralo_min, rino_ave, rino_min, invrgao_sum, inv2rgao_sum, invralo_sum, inv2ralo_sum, invrino_sum, inv2rino_sum = calculate_rij(struct)
        rgaoave.append(rgao_ave)
        rgaomin.append(rgao_min)
        raloave.append(ralo_ave)
        ralomin.append(ralo_min)
        rinoave.append(rino_ave)
        rinomin.append(rino_min)
        invrgaosum.append(invrgao_sum)
        inv2rgaosum.append(inv2rgao_sum)
        invralosum.append(invralo_sum)
        inv2ralosum.append(inv2ralo_sum)
        invrinosum.append(invrino_sum)
        inv2rinosum.append(inv2rino_sum)
    df['r_gao_ave'] = pd.DataFrame(rgaoave)
    df['r_gao_min'] = pd.DataFrame(rgaomin)
    df['r_alo_ave'] = pd.DataFrame(raloave)
    df['r_alo_min'] = pd.DataFrame(ralomin)
    df['r_ino_ave'] = pd.DataFrame(rinoave)
    df['r_ino_min'] = pd.DataFrame(rinomin)
    df['inv_r_gao_sum'] = pd.DataFrame(invrgaosum)
    df['inv2_r_gao_sum'] = pd.DataFrame(inv2rgaosum)
    df['inv_r_alo_sum'] = pd.DataFrame(invralosum)
    df['inv2_r_alo_sum'] = pd.DataFrame(inv2ralosum)
    df['inv_r_ino_sum'] = pd.DataFrame(invrinosum)
    df['inv2_r_ino_sum'] = pd.DataFrame(inv2rinosum)
    return df['r_gao_ave'], df['r_gao_min'], df['r_alo_ave'], df['r_alo_min'], df['r_ino_ave'], df['r_ino_min'], df['inv_r_gao_sum'], df['inv2_r_gao_sum'], df['inv_r_alo_sum'], df['inv2_r_alo_sum'], df['inv_r_ino_sum'], df['inv2_r_ino_sum']

def calculate_rij_mean(struct,rcut=3.0):
    r = np.empty(0)
    if not struct.loc[struct['atom'] == 'O'].empty:
        for _,O in struct.loc[struct['atom'] == 'O'].iterrows():
            if not struct.loc[struct['atom'] == 'Ga'].empty:
                for _,Ga in struct.loc[struct['atom'] == 'Ga'].iterrows():
                    rij = np.linalg.norm([Ga.x-O.x, Ga.y-O.y, Ga.z-O.z])
                    if rij < rcut:
                        r = np.append(r,rij)
            if not struct.loc[struct['atom'] == 'Al'].empty:
                for _,Al in struct.loc[struct['atom'] == 'Al'].iterrows():
                    rij = np.linalg.norm([Al.x-O.x, Al.y-O.y, Al.z-O.z])
                    if rij < rcut:
                        r = np.append(r,rij)
            if not struct.loc[struct['atom'] == 'In'].empty:
                for _,In in struct.loc[struct['atom'] == 'In'].iterrows():
                    rij = np.linalg.norm([In.x-O.x, In.y-O.y, In.z-O.z])
                    if rij < rcut:
                        r = np.append(r,rij)
    r_mean = r.mean()
    return r_mean

def read_struct_get_mean(df,df_name):
    r_mean = []
    for i in range(len(df)):
        filename = str(str(df_name)+'/'+str(i+1)+'/geometry.xyz')
        struct = pd.read_csv(filename, sep=' ',skiprows=6, names = ['-','x','y','z','atom'])
        rmean = calculate_rij_mean(struct)
        r_mean.append(rmean)
    df['r_mean'] = pd.DataFrame(r_mean)
    return df['r_mean']
