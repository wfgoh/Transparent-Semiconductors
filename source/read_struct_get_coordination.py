'''
This function read_struct_get_coord(df,df_name) reads structure file where df=train/test and output coordination number for metal-oxides.
'''

def count_coord(struct,rcut=2.5):
    ngao_tot = np.empty(0)
    nalo_tot = np.empty(0)
    nino_tot = np.empty(0)
    n_tot = np.empty(0)
    if not struct.loc[struct['atom'] == 'O'].empty:
        for _,O in struct.loc[struct['atom'] == 'O'].iterrows():
            ngao = 0
            nalo = 0
            nino = 0
            ntot = 0
            if not struct.loc[struct['atom'] == 'Ga'].empty:
                for _,Ga in struct.loc[struct['atom'] == 'Ga'].iterrows():
                    rij = np.linalg.norm([Ga.x-O.x, Ga.y-O.y, Ga.z-O.z])
                    if rij < rcut:
                        ngao+=1
                        ntot+=1
            ngao_tot = np.append(ngao_tot,ngao)
            if not struct.loc[struct['atom'] == 'Al'].empty:
                for _,Al in struct.loc[struct['atom'] == 'Al'].iterrows():
                    rij = np.linalg.norm([Al.x-O.x, Al.y-O.y, Al.z-O.z])
                    if rij < rcut:
                        nalo+=1
                        ntot+=1
            nalo_tot = np.append(nalo_tot,nalo)
            if not struct.loc[struct['atom'] == 'In'].empty:
                for _,In in struct.loc[struct['atom'] == 'In'].iterrows():
                    rij = np.linalg.norm([In.x-O.x, In.y-O.y, In.z-O.z])
                    if rij < rcut:
                        nino+=1
                        ntot+=1
            nino_tot = np.append(nino_tot,nino)
            n_tot = np.append(n_tot,ntot)
    ngao_mean = ngao_tot.mean()
    nalo_mean = nalo_tot.mean()
    nino_mean = nino_tot.mean()
    n_mean = n_tot.mean()
    return ngao_mean, nalo_mean, nino_mean, n_mean
    
def read_struct_get_coord(df,df_name):
    ngao_mean = []
    nalo_mean = []
    nino_mean = []
    n_mean = []
    for i in range(len(df)):
        filename = str(str(df_name)+'/'+str(i+1)+'/geometry.xyz')
        struct = pd.read_csv(filename, sep=' ',skiprows=6, names = ['-','x','y','z','atom'])
        ngaomean, nalomean, ninomean, nmean = count_coord(struct)
        ngao_mean.append(ngaomean)
        nalo_mean.append(nalomean)
        nino_mean.append(ninomean)
        n_mean.append(nmean)
    df['ngao_mean'] = pd.DataFrame(ngao_mean)
    df['nalo_mean'] = pd.DataFrame(nalo_mean)
    df['nino_mean'] = pd.DataFrame(nino_mean)
    df['n_mean'] = pd.DataFrame(n_mean)
    return df['ngao_mean'], df['nalo_mean'], df['nino_mean'], df['n_mean']
