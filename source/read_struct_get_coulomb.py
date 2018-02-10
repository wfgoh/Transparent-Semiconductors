'''
This function read_struct_get_r_coul(df,df_name) reads structure file where df=train/test and output r_ij and coulomb_ij.
'''

def calc_r_and_coul(struct,rcut1=4.0,rcut2=2.5):
    charge = {'Ga': 3, 'Al':3, 'In':3, 'O':-2}
    
    coul = np.empty(0)
    r = np.empty(0)
    cgaga = np.empty(0)
    rgaga = np.empty(0)
    calal = np.empty(0)
    ralal = np.empty(0)
    cinin = np.empty(0)
    rinin = np.empty(0)
    coo = np.empty(0)
    roo = np.empty(0)
    cgaal = np.empty(0)
    rgaal = np.empty(0)
    cgain = np.empty(0)
    rgain = np.empty(0)
    cgao = np.empty(0)
    rgao = np.empty(0)
    calin = np.empty(0)
    ralin = np.empty(0)
    calo = np.empty(0)
    ralo = np.empty(0)    
    cino = np.empty(0)
    rino = np.empty(0)   
    
    for i,atomi in struct.iterrows():
        for j,atomj in struct.iterrows():
            if j>i:
                rij = np.linalg.norm([atomi.x-atomj.x, atomi.y-atomj.y, atomi.z-atomj.z])
                coulij = (charge[atomi.atom]*charge[atomj.atom])/(np.square(rij))
                if rij < rcut1:
                    if rij < rcut2:
                        coul = np.append(coul,coulij)
                        r = np.append(r,rij)

                    # seperate coul and r according to atomi,j
                    if atomi.atom == 'Ga' and atomj.atom == 'Ga' :
                        cgaga = np.append(cgaga,coulij)
                        rgaga = np.append(rgaga,rij)
                    if atomi.atom == 'Al' and atomj.atom == 'Al' :
                        calal = np.append(calal,coulij)
                        ralal = np.append(ralal,rij)                        
                    if atomi.atom == 'In' and atomj.atom == 'In' :
                        cinin = np.append(cinin,coulij)
                        rinin = np.append(rinin,rij)                        
                    if atomi.atom == 'O' and atomj.atom == 'O' :
                        coo = np.append(coo,coulij)
                        roo = np.append(roo,rij)                        
                    if (atomi.atom == 'Ga' and atomj.atom == 'Al') or (atomj.atom == 'Ga' and atomi.atom == 'Al') :
                        cgaal = np.append(cgaal,coulij)
                        rgaal = np.append(rgaal,rij)
                    if (atomi.atom == 'Ga' and atomj.atom == 'In') or (atomj.atom == 'Ga' and atomi.atom == 'In') :
                        cgain = np.append(cgain,coulij)
                        rgain = np.append(rgain,rij)                        
                    if (atomi.atom == 'Ga' and atomj.atom == 'O') or (atomj.atom == 'Ga' and atomi.atom == 'O') :
                        if rij < rcut2:
                            cgao = np.append(cgao,coulij)
                            rgao = np.append(rgao,rij)                            
                    if (atomi.atom == 'Al' and atomj.atom == 'In') or (atomj.atom == 'Al' and atomi.atom == 'In') :
                        calin = np.append(calin,coulij)
                        ralin = np.append(ralin,rij)                            
                    if (atomi.atom == 'Al' and atomj.atom == 'O') or (atomj.atom == 'Al' and atomi.atom == 'O') :
                        if rij < rcut2:
                            calo = np.append(calo,coulij)
                            ralo = np.append(ralo,rij)                            
                    if (atomi.atom == 'In' and atomj.atom == 'O') or (atomj.atom == 'In' and atomi.atom == 'O') :
                        if rij < rcut2:
                            cino = np.append(cino,coulij)
                            rino = np.append(rino,rij)                            
                                                 
    coul_sum = coul.sum()
    r_mean = r.mean()
    cgaga_sum = cgaga.sum()
    rgaga_mean = rgaga.mean() if rgaga.size !=0 else np.inf
    calal_sum = calal.sum()
    ralal_mean = ralal.mean() if ralal.size !=0 else np.inf
    cinin_sum = cinin.sum()
    rinin_mean = rinin.mean() if rinin.size !=0 else np.inf
    coo_sum = coo.sum()
    roo_mean = roo.mean() if roo.size !=0 else np.inf
    cgaal_sum = cgaal.sum()
    rgaal_mean = rgaal.mean() if rgaal.size !=0 else np.inf
    cgain_sum = cgain.sum()
    rgain_mean = rgain.mean() if rgain.size !=0 else np.inf
    cgao_sum = cgao.sum()
    rgao_mean = rgao.mean() if rgao.size !=0 else np.inf
    calin_sum = calin.sum()
    ralin_mean = ralin.mean() if ralin.size !=0 else np.inf
    calo_sum = calo.sum()
    ralo_mean = ralo.mean() if ralo.size !=0 else np.inf
    cino_sum = cino.sum()
    rino_mean = rino.mean() if rino.size !=0 else np.inf
    
    return coul_sum, r_mean, cgaga_sum, rgaga_mean, calal_sum, ralal_mean, cinin_sum, rinin_mean, coo_sum, roo_mean, cgaal_sum, rgaal_mean, cgain_sum, rgain_mean, cgao_sum, rgao_mean, calin_sum, ralin_mean, calo_sum, ralo_mean, cino_sum, rino_mean
    
def read_struct_get_r_coul(df,df_name):
    coul = []
    r = []
    cgaga = []
    rgaga = []
    calal = []
    ralal = []
    cinin = []
    rinin = []
    coo = []
    roo = []
    cgaal = []
    rgaal = []
    cgain = []
    rgain = []
    cgao = []
    rgao = []
    calin = []
    ralin = []
    calo = []
    ralo = []   
    cino = []
    rino = [] 
    for i in range(len(df)):
        filename = str(str(df_name)+'/'+str(i+1)+'/geometry.xyz')
        struct = pd.read_csv(filename, sep=' ',skiprows=6, names = ['-','x','y','z','atom'])
        coul_sum, r_mean, cgaga_sum, rgaga_mean, calal_sum, ralal_mean, cinin_sum, rinin_mean, coo_sum, roo_mean, cgaal_sum, rgaal_mean, cgain_sum, rgain_mean, cgao_sum, rgao_mean, calin_sum, ralin_mean, calo_sum, ralo_mean, cino_sum, rino_mean = calc_r_and_coul(struct)

        coul.append(coul_sum)
        r.append(r_mean)
        cgaga.append(cgaga_sum)
        rgaga.append(rgaga_mean)
        calal.append(calal_sum)
        ralal.append(ralal_mean)        
        cinin.append(cinin_sum)
        rinin.append(rinin_mean)        
        coo.append(coo_sum)
        roo.append(roo_mean)        
        cgaal.append(cgaal_sum)
        rgaal.append(rgaal_mean)
        cgain.append(cgain_sum)
        rgain.append(rgain_mean)        
        cgao.append(cgao_sum)
        rgao.append(rgao_mean)        
        calin.append(calin_sum)
        ralin.append(ralin_mean)            
        calo.append(calo_sum)
        ralo.append(ralo_mean)
        cino.append(cino_sum)
        rino.append(rino_mean)        

    df['coul'] = pd.DataFrame(coul)
    df['rmean'] = pd.DataFrame(r)
    df['cgaga'] = pd.DataFrame(cgaga)
    df['rgaga'] = pd.DataFrame(rgaga)
    df['calal'] = pd.DataFrame(calal)
    df['ralal'] = pd.DataFrame(ralal)
    df['cinin'] = pd.DataFrame(cinin)
    df['rinin'] = pd.DataFrame(rinin)    
    df['coo'] = pd.DataFrame(coo)
    df['roo'] = pd.DataFrame(roo)    
    df['cgaal'] = pd.DataFrame(cgaal)
    df['rgaal'] = pd.DataFrame(rgaal)
    df['cgain'] = pd.DataFrame(cgain)
    df['rgain'] = pd.DataFrame(rgain)
    df['cgao'] = pd.DataFrame(cgao)
    df['rgao'] = pd.DataFrame(rgao)    
    df['calin'] = pd.DataFrame(calin)
    df['ralin'] = pd.DataFrame(ralin)        
    df['calo'] = pd.DataFrame(calo)
    df['ralo'] = pd.DataFrame(ralo)
    df['cino'] = pd.DataFrame(cino)
    df['rino'] = pd.DataFrame(rino)
 
    return df['coul'], df['rmean'], df['cgaga'], df['rgaga'], df['calal'], df['ralal'], df['cinin'], df['rinin'], df['coo'], df['roo'], df['cgaal'], df['rgaal'], df['cgain'], df['rgain'], df['cgao'], df['rgao'], df['calin'], df['ralin'], df['calo'], df['ralo'], df['cino'], df['rino']
