import pandas as pd
import numpy as np
import read_struct_get_coordination, read_struct_get_r_coulomb, read_struct_get_ewald, read_struct_get_bond, number_mass_volume

train = pd.read_csv('train.csv')
test = pd.read_csv('test.csv')

train['vol'] = get_vol(train.lattice_vector_1_ang,train.lattice_vector_2_ang,train.lattice_vector_3_ang,train.lattice_angle_alpha_degree,train.lattice_angle_beta_degree,train.lattice_angle_gamma_degree)
test['vol'] = get_vol(test.lattice_vector_1_ang,test.lattice_vector_2_ang,test.lattice_vector_3_ang,test.lattice_angle_alpha_degree,test.lattice_angle_beta_degree,test.lattice_angle_gamma_degree)

train.n_al, train.n_ga, train.n_in, train.n_o = get_n_atom(train,'train')
test.n_al, test.n_ga, test.n_in, test.n_o = get_n_atom(test,'test')

train['mass_al'], train['mass_ga'], train['mass_in'], train['mass_o'] = get_mass(train)
test['mass_al'], test['mass_ga'], test['mass_in'], test['mass_o'] = get_mass(test)

train['density'] = (train['mass_al']+train['mass_ga']+train['mass_in']+train['mass_o'])/train['vol']
test['density'] = (test['mass_al']+test['mass_ga']+test['mass_in']+test['mass_o'])/test['vol']

train['mass_tot'] = train['mass_al']+train['mass_ga']+train['mass_in']+train['mass_o']
test['mass_tot'] = test['mass_al']+test['mass_ga']+test['mass_in']+test['mass_o']

train['z_al'] = train.n_al*periodic_table.Element['Al'].Z
train['z_ga'] = train.n_ga*periodic_table.Element['Ga'].Z
train['z_in'] = train.n_in*periodic_table.Element['In'].Z
train['z_o'] = train.n_o*periodic_table.Element['O'].Z
train['z_tot'] = train.z_al+train.z_ga+train.z_in+train.z_o
test['z_al'] = periodic_table.Element['Al'].Z
test['z_ga'] = periodic_table.Element['Ga'].Z
test['z_in'] = periodic_table.Element['In'].Z
test['z_o'] = periodic_table.Element['O'].Z
test['z_tot'] = test.z_al+test.z_ga+test.z_in+test.z_o

train.r_ene, train.k_ene, train.p_ene, train.eta = read_struct_get_ewald(train,'train')
test.r_ene, test.k_ene, test.p_ene, test.eta = read_struct_get_ewald(test,'test')

train.r_gao_ave, train.r_gao_min, train.r_alo_ave, train.r_alo_min, train.r_ino_ave, train.r_ino_min, train.inv_r_gao_sum, train.inv2_r_gao_sum, train.inv_r_alo_sum, train.inv2_r_alo_sum, train.inv_r_ino_sum, train.inv2_r_ino_sum = read_struct_get_bond(train,'train')
test.r_gao_ave, test.r_gao_min, test.r_alo_ave, test.r_alo_min, test.r_ino_ave, test.r_ino_min, test.inv_r_gao_sum, test.inv2_r_gao_sum, test.inv_r_alo_sum, test.inv2_r_alo_sum, test.inv_r_ino_sum, test.inv2_r_ino_sum = read_struct_get_bond(test,'test')

train['inv2_r_sum'] = train.inv2_r_alo_sum + train.inv2_r_gao_sum + train.inv2_r_ino_sum
test['inv2_r_sum'] = test.inv2_r_alo_sum + test.inv2_r_gao_sum + test.inv2_r_ino_sum

train['inv_r_sum'] = train.inv_r_alo_sum + train.inv_r_gao_sum + train.inv_r_ino_sum
test['inv_r_sum'] = test.inv_r_alo_sum + test.inv_r_gao_sum + test.inv_r_ino_sum

train['r_min'] = train.loc[:, ['r_gao_min', 'r_alo_min', 'r_ino_min']].min(axis=1)
test['r_min'] = test.loc[:, ['r_gao_min', 'r_alo_min', 'r_ino_min']].min(axis=1)

train['r_mean'] = read_struct_get_mean(train,'train')
test['r_mean'] = read_struct_get_mean(test,'test')

train['ngao_mean'], train['nalo_mean'], train['nino_mean'], train['n_mean'] = read_struct_get_coord(train,'train')
test['ngao_mean'], test['nalo_mean'], test['nino_mean'], test['n_mean'] = read_struct_get_coord(test,'test')

train.coul, train.rmean, train.cgaga, train.rgaga, train.calal, train.ralal, train.cinin, train.rinin, train.coo, train.roo, train.cgaal, train.rgaal, train.cgain, train.rgain, train.cgao, train.rgao, train.calin, train.ralin, train.calo, train.ralo, train.cino, train.rino = read_struct_get_r_coul(train,'train')
test.coul, test.rmean, test.cgaga, test.rgaga, test.calal, test.ralal, test.cinin, test.rinin, test.coo, test.roo, test.cgaal, test.rgaal, test.cgain, test.rgain, test.cgao, test.rgao, test.calin, test.ralin, test.calo, test.ralo, test.cino, test.rino = read_struct_get_r_coul(test,'test')

train.to_csv('train5.csv',inplace=False)
test.to_csv('test5.csv', inplace=False)
