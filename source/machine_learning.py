from pymatgen import Structure, Lattice
from pymatgen.core import periodic_table
from pymatgen.analysis import ewald
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import ShuffleSplit, cross_val_score, train_test_split

train = pd.read_csv('train5.csv')
test = pd.read_csv('test5.csv')

def ml_regression(train,test,x_feature,y_feature,model,feature_importance=True):
    
    x_train = train[x_feature].values
    y_train = train[y_feature].values
    x_test = test[x_feature].values
    
    print('Defining RFR model...\n')

    model.fit(x_train,y_train)
    score = model.score(x_train,y_train)
    print('Training Score = ', score)

    # CV score
    cv = ShuffleSplit(len(y_train), n_iter=10, test_size=0.1, random_state=0)
    cvscore = cross_val_score(model, x_train, y_train, cv=cv, scoring='r2')
    print('CV Score = ', cvscore)
    print('Mean = ', np.mean(cvscore))
    
    predict = model.predict(x_test)
        
    if feature_importance:
        importances = pd.DataFrame({'feature':x_feature,'importance':np.round(model.feature_importances_,5)})
        importances = importances.sort_values('importance',ascending=False).set_index('feature')
        print(importances)
        importances.plot.bar()
        plt.show()
        
    return model, predict

def submission(predict_bandgap,predict_formation):
    submission = pd.read_csv('sample_submission.csv')
    submission['bandgap_energy_ev'] = predict_bandgap
    submission['formation_energy_ev_natom'] = predict_formation
    submission.to_csv('submission.csv',index=False)

def split_train(train,x_feature,y_feature,model):
    train_2, test_2 = train_test_split(train ,train_size = 0.8)
    x_train_2 = train_2[x_feature].values
    y_train_2 = train_2[y_feature].values
    x_test_2 = test_2[x_feature].values
    y_test_2 = test_2[y_feature].values
    model.fit(x_train_2,y_train_2)
    predict_2 = model.predict(x_test_2)
    plt.scatter(y_test_2,predict_2)
    plt.plot([0,train_2[y_feature].values.max()],[0,train_2[y_feature].values.max()],'--')
    plt.xlabel('Actual')
    plt.ylabel('Predict')
    return plt.show()

features_bandgap = ['mass_al','mass_ga','mass_in','mass_tot','vol','spacegroup', 'percent_atom_al', 'percent_atom_in','density','z_al','z_in','z_tot','inv_r_alo_sum','calo','inv_r_ino_sum','cino','inv2_r_sum','r_min','rmean','nalo_mean','nino_mean','n_mean',
                    'coul','calal','cinin','cgaal','cgain','calin']
features_formation = ['mass_al','mass_ga','mass_in','vol','spacegroup', 'percent_atom_ga','z_ga','inv_r_gao_sum','cgao','r_min','ngao_mean','n_mean',
                      'cgaga','calal','cinin','cgaal','cgain','bandgap_energy_ev']
print('Bandgap features = ', features_bandgap)
print('Formation features = ', features_formation)

print('Bandap Energy:')
model_bandgap, predict_bandgap = ml_regression(train,test,features_bandgap,'bandgap_energy_ev',RandomForestRegressor(n_estimators=1000))
#model_bandgap, predict_bandgap = ml_regression(train,test,features_bandgap,'bandgap_energy_ev',GradientBoostingRegressor(n_estimators=100))
test['bandgap_energy_ev'] = pd.DataFrame(predict_bandgap)
print('Formation Energy')
#model_formation, predict_formation = ml_regression(train,test,features_formation,'formation_energy_ev_natom',RandomForestRegressor(n_estimators=100))
model_formation, predict_formation = ml_regression(train,test,features_formation,'formation_energy_ev_natom',GradientBoostingRegressor(n_estimators=500))

split_train(train,features_bandgap,'bandgap_energy_ev',RandomForestRegressor(n_estimators=1000))
split_train(train,features_formation,'formation_energy_ev_natom',GradientBoostingRegressor(n_estimators=500))

#submission(predict_bandgap,predict_formation)
