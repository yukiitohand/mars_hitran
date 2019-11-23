from hapi import *
import abscoef as abc
import scipy.io as sio

fname_List = ['/Users/yukiitoh/src/matlab/mcd_crism/FFC00017E04_01_IF254L_TRR3_l100t200.mat',
              '/Users/yukiitoh/src/matlab/mcd_crism/FRT00017D33_07_IF165L_TRR3.mat']

# fname = '/Users/yukiitoh/src/matlab/mcd_crism/FFC000061C4_01_IF254L_TRR3_l100t200.mat'
########################################################################################################################
# set up parameters
########################################################################################################################
# general parameters
mol_names = ['h2o']
# mol_name = 'co2'  # {'co2','co','h2o'}
# diluent_mode_list = [2]
# co2_diluent_mode = 1 # {1,2,3}
diluent_mode_co2 = 2
opt_h2o = 2
# db_dir = 'data/data_' + mol_name
# db_name = 'CO2_2000t12500'  # {'CO2_2000t12500','CO_3500t5000','H2O_2000t12500',}
# key_mix_rat = mol_name + '_vol_mix_rat'
rowIDs = None # [277819], None
n_correction_factor = 1.0
gamma0_cor = 1.0
# rowIDs = None


# voigt parameters
nu_strt = 2200  # {2200 for {CO2,H2O}; 3800 for {CO}}
# nu_strt = 4800
nu_end  = 12000  # {12000 for {CO2,H2O}; 3800 for {CO}}
# nu_end = 5200
WavenumberStep = 0.0001 #1/cm
WavenumberWing = 0
WavenumberWingHW = 10

abstol_ctr = 1e-5
abstol_wing = 1e-4

shape_profile = 'voigt' # {'voigt','rautian','rautian2'}
beta_list = [5]
beta = 5

# Diluent_dict = {'co2':{'self':1}, 'co':... (need to determine later), 'h2o':{'air':1}}
# Diluent = Diluent_dict[mol_name]


# fetch database before you begin
# fetch_by_ids('CO2_4800t5200',[7,8,9,10,11,12,13,14,121,15,120,122],4800,5200,
#              ParameterGroups=['LineMixing','Voigt_self','Voigt_CO2','SDVoigt','Voigt_air','Galatry'])
# fetch_by_ids('CO_3500t5000',[26,27,28,29,30,31],3500,5000,
#              ParameterGroups=['LineMixing','Voigt_self','Voigt_CO2','SDVoigt','Voigt_air'])
# fetch_by_ids('H2O_2000t12500',[1,2,3,4,5,6,129],2000,12500,
#              ParameterGroups=['LineMixing','Voigt_self','Voigt_CO2','SDVoigt','Voigt_air'])

# db_begin(db_dir)

for fname in fname_List:
    data = sio.loadmat(fname)
    # 'alt_surf','alt_surf_mid','alt_areoid','L_mid','T_mid','p_Pa_mid',...
    #     'co2_r_mid','co_r_mid','h2ovapor_r_mid','ina_areoid','ema_areoid'
    alt_brd_List = data['alt_areoid'].reshape(-1)
    alt_mid_List = data['alt_surf_mid'].reshape(-1)
    T_List = data['T_mid'].reshape(-1)
    p_Pa_List = data['p_Pa_mid'].reshape(-1)
    co2_r_List = data['co2_r_mid'].reshape(-1)
    co_r_List = data['co_r_mid'].reshape(-1)
    h2ovapor_r_List = data['h2ovapor_r_mid'].reshape(-1)
    ina_areoid = data['ina_areoid'].reshape(-1)
    ema_areoid = data['ema_areoid'].reshape(-1)

    data_ar = np.zeros(alt_mid_List.shape[0],
                       dtype=[('Altitude_surf', 'float'), ('Pressure_Pa', 'float'), ('Temperature', 'float'),
                              ('co2_vol_mix_rat', 'float'), ('co_vol_mix_rat', 'float'),
                              ('h2o_vol_mix_rat', 'float')
                              ])

    data_ar['Altitude_surf'] = alt_mid_List
    data_ar['Pressure_Pa'] = p_Pa_List
    data_ar['Temperature'] = T_List
    data_ar['co2_vol_mix_rat'] = co2_r_List
    data_ar['co_vol_mix_rat'] = co_r_List
    data_ar['h2o_vol_mix_rat'] = h2ovapor_r_List

    ###########################################################################
    # Absorption version 7
    p_atm_ar = abc.pascal2atm(data_ar['Pressure_Pa'])
    theta_in = ina_areoid
    theta_em = ema_areoid
    dsts_in = abc.get_path_length4spherical_layers(alt_brd_List, theta_in)
    dsts_em = abc.get_path_length4spherical_layers(alt_brd_List, theta_em)

    for mol_name in mol_names:
        key_mix_rat = mol_name + '_vol_mix_rat'
        db_dir = 'data/data_' + mol_name
        if mol_name == 'co2':
            mat_suffix = f'_abscoef_{mol_name}_{shape_profile}_w0_HW10_v7_diluent_{diluent_mode_co2}_{n_correction_factor}'
            db_name = 'CO2_2000t12500'
            db_dir = 'data/data_' + mol_name
        elif mol_name == 'h2o':
            if opt_h2o == 1:
                mat_suffix = f'_abscoef_{mol_name}_{shape_profile}_w0_HW10_v7_{opt_h2o}_{n_correction_factor}'
                db_name = 'H2O_2000t12500'
                db_dir = 'data/data_' + mol_name
            elif opt_h2o == 2:
                mat_suffix = f'_abscoef_{mol_name}_{shape_profile}_w0_HW10_v7_{opt_h2o}_{n_correction_factor}'
                db_name = 'Water_for_Mars_H12'
                db_dir = 'data/data_h2o_mars'
        elif mol_name == 'co':
            db_name = 'CO_3500t5000'
            mat_suffix = f'_abscoef_{mol_name}_{shape_profile}_w0_HW10_v7_{n_correction_factor}'
            db_dir = 'data/data_' + mol_name

        db_begin(db_dir)

        print(mat_suffix)
        fctr = 0
        for i in arange(0, len(data_ar)):
            p_atm = p_atm_ar[i]
            T = data_ar['Temperature'][i]
            factor_pL = volumeConcentration(p_atm, T)
            L = (dsts_in[i] + dsts_em[i])*100
            fctr += (factor_pL * data_ar[key_mix_rat][i] * L)

        for i in arange(0,len(data_ar)):
            t = abc.time.time()
            p_atm = p_atm_ar[i]
            p_self = p_atm * data_ar[key_mix_rat][i]
            T = data_ar['Temperature'][i]
            L = (dsts_in[i] + dsts_em[i])*100

            if mol_name == 'co2':
                if diluent_mode_co2 == 1:
                    Diluent = {'self': 1}
                elif diluent_mode_co2 == 2:
                    Diluent = {'self': data_ar[key_mix_rat][i], 'air': 1 - data_ar[key_mix_rat][i]}
                elif diluent_mode_co2 == 3:
                    Diluent = {'self': data_ar[key_mix_rat][i]}
            elif mol_name == 'h2o':
                Diluent = {'air': 1}
            elif mol_name == 'co':
                Diluent = {'co2': 1-data_ar[key_mix_rat][i], 'self': data_ar[key_mix_rat][i]}

            nu1, abscoef1, gamma0_list, sgmd_list, skip_list, wing_expand_list,fac_list \
                = abc.absorptionCoefficient_Voigt_yuki(SourceTables=[db_name],
                                                        Environment={'p':p_atm, 'T':T},
                                                        WavenumberRange=[nu_strt,nu_end],
                                                        WavenumberStep = WavenumberStep,
                                                        WavenumberWing = WavenumberWing,
                                                        WavenumberWingHW = WavenumberWingHW,
                                                        HITRAN_units=False,
                                                        Diluent=Diluent,
                                                        fctr_tol = fctr,
                                                        abstol_ctr=abstol_ctr,
                                                        abstol_wing=abstol_wing,
                                                        shape_profile=shape_profile, beta=beta,
                                                        RowIDs=rowIDs, n_correction_factor=n_correction_factor,
                                                       gamma0_cor=gamma0_cor)

            if i==0:
                nu = nu1
                abscoeff = abscoef1 * data_ar[key_mix_rat][i] * L
            abscoeff = abscoeff + abscoef1 * data_ar[key_mix_rat][i] * L
            elapsed = abc.time.time() - t
            print(elapsed,'[sec]')
        fname_split = fname.split('.')
        mat_fname = fname_split[0] + mat_suffix + '.mat'
        sio.savemat(mat_fname,{'nu':nu,'abscoef':abscoeff})