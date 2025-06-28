# %%
from hapi import *
import abscoef as abc
import scipy.io as sio
import argparse


def mcd_crism_get_atmabsorption(mol_names, fname_List):

    # %%
    # fname_List = ['mcd_output/FRT0000B3CC_07_IF165L_TRR3_mean.mat']
    # fname_List = ['mcd_output/HRS00004025_07_IF174L_TRR3.mat']

    # fname_List = ['mcd_output/HRL000040FF_07_IF183L_TRR3.mat']
    # fname_List = ['/Users/yukiitoh/src/matlab/mcd_crism/FRT0000AD16_07_IF165L_TRR3_mean.mat',
    #               '/Users/yukiitoh/src/matlab/mcd_crism/FRT0000B4B5_07_IF165L_TRR3_mean.mat']
    # fname_List = ['/Users/itohy1/src/matlab/mcd_crism/FRT000030B0_07_IF168L_TRR3_mean.mat']
    # fname_List = ['/Users/yukiitoh/src/matlab/mcd_crism/FRT00002879_07_IF168L_TRR3_mean.mat']
    # fname_List = ['/Users/itohy1/src/matlab/mcd_data/HRL000040A2_07_IF182L_TRR3_mean.mat']
    # fname_List = ['/Users/itohy1/src/matlab/mcd_crism/FRT0000C286_07_IF168L_TRR3_mean.mat']
    # fname_List = [
    #     '/Users/itohy1/src/matlab/mcd_crism/FRT0000C36F_07_IF168L_TRR3_mean.mat',
    #     ]
    # fname_List = [
    #     '/Users/itohy1/src/matlab/mcd_crism/HRL0000CFB5_07_IF185L_TRR3_mean.mat',
    #     '/Users/itohy1/src/matlab/mcd_crism/HRL0000C3EE_07_IF185L_TRR3_mean.mat',
    #     '/Users/itohy1/src/matlab/mcd_crism/HRL0000C67E_07_IF185L_TRR3_mean.mat',
    #     '/Users/itohy1/src/matlab/mcd_crism/HRS0000C628_07_IF177L_TRR3_mean.mat'
    #     ]


    # './mcd_output/FFC00017E04_01_IF254L_TRR3_l100t200.mat'
    # fname = '/Users/yukiitoh/src/matlab/mcd_crism/FFC000061C4_01_IF254L_TRR3_l100t200.mat'
    ########################################################################################################################
    # set up parameters
    ########################################################################################################################
    # general parameters

    # mol_names = ['co2']
    # mol_name = 'co2'  # {'co2','co','h2o'}
    # diluent_mode_list = [2]
    # co2_diluent_mode = 1 # {1,2,3}
    diluent_mode_co2 = 2
    diluent_mode_h2o = 1
    opt_h2o = 1
    # db_dir = 'data/data_' + mol_name
    # db_name = 'CO2_2000t12500'  # {'CO2_2000t12500','CO_3500t5000','H2O_2000t12500',}
    # key_mix_rat = mol_name + '_vol_mix_rat'
    rowIDs = None # [277819], None
    n_correction_factor = 1.0
    gamma0_cor = 1.0
    # rowIDs = None


    # voigt parameters
    # nu_strt = 2200.  # {2200 for {CO2,H2O}; 3800 for {CO}}
    # nu_strt = 12000
    # nu_end  = 12000.  # {12000 for {CO2,H2O}; 4500 for {CO}}
    # nu_end = 5200
    WavenumberStep = 0.0001 #1/cm
    WavenumberWing = 0
    WavenumberWingHW = 10.

    abstol_ctr = 1e-5
    abstol_wing = 1e-4

    shape_profile = 'voigt' # {'voigt','rautian','rautian2','Lorentzian'}
    beta_list = [5]
    beta = 5

    ema_mode = 'ctr'  # {'ctr','mean'}

    # Diluent_dict = {'co2':{'self':1}, 'co':... (need to determine later), 'h2o':{'air':1}}
    # Diluent = Diluent_dict[mol_name]


    # fetch database before you begin
    # fetch_by_ids('CO2_2000t12500',[7,8,9,10,11,12,13,14,121,15,120,122],2000,12500,
    #              ParameterGroups=['160-char', 'Standard','Voigt_self','Voigt_CO2','SDVoigt','Voigt_air','Galatry'])
    # fetch_by_ids('CO_3500t5000',[26,27,28,29,30,31],3500,5000,
    #              ParameterGroups=['160-char', 'Standard', 'Voigt_self', 'Voigt_CO2', 'SDVoigt', 'Voigt_air'])
    # fetch_by_ids('H2O_2000t12500',[1,2,3,4,5,6,129],2000,12500,
    #              ParameterGroups=['160-char', 'Standard'])

    # db_begin(db_dir)

    db_opt = 2020

    if db_opt == 2016:
        db_pdir = './data'
    elif db_opt == 2020:
        db_pdir = './data2020'


    for fname in fname_List:
        data = sio.loadmat(fname)
        fname_split = fname.split('.')
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
        if ema_mode == 'mean':
            theta_em = ema_areoid
        elif ema_mode == 'ctr':
            theta_em = data['ema_areoid_ctr'].reshape(-1)
        dsts_in = abc.get_path_length4spherical_layers(alt_brd_List, theta_in)
        dsts_em = abc.get_path_length4spherical_layers(alt_brd_List, theta_em)

        for mol_name in mol_names:
            if mol_name == 'co':
                nu_strt = 3800
                nu_end = 4500
            elif mol_name in ['h2o','co2']:
                nu_strt = 2200
                nu_end = 12000

            key_mix_rat = mol_name + '_vol_mix_rat'
            db_dir = os.path.join(db_pdir,'data_' + mol_name)
            if mol_name == 'co2':
                mat_suffix = f'_abscoef_{db_opt}_{mol_name}_{shape_profile}_w0_HW10_v7_diluent_{diluent_mode_co2}_{n_correction_factor}_{ema_mode}'
                db_name = 'CO2_2000t12500'
                db_dir = os.path.join(db_pdir,'data_' + mol_name)
            elif mol_name == 'h2o':
                if opt_h2o == 1:
                    mat_suffix = f'_abscoef_{db_opt}_{mol_name}_{shape_profile}_w0_HW10_v7_{opt_h2o}_d{diluent_mode_h2o:1d}_{n_correction_factor}_{ema_mode}'
                    db_name = 'H2O_2000t12500'
                    db_dir = os.path.join(db_pdir,'data_' + mol_name)
                elif opt_h2o == 2:
                    mat_suffix = f'_abscoef_{db_opt}_{mol_name}_{shape_profile}_w0_HW10_v7_{opt_h2o}_d{diluent_mode_h2o:1d}_{n_correction_factor}_{ema_mode}'
                    db_name = 'Water_for_Mars_H12'
                    db_dir = 'data/data_h2o_mars'
            elif mol_name == 'co':
                db_name = 'CO_3500t5000'
                mat_suffix = f'_abscoef_{db_opt}_{mol_name}_{shape_profile}_w0_HW10_v7_{n_correction_factor}_{ema_mode}'
                db_dir = os.path.join(db_pdir,'data_' + mol_name)

            db_begin(db_dir)

            print(mat_suffix)
            fctr = 0
            for i in arange(0, len(data_ar)):
                p_atm = p_atm_ar[i]
                T = data_ar['Temperature'][i]
                factor_pL = volumeConcentration(p_atm, T)
                L = (dsts_in[i] + dsts_em[i])*100
                fctr += (factor_pL * data_ar[key_mix_rat][i] * L)

            # inverse the order
            is_first_idx = True
            for i in arange(len(data_ar)-1, -1, -1):
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
                    if diluent_mode_h2o == 1:
                        Diluent = {'air': 1.}
                    elif diluent_mode_h2o == 2:
                        Diluent = {'self': 1.}
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

                if is_first_idx:
                    nu = nu1
                    abscoeff = abscoef1 * data_ar[key_mix_rat][i] * L
                    is_first_idx=False
                abscoeff = abscoeff + abscoef1 * data_ar[key_mix_rat][i] * L
                elapsed = abc.time.time() - t
                print(elapsed,'[sec]')

                # if p_atm > 0.0050 or i==0:
                if i == 0:
                    mat_fname = fname_split[0] + mat_suffix + f'_idx{i}' + '.mat'
                    sio.savemat(mat_fname,{'nu':nu, 'abscoef':abscoeff})


if __name__ == '__main__':
    
    mol_names = ['co','h2o','co2']

    # parser = argparse.ArgumentParser()
    # parser.add_argument('filename')
    # args = parser.parse_args()
    # print(args.filename)
    # mcd_crism_get_atmabsorption(mol_names, [args.filename])
    # filename_list = [
    #     '/Volumes/LaCie5TB/out/mcd/HRS00002F21_07_IF177L_TRR3_mean.mat',
    #     '/Volumes/LaCie5TB/out/mcd/FRT0000CC76_07_IF168L_TRR3_mean.mat',
    #     '/Volumes/LaCie5TB/out/mcd/FRT00002854_07_IF168L_TRR3_mean.mat',
    # ]
    filename_list = [
        #"/Volumes/LaCie5TB/out/mcd/FRT0000CBA8_07_IF168L_TRR3_mean_MCDv6_1_scena1.mat",
        #'/Volumes/LaCie5TB/out/mcd/FRT0000C0F6_07_IF168L_TRR3_mean_MCDv6_1_scena1.mat',
        '/Volumes/LaCie5TB/out/mcd/FRT00002854_07_IF168L_TRR3_mean_MCDv6_1_scena1.mat',
        '/Volumes/LaCie5TB/out/mcd/HRS0000C1CD_07_IF177L_TRR3_mean_MCDv6_1_scena1.mat',
        '/Volumes/LaCie5TB/out/mcd/HRS00002F21_07_IF177L_TRR3_mean_MCDv6_1_scena1.mat',
        '/Volumes/LaCie5TB/out/mcd/FRT0000CC76_07_IF168L_TRR3_mean_MCDv6_1_scena1.mat',
        '/Volumes/LaCie5TB/out/mcd/FRT0000C0B5_07_IF168L_TRR3_mean_MCDv6_1_scena1.mat',
        '/Volumes/LaCie5TB/out/mcd/FRT0000C36F_07_IF168L_TRR3_mean_MCDv6_1_scena1.mat',
        '/Volumes/LaCie5TB/out/mcd/FRT0000C52D_07_IF168L_TRR3_mean_MCDv6_1_scena1.mat',
    ]
    filename_list = [
        '/Volumes/LaCie5TB/out/mcd/FRT00003452_07_IF168L_TRR3_mean_MCDv6_1_scena1.mat'
    ]
    mcd_crism_get_atmabsorption(mol_names, filename_list)
