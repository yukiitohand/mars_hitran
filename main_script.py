"""
/#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

from abscoef import *
import scipy.io as sio


if __name__ == '__main__':
    vert_prof_extended = extend_vert_profile(vert_profile_61C4_btm)

    #    db_begin('data_co')
    #    #fetch_by_ids('CO_3500t5000',[26,27,28,29,30,31],3500,5000,
    #    #          ParameterGroups=['LineMixing','Voigt_self','Voigt_CO2','SDVoigt','Voigt_air'])
    #
    #    nu_strt = 3800
    #    WavenumberStep = 0.0001 #1/cm
    #
    #    nu_end =4500;
    #
    #    #CO_data = hapi_supple.loadTable('data_co/CO_3500t5000')[0]
    #
    #    mid_alt = (vert_profile_61C4_btm['Altitude'][1:] + vert_profile_61C4_btm['Altitude'][0:len(vert_profile_61C4_btm)-1])/2
    #
    #
    #    for i in arange(0,10):
    #        p = P_atm_61C4_btm[i]
    #        p_self = p*vert_profile_61C4_btm['CO_vol_mix_rat'][i]
    #        T = vert_profile_61C4_btm['Temperature'][i]
    #        if i==0:
    #            L = (mid_alt[i]-vert_profile_61C4_btm['Altitude'][i])*100*2
    #        else:
    #            L = (mid_alt[i]-mid_alt[i-1])*100*2
    #        t = time.time()
    #        nu1,abscoef1 = absorptionCoefficient_Voigt_yuki(SourceTables=['CO_3500t5000'],
    #                                               Environment={'p':p,'T':T},
    #                                               WavenumberRange=[nu_strt,nu_end],
    #                                               WavenumberStep = WavenumberStep,
    #                                               WavenumberWing = 0,
    #                                               WavenumberWingHW = 100,
    #                                               HITRAN_units=False,
    #                                               Diluent={'co2':1-vert_profile_61C4['CO_vol_mix_rat'][i],
    #                                                        'self':vert_profile_61C4['CO_vol_mix_rat'][i]})
    #        print(time.time()-t,'[sec]')
    #        if i==0:
    #            nu = nu1
    #            abscoef = abscoef1*vert_profile_61C4_btm['CO_vol_mix_rat'][i]*L
    #        else:
    #            abscoef = abscoef+abscoef1*vert_profile_61C4_btm['CO_vol_mix_rat'][i]*L
    #
    #    mat_fname = 'abscoef_CO_OlympusMons_61C4_W0_hw100_v2.mat'
    #    sio.savemat(mat_fname,{'nu':nu,'abscoef':abscoef})

    ###########################################################################
    #    # CO absorption version 6
    #    db_begin('data_co')
    #    #fetch_by_ids('CO_3500t5000',[26,27,28,29,30,31],3500,5000,
    #    #          ParameterGroups=['LineMixing','Voigt_self','Voigt_CO2','SDVoigt','Voigt_air'])
    #
    #    nu_strt = 3800
    #    WavenumberStep = 0.0001 #1/cm
    #
    #    nu_end =4500;
    #
    #    vert_profile_61C4_btm = vert_prof_extended
    #    P_atm_61C4_btm = pascal2atm(vert_profile_61C4_btm['Pressure'])
    #
    #    mid_alt = (vert_profile_61C4_btm['Altitude'][1:] + vert_profile_61C4_btm['Altitude'][0:-1])/2
    #
    #    mid_alt = np.concatenate(([vert_profile_61C4_btm['Altitude'][0]],mid_alt,[vert_profile_61C4_btm['Altitude'][-1]]))
    #    theta = 57.6
    #    dsts = get_path_length4spherical_layers(mid_alt,theta)
    #
    #    fctr = 0
    #    for i in arange(0,85):
    #        p = P_atm_61C4_btm[i]
    #        T = vert_profile_61C4_btm['Temperature'][i]
    #        factor_pL = volumeConcentration(p,T)
    #        # if i==0:
    #        #     L = (mid_alt[i]-vert_profile_61C4_btm['Altitude'][i])*100*(1+1/np.cos(np.deg2rad(57.6)))
    #        # else:
    #        #     L = (mid_alt[i]-mid_alt[i-1])*100*(1+1/np.cos(np.deg2rad(57.6)))
    #        L = (mid_alt[i+1]-mid_alt[i])*100 + dsts[i]*100
    #        fctr += (factor_pL * vert_profile_61C4_btm['CO_vol_mix_rat'][i] * L)
    #
    #
    #    for i in arange(0,85):
    #        t = time.time()
    #        p = P_atm_61C4_btm[i]
    #        p_self = p*vert_profile_61C4_btm['CO_vol_mix_rat'][i]
    #        T = vert_profile_61C4_btm['Temperature'][i]
    #        # if i==0:
    #        #     L = (mid_alt[i]-vert_profile_61C4_btm['Altitude'][i])*100*(1+1/np.cos(np.deg2rad(57.6)))
    #        # else:
    #        #     L = (mid_alt[i]-mid_alt[i-1])*100*(1+1/np.cos(np.deg2rad(57.6)))
    #        L = (mid_alt[i+1]-mid_alt[i])*100 + dsts[i]*100
    #        nu1,abscoef1,skip_list,wing_exapand_list \
    #        = absorptionCoefficient_Voigt_yuki(SourceTables=['CO_3500t5000'],
    #                                              Environment={'p':p,'T':T},
    #                                              WavenumberRange=[nu_strt,nu_end],
    #                                              WavenumberStep = WavenumberStep,
    #                                              WavenumberWing = 0,
    #                                              WavenumberWingHW = 10,
    #                                              HITRAN_units=False,
    #                                              Diluent={'co2':1-vert_profile_61C4_btm['CO_vol_mix_rat'][i],
    #                                                       'self':vert_profile_61C4_btm['CO_vol_mix_rat'][i]},
    #                                              fctr_tol = fctr,
    #                                              abstol_ctr = 1e-5,
    #                                              abstol_wing = 1e-4)
    #
    #
    #        if i==0:
    #            nu = nu1
    #            abscoef = abscoef1*vert_profile_61C4_btm['CO_vol_mix_rat'][i]*L
    #        else:
    #            abscoef = abscoef+abscoef1*vert_profile_61C4_btm['CO_vol_mix_rat'][i]*L
    #
    #        elapsed = time.time() - t
    #        print(elapsed,'[sec]')
    #
    #    mat_fname = 'abscoef_CO_OlympusMons_61C4_w0_HW10_v6.mat'
    #    sio.savemat(mat_fname,{'nu':nu,'abscoef':abscoef})
    #

    ###############################################################################
    #    """
    #    next we perform the computation of CO2 absorption spectra
    #    """
    #    db_begin('data_co2')
    #    #fetch_by_ids('CO2_2000t12500',[7,8,9,10,11,12,13,14,121,15,120,122],2000,12500,
    #    #         ParameterGroups=['LineMixing','Voigt_self','Voigt_CO2','SDVoigt','Voigt_air'])
    #
    #    nu_strt = 2200
    #    WavenumberStep = 0.0001 #1/cm
    #
    #    nu_end =12000;
    #
    #    vert_profile_61C4_btm = vert_prof_extended
    #    P_atm_61C4_btm = pascal2atm(vert_profile_61C4_btm['Pressure'])
    #
    #    mid_alt = (vert_profile_61C4_btm['Altitude'][1:] + vert_profile_61C4_btm['Altitude'][0:-1])/2
    #
    #    mid_alt = np.concatenate(([vert_profile_61C4_btm['Altitude'][0]],mid_alt,[vert_profile_61C4_btm['Altitude'][-1]]))
    #    theta = 57.6
    #    dsts = get_path_length4spherical_layers(mid_alt,theta)
    #
    #    fctr = 0
    #    for i in arange(0,85):
    #        p = P_atm_61C4_btm[i]
    #        T = vert_profile_61C4_btm['Temperature'][i]
    #        factor_pL = volumeConcentration(p,T)
    #        # if i==0:
    #        #     L = (mid_alt[i]-vert_profile_61C4_btm['Altitude'][i])*100*(1+1/np.cos(np.deg2rad(57.6)))
    #        # else:
    #        #     L = (mid_alt[i]-mid_alt[i-1])*100*(1+1/np.cos(np.deg2rad(57.6)))
    #        L = (mid_alt[i+1]-mid_alt[i])*100 + dsts[i]*100
    #        fctr += (factor_pL * vert_profile_61C4_btm['CO2_vol_mix_rat'][i] * L)
    #
    #
    #    for i in arange(0,85):
    #        t = time.time()
    #        p = P_atm_61C4_btm[i]
    #        p_self = p*vert_profile_61C4_btm['CO2_vol_mix_rat'][i]
    #        T = vert_profile_61C4_btm['Temperature'][i]
    #        # if i==0:
    #        #     L = (mid_alt[i]-vert_profile_61C4_btm['Altitude'][i])*100*(1+1/np.cos(np.deg2rad(57.6)))
    #        # else:
    #        #     L = (mid_alt[i]-mid_alt[i-1])*100*(1+1/np.cos(np.deg2rad(57.6)))
    #        L = (mid_alt[i+1]-mid_alt[i])*100 + dsts[i]*100
    #        nu1,abscoef1,skip_list,wing_exapand_list \
    #        = absorptionCoefficient_Voigt_yuki(SourceTables=['CO2_2000t12500'],
    #                                              Environment={'p':p,'T':T},
    #                                              WavenumberRange=[nu_strt,nu_end],
    #                                              WavenumberStep = WavenumberStep,
    #                                              WavenumberWing = 0,
    #                                              WavenumberWingHW = 10,
    #                                              HITRAN_units=False,
    #                                              Diluent={'self':1},
    #                                              fctr_tol = fctr,
    #                                              abstol_ctr = 1e-5,
    #                                              abstol_wing = 1e-4)
    #        if i==0:
    #            nu = nu1
    #            abscoef = abscoef1*vert_profile_61C4_btm['CO2_vol_mix_rat'][i]*L
    #        else:
    #            abscoef = abscoef+abscoef1*vert_profile_61C4_btm['CO2_vol_mix_rat'][i]*L
    #
    #        elapsed = time.time() - t
    #        print(elapsed,'[sec]')
    #
    #    mat_fname = 'abscoef_CO2_OlympusMons_61C4_w0_HW10_v6.mat'
    #    sio.savemat(mat_fname,{'nu':nu,'abscoef':abscoef})

    ###############################################################################
    #    """
    #    Now converting H2O
    #    """
    #    db_begin('data_h2o')
    #    # fetch_by_ids('H2O_2000t12500',[1,2,3,4,5,6,129],2000,12500,
    #    #          ParameterGroups=['LineMixing','Voigt_self','Voigt_CO2','SDVoigt','Voigt_air'])
    #
    #    nu_strt = 2200
    #    WavenumberStep = 0.0001 #1/cm
    #
    #    nu_end =12000;
    #
    #    mid_alt = (vert_profile_61C4_btm['Altitude'][1:] + vert_profile_61C4_btm['Altitude'][0:len(vert_profile_61C4_btm)-1])/2
    #
    #    for i in arange(0,12):
    #        t = time.time()
    #        p = P_atm_61C4_btm[i]
    #        p_self = p*vert_profile_61C4_btm['H2O_vol_mix_rat'][i]
    #        T = vert_profile_61C4_btm['Temperature'][i]
    #        if i==0:
    #            L = (mid_alt[i]-vert_profile_61C4_btm['Altitude'][i])*100*2
    #        else:
    #            L = (mid_alt[i]-mid_alt[i-1])*100*2
    #        nu1,abscoef1 = absorptionCoefficient_Voigt_yuki(SourceTables=['H2O_2000t12500'],
    #                                              Environment={'p':p,'T':T},
    #                                              WavenumberRange=[nu_strt,nu_end],
    #                                              WavenumberStep = WavenumberStep,
    #                                              WavenumberWing = 0,
    #                                              WavenumberWingHW = 100,
    #                                              HITRAN_units=False,
    #                                              Diluent={'air':1})
    #        if i==0:
    #            nu = nu1
    #            abscoef = abscoef1*vert_profile_61C4_btm['H2O_vol_mix_rat'][i]*L
    #        else:
    #            abscoef = abscoef+abscoef1*vert_profile_61C4_btm['H2O_vol_mix_rat'][i]*L
    #
    #        elapsed = time.time() - t
    #        print(elapsed,'[sec]')
    #
    #    mat_fname = 'abscoef_H2O_OlympusMons_61C4_w0_HW100_v2.mat'
    #    sio.savemat(mat_fname,{'nu':nu,'abscoef':abscoef})
    #

    ###########################################################################
    #    # H2O absorption version 6
    #    db_begin('data_h2o')
    #    # fetch_by_ids('H2O_2000t12500',[1,2,3,4,5,6,129],2000,12500,
    #    #          ParameterGroups=['LineMixing','Voigt_self','Voigt_CO2','SDVoigt','Voigt_air'])
    #
    #    nu_strt = 2200
    #    WavenumberStep = 0.0001 #1/cm
    #
    #    nu_end =12000;
    #
    #    vert_profile_61C4_btm = vert_prof_extended
    #    P_atm_61C4_btm = pascal2atm(vert_profile_61C4_btm['Pressure'])
    #    mid_alt = (vert_profile_61C4_btm['Altitude'][1:] + vert_profile_61C4_btm['Altitude'][0:-1])/2
    #    mid_alt = np.concatenate(([vert_profile_61C4_btm['Altitude'][0]],mid_alt,[vert_profile_61C4_btm['Altitude'][-1]]))
    #    theta = 57.6
    #    dsts = get_path_length4spherical_layers(mid_alt,theta)
    #
    #    fctr = 0
    #    for i in arange(0,85):
    #        p = P_atm_61C4_btm[i]
    #        T = vert_profile_61C4_btm['Temperature'][i]
    #        factor_pL = volumeConcentration(p,T)
    #        # if i==0:
    #        #     L = (mid_alt[i]-vert_profile_61C4_btm['Altitude'][i])*100*(1+1/np.cos(np.deg2rad(57.6)))
    #        # else:
    #        #     L = (mid_alt[i]-mid_alt[i-1])*100*(1+1/np.cos(np.deg2rad(57.6)))
    #        L = (mid_alt[i+1]-mid_alt[i])*100 + dsts[i]*100
    #        fctr += (factor_pL * vert_profile_61C4_btm['H2O_vol_mix_rat'][i] * L)
    #
    #
    #    for i in arange(0,85):
    #        t = time.time()
    #        p = P_atm_61C4_btm[i]
    #        p_self = p*vert_profile_61C4_btm['H2O_vol_mix_rat'][i]
    #        T = vert_profile_61C4_btm['Temperature'][i]
    #        # if i==0:
    #        #     L = (mid_alt[i]-vert_profile_61C4_btm['Altitude'][i])*100*(1+1/np.cos(np.deg2rad(57.6)))
    #        # else:
    #        #     L = (mid_alt[i]-mid_alt[i-1])*100*(1+1/np.cos(np.deg2rad(57.6)))
    #        L = (mid_alt[i+1]-mid_alt[i])*100 + dsts[i]*100
    #        nu1,abscoef1,skip_list,wing_exapand_list \
    #        = absorptionCoefficient_Voigt_yuki(SourceTables=['H2O_2000t12500'],
    #                                              Environment={'p':p,'T':T},
    #                                              WavenumberRange=[nu_strt,nu_end],
    #                                              WavenumberStep = WavenumberStep,
    #                                              WavenumberWing = 0,
    #                                              WavenumberWingHW = 10,
    #                                              HITRAN_units=False,
    #                                              Diluent={'air':1},
    #                                              fctr_tol = fctr,
    #                                              abstol_ctr = 1e-5,
    #                                              abstol_wing = 1e-4)
    #
    #        if i==0:
    #            nu = nu1
    #            abscoef = abscoef1*vert_profile_61C4_btm['H2O_vol_mix_rat'][i]*L
    #        else:
    #            abscoef = abscoef+abscoef1*vert_profile_61C4_btm['H2O_vol_mix_rat'][i]*L
    #
    #        elapsed = time.time() - t
    #        print(elapsed,'[sec]')

    # mat_fname = 'abscoef_H2O_OlympusMons_61C4_w0_HW10_v6.mat'
    # sio.savemat(mat_fname,{'nu':nu,'abscoef':abscoef})

    ###########################################################################
    # N2 absorption version 6
    # db_begin('data_n2')
    # fetch_by_ids('N2_2000t12500',[69,118],2000,12500,
    #          ParameterGroups=['LineMixing','Voigt_self','Voigt_CO2','SDVoigt','Voigt_air'])
    # db_begin('data_h2')
    # db_name = 'H2_2000t12500'
    # fetch_by_ids(db_name,[103,115],2000,12500,
    #           ParameterGroups=['LineMixing','Voigt_self','Voigt_CO2','SDVoigt','Voigt_air'])
    db_begin('data_h2o_mars')
    db_name = 'Water_for_Mars_H12'
    key_mix_rat = 'H2O_vol_mix_rat'
    mat_fname = 'abscoef_H2O_mars_61C4_w0_HW10_v6_nohdo.mat'

    nu_strt = 2200
    WavenumberStep = 0.0001  # 1/cm

    # nu_end =12000;
    nu_end = 5000

    theta_in = 57.6
    theta_em = 0


    vert_profile = vert_prof_extended
    P_atm = pascal2atm(vert_profile['Pressure'])
    mid_alt = (vert_profile['Altitude'][1:] + vert_profile['Altitude'][0:-1]) / 2
    mid_alt = np.concatenate(([vert_profile['Altitude'][0]], mid_alt, [vert_profile['Altitude'][-1]]))

    dsts_in = get_path_length4spherical_layers(mid_alt,theta_in)
    dsts_em = get_path_length4spherical_layers(mid_alt,theta_em)

   fctr = 0
   for i in arange(0,85):
       p = P_atm[i]
       T = vert_profile['Temperature'][i]
       factor_pL = volumeConcentration(p,T)
       #if i==0:
       #    L = (mid_alt[i]-vert_profile_61C4_btm['Altitude'][i])*100*(1+1/np.cos(np.deg2rad(57.6)))
       #else:
       #    L = (mid_alt[i]-mid_alt[i-1])*100*(1+1/np.cos(np.deg2rad(57.6)))
       #L = (dsts_in[i] + dsts_em[i])*100
       #fctr += (factor_pL * vert_profile[key_mix_rat][i] * L)


#
#
#    for i in arange(0,85):
#        t = time.time()
#        p = P_atm[i]
#        p_self = p*vert_profile[key_mix_rat][i]
#        T = vert_profile['Temperature'][i]
#        # if i==0:
#        #     L = (mid_alt[i]-vert_profile_61C4_btm['Altitude'][i])*100*(1+1/np.cos(np.deg2rad(57.6)))
#        # else:
#        #     L = (mid_alt[i]-mid_alt[i-1])*100*(1+1/np.cos(np.deg2rad(57.6)))
#        L = (dsts_in[i] + dsts_em[i])*100
#        nu1,abscoef1,skip_list,wing_exapand_list \
#        = absorptionCoefficient_Voigt_yuki(SourceTables=[db_name],Components=[(1,1),(1,2),(1,3),(1,4),(1,5),(1,6)],
#                                              Environment={'p':p,'T':T},
#                                              WavenumberRange=[nu_strt,nu_end],
#                                              WavenumberStep = WavenumberStep,
#                                              WavenumberWing = 0,
#                                              WavenumberWingHW = 10,
#                                              HITRAN_units=False,
#                                              Diluent={'air':1},
#                                              fctr_tol = fctr,
#                                              abstol_ctr = 1e-5,
#                                              abstol_wing = 1e-4)
#
#        if i==0:
#            nu = nu1
#            abscoef = abscoef1*vert_profile[key_mix_rat][i]*L
#        else:
#            abscoef = abscoef+abscoef1*vert_profile[key_mix_rat][i]*L
#
#        elapsed = time.time() - t
#        print(elapsed,'[sec]')
#
#
#    sio.savemat(mat_fname,{'nu':nu,'abscoef':abscoef})