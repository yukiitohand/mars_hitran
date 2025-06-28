#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 18 16:52:52 2018

@author: yukiitoh
"""
import sys

# sys.path.insert(0, '/Users/yukiitoh/src/python/voigt_cy')
from hapi import *
# import hapi_supple
import scipy.io as sio
import matplotlib.pyplot as plt
import numpy as np
import time
import pickle
import time
from scipy.special import wofz
# import voigt_cy.voigt_cy as voigt_cy

from mars_vertical_profile import *

# define precision
__ComplexType__ = complex128
__IntegerType__ = int64
__FloatType__ = float64


def get_DB_value_robust(TableName, vartypename, species, RowID, subst_species='air', missing_substitute=0.):
    entry_name = vartypename + '_' + species.lower()
    sub_etrynm = vartypename + '_' + subst_species.lower()
    if entry_name in LOCAL_TABLE_CACHE[TableName]['data'].keys():
        if not np.ma.is_masked(LOCAL_TABLE_CACHE[TableName]['data'][entry_name][RowID]):
            varval = LOCAL_TABLE_CACHE[TableName]['data'][entry_name][RowID]
            varname = entry_name
        else:
            if sub_etrynm in LOCAL_TABLE_CACHE[TableName]['data'].keys():
                if not np.ma.is_masked(LOCAL_TABLE_CACHE[TableName]['data'][sub_etrynm][RowID]):
                    varval = LOCAL_TABLE_CACHE[TableName]['data'][sub_etrynm][RowID]
                    varname = sub_etrynm
                else:
                    varval = missing_substitute
                    varname = ''
            else:
                varval = missing_substitute
                varname = ''
    else:
        if sub_etrynm in LOCAL_TABLE_CACHE[TableName]['data'].keys():
            if not np.ma.is_masked(LOCAL_TABLE_CACHE[TableName]['data'][sub_etrynm][RowID]):
                varval = LOCAL_TABLE_CACHE[TableName]['data'][sub_etrynm][RowID]
                varname = sub_etrynm
            else:
                varval = missing_substitute
                varname = ''
        else:
            varval = missing_substitute
            varname = ''
    
    return varval, varname

def get_path_length4spherical_layers(layer_boundary_heights, theta, rMars=3396190.0):
    # modeled as spherical layers.
    # Solve a quadratic equation
    if theta == 0:
        dsts = layer_boundary_heights[2:-1] - layer_boundary_heights[1:-2]
    else:
        tan_theta = np.tan(np.deg2rad(theta))
        a = 1 + 1 / tan_theta ** 2
        b = 2 * rMars / tan_theta
        c = -(2 * rMars + layer_boundary_heights) * layer_boundary_heights

        xhat = (-b / 2 + sqrt((b / 2) ** 2 - a * c)) / a
        yhat = xhat / tan_theta + rMars
        dsts = sqrt((xhat[1:] - xhat[0:-1]) ** 2 + (yhat[1:] - yhat[0:-1]) ** 2)

    return dsts

def V(x, sigma, gamma):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    """
    # sigma = alpha / np.sqrt(2 * np.log(2))

    return np.real(wofz((x + 1j*gamma)/sigma/np.sqrt(2))) / sigma\
                                                           /np.sqrt(2*np.pi)

def PROFILE_VOIGT_w(LineCenter, SigmaD, Gamma0, nu_grid):
    return np.real(wofz((nu_grid - (LineCenter - 1j * Gamma0)) / (SigmaD * np.sqrt(2)))) / (SigmaD * np.sqrt(2 * np.pi))


# def PROFILE_VOIGT_cy(LineCenter, SigmaD, Gamma0, nu_grid):
#     return np.real(voigt_cy.wofz(np.complex128((nu_grid - (LineCenter - 1j * Gamma0)) / (SigmaD * np.sqrt(2))) )) / (
#                 SigmaD * np.sqrt(2 * np.pi))

# def PROFILE_RAUTIAN_cy(LineCenter, SigmaD, Gamma0, nu_grid, beta):
#     beta_p = beta * Gamma0 / (SigmaD * np.sqrt(2))
#     # print(beta_p)
#     wz = voigt_cy.wofz((nu_grid - (LineCenter - 1j * (Gamma0*(1+beta)))) / (SigmaD * np.sqrt(2)))
#     return np.real( wz / (1-np.sqrt(pi)*beta_p*wz) ) / (SigmaD * np.sqrt(2 * np.pi))

# def PROFILE_RAUTIAN_cy2(LineCenter, SigmaD, Gamma0, nu_grid, beta):
#     beta_p = beta / (SigmaD * np.sqrt(2))
#     wz = voigt_cy.wofz((nu_grid - (LineCenter - 1j * (Gamma0+beta))) / (SigmaD * np.sqrt(2)))
#     return np.real( wz / (1-np.sqrt(pi)*beta_p*wz) ) / (SigmaD * np.sqrt(2 * np.pi))


def absorptionCoefficient_Voigt_yuki(Components=None, SourceTables=None, partitionFunction=PYTIPS2017,
                                     Environment=None, OmegaRange=None, OmegaStep=None, OmegaWing=None,
                                     IntensityThreshold=DefaultIntensityThreshold,
                                     OmegaWingHW=DefaultOmegaWingHW,
                                     GammaL='gamma_air', HITRAN_units=True, LineShift=True,
                                     File=None, Format=None, OmegaGrid=None,
                                     WavenumberRange=None, WavenumberStep=None, WavenumberWing=None,
                                     WavenumberWingHW=None, WavenumberGrid=None,
                                     Diluent={}, EnvDependences=None,
                                     fctr_tol=1, abstol_ctr=0., abstol_wing=0.,
                                     shape_profile='voigt', beta=0, RowIDs=None,
                                     n_correction_factor=1.,gamma0_cor=1.):
    """
    INPUT PARAMETERS:
        Components:  list of tuples [(M,I,D)], where
                        M - HITRAN molecule number,
                        I - HITRAN isotopologue number,
                        D - relative abundance (optional)
        SourceTables:  list of tables from which to calculate cross-section   (optional)
        partitionFunction:  pointer to partition function (default is PYTIPS) (optional)
        Environment:  dictionary containing thermodynamic parameters.
                        'p' - pressure in atmospheres,
                        'T' - temperature in Kelvin
                        Default={'p':1.,'T':296.}
        WavenumberRange:  wavenumber range to consider.
        WavenumberStep:   wavenumber step to consider.
        WavenumberWing:   absolute wing for calculating a lineshape (in cm-1)
        WavenumberWingHW:  relative wing for calculating a lineshape (in halfwidths)
        IntensityThreshold:  threshold for intensities
        GammaL:  specifies broadening parameter ('gamma_air' or 'gamma_self')
        HITRAN_units:  use cm2/molecule (True) or cm-1 (False) for absorption coefficient
        File:   write output to file (if specified)
        Format:  c-format of file output (accounts for significant digits in WavenumberStep)
    OUTPUT PARAMETERS:
        Wavenum: wavenumber grid with respect to parameters WavenumberRange and WavenumberStep
        Xsect: absorption coefficient calculated on the grid
    ---
    DESCRIPTION:
        Calculate absorption coefficient using Voigt profile.
        Absorption coefficient is calculated at arbitrary temperature and pressure.
        User can vary a wide range of parameters to control a process of calculation.
        The choise of these parameters depends on properties of a particular linelist.
        Default values are a sort of guess which gives a decent precision (on average)
        for a reasonable amount of cpu time. To increase calculation accuracy,
        user should use a trial and error method.
    ---
    EXAMPLE OF USAGE:
        nu,coef = absorptionCoefficient_Voigt(((2,1),),'co2',WavenumberStep=0.01,
                                              HITRAN_units=False,GammaL='gamma_self')
    ---
    """

    # Paremeters OmegaRange,OmegaStep,OmegaWing,OmegaWingHW, and OmegaGrid
    # are deprecated and given for backward compatibility with the older versions.
    if WavenumberRange:  OmegaRange = WavenumberRange
    if WavenumberStep:   OmegaStep = WavenumberStep
    if WavenumberWing:   OmegaWing = WavenumberWing
    if WavenumberWingHW: OmegaWingHW = WavenumberWingHW
    if WavenumberGrid:   OmegaGrid = WavenumberGrid

    # "bug" with 1-element list
    Components = listOfTuples(Components)
    SourceTables = listOfTuples(SourceTables)

    # determine final input values
    Components, SourceTables, Environment, OmegaRange, OmegaStep, OmegaWing, \
    IntensityThreshold, Format = \
        getDefaultValuesForXsect(Components, SourceTables, Environment, OmegaRange,
                                 OmegaStep, OmegaWing, IntensityThreshold, Format)

    # warn user about too large omega step
    if OmegaStep > 0.1: warn('Big wavenumber step: possible accuracy decline')

    # get uniform linespace for cross-section
    # number_of_points = (OmegaRange[1]-OmegaRange[0])/OmegaStep + 1
    # Omegas = linspace(OmegaRange[0],OmegaRange[1],number_of_points)
    if OmegaGrid is not None:
        Omegas = npsort(OmegaGrid)
    else:
        # Omegas = arange(OmegaRange[0],OmegaRange[1],OmegaStep)
        Omegas = arange_(OmegaRange[0], OmegaRange[1], OmegaStep)  # fix
    number_of_points = len(Omegas)
    Xsect = zeros(number_of_points)

    # reference temperature and pressure
    Tref = np.float64(296.)  # K
    pref = np.float64(1.)  # atm

    # actual temperature and pressure
    T = Environment['T']  # K
    p = Environment['p']  # atm

    # create dictionary from Components
    ABUNDANCES = {}
    NATURAL_ABUNDANCES = {}
    for Component in Components:
        M = Component[0]
        I = Component[1]
        if len(Component) >= 3:
            ni = Component[2]
        else:
            try:
                ni = ISO[(M, I)][ISO_INDEX['abundance']]
            except KeyError:
                raise Exception('cannot find component M,I = %d,%d.' % (M, I))
        ABUNDANCES[(M, I)] = ni
        NATURAL_ABUNDANCES[(M, I)] = ISO[(M, I)][ISO_INDEX['abundance']]

    # precalculation of volume concentration
    if HITRAN_units:
        factor = __FloatType__(1.0)
    else:
        factor = volumeConcentration(p, T)

        # setup the default empty environment dependence function
    if not EnvDependences:
        EnvDependences = lambda ENV, LINE: {}
    Env = Environment.copy()
    Env['Tref'] = Tref
    Env['pref'] = pref

    # setup the Diluent variable
    GammaL = GammaL.lower()
    if not Diluent:
        if GammaL == 'gamma_air':
            Diluent = {'air': 1.}
        elif GammaL == 'gamma_self':
            Diluent = {'self': 1.}
        else:
            raise Exception('Unknown GammaL value: %s' % GammaL)

    # Simple check
    print(Diluent)  # Added print statement # CHANGED RJH 23MAR18  # Simple check
    for key in Diluent:
        val = Diluent[key]
        if val < 0 or val > 1:  # if val < 0 and val > 1:# CHANGED RJH 23MAR18
            raise Exception('Diluent fraction must be in [0,1]')

    # SourceTables contain multiple tables
    for TableName in SourceTables:

        # get the number of rows
        if not RowIDs:
            nline = LOCAL_TABLE_CACHE[TableName]['header']['number_of_rows']
            line_list = np.arange(0, nline)
        else:
            nline = len(RowIDs)
            line_list = RowIDs

        skip_list = []
        wing_expand_list = np.ones(nline, )
        Gamma0_list = np.zeros(nline, )
        SgmD_List = np.zeros(nline, )
        facList = np.zeros(nline, )

        # get parameter names for each table
        parnames = LOCAL_TABLE_CACHE[TableName]['data'].keys()

        # loop through line centers (single stream)
        for rowid_i, RowID in enumerate(line_list):

            # Get the custom environment dependences
            Line = {}
            for parname in parnames:
                Line[parname] = LOCAL_TABLE_CACHE[TableName]['data'][parname][RowID]
            CustomEnvDependences = EnvDependences(Env, Line)

            # get basic line parameters (lower level)
            LineCenterDB = LOCAL_TABLE_CACHE[TableName]['data']['nu'][RowID]
            LineIntensityDB = LOCAL_TABLE_CACHE[TableName]['data']['sw'][RowID]
            LowerStateEnergyDB = LOCAL_TABLE_CACHE[TableName]['data']['elower'][RowID]
            MoleculeNumberDB = LOCAL_TABLE_CACHE[TableName]['data']['molec_id'][RowID]
            IsoNumberDB = LOCAL_TABLE_CACHE[TableName]['data']['local_iso_id'][RowID]

            # filter by molecule and isotopologue
            if (MoleculeNumberDB, IsoNumberDB) not in ABUNDANCES: continue

            # partition functions for T and Tref
            SigmaT = partitionFunction(MoleculeNumberDB, IsoNumberDB, T)
            SigmaTref = partitionFunction(MoleculeNumberDB, IsoNumberDB, Tref)

            # get all environment dependences from voigt parameters

            #   intensity
            if 'sw' in CustomEnvDependences:
                LineIntensity = CustomEnvDependences['sw']
            else:
                LineIntensity = EnvironmentDependency_Intensity(LineIntensityDB, T, Tref, SigmaT, SigmaTref,
                                                                LowerStateEnergyDB, LineCenterDB)

            #   FILTER by LineIntensity: compare it with IntencityThreshold
            if LineIntensity < IntensityThreshold: continue

            #   doppler broadening coefficient (GammaD)
            cMassMol = 1.66053873e-27  # hapi
            m = molecularMass(MoleculeNumberDB, IsoNumberDB) * cMassMol * 1000
            GammaD = sqrt(2 * cBolts * T * log(2) / m / cc ** 2) * LineCenterDB

            #   pressure broadening coefficient
            Gamma0 = 0.
            Shift0 = 0.
            for species in Diluent:
                species_lower = species  # species_lower = species.lower() # CHANGED RJH 23MAR18

                abun = Diluent[species]

                Gamma0DB, gamma_name = get_DB_value_robust(TableName, 'gamma', species_lower, RowID, subst_species='air', missing_substitute=0.)
                # gamma_name = 'gamma_' + species_lower

                TempRatioPowerDB, n_name = get_DB_value_robust(TableName, 'n', species, RowID, subst_species='air', missing_substitute=0.)
                # TempRatioPowerDB = TempRatioPowerDB * n_correction_factor

                # Add to the final Gamma0
                Gamma0 += abun * CustomEnvDependences.get(gamma_name,  # default ->
                                                          EnvironmentDependency_Gamma0(Gamma0DB, T, Tref, p, pref,
                                                                                       TempRatioPowerDB))

                Shift0DB, delta_name = get_DB_value_robust(TableName, 'delta', species_lower, RowID, subst_species='air', missing_substitute=0.)
                # delta_name = 'delta_' + species_lower
                deltap, deltap_name = get_DB_value_robust(TableName, 'deltap', species_lower, RowID, subst_species='air', missing_substitute=0.)

                Shift0 += abun * CustomEnvDependences.get(delta_name,  # default ->
                                                          ((Shift0DB + deltap * (T - Tref)) * p / pref))

            #   get final wing of the line according to Gamma0, OmegaWingHW and OmegaWing
            OmegaWingF = max(OmegaWing, OmegaWingHW * Gamma0, OmegaWingHW * GammaD)
            SgmD = GammaD / np.sqrt(2 * np.log(2))
            abu_rat = ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)] / NATURAL_ABUNDANCES[(MoleculeNumberDB, IsoNumberDB)]

            # if the value of the center is smaller than a certain value, if will skip
            if abstol_ctr > fctr_tol * LineIntensity * abu_rat * PROFILE_VOIGT_w(0., SgmD, Gamma0, 0.):
                skip_list.append(rowid_i)
                continue

            # if the value at the end of the wing is greater than the tolerance, then expand the wing
            # OmegaWingFF = OmegaWingF
            while abstol_wing < fctr_tol * LineIntensity * abu_rat * PROFILE_VOIGT_w(0., SgmD, Gamma0, OmegaWingF):
                wing_expand_list[rowid_i] += 1
                OmegaWingF *= 10
            # print(RowID,wing_expand_list[RowID])

            # no speed up with this computations...
            # BoundIndexLower = int(round( min(max((LineCenterDB-OmegaWingF-OmegaRange[0])/OmegaStep,0),number_of_points-1)))
            # BoundIndexUpper = int(round( min(max((LineCenterDB+OmegaWingF-OmegaRange[0])/OmegaStep,0),number_of_points-1)))

            # print(BoundIndexLower,BoundIndexUpper)

            #
            BoundIndexLower = bisect(Omegas, LineCenterDB - OmegaWingF)
            BoundIndexUpper = bisect(Omegas, LineCenterDB + OmegaWingF)

            Gamma0 = Gamma0 * gamma0_cor

            Gamma0_list[rowid_i] = Gamma0
            SgmD_List[rowid_i] = SgmD
            facList[rowid_i] = factor * abu_rat * LineIntensity

            if not (BoundIndexUpper == 0 or BoundIndexLower == number_of_points - 1):
                # print(RowID, BoundIndexLower,BoundIndexUpper,LineCenterDB,'LineIntensity',fctr_tol*LineIntensity,'GammaD',GammaD,'Gamma0',Gamma0,'Shift0',Shift0)
                # t = time.time()
                if shape_profile == 'voigt':
                    lineshape_vals = PROFILE_VOIGT_w(LineCenterDB + Shift0, SgmD, Gamma0,
                                                    Omegas[BoundIndexLower:BoundIndexUpper])
                # elif shape_profile == 'rautian':
                #     # beta_p = beta*Gamma0/(SgmD*np.sqrt(2))
                #     lineshape_vals = PROFILE_RAUTIAN_cy(LineCenterDB + Shift0, SgmD, Gamma0,
                #                                       Omegas[BoundIndexLower:BoundIndexUpper], beta)
                # elif shape_profile == 'rautian2':
                #     beta_p = beta*(p/pref)
                #     lineshape_vals = PROFILE_RAUTIAN_cy2(LineCenterDB + Shift0, SgmD, Gamma0,
                #                                       Omegas[BoundIndexLower:BoundIndexUpper], beta_p)
                elif shape_profile == "doppler":
                    lineshape_vals = PROFILE_DOPPLER(LineCenterDB + Shift0, GammaD,
                                                     Omegas[BoundIndexLower:BoundIndexUpper])
                elif shape_profile == "Lorentzian":
                    lineshape_vals = PROFILE_LORENTZ(LineCenterDB + Shift0, Gamma0, 0., 
                                                     Omegas[BoundIndexLower:BoundIndexUpper])
                # lineshape_vals = voigt_cy.VOIGT_PROFILE(LineCenterDB+Shift0,GammaD/np.sqrt(2*np.log(2)),Gamma0,Omegas[BoundIndexLower:BoundIndexUpper])
                # elapsed = time.time() - t
                Xsect[BoundIndexLower:BoundIndexUpper] += (factor * abu_rat * LineIntensity) * lineshape_vals

            # print(elapsed,'[sec]')

    if File: save_to_file(File, Format, Omegas, Xsect)
    return Omegas, Xsect, Gamma0_list, SgmD_List, skip_list, wing_expand_list, facList


def extend_vert_profile(vert_profile):
    alt = vert_profile['Altitude']
    alt_new = []
    for i in range(len(alt) - 1):
        alt_new_tmp = np.linspace(alt[i], alt[i + 1], 5, endpoint=False)
        alt_new.extend(alt_new_tmp)
    alt_new.append(alt[-1])
    alt_new = np.array(alt_new)

    vert_prof_extended = np.array([])
    dtypes_new = []
    names = vert_profile.dtype.names
    dtypes = vert_profile.dtype
    for i in range(len(names)):
        name = names[i]
        dtype = dtypes[i]
        if name == 'Altitude':
            y = alt_new
        else:
            y = np.interp(alt_new, vert_profile['Altitude'], vert_profile[name])
        if i == 0:
            vert_prof_extended = y.reshape(len(alt_new), 1)
            dtypes_new = [(name, '<f8')]
        else:
            vert_prof_extended = np.hstack((vert_prof_extended, y.reshape(len(alt_new), 1)))
            dtypes_new.append((name, '<f8'))

    vv = []
    for i in range(vert_prof_extended.shape[0]):
        v = tuple(vert_prof_extended[i, :])
        vv.append(v)

    vert_prof_extended = np.array(vv, dtype=dtypes_new)

    return vert_prof_extended
