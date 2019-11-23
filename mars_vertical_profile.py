#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 17:29:36 2018

@author: yukiitoh
"""
import numpy as np

def pascal2atm(p_pascal):
    p_atm = p_pascal/(1.01325*10**5)
    return p_atm


vert_profile_61C4 = np.array(
   [(0.00000e+00, 2.56986e+02, 1.05284e+02, 2.40943e-03, 9.64510e-01, 1.68734e-02, 1.65242e-02, 6.64719e-04, 1.55559e-04, 2.98487e-09, 7.18400e-10, 1.25432e-03, 1.02501e-11, 1.46003e-05),
    (7.35294e+03, 1.97821e+02, 5.28003e+01, 1.40903e-03, 9.64971e-01, 1.66578e-02, 1.63132e-02, 6.55715e-04, 1.46139e-04, 1.16041e-09, 3.28525e-09, 1.23806e-03, 3.01215e-10, 1.44159e-05),
    (1.47059e+04, 1.79775e+02, 2.47562e+01, 7.29854e-04, 9.65643e-01, 1.63402e-02, 1.60022e-02, 6.42556e-04, 1.34026e-04, 1.00080e-09, 9.75659e-09, 1.21415e-03, 1.71713e-09, 1.41693e-05),
    (2.20588e+04, 1.72587e+02, 1.09962e+01, 3.37550e-04, 9.66061e-01, 1.61577e-02, 1.58233e-02, 6.33077e-04, 1.08125e-04, 1.14487e-09, 5.15612e-08, 1.19960e-03, 1.25666e-08, 1.41537e-05),
    (2.94118e+04, 1.66711e+02, 4.80277e+00, 1.51118e-04, 9.66361e-01, 1.60320e-02, 1.57001e-02, 6.26865e-04, 7.16881e-05, 2.06016e-09, 8.26095e-07, 1.18968e-03, 1.06122e-07, 1.43577e-05),
    (3.67647e+04, 1.58873e+02, 2.03047e+00, 6.70307e-05, 9.66355e-01, 1.60432e-02, 1.57108e-02, 6.40882e-04, 2.90609e-05, 3.82951e-09, 5.12933e-06, 1.19677e-03, 4.65390e-07, 1.48967e-05),
    (4.41176e+04, 1.58015e+02, 8.41618e-01, 2.80078e-05, 9.65947e-01, 1.61829e-02, 1.58460e-02, 7.20964e-04, 2.79428e-05, 1.87872e-09, 2.32986e-05, 1.23037e-03, 2.73122e-06, 1.58516e-05),
    (5.14706e+04, 1.53171e+02, 3.45034e-01, 1.18233e-05, 9.65346e-01, 1.63464e-02, 1.60036e-02, 8.81704e-04, 2.81519e-05, 1.54783e-09, 1.01305e-04, 1.26623e-03, 7.16775e-06, 1.65186e-05),
    (5.88235e+04, 1.44478e+02, 1.36303e-01, 4.94108e-06, 9.64785e-01, 1.64690e-02, 1.61217e-02, 1.05163e-03, 2.68477e-05, 1.21373e-09, 2.29148e-04, 1.28751e-03, 1.15172e-05, 1.69471e-05),
    (6.61765e+04, 1.36244e+02, 5.13153e-02, 1.97754e-06, 9.64226e-01, 1.65811e-02, 1.62299e-02, 1.22213e-03, 2.57494e-05, 6.50389e-10, 3.84955e-04, 1.30045e-03, 1.50206e-05, 1.72954e-05),
    (7.35294e+04, 1.32164e+02, 1.84798e-02, 7.33836e-07, 9.63434e-01, 1.67300e-02, 1.63734e-02, 1.47293e-03, 2.47248e-05, 2.14705e-10, 6.27526e-04, 1.31236e-03, 1.82849e-05, 1.76565e-05),
    (8.08824e+04, 1.29078e+02, 6.52559e-03, 2.65271e-07, 9.62829e-01, 1.67015e-02, 1.63447e-02, 1.79707e-03, 2.38964e-05, 5.11414e-11, 9.75901e-04, 1.30052e-03, 2.18056e-05, 1.78494e-05),
    (8.82353e+04, 1.17742e+02, 2.17928e-03, 9.65455e-08, 9.61310e-01, 1.68066e-02, 1.64470e-02, 2.46075e-03, 2.22977e-05, 1.49995e-11, 1.65358e-03, 1.29455e-03, 2.51288e-05, 1.80107e-05),
    (9.55882e+04, 9.90449e+01, 6.30687e-04, 3.31957e-08, 9.50020e-01, 2.09129e-02, 2.02272e-02, 4.27755e-03, 2.50255e-05, 6.21147e-12, 3.15826e-03, 1.60023e-03, 3.08995e-05, 1.99169e-05),
    (1.02941e+05, 1.07493e+02, 1.66017e-04, 8.09489e-09, 9.11509e-01, 3.52187e-02, 3.29216e-02, 1.05853e-02, 3.64303e-05, 1.75655e-12, 7.04603e-03, 3.12255e-03, 4.25809e-05, 2.41604e-05),
    (1.10294e+05, 1.46923e+02, 5.94658e-05, 2.13684e-09, 8.84502e-01, 3.92513e-02, 4.05757e-02, 1.88568e-02, 4.63680e-05, 2.74812e-13, 1.23885e-02, 4.70386e-03, 6.66132e-05, 3.43425e-05),
    (1.17647e+05, 1.75896e+02, 2.65994e-05, 7.80613e-10, 8.63457e-01, 3.91333e-02, 4.70319e-02, 2.60265e-02, 5.91569e-05, 6.15810e-14, 1.86288e-02, 5.85160e-03, 1.08480e-04, 5.44519e-05),
    (1.25000e+05, 1.92717e+02, 1.33308e-05, 3.50925e-10, 8.39559e-01, 3.88661e-02, 5.49359e-02, 3.29302e-02, 7.80306e-05, 2.22199e-14, 2.67791e-02, 6.85449e-03, 1.83341e-04, 9.17223e-05),
    (1.32353e+05, 2.02045e+02, 7.08473e-06, 1.74882e-10, 8.10416e-01, 3.87007e-02, 6.45444e-02, 4.01846e-02, 1.04761e-04, 1.08808e-14, 3.78919e-02, 7.83411e-03, 3.16062e-04, 1.57992e-04),
    (1.39706e+05, 2.07222e+02, 3.90390e-06, 9.21111e-11, 7.75076e-01, 3.84271e-02, 7.56827e-02, 4.80102e-02, 1.41005e-04, 6.46596e-15, 5.29930e-02, 8.81248e-03, 5.44914e-04, 2.71389e-04),
    (1.47059e+05, 2.10090e+02, 2.21026e-06, 5.02042e-11, 7.33070e-01, 3.78866e-02, 8.80395e-02, 5.63838e-02, 1.88103e-04, 4.37624e-15, 7.29630e-02, 9.77393e-03, 9.21820e-04, 4.56112e-04),
    (1.54412e+05, 2.11664e+02, 1.28055e-06, 2.79977e-11, 6.84372e-01, 3.69885e-02, 1.01247e-01, 6.51651e-02, 2.46723e-04, 3.22023e-15, 9.83549e-02, 1.06924e-02, 1.49906e-03, 7.36130e-04),
    (1.61765e+05, 2.12417e+02, 7.64504e-07, 1.60886e-11, 6.26891e-01, 3.54936e-02, 1.14773e-01, 7.40472e-02, 3.21900e-04, 2.54717e-15, 1.32064e-01, 1.14833e-02, 2.48593e-03, 1.20778e-03),
    (1.69118e+05, 2.12850e+02, 4.69948e-07, 9.46899e-12, 5.63148e-01, 3.34174e-02, 1.27252e-01, 8.22022e-02, 4.10860e-04, 2.04299e-15, 1.73532e-01, 1.20601e-02, 4.07627e-03, 1.95498e-03),
    (1.76471e+05, 2.13122e+02, 2.95622e-07, 5.64209e-12, 4.96152e-01, 3.09399e-02, 1.38300e-01, 8.94037e-02, 5.08104e-04, 1.62931e-15, 2.20231e-01, 1.24404e-02, 6.25611e-03, 2.96635e-03),
    (1.83824e+05, 2.13259e+02, 1.92460e-07, 3.46689e-12, 4.26271e-01, 2.79320e-02, 1.46519e-01, 9.47466e-02, 6.14671e-04, 1.30080e-15, 2.73682e-01, 1.24960e-02, 9.44412e-03, 4.42488e-03),
    (1.91176e+05, 2.13343e+02, 1.29822e-07, 2.19795e-12, 3.59024e-01, 2.46687e-02, 1.50579e-01, 9.73854e-02, 7.20025e-04, 1.02858e-15, 3.29617e-01, 1.22026e-02, 1.40833e-02, 6.50405e-03),
    (1.98529e+05, 2.13399e+02, 8.93098e-08, 1.40314e-12, 2.93238e-01, 2.12639e-02, 1.52333e-01, 9.85249e-02, 8.24706e-04, 7.87572e-16, 3.86927e-01, 1.17158e-02, 1.95271e-02, 8.92735e-03),
    (2.05882e+05, 2.13431e+02, 6.43576e-08, 9.47486e-13, 2.38027e-01, 1.80235e-02, 1.48224e-01, 9.58656e-02, 9.09790e-04, 6.10586e-16, 4.38949e-01, 1.08706e-02, 2.76498e-02, 1.24413e-02),
    (2.13235e+05, 2.13455e+02, 4.70170e-08, 6.40567e-13, 1.85941e-01, 1.48318e-02, 1.42382e-01, 9.20840e-02, 9.89083e-04, 4.52515e-16, 4.89409e-01, 9.91956e-03, 3.65640e-02, 1.62774e-02),
    (2.20588e+05, 2.13476e+02, 3.53820e-08, 4.48360e-13, 1.45111e-01, 1.21481e-02, 1.34050e-01, 8.66934e-02, 1.04558e-03, 3.36218e-16, 5.30489e-01, 8.91753e-03, 4.81143e-02, 2.11295e-02),
    (2.27941e+05, 2.13495e+02, 2.73010e-08, 3.26094e-13, 1.13782e-01, 9.89322e-03, 1.23616e-01, 7.99446e-02, 1.08285e-03, 2.55185e-16, 5.63653e-01, 7.87248e-03, 6.18898e-02, 2.68392e-02),
    (2.35294e+05, 2.13513e+02, 2.10655e-08, 2.31753e-13, 8.24532e-02, 7.63839e-03, 1.13181e-01, 7.31958e-02, 1.12011e-03, 1.74152e-16, 5.96816e-01, 6.82742e-03, 7.56654e-02, 3.25490e-02),
    (2.42647e+05, 2.13516e+02, 1.63670e-08, 3.31881e-13, 5.06358e-02, 5.08506e-03, 9.50946e-02, 6.14991e-02, 1.13389e-03, 9.86283e-17, 6.26794e-01, 5.31462e-03, 1.04023e-01, 4.39642e-02),
    (2.50000e+05, 2.13516e+02, 1.30380e-08, 2.44444e-13, 3.50673e-02, 3.71574e-03, 8.16240e-02, 5.27875e-02, 1.11299e-03, 6.47353e-17, 6.31976e-01, 4.32345e-03, 1.28262e-01, 5.34860e-02)],
    dtype = [('Altitude','float'),('Temperature','float'),('Pressure','float'),('Density','float'),
    ('CO2_vol_mix_rat','float'),('Ar_vol_mix_rat','float'),('N2_vol_mix_rat','float'),('CO_vol_mix_rat','float'),('H2O_vol_mix_rat','float'),('O3_vol_mix_rat','float'),
    ('O_vol_mix_rat','float'),('O2_vol_mix_rat','float'),('H_vol_mix_rat','float'),('H2_vol_mix_rat','float')])
P_atm_61C4 = pascal2atm(vert_profile_61C4['Pressure'])

vert_profile_61C4_btm = np.array([
    (0.00000e+00, 6.01933e+02, 2.65548e+02, 9.64293e-01, 6.68940e-04, 1.41776e-04),
    (7.35294e+03, 3.31305e+02, 2.31962e+02, 9.64342e-01, 6.68045e-04, 1.43454e-04),
    (1.47059e+04, 1.75956e+02, 2.20814e+02, 9.64463e-01, 6.65883e-04, 1.36722e-04),
    (2.20588e+04, 9.07147e+01, 2.08574e+02, 9.64922e-01, 6.56963e-04, 1.43454e-04),
    (2.94118e+04, 4.46006e+01, 1.90503e+02, 9.65685e-01, 6.41866e-04, 1.53239e-04),
    (3.67647e+04, 2.05021e+01, 1.76309e+02, 9.65971e-01, 6.35561e-04, 1.20998e-04),
    (4.41176e+04, 9.11768e+00, 1.71221e+02, 9.66165e-01, 6.29968e-04, 1.01298e-04),
    (5.14706e+04, 3.95336e+00, 1.64154e+02, 9.66384e-01, 6.27706e-04, 5.83864e-05),
    (5.88235e+04, 1.66149e+00, 1.58129e+02, 9.66317e-01, 6.48253e-04, 2.82376e-05),
    (6.61765e+04, 6.90412e-01, 1.58398e+02, 9.65829e-01, 7.51919e-04, 2.81796e-05),
    (7.35294e+04, 2.83622e-01, 1.52494e+02, 9.65252e-01, 9.04935e-04, 2.79219e-05),
    (8.08824e+04, 1.11534e-01, 1.42553e+02, 9.64741e-01, 1.05848e-03, 2.68271e-05),
    (8.82353e+04, 4.13790e-02, 1.34764e+02, 9.64146e-01, 1.24256e-03, 2.57729e-05),
    (9.55882e+04, 1.48208e-02, 1.31758e+02, 9.63277e-01, 1.53190e-03, 2.45610e-05),
    (1.02941e+05, 5.21351e-03, 1.28147e+02, 9.62704e-01, 1.90500e-03, 2.35769e-05),
    (1.10294e+05, 1.72120e-03, 1.15004e+02, 9.60728e-01, 2.67188e-03, 2.20267e-05),
    (1.17647e+05, 4.79134e-04, 9.67730e+01, 9.44034e-01, 5.11492e-03, 2.67020e-05),
    (1.25000e+05, 1.31343e-04, 1.13907e+02, 9.06688e-01, 1.21573e-02, 3.73729e-05),
    (1.32353e+05, 4.81570e-05, 1.50720e+02, 8.81469e-01, 2.02123e-02, 4.68125e-05),
    (1.39706e+05, 2.17821e-05, 1.75755e+02, 8.60248e-01, 2.73050e-02, 5.98932e-05),
    (1.47059e+05, 1.08798e-05, 1.89602e+02, 8.35292e-01, 3.43489e-02, 7.94912e-05),
    (1.54412e+05, 5.72583e-06, 1.96863e+02, 8.04432e-01, 4.19164e-02, 1.07698e-04),
    (1.61765e+05, 3.11770e-06, 2.00664e+02, 7.66333e-01, 5.02284e-02, 1.46942e-04),
    (1.69118e+05, 1.74390e-06, 2.02732e+02, 7.20445e-01, 5.91717e-02, 1.99006e-04),
    (1.76471e+05, 1.00066e-06, 2.03864e+02, 6.66813e-01, 6.84219e-02, 2.64911e-04),
    (1.83824e+05, 5.88688e-07, 2.04481e+02, 6.06268e-01, 7.75767e-02, 3.44187e-04),
    (1.91176e+05, 3.56562e-07, 2.04792e+02, 5.38895e-01, 8.61330e-02, 4.37931e-04),
    (1.98529e+05, 2.24989e-07, 2.04947e+02, 4.67142e-01, 9.26095e-02, 5.44418e-04),
    (2.05882e+05, 1.45480e-07, 2.05052e+02, 3.95573e-01, 9.74107e-02, 6.53322e-04),
    (2.13235e+05, 9.77250e-08, 2.05105e+02, 3.25914e-01, 9.93472e-02, 7.62992e-04),
    (2.20588e+05, 6.79467e-08, 2.05142e+02, 2.63446e-01, 9.83943e-02, 8.60537e-04),
    (2.27941e+05, 4.85422e-08, 2.05165e+02, 2.06645e-01, 9.51644e-02, 9.48527e-04),
    (2.35294e+05, 3.56601e-08, 2.05187e+02, 1.60009e-01, 9.01934e-02, 1.01566e-03),
    (2.42647e+05, 2.70077e-08, 2.05205e+02, 1.24483e-01, 8.33192e-02, 1.05998e-03),
    (2.50000e+05, 2.04547e-08, 2.05224e+02, 8.89580e-02, 7.64449e-02, 1.10430e-03)
    ],
    dtype = [('Altitude','float'),('Pressure','float'),('Temperature','float'),
             ('CO2_vol_mix_rat','float'),('CO_vol_mix_rat','float'),
             ('H2O_vol_mix_rat','float')
            ]
    )
P_atm_61C4_btm = pascal2atm(vert_profile_61C4_btm['Pressure'])

vert_profile_40FF = np.array([
    (0.00000e+00, 6.21048e+02, 2.69101e+02, 9.60146e-01, 7.35076e-04, 2.76836e-04, 1.85149e-02, 1.40210e-03, 1.89063e-02, 5.52141e-09, 2.01460e-10, 1.59496e-05, 1.20607e-12),
    (7.35294e+03, 3.36845e+02, 2.26209e+02, 9.59816e-01, 7.43980e-04, 3.04902e-04, 1.86531e-02, 1.41397e-03, 1.90474e-02, 3.40265e-09, 2.88484e-10, 1.60749e-05, 3.85194e-12),
    (1.47059e+04, 1.74681e+02, 2.09310e+02, 9.58879e-01, 7.67944e-04, 3.09075e-04, 1.90854e-02, 1.44965e-03, 1.94888e-02, 2.02541e-09, 4.93963e-10, 1.64618e-05, 1.60854e-11),
    (2.20588e+04, 8.55645e+01, 1.91191e+02, 9.58602e-01, 7.76111e-04, 1.80737e-04, 1.92745e-02, 1.46421e-03, 1.96819e-02, 1.69984e-09, 1.28899e-09, 1.66840e-05, 6.86200e-11),
    (2.94118e+04, 3.96924e+01, 1.79387e+02, 9.58529e-01, 7.78719e-04, 9.60233e-05, 1.93486e-02, 1.46974e-03, 1.97576e-02, 1.91712e-09, 5.44275e-09, 1.68056e-05, 4.16086e-10),
    (3.67647e+04, 1.78386e+01, 1.74448e+02, 9.58522e-01, 7.78625e-04, 8.39389e-05, 1.93573e-02, 1.47041e-03, 1.97665e-02, 1.68353e-09, 2.08560e-08, 1.68758e-05, 3.00439e-09),
    (4.41176e+04, 7.88279e+00, 1.70455e+02, 9.58434e-01, 7.81590e-04, 6.70408e-05, 1.94052e-02, 1.47542e-03, 1.98154e-02, 2.51699e-09, 2.32530e-07, 1.70822e-05, 3.90702e-08),
    (5.14706e+04, 3.41268e+00, 1.63639e+02, 9.58329e-01, 7.94092e-04, 3.03169e-05, 1.94623e-02, 1.48562e-03, 1.98741e-02, 7.28107e-09, 2.75593e-06, 1.73255e-05, 1.62173e-07),
    (5.88235e+04, 1.43499e+00, 1.58951e+02, 9.58110e-01, 8.36468e-04, 2.49986e-05, 1.95362e-02, 1.50669e-03, 1.99509e-02, 4.75406e-09, 1.25848e-05, 1.76744e-05, 9.01963e-07),
    (6.61765e+04, 6.03859e-01, 1.62174e+02, 9.57651e-01, 9.67530e-04, 2.19364e-05, 1.96534e-02, 1.53575e-03, 2.00743e-02, 2.96398e-09, 7.09984e-05, 1.79454e-05, 3.53279e-06),
    (7.35294e+04, 2.56723e-01, 1.60746e+02, 9.57206e-01, 1.14338e-03, 1.97222e-05, 1.97146e-02, 1.54868e-03, 2.01404e-02, 2.10361e-09, 1.99659e-04, 1.80255e-05, 6.69890e-06),
    (8.08824e+04, 1.07056e-01, 1.53889e+02, 9.56575e-01, 1.43876e-03, 1.81258e-05, 1.97419e-02, 1.55030e-03, 2.01743e-02, 1.84235e-09, 4.70827e-04, 1.80239e-05, 9.60598e-06),
    (8.82353e+04, 4.23404e-02, 1.42307e+02, 9.55665e-01, 1.77244e-03, 1.68294e-05, 1.98679e-02, 1.55652e-03, 2.03114e-02, 1.03644e-09, 7.84340e-04, 1.80768e-05, 1.20419e-05),
    (9.55882e+04, 1.57318e-02, 1.34372e+02, 9.54628e-01, 2.09871e-03, 1.56245e-05, 2.00565e-02, 1.56361e-03, 2.05113e-02, 3.36128e-10, 1.10763e-03, 1.81781e-05, 1.43886e-05),
    (1.02941e+05, 5.59833e-03, 1.28383e+02, 9.53540e-01, 2.55469e-03, 1.37669e-05, 2.01351e-02, 1.55559e-03, 2.06004e-02, 7.39030e-11, 1.58287e-03, 1.81555e-05, 1.66593e-05),
    (1.10294e+05, 1.84572e-03, 1.15077e+02, 9.50249e-01, 3.58085e-03, 1.25431e-05, 2.07471e-02, 1.57893e-03, 2.12654e-02, 2.26419e-11, 2.58178e-03, 1.84549e-05, 1.93855e-05),
    (1.17647e+05, 5.18964e-04, 9.79039e+01, 9.30194e-01, 6.82181e-03, 1.43157e-05, 2.74003e-02, 2.12236e-03, 2.87457e-02, 9.31841e-12, 5.00024e-03, 2.04371e-05, 2.39982e-05),
    (1.25000e+05, 1.44670e-04, 1.15334e+02, 8.88070e-01, 1.54712e-02, 1.87495e-05, 4.06315e-02, 3.82136e-03, 4.22030e-02, 1.99184e-12, 1.01243e-02, 2.47431e-05, 3.33158e-05),
    (1.32353e+05, 5.41350e-05, 1.51028e+02, 8.59558e-01, 2.47172e-02, 2.29547e-05, 4.92498e-02, 5.34307e-03, 4.50373e-02, 2.97034e-13, 1.62455e-02, 3.47565e-05, 5.07310e-05),
    (1.39706e+05, 2.46139e-05, 1.75437e+02, 8.34643e-01, 3.28628e-02, 2.89060e-05, 5.75992e-02, 6.50415e-03, 4.49873e-02, 6.83795e-14, 2.34655e-02, 5.41685e-05, 8.11013e-05),
    (1.47059e+05, 1.23633e-05, 1.89185e+02, 8.05873e-01, 4.09523e-02, 3.78260e-05, 6.76975e-02, 7.55532e-03, 4.48367e-02, 2.56850e-14, 3.30655e-02, 9.02275e-05, 1.35999e-04),
    (1.54412e+05, 6.55068e-06, 1.96635e+02, 7.71076e-01, 4.95542e-02, 5.05593e-05, 7.96859e-02, 8.58322e-03, 4.45924e-02, 1.28720e-14, 4.63279e-02, 1.55055e-04, 2.34749e-04),
    (1.61765e+05, 3.59465e-06, 2.00691e+02, 7.29343e-01, 5.87921e-02, 6.79081e-05, 9.32708e-02, 9.58947e-03, 4.40172e-02, 7.70129e-15, 6.44862e-02, 2.67593e-04, 4.07865e-04),
    (1.69118e+05, 2.02896e-06, 2.02907e+02, 6.80415e-01, 6.85028e-02, 9.04794e-05, 1.07942e-01, 1.05413e-02, 4.29591e-02, 5.18508e-15, 8.85737e-02, 4.54274e-04, 6.98702e-04),
    (1.76471e+05, 1.17477e-06, 2.04111e+02, 6.24701e-01, 7.83649e-02, 1.18461e-04, 1.23062e-01, 1.13962e-02, 4.13501e-02, 3.76403e-15, 1.19147e-01, 7.44303e-04, 1.15628e-03),
    (1.83824e+05, 6.98455e-07, 2.04731e+02, 5.62357e-01, 8.80008e-02, 1.52281e-04, 1.37973e-01, 1.21004e-02, 3.91145e-02, 2.87341e-15, 1.57064e-01, 1.17720e-03, 1.84696e-03),
    (1.91176e+05, 4.31643e-07, 2.05025e+02, 4.94127e-01, 9.60858e-02, 1.93074e-04, 1.50569e-01, 1.24768e-02, 3.60514e-02, 2.26992e-15, 2.04873e-01, 1.91912e-03, 3.05475e-03),
    (1.98529e+05, 2.73717e-07, 2.05220e+02, 4.25086e-01, 1.02465e-01, 2.36562e-04, 1.60539e-01, 1.25995e-02, 3.25976e-02, 1.76825e-15, 2.57546e-01, 2.94122e-03, 4.74433e-03),
    (2.05882e+05, 1.78299e-07, 2.05331e+02, 3.55693e-01, 1.06959e-01, 2.82552e-04, 1.67599e-01, 1.24520e-02, 2.87622e-02, 1.35855e-15, 3.15016e-01, 4.26727e-03, 6.95806e-03),
    (2.13235e+05, 1.21486e-07, 2.05389e+02, 2.92728e-01, 1.07417e-01, 3.26898e-04, 1.68366e-01, 1.18668e-02, 2.47749e-02, 1.05179e-15, 3.74542e-01, 6.31449e-03, 1.04627e-02),
    (2.20588e+05, 8.35907e-08, 2.05436e+02, 2.31164e-01, 1.06997e-01, 3.70885e-04, 1.67763e-01, 1.11862e-02, 2.07545e-02, 7.67458e-16, 4.34516e-01, 8.51876e-03, 1.42485e-02),
    (2.27941e+05, 6.06342e-08, 2.05460e+02, 1.84692e-01, 1.01903e-01, 4.04373e-04, 1.59822e-01, 1.01559e-02, 1.72452e-02, 5.87957e-16, 4.87356e-01, 1.19854e-02, 2.03985e-02),
    (2.35294e+05, 4.41075e-08, 2.05482e+02, 1.39034e-01, 9.65577e-02, 4.37295e-04, 1.51486e-01, 9.10666e-03, 1.37634e-02, 4.14110e-16, 5.39811e-01, 1.55201e-02, 2.66761e-02),
    (2.42647e+05, 3.33937e-08, 2.05502e+02, 1.07948e-01, 8.90235e-02, 4.54138e-04, 1.39693e-01, 8.03218e-03, 1.11404e-02, 3.07964e-16, 5.75346e-01, 2.02212e-02, 3.53375e-02),
    (2.50000e+05, 2.56229e-08, 2.05521e+02, 8.17423e-02, 8.07565e-02, 4.65598e-04, 1.26743e-01, 6.94922e-03, 8.80490e-03, 2.24489e-16, 6.05215e-01, 2.53130e-02, 4.47972e-02)
    ],
    dtype = [('Altitude','float'),('Pressure','float'),('Temperature','float'),
             ('CO2_vol_mix_rat','float'),('CO_vol_mix_rat','float'),
             ('H2O_vol_mix_rat','float'),('N2_vol_mix_rat','float'),
             ('O2_vol_mix_rat','float'),('Ar_vol_mix_rat','float'),
             ('O3_vol_mix_rat','float'),('O_vol_mix_rat','float'),
             ('H2_vol_mix_rat','float'),('H_vol_mix_rat','float'),
            ]
    )
P_atm_40FF = pascal2atm(vert_profile_40FF['Pressure'])