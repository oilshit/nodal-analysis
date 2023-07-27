from math import pow, pi

import matplotlib.pyplot as plt

import conversion as cvt
import calculation as clc

DAY_TO_SECONDS = 86400
STANDARD_PRESSURE = 14.7
STANDARD_TEMPERATURE = cvt.fahrenheit_to_rankine(60)

ipr_q = [
    23.5,
    1258,
    2492.6,
    3727.1,
    4961.7,
    6196.2,
    7430.8,
    8665.3,
    9899.9,
    11134.4,
    12369,
    13603.5,
    14838.1,
    16072.6,
    17307.2,
    18541.7,
    19776.3,
    21010.8,
    22245.4,
    23479.9
]

ipr_p = [
    3997.07,
    3843.15,
    3689.24,
    3535.32,
    3381.4,
    3227.48,
    3073.57,
    2919.65,
    2765.73,
    2611.81,
    2457.9,
    2301.9,
    2136.46,
    1958.94,
    1766.23,
    1553.6,
    1313.2,
    1030.08,
    667.43,
    21.43
]

tpr_p = []

# wellhead properties
p_wellhead = 250          # psia
t_wellhead = cvt.fahrenheit_to_rankine(159.4)

# tubing properties
i_d = cvt.inch_to_feet(3.992)
roughness = cvt.inch_to_feet(0.0018)

area = 0.25 * pi * pow(i_d, 2)

for q_oil in ipr_q:
    # rate properties
    # q_oil = 11489.7           # in stb/d
    q_gas = 4.596             # in MMSCFD

    # density
    # rho_oil = 52.829
    rho_oil = 52.829
    rho_gas = 0.89886

    # viscosity properties
    mu_oil = 3.3726
    mu_gas = 0.012433

    # Z-factor
    z = 0.96312

    # interfacial tension
    sigma = 14.9966

    # superficial velocity
    # u_sl = 8.982
    u_sl = q_oil / area * 5.615 / 86400
    # u_sg = 34.182
    u_sg = q_gas / area * z * t_wellhead / STANDARD_TEMPERATURE * STANDARD_PRESSURE / p_wellhead

    # mixture superficial velocity
    u_m = u_sg + u_sl

    lambda_g = u_sg / u_m
    l_b = clc.l_b(u_m, i_d, output=False)

    # Hagedorn-Brown constants
    n_vl, n_vg, n_d, n_l = [
        props[1] for props in dict.items(
            clc.hagedorn_brown_number(
            u_sl=u_sl,
            u_sg=u_sg,
            sigma=sigma,
            rho_l=rho_oil,
            i_d=i_d,
            mu_l=mu_oil,
            output=False)
        )
    ]

    print(n_vl, n_vg, n_d, n_l)

    # graph correlation
    c_nl = clc.first_graph(n_l)
    yl_per_phi = clc.second_graph(n_vl, n_vg, p_wellhead, c_nl, n_d)
    phi = clc.third_graph(n_vg, n_l, n_d)

    y_l = phi * yl_per_phi
    print("graph:", c_nl, yl_per_phi, phi, y_l)

    # mixture mass flow rate
    mt = area * (u_sl * rho_oil + u_sg * rho_gas) * DAY_TO_SECONDS
    # print("mt:", mt)

    # Reynolds number
    n_re = 2.2E-2 * mt / (i_d * pow(mu_oil, y_l) * pow(mu_gas, (1 - y_l)))
    # print("Reynolds Number =", n_re)

    # FRICTION FACTOR
    # ====================================
    f_f = clc.fanning_friction_factor(n_re, roughness)
    # print("Fanning Friction Factor =", f_f)
    # print()
    # ====================================

    # INSITU AVERAGE DENSITY
    # ====================================
    rho_avg = y_l * rho_oil + (1 - y_l) * rho_gas
    print("Insitu average density (lbm/ft^3) =", rho_avg)
    # print()
    # ====================================

    # INSITU AVERAGE DENSITY
    # ====================================
    dp_dz = (1 / 144) * (
        rho_avg + (
            f_f * pow(mt, 2) / (7.143E10 * pow(i_d, 5) * rho_avg)
        )
    )

    print("Pressure gradient at top of tubing (wellhead, psi/ft) =", dp_dz)
    # ====================================

    depth = 8000
    p_wf = p_wellhead + (dp_dz * depth)

    tpr_p.append(p_wf)

print(tpr_p)

# print(q_oil, p_wf)

plt.plot(ipr_q, ipr_p)
plt.plot(ipr_q, tpr_p)

plt.xlabel("Flow rate (stb/d)")
plt.ylabel("Pressure (psia)")

plt.grid()
plt.xlim(0, max(ipr_q) + 1000)
plt.ylim(0, max(ipr_p) + 300)

plt.show()