from math import pow, log10

STANDARD_PRESSURE = 14.7

def l_b(u_m, i_d, output=False):
    l_b = 1.071 - 0.2218 * (pow(u_m, 2) / i_d)

    l_b = 0.13 if l_b < 0.13 else l_b

    if (output):
        print(l_b)

    return l_b

def hagedorn_brown_number(u_sl, u_sg, sigma, rho_l, i_d, mu_l, output=False):
    n_vl = 1.938 * u_sl * pow(rho_l / sigma, 0.25)
    n_vg = 1.938 * u_sg * pow(rho_l / sigma, 0.25)
    n_d = 120.872 * i_d * pow(rho_l / sigma, 0.5)
    n_l = 0.1572 * mu_l * pow(1 / (rho_l * pow(sigma, 3)), 0.25)

    result = {
        "n_vl": n_vl,
        "n_vg": n_vg,
        "n_d": n_d,
        "n_l": n_l,
    }

    if (output):
        print(result)

    return result

def first_graph(n_L):
    """
    Relation of graph between c_NL and n_L
    """

    X = log10(n_L + 3)

    Y1 = -2.6951 + 0.15841 * X - 0.551 * pow(X, 2)
    Y2 = 0.54785 * pow(X, 3) - 0.12195 * pow(X, 4)
    Y = Y1 + Y2

    c_NL = pow(10, Y)

    return c_NL


def second_graph(n_vl, n_vg, p, c_NL, n_D):
    A = n_vl / pow(n_vg, 0.575)
    B = pow(p / STANDARD_PRESSURE, 0.1)
    C = c_NL / n_D

    X = A * B * C

    log_X = log10(X) + 6

    Y1 = -0.10307 + 0.61777 * log_X
    Y2 = -0.63295 * pow(log_X, 2) + 0.29598 * pow(log_X, 3)
    Y3 = -0.0401 * pow(log_X, 4)

    Y = Y1 + Y2 + Y3

    yl_per_phi = Y

    return yl_per_phi

def third_graph(n_vg, n_L, n_D):
    """
    Relation of third graph
    phi versus A function
    """

    A = n_vg * pow(n_L, 0.38) / pow(n_D, 2.14)

    if (A <= 0.01):
        phi = 1
    
    else:
        phi1 = 0.91163 - 4.82176 * A + 1_232.25 * pow(A, 2)
        phi2 = -22_253.6 * pow(A, 3) + 116_174.3 * pow(A, 4)

        phi = phi1 = phi2
    
    return phi

def fanning_friction_factor(n_re, roughness):
    A = roughness / 3.7605
    B = -5.0452 / n_re
    
    C = pow(roughness, 1.1098) / 2.8257
    D = pow(7.149 / n_re, 0.8981)
    E = log10(C + D)

    F = A + B * E

    f_F = pow(1 / (4 * log10(F)), 2)

    return f_F