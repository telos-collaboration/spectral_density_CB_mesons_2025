import fitting
import bootstrap
import numpy as np
import itertools
import matplotlib.pyplot as plt
from scipy import linalg

import gvar as gv

marker = itertools.cycle(("o", "v", "*", "s", "p", "x"))
sperater = np.nditer(np.linspace(0, 0.5, 12))


def meson_mass(C_tmp, TI, TF):
    # load the ensamble info
    GLB_T = np.shape(C_tmp)[1]

    C_fold = bootstrap.fold_correlators(C_tmp)
    C_boot = bootstrap.Correlator_resample(C_fold)

    E_fit, E_fit_err, X2 = fitting.fit_cosh_std(C_boot, TI, TF, GLB_T)

    # E_fit , E_fit_err, X2 = fitting.X2_single_state_fit(C_boot, TI, TF+1 )

    # print(fitting.print_non_zero(E_fit[-1],  E_fit_err)+'--- Interval=['+str(TI)+','+str(TF)+'] Xsqr= '+ str(round(X2,2)))
    print()

    return E_fit, E_fit_err, round(X2, 2)


def meson_mass_boo(C_tmp, TI, TF):
    # load the ensamble info
    GLB_T = np.shape(C_tmp)[1]

    C_fold = bootstrap.fold_correlators(C_tmp)
    C_boot = bootstrap.Correlator_resample(C_fold)

    E_fit, E_fit_err, X2 = fitting.fit_cosh_booerr(C_boot, TI, TF, GLB_T)

    # E_fit , E_fit_err, X2 = fitting.X2_single_state_fit(C_boot, TI, TF+1 )

    # print(fitting.print_non_zero(E_fit[-1],  E_fit_err)+'--- Interval=['+str(TI)+','+str(TF)+'] Xsqr= '+ str(round(X2,2)))
    print()

    return E_fit, E_fit_err, round(X2, 2)


def meson_mass_X2(C_tmp, TI, TF):
    # load the ensamble info
    GLB_T = np.shape(C_tmp)[1]

    C_fold = bootstrap.fold_correlators(C_tmp)
    C_boot = bootstrap.Correlator_resample(C_fold)

    E_fit, E_fit_err, X2 = fitting.X2_single_state_fit(C_boot, TI, TF + 1)

    # print(fitting.print_non_zero(E_fit[-1],  E_fit_err)+'--- Interval=['+str(TI)+','+str(TF)+'] Xsqr= '+ str(round(X2,2)))
    print()

    return E_fit[-1], E_fit_err, round(X2, 2)


def bin_proj_Baryon(corr_e, corr_o):
    def flip_boundary(C, t):
        if t == 2:
            return C
        else:
            C[:, -(t - 2) :] = -C[:, -(t - 2) :]
            return C

    size = np.shape(corr_e)

    C_e = flip_boundary(corr_e, 0)  # source location at t=0
    C_o = flip_boundary(corr_o, 0)

    C_Ebin = np.zeros(shape=(size[0], size[1] + 1))
    C_Obin = np.zeros(shape=(size[0], size[1] + 1))

    C_etmp = np.zeros(shape=(size[0], size[1] + 1))
    C_etmp[:, 0 : size[1]] = C_e
    C_etmp[:, size[1]] = C_e[:, 0]

    C_otmp = np.zeros(shape=(size[0], size[1] + 1))
    C_otmp[:, 0 : size[1]] = C_o
    C_otmp[:, size[1]] = C_o[:, 0]

    C_Ebin += (C_etmp - np.flip(C_otmp, axis=1)) / 2
    C_Obin += (C_otmp - np.flip(C_etmp, axis=1)) / 2

    bootstrap.num_confs = np.shape(corr_e)[0]

    return bootstrap.Correlator_resample(C_Ebin), -bootstrap.Correlator_resample(C_Obin)


def baryon_mass(C_tmp, TI, TF):
    # load the ensamble info
    # GLB_T = np.shape(C_tmp)[1]
    # num_config = np.shape(C_tmp)[0]

    E_fit, E_fit_err, X2 = fitting.X2_single_exp_fit(C_tmp, TI, TF)
    print(
        fitting.print_non_zero(E_fit[-1], E_fit_err)
        + "--- Interval=["
        + str(TI)
        + ","
        + str(TF)
        + "] Xsqr= "
        + str(round(X2 / (TF - TI - 2), 2))
    )

    return E_fit, E_fit_err, round(X2 / (TF - TI - 2), 2)


def Analysis_lnC(C_resample, ti, tf, measurement):
    num_sample = np.shape(C_resample)[0]
    GLB_T = np.shape(C_resample)[1]

    Mass_channel_tmp = np.zeros(shape=(num_sample, GLB_T))
    T_dot = np.arange(ti, tf, 1, dtype=int)

    mass_tmp = []
    err_tmp = []
    for t in T_dot:
        for N in range(num_sample):
            c_t = C_resample[N, t]
            Mass_channel_tmp[N, t] = np.log(c_t)

        Mass_err = bootstrap.bootstrap_error(
            Mass_channel_tmp[0:-1, t], Mass_channel_tmp[-1, t]
        )

        mass_tmp.append(Mass_channel_tmp[-1, t])
        err_tmp.append(Mass_err)

    sprnext = next(sperater)
    plt.errorbar(
        T_dot + sprnext,
        mass_tmp,
        err_tmp,
        linestyle="",
        marker=next(marker),
        alpha=0.6,
        label=measurement,
    )

    return Mass_channel_tmp


def Analysis_Mass_eff_simple(C_resample, ti, tf, dt, measurement):
    num_sample = np.shape(C_resample)[0]
    GLB_T = np.shape(C_resample)[1]

    Mass_channel_tmp = np.zeros(shape=(num_sample, GLB_T))
    T_dot = np.arange(ti, tf, 1, dtype=int)

    mass_tmp = []
    err_tmp = []

    Mass_channel_tmp = -np.log(np.roll(C_resample, -1, axis=1) / C_resample)

    for t in T_dot:
        Mass_err = bootstrap.bootstrap_error(
            Mass_channel_tmp[0:-1, t], Mass_channel_tmp[-1, t]
        )
        mass_tmp.append(Mass_channel_tmp[-1, t])
        err_tmp.append(Mass_err)

    mass_tmp = np.array(mass_tmp)
    err_tmp = np.array(err_tmp)

    sprnext = next(sperater)
    mk = next(marker)

    select = (abs(err_tmp / mass_tmp)) < 0.4

    plt.errorbar(
        T_dot[select] + sprnext,
        mass_tmp[select],
        err_tmp[select],
        linestyle="",
        marker=mk,
        alpha=0.6,
        label=measurement,
    )

    plt.errorbar(
        T_dot[~select] + sprnext,
        mass_tmp[~select],
        err_tmp[~select],
        linestyle="",
        marker=mk,
        color=plt.gca().lines[-1].get_color(),
        alpha=0.1,
    )

    return Mass_channel_tmp

def Analysis_Mass_eff_simple2(C_resample, ti, tf, dt, measurement, n, E_string):
    num_sample = np.shape(C_resample)[0]
    GLB_T = np.shape(C_resample)[1]

    Mass_channel_tmp = np.zeros(shape=(num_sample, GLB_T))
    T_dot = np.arange(ti, tf, 1, dtype=int)

    mass_tmp = []
    err_tmp = []

    Mass_channel_tmp = -np.log(np.roll(C_resample, -1, axis=1) / C_resample)

    for t in T_dot:
        Mass_err = bootstrap.bootstrap_error(
            Mass_channel_tmp[0:-1, t], Mass_channel_tmp[-1, t]
        )
        mass_tmp.append(Mass_channel_tmp[-1, t])
        err_tmp.append(Mass_err)

    mass_tmp = np.array(mass_tmp)
    err_tmp = np.array(err_tmp)

    sprnext = next(sperater)
    mk = next(marker)

    select = (abs(err_tmp / mass_tmp)) < 0.10

    plt.errorbar(
        T_dot[select] + sprnext,
        mass_tmp[select],
        1.5*err_tmp[select],
        linestyle="",
        marker=mk,
        alpha=0.6,
        label=f"$aE_{n}$ fit: "+E_string,
    )

    plt.errorbar(
        T_dot[~select] + sprnext,
        mass_tmp[~select],
        1.5*err_tmp[~select],
        linestyle="",
        marker=mk,
        color=plt.gca().lines[-1].get_color(),
        alpha=0.1,
    )
    plt.xlabel('$t / a$', fontsize=15)
    plt.ylabel('$am_{\\rm{eff}}$', fontsize=15)

    return Mass_channel_tmp


def Analysis_Mass_eff_cosh(C_resample, ti, tf, measurement):
    num_sample = np.shape(C_resample)[0]
    GLB_T = np.shape(C_resample)[1]

    Mass_channel_tmp = np.zeros(shape=(num_sample, GLB_T))
    T_dot = np.arange(ti, tf, 1, dtype=int)

    mass_tmp = []
    err_tmp = []
    Mass_channel_tmp = np.arccosh(
        (np.roll(C_resample, 1, axis=1) + np.roll(C_resample, -1, axis=1))
        / C_resample
        / 2
    )
    for t in T_dot:
        Mass_err = bootstrap.bootstrap_error(
            Mass_channel_tmp[0:-1, t], Mass_channel_tmp[-1, t]
        )
        mass_tmp.append(Mass_channel_tmp[-1, t])
        err_tmp.append(Mass_err)

    mass_tmp = np.array(mass_tmp)
    err_tmp = np.array(err_tmp)

    sprnext = next(sperater)
    mk = next(marker)

    select = abs(err_tmp / mass_tmp) < 0.4

    plt.errorbar(
        T_dot[select] + sprnext,
        mass_tmp[select],
        1.5*err_tmp[select],
        linestyle="",
        elinewidth=2.0,
        capsize=2,
        marker=mk,
        alpha=0.6,
        label=measurement,
    )

    plt.errorbar(
        T_dot[~select] + sprnext,
        mass_tmp[~select],
        1.5*err_tmp[~select],
        linestyle="",
        elinewidth=2.0,
        capsize=2,
        marker=mk,
        color=plt.gca().lines[-1].get_color(),
        alpha=0.1,
    )
    return Mass_channel_tmp


def Analysis_Mass_eff_cosh2(C_resample, ti, tf, measurement, n, E_string):
    num_sample = np.shape(C_resample)[0]
    GLB_T = np.shape(C_resample)[1]

    Mass_channel_tmp = np.zeros(shape=(num_sample, GLB_T))
    T_dot = np.arange(ti, tf, 1, dtype=int)

    mass_tmp = []
    err_tmp = []
    Mass_channel_tmp = np.arccosh(
        (np.roll(C_resample, 1, axis=1) + np.roll(C_resample, -1, axis=1))
        / C_resample
        / 2
    )
    for t in T_dot:
        Mass_err = bootstrap.bootstrap_error(
            Mass_channel_tmp[0:-1, t], Mass_channel_tmp[-1, t]
        )
        mass_tmp.append(Mass_channel_tmp[-1, t])
        err_tmp.append(Mass_err)

    mass_tmp = np.array(mass_tmp)
    err_tmp = np.array(err_tmp)

    sprnext = next(sperater)
    mk = next(marker)

    select = abs(err_tmp / mass_tmp) < 0.08

    plt.errorbar(
        T_dot[select] + sprnext,
        mass_tmp[select],
        1.5*err_tmp[select],
        linestyle="",
        elinewidth=2.0,
        capsize=2,
        marker=mk,
        alpha=0.6,
        label=f"$aE_{n}$ fit: "+E_string,
    )
    
    plt.errorbar(
        T_dot[~select] + sprnext,
        mass_tmp[~select],
        1.5*err_tmp[~select],
        linestyle="",
        elinewidth=2.0,
        capsize=2,
        marker=mk,
        color=plt.gca().lines[-1].get_color(),
        alpha=0.1,
    )
    
    
    plt.gca().set_title(None)
    plt.ylim(0.0, 2.0)
    plt.xlabel('$t / a$', fontsize=15)
    plt.ylabel('$am_{\\rm{eff}}$', fontsize=15)
    plt.ylim(0.0, 2.0)
    return Mass_channel_tmp


def GEVP_fixT(Cmat, t0, ti, tf):
    Mshape = Cmat.shape

    Lambda_n = np.zeros(shape=(Mshape[0], Mshape[1], Mshape[2]))
    Vector_n = np.zeros(shape=Mshape, dtype=complex)

    T_dot = np.arange(ti, tf, 1, dtype=int)
    for t in T_dot:
        for N in range(Mshape[0]):
            value, vector = linalg.eig(
                Cmat[N, t], Cmat[N, t0], overwrite_a=True, overwrite_b=True
            )

            Vector_n[N, t] = vector

            for n in range(Mshape[2]):
                Lambda_n[N, t, n] = np.real(value[n])

    Lambda_n.sort(axis=2)

    return Lambda_n[:, :, ::-1], Vector_n
