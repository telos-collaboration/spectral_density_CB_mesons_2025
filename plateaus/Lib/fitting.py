import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy.optimize import minimize
from scipy.optimize import approx_fprime
from scipy.stats import chisquare

import corrfitter as cf
import gvar as gv

from bootstrap import bootstrap_error

######## print ########


def find_non_zero(x):
    if x < 1:
        A = "%e" % x
        return int(A.partition("-")[2]) - 1
    else:
        print("error > 1")
        return 0


def print_non_zero(v, err):
    dig = find_non_zero(err) + 2
    buff = 0

    if dig == 0:
        if v < 0:
            buff = 1
        return str(v)[0 : dig + 3 + buff] + "(" + str(err)[dig : dig + 3] + ")"
    else:
        if v < 0:
            buff = 1
        return str(v)[0 : dig + 2 + buff] + "(" + str(err)[dig : dig + 2] + ")"


def printX(x):
    return np.array2string(x, formatter={"float_kind": lambda x: "%.12f" % x})


######## mass extraction with exp. function ########


def X2_single_state_fit(C, ti, tf):
    num_sample = np.shape(C)[0]
    # print(num_sample)
    T = np.shape(C)[1]

    def func(t, a, M):
        return a * a * M * (np.exp(-M * t) + np.exp(-M * (T - t))) / 2

    def Cov(ti, tj):
        return np.mean((C[:, ti] - np.mean(C[:, ti])) * (C[:, tj] - np.mean(C[:, tj])))
        # return np.mean( ( C[:, ti]  *  C[:, tj]) - np.mean(C[:, ti])*np.mean(C[:, tj])  )

    size = tf - ti
    M = np.mat(np.zeros(shape=(size, size)))
    for a in np.arange(ti, tf, 1):
        for b in np.arange(ti, tf, 1):
            M[a - ti, b - ti] = Cov(a, b)

    M = np.mat(np.cov(C[0:-1, ti:tf].T))

    M_I = M.I

    def X2_boot_const(X):
        def Cf_vector(ti, tf):
            V = np.zeros(tf - ti)
            for time in range(tf - ti):
                V[time] = C[N, time + ti] - func(time + ti, a, E)

            return V

        (a, E) = X

        V = np.mat(Cf_vector(ti, tf))

        return V * M_I * V.T

    def nor_X2(a, E1):
        def Cf_vector_nor(ti, tf):
            t_dot = np.arange(ti, tf - ti)
            V = np.zeros(tf - ti)
            for time in range(tf - ti):
                V[time] = C[N, time + ti] - func(time + ti, a, E1)

            return V

        V = np.mat(Cf_vector_nor(ti, tf))

        # print(r'Xsqr/d.o.f.='+ str((V*M_I*V.T)[0,0] / (tf-ti-2)))
        return (V * M_I * V.T)[0, 0]

    E1 = []
    a1 = []

    t_slice = np.arange(ti, tf)

    x0, pcov = curve_fit(func, t_slice, C[-1, ti:tf])
    print("initial: ", x0)

    print("one-state fitting time region ", ti, "to", tf - 1, ": ")

    for N in range(num_sample):
        # print N

        res = minimize(X2_boot_const, x0, method="Nelder-Mead", tol=10**-8)

        E1.append(res.x[1])
        a1.append(res.x[0])
        # XX.append(nor_X2(res.x[0], res.x[1]))

        # print res

    X2 = nor_X2(a1[-1], E1[-1])
    E_err = bootstrap_error(E1[0:-1], E1[-1])

    print(E1[-1], E_err, round(X2 / (tf - ti - 2), 2))

    return np.array(E1), E_err, X2 / (tf - ti - 2)


def X2_single_exp_fit(C, ti, tf):
    num_sample = np.shape(C)[0]
    T = np.shape(C)[1] - 1

    def func(t, a, M):
        return a * M * np.exp(-M * t)

    def Cov(ti, tj):
        return np.mean((C[:, ti] - np.mean(C[:, ti])) * (C[:, tj] - np.mean(C[:, tj])))

    size = tf - ti
    M = np.mat(np.zeros(shape=(size, size)))
    for a in np.arange(ti, tf, 1):
        for b in np.arange(ti, tf, 1):
            M[a - ti, b - ti] = Cov(a, b)

    M_I = M.I

    def X2_boot_const(X):
        def Cf_vector(ti, tf):
            V = np.zeros(tf - ti)
            for time in range(tf - ti):
                V[time] = C[N, time + ti] - func(time + ti, a, E)

            return V

        (a, E) = X

        V = np.mat(Cf_vector(ti, tf))

        return V * M_I * V.T

    def nor_X2(a, E1):
        def Cf_vector_nor(ti, tf):
            t_dot = np.arange(ti, tf - ti)
            V = np.zeros(tf - ti)
            for time in range(tf - ti):
                V[time] = C[N, time + ti] - func(time + ti, a, E1)

            return V

        # print 'A=', a

        V = np.mat(Cf_vector_nor(ti, tf))

        print(r"Xsqr/d.o.f.=" + str((V * M_I * V.T)[0, 0] / float(tf - ti - 2)))

        return (V * M_I * V.T)[0, 0]

    E1 = []
    a1 = []

    t_slice = np.arange(ti, tf)

    x0, pcov = curve_fit(func, t_slice, C[-1, ti:tf])
    print("x0=", x0)

    print("A*m*exp(-mt) fitting time region ", ti, "to", tf - 1, ": ", end="")

    for N in range(num_sample):
        res = minimize(X2_boot_const, x0, method="Nelder-Mead", tol=10**-8)

        E1.append(res.x[1])
        a1.append(res.x[0])

    X2 = nor_X2(a1[-1], E1[-1])
    E_err = bootstrap_error(E1[0:-1], E1[-1])

    return np.array(E1), E_err, X2


def make_models(T, tmin, tmax, tp):
    """Create corrfitter model for G(t)."""
    return [cf.Corr2(datatag="Gab", tp=tp, tmin=tmin, tmax=tmax, a="a", b="a", dE="dE")]


def make_prior(N):
    prior = gv.BufferDict()
    # setting the sdev of the prioir to infinity amounts to turning off the prior contribution to chi2
    prior["log(a)"] = gv.log(gv.gvar(N * [0.1], N * [np.Inf]))
    prior["log(dE)"] = gv.log(gv.gvar(N * [0.1], N * [np.Inf]))
    return prior


def first_fit_parameters(fit):
    p = fit.p
    E = np.cumsum(p["dE"])
    a = p["a"]
    chi2 = fit.chi2
    dof = fit.dof
    return E, a, chi2, dof


def bootstrap_fit(fitter, dset, T, tmin, tmax, tp):
    gv.ranseed(12)
    # dset = gv.dataset.Dataset(ifile)
    pdatalist = (
        cf.process_dataset(ds, make_models(T, tmin, tmax, tp))
        for ds in gv.dataset.bootstrap_iter(dset, n=10)
    )
    bs = gv.dataset.Dataset()
    for bsfit in fitter.bootstrapped_fit_iter(pdatalist=pdatalist):
        bs.append(
            E=np.cumsum(bsfit.pmean["dE"][:2]),
            a=bsfit.pmean["a"][:2],
        )
    bs = gv.dataset.avg_data(bs, bstrap=True)
    E = bs["E"]
    a = bs["a"]
    print("{:2}  {:15}  {:15}".format("E", E[0], E[1]))
    print("{:2}  {:15}  {:15}".format("a", a[0], a[1]))


def fit_correlator_without_bootstrap(
    data_corr, T, tmin, tmax, Nmax, tp, plotting=False, printing=False
):
    T = abs(T)

    fitter = cf.CorrFitter(models=make_models(T, tmin, tmax, tp))
    p0 = None
    for N in range(1, Nmax + 1):
        prior = make_prior(N)
        fit = fitter.lsqfit(data=data_corr, prior=prior, p0=p0)
        p0 = fit.pmean

        if printing:
            print("nterm =", N, 30 * "=")
            print(fit)

    E, a, chi2, dof = first_fit_parameters(fit)
    if plotting:
        # fit.show_plots(view='ratio')
        fit.show_plots(view="log")

    return E, a, chi2, dof


def fit_exp_std(C_boot, TI, TF, GLB_T):
    # This function fits the correlators with a cosh function
    # the error is estimated by standard deviation with a covariance matrix

    # load the ensamble info
    num_sample = C_boot.shape[0]

    print("A* ( exp[-mt] ) fitting time region ", TI, "to", TF, ": ")

    # std = np.std(C_boot, axis=0)
    cov = np.cov(C_boot[0:-1].T)

    if np.isnan(cov.sum()):
        print("cov contain nan")

    if np.isnan(cov.T.sum()):
        print("cov contain nan")

    CORR = dict(Gab=gv.gvar(C_boot[-1], cov))

    E, a, chi2, dof = fit_correlator_without_bootstrap(
        CORR, 0, TI, TF, 1, None, plotting=False, printing=True
    )

    print(gv.mean(E[0]), gv.sdev(E[0]), round(chi2 / dof, 2))

    return gv.mean(E[0]), gv.sdev(E[0]), chi2 / dof


def fit_cosh_std(C_boot, TI, TF, GLB_T):
    # This function fits the correlators with a cosh function
    # the error is estimated by standard deviation with a covariance matrix

    # load the ensamble info
    num_sample = C_boot.shape[0]

    print("A* ( exp[-mt] + exp[-m(T-t)] ) fitting time region ", TI, "to", TF, ": ")

    std = np.std(C_boot, axis=0)
    cov = np.cov(C_boot[0:-1].T)

    if np.isnan(cov.sum()):
        print("cov contain nan")

    if np.isnan(cov.T.sum()):
        print("cov contain nan")

    CORR = dict(Gab=gv.gvar(C_boot[-1], cov))

    E, a, chi2, dof = fit_correlator_without_bootstrap(
        CORR, 0, TI, TF, 1, GLB_T, plotting=False, printing=True
    )

    print(gv.mean(E[0]), gv.sdev(E[0]), round(chi2 / dof, 2))

    return gv.mean(E[0]), gv.sdev(E[0]), chi2 / dof


def fit_cosh_booerr(C_boot, TI, TF, GLB_T):
    # This function fits the correlators with a cosh function
    # the error is estimated by the bootstrap distribution width ( 68% )

    # load the ensamble info
    num_sample = C_boot.shape[0]

    print("A* ( exp[-mt] + exp[-m(T-t)] ) fitting time region ", TI, "to", TF, ": ")

    E_sample = np.zeros(num_sample)
    a_sample = np.zeros(num_sample)
    chi2_dof = np.zeros(num_sample)

    std = np.std(C_boot, axis=0)
    cov = np.cov(C_boot.T)

    for n in range(num_sample):
        CORR = dict(Gab=gv.gvar(C_boot[n], cov))

        pp = False

        if n == num_sample - 1:
            pp = True

        E, a, chi2, dof = fit_correlator_without_bootstrap(
            CORR, 0, TI, TF, 1, GLB_T, plotting=False, printing=pp
        )
        E_sample[n] = gv.mean(E[0])
        a_sample[n] = gv.mean(a[0])
        chi2_dof[n] = chi2 / dof

    err = bootstrap_error(E_sample[0:-1], E_sample[-1])
    print(E_sample[-1], err, round(chi2_dof[-1], 2), gv.sdev(E[0]))

    return E_sample[-1], err, chi2_dof[-1]
