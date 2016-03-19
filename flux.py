import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.analytic_functions import blackbody_lambda


def fluxcalculator(temp):
    """Given a temperature this functions plots a graph between the flux emitted and a set of wavelengths"""      

    temp = temp * u.K
    wavemax = (const.b_wien / temp).to(u.AA)  # Wien's displacement law
    waveset = np.logspace(
    0, np.log10(wavemax.value + 10*wavemax.value), num=100) * u.AA
    with np.errstate(all='ignore'):
        flux = blackbody_lambda(waveset, temp)

    #Plotting
    
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot(waveset.value, flux.value)
    ax.axvline(wavemax.value, ls='--')
    ax.get_yaxis().get_major_formatter().set_powerlimits((0, 1))
    ax.set_xlabel(r'$\lambda$ ({0})'.format(waveset.unit))
    ax.set_ylabel(r'$B_{\lambda}(T)$' + ' ' + '({0})'.format(flux.unit))
    ax.set_title('Blackbody, T = {0}'.format(temp))
    plt.show(block=True)

a = fluxcalculator(6000)
