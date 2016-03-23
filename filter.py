from astropy.analytic_functions import blackbody_lambda
from astropy import units as u
from astropy import constants as const
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from numpy import trapz
temperature = 10000.0 * u.K
from math import log10
from math import ceil

R = const.R_sun / 2.0
r = 10.0 * const.pc
SI = u.J / (u.second * (u.meter ** 2) * u.meter * u.sr)


# Calculating flux through different filters
def flux_filter(central_wavelength, bandwidth, temperature):

	domain = (np.logspace(np.log10(central_wavelength - bandwidth / 2), np.log10(central_wavelength + bandwidth / 2), num=100) * u.AA).to(u.meter)
	flux = blackbody_lambda(domain, temperature)
        flux = flux.to(SI)
        with np.errstate(all='ignore'):
	    flux_earth = flux * (R / r) ** 2 #flux recieved on earth
	area = trapz(flux_earth.value, domain) * np.pi * (u.W / u.meter ** 2)

	return area

# Function for calculating magnitude 
def magnitude(f, fr):
	magnitude = - 2.5 * log10(f.value / fr.value)
	return magnitude

#Ultraviolet filter

central_wavelength = 3735.0 * u.AA
bandwidth = 485.0 * u.AA
flux_earth = 1.81e-9 * (u.W / u.meter ** 2)

flux = flux_filter(central_wavelength.value, bandwidth.value, temperature)
magnitude = magnitude(flux, flux_earth)

print "The magnitude for ultraviolet filter is %s" %round(magnitude, 2)

#Blue filter

central_wavelength = 4443.0 * u.AA
bandwidth = 831.0 * u.AA
flux_earth = 4.26e-9 * (u.W / u.meter ** 2)

flux = flux_filter(central_wavelength.value, bandwidth.value, temperature)
magnitude = magnitude(flux, flux_earth)

print "The magnitude for Blue filter is %s" %round(magnitude, 2)

#Violet filter

central_wavelength = 5483.0 * u.AA
bandwidth = 827.0 * u.AA
flux_earth = 3.64e-9 * (u.W / u.meter ** 2)

flux = flux_filter(central_wavelength.value, bandwidth.value, temperature)
magnitude = magnitude(flux, flux_earth)

print "The magnitude for Violet filter is %s" %round(magnitude, 2)

#Red filter

central_wavelength = 6855.0 * u.AA
bandwidth = 1742.0 * u.AA
flux_earth = 3.08e-9 * (u.W / u.meter ** 2)

flux = flux_filter(central_wavelength.value, bandwidth.value, temperature)
magnitude = magnitude(flux, flux_earth)

print "The magnitude for Red filter is %s" %round(magnitude, 2)
