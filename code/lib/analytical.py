import scipy.special
import numpy as np

def psd_alphatau(x, alpha, tau):
    return (np.exp(-((1.0/2)) * alpha * (x + np.abs(x))) * (scipy.special.erfc((-alpha * tau + np.abs(x))/(
    2 * np.sqrt(tau))) - 
   np.exp((x * alpha)/np.sign(x))*
     scipy.special.erfc((alpha * tau + np.abs(x))/(
     2 * np.sqrt(tau)))))/(2 * alpha * tau)

def psd_whitenoise(x, f0, T, sigma):
    return (np.exp((f0*np.abs(x)*(np.sign(x)-1))/(2 *sigma**2))*
        (scipy.special.erfc((np.abs(x) - f0 * T)/(2 * np.sqrt(T) * sigma)) - 
 np.exp((f0 * np.abs(x))/sigma**2) * scipy.special.erfc((f0 * T + np.abs(x))/(2 * np.sqrt(T) * sigma)))/(2 * f0 * T))

def csd_whitenoise(x, f0, T, sigma):
    xl = x
    return (f0**2 * T + (
 2 * np.exp(-((-f0 * T + xl)**2/(4 * T * sigma**2))) * f0 * np.sqrt(
  T) * sigma)/np.sqrt(np.pi) - sigma**2 + (f0**2 * T - sigma**2) * scipy.special.erf((
   f0 * T - xl)/(2 * np.sqrt(T) * sigma)) -
 f0 * xl * scipy.special.erfc((-f0 * T + xl)/(2 * np.sqrt(T) * sigma)) +
 np.exp((f0 * xl)/sigma**2) * sigma**2 * scipy.special.erfc((f0 * T + xl)/(
   2 * np.sqrt(T) * sigma)))/(2 * f0**2 * T)

def csd_alphatau(x, alpha, tau):
    xl = x
    
    xpos = (-((2 * np.exp(-((xl + alpha * tau)**2/(
    4 * tau))) * alpha * np.sqrt(tau))/np.sqrt(np.pi)) + 
 np.exp(-xl * alpha)*
   scipy.special.erfc((xl - alpha * tau)/(
   2 * np.sqrt(tau))) + (-1 + alpha * (xl + alpha * tau)) * scipy.special.erfc((
   xl + alpha * tau)/(2 * np.sqrt(tau))))/(2 * alpha**2 * tau)
    
    xneg = (1 + alpha * (-xl - (
    2 * np.exp(-((xl + alpha * tau)**2/(4 * tau))) * np.sqrt(tau))/
    np.sqrt(np.pi) + alpha * tau) - (-1 + alpha * (xl + alpha*tau))
            * scipy.special.erf((xl + alpha * tau)/(2 * np.sqrt(tau))) - 
 np.exp(-xl * alpha)*
   scipy.special.erfc((-xl + alpha * tau)/(
   2 * np.sqrt(tau))))/(2 * alpha**2 * tau)

    return np.where(x>0, xpos, xneg)


def psd_early(x, f0, T, sigma, tp):
    return (-1 + scipy.special.erf((f0 * T - x)/(2 * np.sqrt(T) * sigma)) + 
 np.exp((f0 * x)/sigma**2) * (scipy.special.erf((f0 * T + x)/(2 * np.sqrt(T) * sigma)) - 
    scipy.special.erf((f0 * (T - tp) + x)/(2 * np.sqrt(T - tp) * sigma))) + 
 scipy.special.erfc((f0 * (T - tp) - x)/(2 * np.sqrt(T - tp) * sigma)))/(2 * f0)

