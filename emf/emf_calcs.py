from .. import np

def phasors_to_magnitudes(*args):
    """Convert vectors of complex x and y phasors into real quantities, namely the amplitude of the field in the x and y directions, the product (the hypotenuse of the amplitudes), and the maximum field. Results of E_field and B_field can be passed directly to this function for conversion from phasor form into usable form.
    args:
        Ph_x - complex 1D numpy array, phasor horizontal components
        Ph_y - complex 1D numpy array, phasor vertical components
    returns:
        mag_x - 1D numpy array, maximum horizontal field
        mag_y - 1D numpy array, maximum vertical field
        prod - 1D numpy array, sqrt(mag_x**2 + mag_y**2)
        maximum - 1D numpy array, maximum field at any time"""

    #amplitude along each component, storing squared magnitudes for later
    mag_x_sq = np.real(Ph_x)**2 + np.imag(Ph_x)**2
    mag_x = np.sqrt(mag_x_sq)
    mag_y_sq = np.real(Ph_y)**2 + np.imag(Ph_y)**2
    mag_y = np.sqrt(mag_y_sq)

    #phase angle of each component
    phase_x = np.arctan2(np.imag(Ph_x), np.real(Ph_x))
    phase_y = np.arctan2(np.imag(Ph_y), np.real(Ph_y))

    #"product"
    prod = np.sqrt(mag_x**2 + mag_y**2)

    #maximum resultant value found by setting the time derivative of the
    #squared resultant magnitude to zero (Appendix 8.1 EPRI's "Big Red Book")
    num = mag_x_sq*np.sin(2*phase_x) + mag_y_sq*np.sin(2*phase_y)
    den = mag_x_sq*np.cos(2*phase_x) + mag_y_sq*np.cos(2*phase_y)
    t1 = (0.5)*np.arctan2(-num, den)
    t2 = t1 + np.pi/2
    term1 = mag_x_sq*(np.cos(t1 + phase_x))**2
    term2 = mag_y_sq*(np.cos(t1 + phase_y))**2
    ax_mag1 = np.sqrt(term1 + term2)
    term1 = mag_x_sq*(np.cos(t2 + phase_x))**2
    term2 = mag_y_sq*(np.cos(t2 + phase_y))**2
    ax_mag2 = np.sqrt(term1 + term2)

    #pick out the semi-major axis magnitude from the two semi-axis results
    maximum = np.maximum(ax_mag1, ax_mag2)

    #return the 4 output columns
    return(mag_x, mag_y, prod, maximum)
