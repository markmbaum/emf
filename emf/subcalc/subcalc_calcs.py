from .. import os, np, ctypes, itertools, _resource_filename

#laod C functions in the lift.so extension or fall back to the Python module
fn = _resource_filename(__name__, 'lift')
try:
    lift = np.ctypeslib.load_library('lift.so', fn)
except OSError as err:
    _lift_loaded = False
    print('The C extension for emf.subcalc cannot be loaded. System error message:\n\n%s\n\nemf.subcalc will use pure Python for calculations, which is slower.\nYou may need to recompile the C extension before it can be used.\nThe extension is in this directory:\n%s' % (str(err), fn))
    from . import lift_backup as lift
else:
    _lift_loaded = True

#pull pi out of numpy
pi = np.pi

#magnetic permeability constant, in SI units
MU_0 = 4*pi*1e-7
#convenient constant, converted for mG
magnetic_prefactor = 1.0e7*MU_0/(4.*pi)

def B_field_segment(a, b, I, ph, x, y, z):
    """Calculate the magnetic field generated by a current carrying wire segment at an arbitrary set of points in 3 dimensions
    args:
        a - 1D numpy array containing the x,y,z coordinates of the beginning
            of the wire segment (ft)
        b - 1D numpy array containing the x,y,z coordinates of the end
            of the wire segment (ft)
        I - float, the current of the wire segment (Amps)
        ph - float, the phase of the current (degrees)
        x - 1D numpy array, x coordinates of the sample points (ft)
        y - 1D numpy array, y coordinates of the sample points (ft)
        z - 1D numpy array, z coordinates of the sample points (ft)
    returns:
        Ph_x - 1D numpy array of complex numbers, magnetic field phasors in
               the x direction
        Ph_y - 1D numpy array of complex numbers, magnetic field phasors in
               the y direction
        Ph_z - 1D numpy array of complex numbers, magnetic field phasors in
               the z direction"""

    #conversions
    a  = a*0.3048           #convert to meters
    b  = b*0.3048           #convert to meters
    x  = x*0.3048           #convert to meters
    y  = y*0.3048           #convert to meters
    z  = z*0.3048           #convert to meters
    ph = ph*2.0*pi/360.0 #convert to radians

    #pull out the coordinates of the starting and ending point of the wire
    xa, ya, za = a
    xb, yb, zb = b

    #store the differences bewteen each dimension's coordinates
    delta_x = xb - xa
    delta_y = yb - ya
    delta_z = zb - za

    #convert the current and phase into a complex phasor
    I = I*(np.cos(ph) + complex(0,1)*np.sin(ph))

    #calculate fac_x, fac_y, fac_z, and integral with the C extension or the
    #backup Python module
    if(_lift_loaded):

        #use the functions in lift.c to do the heavy computational lifting of
        #these calculations, which are mostly made heavy by the path integral
        N = len(x)
        #variables to be filled in by lift.lift
        fac_x = np.empty((N,), dtype=np.double)
        fac_y = np.empty((N,), dtype=np.double)
        fac_z = np.empty((N,), dtype=np.double)
        integral = np.empty((N,), dtype=np.double)
        #variable of convenience
        c_double_ptr = ctypes.POINTER(ctypes.c_double)
        #call the calculations, converting everything to its C type
        lift.lift(x.ctypes.data_as(c_double_ptr),
                y.ctypes.data_as(c_double_ptr),
                z.ctypes.data_as(c_double_ptr),
                ctypes.c_double(delta_x),
                ctypes.c_double(delta_y),
                ctypes.c_double(delta_z),
                ctypes.c_double(xa),
                ctypes.c_double(ya),
                ctypes.c_double(za),
                ctypes.c_long(N),
                fac_x.ctypes.data_as(c_double_ptr),
                fac_y.ctypes.data_as(c_double_ptr),
                fac_z.ctypes.data_as(c_double_ptr),
                integral.ctypes.data_as(c_double_ptr))

    else:

        fac_x, fac_y, fac_z, integral = lift.lift(x, y, z,
                delta_x, delta_y, delta_z, xa, ya, za, len(x))

    #carry through with the rest of the calculations
    con  = magnetic_prefactor*I*integral
    Ph_x = fac_x*con
    Ph_y = fac_y*con
    Ph_z = fac_z*con

    return(Ph_x, Ph_y, Ph_z)

def grid_segment_results(Ph_x, Ph_y, Ph_z, x, y):
    """Convert the flat arrays returned by B_field_grid into 2D arrays
    args:
        Ph_x - complex 1D numpy array representing phasors in the x direction
        Ph_y - complex 1D numpy array representing phasors in the y direction
        Ph_z - complex 1D numpy array representing phasors in the z direction
        x - 1D numpy array, sorted x coordinates of sample grid (ft)
        y - 1D numpy array, sorted y coordinates of sample grid (ft)
    returns:
        Ph_x_grid - complex 2D numpy array, phasors in the x direction
        Ph_y_grid - complex 2D numpy array, phasors in the y direction
        Ph_z_grid - complex 2D numpy array, phasors in the z direction
        X - 2D array containing x coordinates of the sample grid
        Y - 2D array containing y coordinates of the sample grid"""

    Nx, Ny = len(x), len(y)
    Ph_x_grid = np.flipud(np.reshape(Ph_x, (Nx, Ny)).T)
    Ph_y_grid = np.flipud(np.reshape(Ph_y, (Nx, Ny)).T)
    Ph_z_grid = np.flipud(np.reshape(Ph_z, (Nx, Ny)).T)
    X, Y = np.meshgrid(x, y[::-1], indexing='xy')

    return(Ph_x_grid, Ph_y_grid, Ph_z_grid, X, Y)

def phasors_to_magnitudes(Ph_x, Ph_y, Ph_z):
    """Convert vectors of complex x, y, and z phasors into real quantities, namely the amplitude of the field in the x, y, and z directions, the product/resultant (the hypotenuse of the amplitudes), and the time-dependent maximum field
    args:
        Ph_x - complex 3D numpy array, phasor x components
        Ph_y - complex 3D numpy array, phasor y components
        Ph_z - complex 3D numpy array, phasor z components
    returns:
        mag_x - 3D numpy array, maximum x field
        mag_y - 3D numpy array, maximum y field
        mag_z - 3D numpy array, maximum z field
        prod - 3D numpy array, sqrt(mag_x**2 + mag_y**2)
        maximum - 3D numpy array, maximum field at any time"""

    #amplitude along each component, storing squared magnitudes for later
    mag_x_sq = np.real(Ph_x)**2 + np.imag(Ph_x)**2
    mag_x    = np.sqrt(mag_x_sq)
    mag_y_sq = np.real(Ph_y)**2 + np.imag(Ph_y)**2
    mag_y    = np.sqrt(mag_y_sq)
    mag_z_sq = np.real(Ph_z)**2 + np.imag(Ph_z)**2
    mag_z    = np.sqrt(mag_z_sq)

    #phase angle of each component
    phase_x = np.arctan2(np.imag(Ph_x), np.real(Ph_x))
    phase_y = np.arctan2(np.imag(Ph_y), np.real(Ph_y))
    phase_z = np.arctan2(np.imag(Ph_z), np.real(Ph_z))

    #"product"
    prod = np.sqrt(mag_x**2 + mag_y**2 + mag_z**2)

    #maximum resultant value found by setting the time derivative of the
    #squared resultant magnitude to zero (Appendix 8.1 EPRI's "Big Red Book")
    num     = (mag_x_sq*np.sin(2*phase_x)
                + mag_y_sq*np.sin(2*phase_y)
                + mag_z_sq*np.sin(2*phase_z))
    den     = (mag_x_sq*np.cos(2*phase_x)
                + mag_y_sq*np.cos(2*phase_y)
                + mag_z_sq*np.cos(2*phase_z))
    t1      = (0.5)*np.arctan2(-num, den)
    t2      = t1 + pi/2
    x_term  = mag_x_sq*(np.cos(t1 + phase_x))**2
    y_term  = mag_y_sq*(np.cos(t1 + phase_y))**2
    z_term  = mag_z_sq*(np.cos(t1 + phase_z))**2
    ax_mag1 = np.sqrt(x_term + y_term + z_term)
    x_term  = mag_x_sq*(np.cos(t2 + phase_x))**2
    y_term  = mag_y_sq*(np.cos(t2 + phase_y))**2
    z_term  = mag_z_sq*(np.cos(t2 + phase_z))**2
    ax_mag2 = np.sqrt(x_term + y_term + z_term)

    #pick out the semi-major axis magnitude from the two semi-axis results
    maximum = np.maximum(ax_mag1, ax_mag2)

    #return the 5 output columns
    return(mag_x, mag_y, mag_z, prod, maximum)
