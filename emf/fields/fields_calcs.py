from .. import np

epsilon = 8.854e-12   #electric permeability constant, in SI units
C_Eepsilon = 1./(2.*np.pi*epsilon)  #convenient constant

mu = 4*np.pi*1e-7     #magnetic permeability constant, in SI units
C_B = 1.0e7*mu/(2.*np.pi)   #convenient constant, converted for mG

def E_field(x_cond, y_cond, subconds, d_cond, d_bund, V_cond, p_cond, x, y):
    """Calculate the approximate electric field generated by a group of
    conductors. Each of the inputs labeled '_cond' should be an numpy
    array of parameters, where each index in those arrays describes a
    unique conductor, i.e. the 0th value in each variable is attributed to
    one power line."""

    #conversions and screening out underground lines
    ohd = y_cond > 0.
    x_cond = x_cond[ohd]*0.3048         #convert to meters
    y_cond = y_cond[ohd]*0.3048         #convert to meters
    subconds = subconds[ohd]
    d_cond = d_cond[ohd]*0.0254         #convert to meters
    d_bund = d_bund[ohd]*0.0254         #convert to meters
    V_cond = V_cond[ohd]/np.sqrt(3)     #convert to ground reference from
                                            #line-line reference, leave in kV
    p_cond = p_cond[ohd]*2*np.pi/360.   #convert to radians
    x = x*0.3048                        #convert to meters
    y = y*0.3048                        #convert to meters

    #array length variables
    N = len(p_cond)     #number of conductors
    Z = len(x)          #number of sample points or x,y pairs

    #calculate the effective conductor diameters
    d_cond  = d_bund*((subconds*d_cond/d_bund)**(1./subconds))

    #compute the matrix of potential coefficients
    range_N = range(N)
    P = np.empty((N,N))
    #diagonals
    P[range_N, range_N] = C_Eepsilon*np.log(4*y_cond/d_cond)
    #other elements
    for a in range_N:
        for b in range_N:
            if(a != b):
                n = (x_cond[a] - x_cond[b])**2 + (y_cond[a] + y_cond[b])**2
                d = (x_cond[a] - x_cond[b])**2 + (y_cond[a] - y_cond[b])**2
                P[a,b] = C_Eepsilon*np.log(np.sqrt(n/d))

    #initialize complex voltage phasors
    V = V_cond*(np.cos(p_cond) + complex(0,1)*np.sin(p_cond))

    #compute real and imaginary charge phasors
    Q = np.linalg.solve(P, V)

    #compute components of the electric field phasors at each point, with each
    #column of E_x and E_y storing components due to each conductor, so that the
    #rows represent a spatial point across the ROW or an x,y pair
    Ex = np.empty((N,Z))
    Ey = np.empty((N,Z))
    #first compute the coefficients without the charges
    for a in range(Z):
        #denominators, squared distance between the point and the conductors
        d1 = (x[a] - x_cond)**2 + (y[a] - y_cond)**2
        d2 = (x[a] - x_cond)**2 + (y[a] + y_cond)**2
        #x component numerator, the same for the conductor and its image
        nx = C_Eepsilon*(x[a] - x_cond)
        #y component numerators, different for the conductor and its image
        ny1 = C_Eepsilon*(y[a] - y_cond)
        ny2 = C_Eepsilon*(y[a] + y_cond)
        #evaluate
        Ex[:,a] = nx/d1 - nx/d2
        Ey[:,a] = ny1/d1 - ny2/d2

    #multiply the charges by the field coefficients calculated above
    Q = np.tile(np.reshape(Q, (N,1)), (1,Z))
    Ex = Ex*Q
    Ey = Ey*Q

    #sum the phasors for each sample point, yielding the sum of phasors for each
    #conductor, which are the final phasors for each point
    Ex = np.sum(Ex, axis = 0)
    Ey = np.sum(Ey, axis = 0)

    #return phasors, complex numbers, for the x and y components
    #   - these complex phasors are converted to real valued outputs by the
    #   - phasors_to_magnitudes() function
    return(Ex, Ey)

def B_field(x_cond, y_cond, I_cond, p_cond, x, y):
    """Calculate the approximate magnetic field generated by a group of
    conductors. Each of the variables labeled '_cond' should be an numpy
    array of parameters, where each index in those arrays describes a
    unique conductor, i.e. the 0th value in each variable is attributed to
    one power line."""

    #array length variables
    N = len(p_cond)         #number of conductors
    Z = len(x)              #number of sample points or x,y pairs

    #conversions
    x_cond = x_cond*0.3048          #convert to meters
    y_cond = y_cond*0.3048          #convert to meters
    p_cond = p_cond*2*np.pi/360.    #convert to radians
    x = x*0.3048                    #convert to meters
    y = y*0.3048                    #convert to meters

    #initialize complex current phasors
    I = I_cond*(np.cos(p_cond) + complex(0,1)*np.sin(p_cond))

    #compute magnetic field component phasors for each x,y point
    Bx = np.zeros((Z,), dtype = complex)
    By = np.zeros((Z,), dtype = complex)
    for a in range(Z): #x,y pairs
        for b in range(N): #conductors
            dx = x[a] - x_cond[b]
            dy = y[a] - y_cond[b]
            #magnitude phasor
            B = C_B*I[b]/np.sqrt(dx**2 + dy**2)
            #break it up into components
            theta = np.arctan(abs(dy/dx))
            #x component calculated with sine and y component with cosine
            #because the field is perpendicular to the line to the conductor,
            #potentially contrary to one's first instinct that sine goes
            #with y components and cosine with x components, but it's right
            Bx[a] -= np.sign(dy)*np.sin(theta)*B
            By[a] += np.sign(dx)*np.cos(theta)*B

    #return phasors, complex numbers, for the x and y components
    #   - these complex phasors are converted to real valued outputs by the
    #   - phasors_to_magnitudes() function
    return(Bx, By)

def phasors_to_magnitudes(Ph_x, Ph_y):
    """Convert vectors of complex x and y phasors into real quantities,
    namely the amplitude of the field in the x and y directions, the
    product (the hypotenuse of the amplitudes), and the maximum field.
    args:
        Ph_x - numpy vector or number, phasor x component
        Ph_y - numpy vector or number, phasor y component"""
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
    #return the 4 output colums
    return(mag_x, mag_y, prod, maximum)
