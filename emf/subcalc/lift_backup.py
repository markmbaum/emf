from .. import np

def B_integral(xp, yp, zp, delta_x, delta_y, delta_z, xa, ya, za):

    delta_x_sq = delta_x**2.0
    delta_y_sq = delta_y**2.0
    delta_z_sq = delta_z**2.0

    xa_sq = xa**2.0
    ya_sq = ya**2.0
    za_sq = za**2.0

    xp_sq = xp**2.0
    yp_sq = yp**2.0
    zp_sq = zp**2.0

    #store the numerator of the integrated function, evaluated at 1
    lim_eq1_num = (delta_x_sq + delta_x*(xa - xp)
                + delta_y_sq + delta_y*(ya - yp)
                + delta_z*(delta_z + za - zp))

    #store a portion of the denominator of the integrated function,
    #evaluated at 1
    lim_eq1_sqrt = np.sqrt(
                delta_x_sq + 2*delta_x*xa - 2*delta_x*xp + xa_sq - 2*xa*xp
                + delta_y_sq + 2*delta_y*ya - 2*delta_y*yp + ya_sq - 2*ya*yp
                + delta_z_sq + 2*delta_z*za - 2*delta_z*zp + za_sq - 2*za*zp
                + xp_sq + yp_sq + zp_sq)

    #store the numerator of the integrated function, evaluated at 0*/
    lim_eq0_num = delta_x*(xa - xp) + delta_y*(ya - yp) + delta_z*(za - zp)

    #store a portion of the denominator of the integrated function,
    #evaluated at 0
    lim_eq0_sqrt = np.sqrt(
                xa_sq - 2*xa*xp
                + ya_sq - 2*ya*yp
                + za_sq - 2*za*zp
                + xp_sq + yp_sq + zp_sq)

    #store the other portion of the integrated function's denominator,
    #which is the same regardless of the evaluation value/point
    body = (delta_x_sq*(ya_sq - 2.0*ya*yp + za_sq
        - 2.0*za*zp + yp_sq + zp_sq) - 2.0*xa*(delta_x*(delta_y*ya - delta_y*yp
        + delta_z*za - delta_z*zp) + xp*(delta_y_sq + delta_z_sq))
        + 2.0*delta_x*xp*(delta_y*ya - delta_y*yp + delta_z*za - delta_z*zp)
        + xa_sq*(delta_y_sq + delta_z_sq) + delta_y_sq*za_sq
        - 2.0*delta_y_sq*za*zp + delta_y_sq*xp_sq + delta_y_sq*zp_sq
        - 2.0*delta_y*ya*delta_z*za + 2.0*delta_y*ya*delta_z*zp
        + 2.0*delta_y*delta_z*za*yp - 2.0*delta_y*delta_z*yp*zp
        + ya_sq*delta_z_sq - 2.0*ya*delta_z_sq*yp + delta_z_sq*xp_sq
        + delta_z_sq*yp_sq)

    #put it all together
    result = lim_eq1_num/(body*lim_eq1_sqrt) - lim_eq0_num/(body*lim_eq0_sqrt)

    return(result)

def lift(x, y, z, delta_x, delta_y, delta_z, xa, ya, za, N):

    #allocate return variables
    fac_x = np.zeros((N,), dtype=float)
    fac_y = np.zeros((N,), dtype=float)
    fac_z = np.zeros((N,), dtype=float)
    integral = np.zeros((N,), dtype=float)

    for i in range(N):
        #store the sample point's coordinates
        xp = x[i]
        yp = y[i]
        zp = z[i]

        #store the "fac" values, which are the components of a
        #three dimensional vector that will be scaled by the path
        #integral, the current, and the leading Biot-Savart constants
        fac_x[i] = delta_y*(zp - za) + delta_z*(ya - yp)
        fac_y[i] = delta_z*(xp - xa) + delta_x*(za - zp)
        fac_z[i] = delta_x*(yp - ya) + delta_y*(xa - xp)

        #store the path integral value for this location
        integral[i] = B_integral(xp, yp, zp,
                                delta_x, delta_y, delta_z,
                                xa, ya, za)

    return(fac_x, fac_y, fac_z, integral)
