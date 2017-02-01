double B_longegral (double xp, double yp, double zp,
        double delta_x, double delta_y, double delta_z,
        double xa, double ya, double za);

void lift_grid (double *x, double *y, double *z,
        double delta_x, double delta_y, double delta_z,
        double xa, double ya, double za,
        long Nx, long Ny, long Nz,
        double *fac_x, double *fac_y, double *fac_z,
        double *integral);

void lift_general (double *x, double *y, double *z,
        double delta_x, double delta_y, double delta_z,
        double xa, double ya, double za, long N,
        double *fac_x, double *fac_y, double *fac_z,
        double *integral);
