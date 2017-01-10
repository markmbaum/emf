typedef struct Results {
    double *fac_x;
    double *fac_y;
    double *fac_z;
    double *integral;
} Results;

double B_integral (double xp, double yp, double zp,
        double delta_x, double delta_y, double delta_z,
        double xa, double ya, double za);

Results compute (double *x, double *y, double *z,
        double delta_x, double delta_y, double delta_z,
        double xa, double ya, double za,
        int Nx, int Ny, int Nz);
