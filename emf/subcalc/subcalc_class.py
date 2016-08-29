import numpy as np

from ..emf_class import EMFError

class Footprint(object):

    def get_x(self):
        if(self.close_pts):
            return(np.array(self._x + [self._x[0]]))
    x = property(get_x)

    def get_y(self):
        if(self.close_pts):
            return(np.array(self._y + [self._y[0]]))
    y = property(get_y)

    def __init__(self, name, x, y, of_concern, close_pts):
        #check x and y are the same length
        if(len(x) != len(y)):
            raise(EMFError("""
            Footprints must have the same number of x and y values to form
            spatial coordinates"""))
        #set attributes
        self.name = name
        self._x = x
        self._y = y
        self.of_concern = of_concern
        self.close_pts = close_pts

class Model(object):

    def __init__(self, data, grid_info):
        #dictionary of reference grid information from the SubCalc model
        self.grid_info = grid_info
        #1D results vectors
        self.x = data['x']
        self.y = data['y']
        self.Bmax = data['Bmax']
        #2D reference grid arrays
        self.grid_x, self.grid_y, self.grid_Bmax = self._meshgrid()
        #other reference objects in the model, like substation boundaries,
        #stored in a list of Footprint objects
        self.footprints = []

    def _meshgrid(self):

        #find the number of points in a row
        x = self.x
        y = self.y
        Bmax = self.Bmax
        count = 0
        v = y[count]
        while(y[count+1] == y[count]):
            count += 1
        count += 1

        #convert to ncols and nrows
        L = len(x)
        ncols = count
        nrows = L/ncols

        #fill 2D arrays
        X, Y, BMax = [np.empty((nrows,ncols)) for i in range(3)]
        count = 0
        for i in range(nrows):
            X[i,:] = x[count:count+ncols]
            Y[i,:] = y[count:count+ncols]
            BMax[i,:] = Bmax[count:count+ncols]
            count += ncols

        return(X, Y, BMax)

    def field_at_point(self, x, y):
        """Interpolate in the x and y directions to find an estimated Bmax
        value at an x,y location within the model"""
        #check that the point is in the grid
        if((x > self.grid_info['Maximum X Coordinate']) or
            (x < self.grid_info['Minimum X Coordinate']) or
            (y > self.grid_info['Maximum Y Coordinate']) or
            (y < self.grid_info['Minimum Y Coordinate'])):
            raise(EMFError("""
            x,y coordinates must fall inside the reference grid"""))
        #get indices of two points closest to the x,y along both axes
        xdiff = x - self.grid_x[0,:]
        xidx = np.arange(len(xdiff))[xdiff <= 0.][:2]
        ydiff = y - self.grid_y[:,0]
        yidx = np.arange(len(ydiff))[ydiff >= 0.][:2][::-1]
        #interpolate
        mx = (self.grid_Bmax[yidx[0], xidx[1]]
                - self.grid_Bmax[yidx[0], xidx[0]])/(
                    self.grid_x[yidx[0], xidx[1]]
                        - self.grid_x[yidx[0], xidx[0]])
        my = (self.grid_Bmax[yidx[1], xidx[0]]
                - self.grid_Bmax[yidx[0], xidx[0]])/(
                    self.grid_y[yidx[1], xidx[0]]
                        - self.grid_y[yidx[0], xidx[0]])
        dx = xdiff[xidx[0]]
        dy = ydiff[yidx[0]]
        print mx,dx,my,dy
        B = self.grid_Bmax[yidx[0], xidx[0]] + mx*dx + my*dy

        return(B)
