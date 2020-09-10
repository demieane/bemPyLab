import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt

# solution of Laplace in 2d with neumman boundary conditions
# with bem

def funcSource(x, y, xo, yo):
    Phi = np.log(math.sqrt((x-xo)**2 + (y-yo)**2))
    return Phi

def boundary():
    # the normal and tangent vectors
    # visualization of the domain
    # for a rectangular domain
    pass

def lineParametrization(xa, ya, xb, yb, Npanels):
    # from A to B
    t = np.linspace(0,1,Npanels)
    #print(t)
    xline = xa + t*(xb-xa)
    yline = ya + t*(yb-ya)
    return xline, yline

def boundaryConditions(nx, ny, x, y, xo, yo):
    # Neumann BCs
    dPhidn = nx*(x-xo)/((x-xo)**2-(y-yo)**2)+ny*(y-yo)/((x-xo)**2-(y-yo)**2)
    return dPhidn

def collocationScheme(xi, yi):
    #normal and tangent unit vectors
    #xtest = [1,2]
    sizeN= len(xi)-1
    xcolloc = np.zeros(sizeN);ycolloc = np.zeros(sizeN)
    sin_theta = np.zeros(sizeN);cos_theta = np.zeros(sizeN)
    nx = np.zeros(sizeN);ny = np.zeros(sizeN)
    tx = np.zeros(sizeN);ty = np.zeros(sizeN)

    for ii in range(0,len(xcolloc)):
        xcolloc[ii] = 0.5*(xi[ii+1]+xi[ii])
        ycolloc[ii] = 0.5*(yi[ii+1]+yi[ii])

        lpanel = math.sqrt((xi[ii+1]-xi[ii])**2+(yi[ii+1]-yi[ii])**2) #panel length
        sin_theta[ii] = (yi[ii+1]-yi[ii])/lpanel
        cos_theta[ii] = (xi[ii+1]-xi[ii])/lpanel
        nx[ii] = sin_theta[ii]
        ny[ii] = -cos_theta[ii]
        tx[ii] = -ny[ii]
        ty[ii] = nx[ii]
    return xcolloc, ycolloc, nx, ny, tx, ty

def singularIntegral():
    # induction factors for sources, analytic solution
    pass

def eqs():
    #set up matrtices for the system of equations
    pass

def show():
   return matplotlib.pyplot.show(block=True)

#===========================MAIN CODE BELOW====================================#

#==============CREATE MESH=====================================================#
Np = 3 #number of panels
xo, yo = -2, 0 #source position
# square domain creation
xa, ya = 0, 0; xb, yb = 1, 1; xc, yc = 0, 2; xd, yd = -1, 1
[x1, y1] =lineParametrization(xa, ya, xb, yb, Np)
[x2, y2] =lineParametrization(xb, yb, xc, yc, Np)
[x3, y3] =lineParametrization(xc, yc, xd, yd, Np)
[x4, y4] =lineParametrization(xd, yd, xa, ya, Np)

# plot the discretized domain
plt.plot(xo, yo, 'ks')
plt.plot(x1, y1, x2, y2, x3, y3, x4, y4)
plt.plot(x1, y1, 'bo', x2, y2, 'bo', x3, y3, 'bo', x4, y4, 'bo')
plt.plot(xa, ya, 'r*', xb, yb, 'r*', xc, yc, 'r*', xd, yd, 'r*')
plt.grid(linestyle='-', linewidth=0.5)
##show()

#nodes on the boundary with counter-clockwise numbering
xi =np.concatenate((x1, x2[1:len(x2)], x3[1:len(x3)], x4[1:len(x4)]))
yi =np.concatenate((y1, y2[1:len(y2)], y3[1:len(y3)], y4[1:len(y4)]))
print(xi)
#print(yi)

[xcolloc, ycolloc, nx, ny, tx, ty] = collocationScheme(xi, yi)

plt.plot(xcolloc, ycolloc, 'gs')
plt.quiver(xcolloc, ycolloc, nx, ny)
plt.quiver(xcolloc, ycolloc, tx, ty)
plt.grid(linestyle='-', linewidth=0.5)
show()
