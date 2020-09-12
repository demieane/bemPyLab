import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt

# solution of Laplace in 2d with neumman boundary conditions
# with bem

def funcSource(x, y, xo, yo):
    Phi = np.log(math.sqrt((x-xo)**2 + (y-yo)**2))
    return Phi

def constantStrengthSource(x1, y1, x2, y2, xp, yp):
    # The induced velocity and velocity potential due to SL panel
    # At the centre of the element
    # curvilinear coordinate system
    #uu = 0.0  #self-induced velocity in x-axis
    #vv = -0.5 #self-induced velocity in y-axis
    xm=0.5*(x1+x2);
    ym=0.5*(y1+y2);

    xi1 = 0
    xi2 = math.sqrt((x2-x1)**2+(y2-y1)**2)

    fi = 2*xi2*math.log(0.5*xi2)/4/math.pi #potential
    th0 = math.atan2(y2-y1,x2-x1)
    rcon = math.sqrt((xp-xm)**2+(yp-ym)**2)/xi2

    if (rcon>0.0001):
        xip =  (xp-x1)*math.cos(th0)+(yp-y1)*math.sin(th0)
        yip = -(xp-x1)*math.sin(th0)+(yp-y1)*math.cos(th0)

        r1 = math.sqrt((xp-x1)**2+(yp-y1)**2)
        r2 = math.sqrt((xp-x2)**2+(yp-y2)**2)
        th1 = math.atan2(yip,xip-xi1)
        th2 = math.atan2(yip,xip-xi2)

        #uu = math.log(r1/r2)/2/math.pi
        #vv = (th2-th1)/2/math.pi

        fi = (2*(xip-xi1)*math.log(r1) - 2*(xip-xi2)*math.log(r2) -2*(x2-x1) + 2*yip*(th2-th1))/4/math.pi

    # global coordinate system

    #f1 = uu*math.cos(th0)-vv*math.sin(th0) #u
    #f2 = uu*math.sin(th0)+vv*math.cos(th0) #v
    #f3 = fi; #potential
    f = fi#[f1, f2, f3]
    return f

def constantStrengthDoublet(x1, y1, x2, y2, xp, yp):
    xm=0.5*(x1+x2)
    ym=0.5*(y1+y2)
    xi1=0.0
    xi2=math.sqrt((x2-x1)**2+(y2-y1)**2)
    th0=math.atan2(y2-y1,x2-x1)

    xip =  (xp-x1)*math.cos(th0)+(yp-y1)*math.sin(th0)
    yip = -(xp-x1)*math.sin(th0)+(yp-y1)*math.cos(th0)
    fi=0.0
    rcon=math.sqrt((xp-xm)**2+(yp-ym)**2)/xi2;
    if (rcon>0.0001):
        th1=math.atan2(yip,xip-xi1);
        th2=math.atan2(yip,xip-xi2);
        fi=-(th2-th1)/2/math.pi;
    return fi

def lineParametrization(xa, ya, xb, yb, Npanels):
    # from A to B
    t = np.linspace(0,1,Npanels)
    #print(t)
    xline = xa + t*(xb-xa)
    yline = ya + t*(yb-ya)
    return xline, yline

def funcSourceDerivative(nx, ny, x, y, xo, yo):
    # Neumann BCs
    rr = ((x-xo)**2-(y-yo)**2)
    dPhidn = nx*(x-xo)/rr+ny*(y-yo)/rr
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
        #"""
        nx[ii] = sin_theta[ii]
        ny[ii] = -cos_theta[ii]
        tx[ii] = -ny[ii]
        ty[ii] = nx[ii]
        """
        nx[ii] = -sin_theta[ii]
        ny[ii] = cos_theta[ii]
        tx[ii] = ny[ii]
        ty[ii] = -nx[ii]
        """
    return xcolloc, ycolloc, nx, ny, tx, ty

def show():
   return matplotlib.pyplot.show(block=True)

#===========================MAIN CODE BELOW====================================#

#==============CREATE MESH=====================================================#
Np = 10 #number of panels
xo, yo = 2,0 #source position
# square domain creation
xa, ya = 1, 0
xb, yb = 0, -1
xc, yc = -1, 0
xd, yd = 0, 1
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
#print(xi)
#print(yi)

[xcolloc, ycolloc, nx, ny, tx, ty] = collocationScheme(xi, yi)

plt.plot(xcolloc, ycolloc, 'gs')
plt.quiver(xcolloc, ycolloc, nx, ny)
plt.quiver(xcolloc, ycolloc, tx, ty)
plt.grid(linestyle='-', linewidth=0.5)
show()

# analytic solution
dirichletBC = np.zeros(len(xcolloc))
neumannBC = np.zeros(len(xcolloc))

for jj in range(len(xcolloc)):
    dirichletBC[jj] = funcSource(xcolloc[jj], ycolloc[jj], xo, yo)
    neumannBC[jj] = funcSourceDerivative(nx[jj], ny[jj], xcolloc[jj], ycolloc[jj], xo, yo)

xindex = np.linspace(0, 1, len(xcolloc))
plt.plot(xindex, dirichletBC)
plt.plot(xindex, dirichletBC, 'bo', label='dirichlet data')
plt.plot(xindex, neumannBC)
plt.plot(xindex, neumannBC,'ko', label='neumann data')
strLabel = 'Collocation points counter-clockwise indexing from (' + str(xa) + ',' + str(ya) +')'
plt.xlabel(strLabel);
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
plt.grid(linestyle='-', linewidth=0.5)
show()

"""
# test the constant strength source elements
x1, y1 = 0, 0# = xi[0];y1=yi[0];
x2, y2 = 0, 1# = xi[1];y2=yi[1];

xp, yp = 0, 0.5# = xcolloc[0]; yp = ycolloc[0];
ftest = constantStrengthSource(x1, y1, x2, y2, xp, yp)
print(ftest)
"""
sizeA = len(xcolloc)
Hij = np.zeros((sizeA, sizeA))
Gij = np.zeros((sizeA, sizeA))
dudn = np.zeros((sizeA,1))
#print(len(neumannBC))

for ii in range(0, sizeA): #for each collocation point
    #sumGij = 0
    xp, yp = xcolloc[ii], ycolloc[ii]
    for jj in range(0, sizeA): #we evaluate the effect of source distribution on the panels
        # [u, v, fi]
        x1, y1 = xi[jj], yi[jj]
        x2, y2 = xi[jj+1], yi[jj+1]
        fi_source = constantStrengthSource(x1, y1, x2, y2, xp, yp) #w
        fi_doublet = constantStrengthDoublet(x1, y1, x2, y2, xp, yp)#dw/dn

        Hij[ii][jj] = fi_doublet
        if (ii==jj):
            Hij[ii][jj] = Hij[ii][jj] - 0.5

        Gij[ii][jj] = fi_source

    dudn[ii][0] = neumannBC[ii]

print(Hij)
print(Gij)
b = np.matmul(Gij, dudn)
ubem = np.linalg.solve(Hij,b)

# comparison
xindex = np.linspace(0, 1, len(xcolloc))
plt.plot(xindex, dirichletBC)
plt.plot(xindex, dirichletBC, 'bo', label='dirichlet data')
plt.plot(xindex, ubem, 'r*')
#plt.plot(xindex, neumannBC)
#plt.plot(xindex, neumannBC,'ko', label='neumann data')
strLabel = 'Collocation points counter-clockwise indexing from (' + str(xa) + ',' + str(ya) +')'
plt.xlabel(strLabel);
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
plt.grid(linestyle='-', linewidth=0.5)
show()
