#===============================================================================
import numpy as np
import math
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
#===============================================================================
# Prediction of the wave resistance of a bluff body
#
# Last review at 11-2-2020
# Written by: Anevlavi Dimitra
#===============================================================================

def constantStrengthSource(x1, y1, x2, y2, xp, yp):
    # The induced velocity and velocity potential due to SL panel
    # At the centre of the element
    # curvilinear coordinate system
    # from analytic formulas and sign convention
    uu = 0.;
    vv = -0.5;
    # local coordinate system
    xi1 = 0;
    xi2 = math.sqrt((x2-x1)**2+(y2-y1)**2);
    th0 = math.atan2(y2-y1,x2-x1);

    xm = 0.5*(x1+x2); ym = 0.5*(y1+y2);
    rcon=math.sqrt((xp-xm)**2+(yp-ym)**2)/xi2;

    #flag = (rcon>0.0001); % if flag=0 self-induced quantities
    flag = (rcon>1e-10);

    # Local coordinate system, aligned with panel j with
    # (0,0) positioned at the previous (x1,y1) node.
    xip= (xp-x1)*math.cos(th0)+(yp-y1)*math.sin(th0); #from global to jt
    yip=-(xp-x1)*math.sin(th0)+(yp-y1)*math.cos(th0);

    r1=math.sqrt((xp-x1)**2+(yp-y1)**2);
    r2=math.sqrt((xp-x2)**2+(yp-y2)**2);
    th1=math.atan2(yip,xip-xi1);
    th2=math.atan2(yip,xip-xi2);

    uu = flag*math.log(r1/r2)/2/math.pi + (flag-1)*uu;
    vv = flag*(th2-th1)/2/math.pi + (flag-1)*vv;

    # for source
    fi= (2*(xip-xi1)*math.log(r1)-2*(xip-xi2)*math.log(r2)- 2*(xi2-xi1) + 2*yip*(th2-th1))/(4*math.pi);
    #for sinks
    #fi= -(2*(xip-xi1)*math.log(r1)-2*(xip-xi2)*math.log(r2)- 2*xi2 + 2*yip*(th2-th1))/(4*math.pi);

    # back to the global coordinate system
    u = uu*math.cos(th0)-vv*math.sin(th0); #u
    v = uu*math.sin(th0)+vv*math.cos(th0); #v

    f = [fi, u, v]
    return f

def lineParametrization(xa, ya, xb, yb, Npanels):
    # line segment from point A to point B
    t = np.linspace(0,1,Npanels)
    xline = xa + t*(xb-xa)
    yline = ya + t*(yb-ya)
    return xline, yline

def collocationScheme(xi, yi):
    #normal and tangent unit vectors
    sizeN= len(xi)-1
    xcolloc = np.zeros(sizeN);ycolloc = np.zeros(sizeN)
    sin_theta = np.zeros(sizeN);cos_theta = np.zeros(sizeN)
    nx = np.zeros(sizeN);ny = np.zeros(sizeN)
    tx = np.zeros(sizeN);ty = np.zeros(sizeN)

    for ii in range(0,len(xcolloc)):
        xcolloc[ii] = 0.5*(xi[ii+1]+xi[ii])
        ycolloc[ii] = 0.5*(yi[ii+1]+yi[ii])

        lpanel = np.sqrt((xi[ii+1]-xi[ii])**2+(yi[ii+1]-yi[ii])**2) #panel length
        sin_theta[ii] = (yi[ii+1]-yi[ii])/lpanel
        cos_theta[ii] = (xi[ii+1]-xi[ii])/lpanel

        ## The normal vector outside the domain
        #nx[ii] = sin_theta[ii]
        #ny[ii] = -cos_theta[ii]
        #tx[ii] = -ny[ii]
        #ty[ii] = nx[ii]

        # The normal vector inside the domain
        nx[ii] = sin_theta[ii]
        ny[ii] = -cos_theta[ii]
        tx[ii] = -ny[ii]
        ty[ii] = nx[ii]

    return xcolloc, ycolloc, nx, ny, tx, ty

#===============================================================================
# CREATE MESH FOR SQUARE DOMAIN & BOYNDARY CONDITIONS
#===============================================================================
plotMesh = 1
debuggTool = 1

R = 0.5 #radius of cylinder
h = 2*R #cylinder submergence
Froude = 1.2
gi = 9.81
U = Froude*np.sqrt(gi*h)
kappa = gi/U**2

Np1  = 100 #number of nodes on free-surface
Np2  = 50 #number of nodes on cylinder

# free domain creation
xa, ya, xb, yb = 3, 0, -15, 0
[x1, y1] =lineParametrization(xa, ya, xb, yb, Np1)

# cylinder surface
theta = np.linspace(2*np.pi, 0, Np2)
dll = R*abs(theta[1]-theta[0]) #ds approximation
x2 = R*np.cos(theta)
y2 = R*np.sin(theta) - h

#nodes on the boundary with clockwise numbering
[xcolloc1, ycolloc1, nx1, ny1, tx1, ty1] = collocationScheme(x1, y1)
[xcolloc2, ycolloc2, nx2, ny2, tx2, ty2] = collocationScheme(x2, y2)

#make numpy array
xi = np.array([x1,x2]) #nodal values
yi = np.array([y1,y2]) #nodal values

if debuggTool == 1:
    print("First row of xi -> free surface x-nodes \n", xi[0])
    print("Second row of xi -> cylinder x-nodes \n", xi[1])
    print("\n", yi[0])
    print("\n", yi[1])

xcolloc = np.array([xcolloc1,xcolloc2])
ycolloc = np.array([ycolloc1,ycolloc2])
Nel = len(xcolloc[0][:]) + len(xcolloc[1][:])
Nfs = len(xcolloc[0][:]); Nbody = len(xcolloc[1][:])

if debuggTool == 1: print(Nel)

nx = np.array([nx1,nx2])
ny = np.array([ny1,ny2])
tx = np.array([tx1,tx2])
ty = np.array([ty1,ty2])

# plot the discretized domain
if (plotMesh==1):
    plt.plot(xi[0], yi[0], xi[1], yi[1])
    plt.plot(xi[0], yi[0], 'bo', xi[1], yi[1], 'ro')

    plt.grid(linestyle='-', linewidth=0.5)
    plt.plot(xcolloc[0], ycolloc[0], 'gs', xcolloc[1], ycolloc[1], 'gs')
    plt.quiver(xcolloc[0][0], ycolloc[0][0], nx[0][0], ny[0][0])
    plt.quiver(xcolloc[1][0], ycolloc[1][0], nx[1][0], ny[1][0])
    plt.quiver(xcolloc[0][0], ycolloc[0][0], tx[0][0], ty[0][0])
    plt.quiver(xcolloc[1][0], ycolloc[1][0], tx[1][0], ty[1][0])
    plt.grid(linestyle='-', linewidth=0.5)
    plt.title('Nodes, control points and orientation vectors')
    plt.axis("equal")
    plt.show()

#===============================================================================
# LINEAR SYSTEM OF EQUATIONS
#===============================================================================
Bv = np.zeros((Nel, 1)) #RHS: Bvector
Am = np.zeros(Nel) #LHS: Amatrix

kel = -1
for isec in range(2):
    #print("BOUNDARY:", isec)
    for i in range(len(xcolloc[isec][:])): # it runs through all the collocation points
        kel = kel + 1
        #print("index = ", i, " xi=", xcolloc[isec][i])
        #print("kel", kel)
        if isec == 1:
            Bv[kel][0] = U*nx[1][i]

if debuggTool == 1: print("RHS Matrix =>", Bv)

#===============================================================================
Potm = np.zeros((Nel, Nel))
xalfa = np.zeros((Nel, Nel))
yalfa = np.zeros((Nel, Nel))
un = np.zeros((Nel, Nel))
#print("Nel", Nel)
#print("un", un)
kel = -1
for isec in range(2): #for every boundary [0, 1]
    #print("BOUNDARY isec:", isec)
    for i in range(len(xcolloc[isec][:])): # for every panel on that boundary
        kel = kel + 1
        #print("index = ", i, " xi=", xcolloc[isec][i])
        #print("kel", kel)
        xp = xcolloc[isec][i];  yp = ycolloc[isec][i]
        jel = -1
        for jsec in range(2): #[0, 1]
            #print("BOUNDARY jsec:", jsec)
            for j in range(len(xcolloc[jsec][:])):
                jel = jel + 1
                #print("index = ", i, " xi=", xi[jsec][j])
                #print("jel", jel)
                f = constantStrengthSource(xi[jsec][j], yi[jsec][j], xi[jsec][j+1], yi[jsec][j+1], xp, yp)

                fi = f[0]
                up = f[1]
                vp = f[2]
                Potm[kel][jel] = fi
                xalfa[kel][jel] = up
                yalfa[kel][jel] = vp

                if kel == jel: #for self induced velocities expression on the boundary
                    # save values of matrices
                    xalfa[kel][jel] = -0.5*nx[isec][i] #sinks
                    yalfa[kel][jel]= -0.5*ny[isec][i]
                #print("matrix: ", kel, jel)

                un[kel][jel] = up*nx[isec][i] + vp*ny[isec][i]
                #print("RESULTS", un[kel][jel], xi[jsec][j], yi[jsec][j], xi[jsec][j+1], yi[jsec][j+1], xp, yp)

if debuggTool == 1: print("LHS Matrix =>", un)
Am = un

#===============================================================================
# DAWSON FINITE DIFFERENT SCHEME
#===============================================================================
delx = xcolloc[0][0]-xcolloc[0][1]
print("delx", delx)
for i in range(len(xcolloc[0][:])): # for every panel on the free surface boundary
    ir = len(xcolloc[0][:])-1 -i;
    irDawson = ir + 1
    print("i", i, "ir", ir)
    if irDawson>3:
        xp1 = xcolloc[0][ir]
        xp2 = xcolloc[0][ir-1]
        xp3 = xcolloc[0][ir-2]
        xp4 = xcolloc[0][ir-3]
    if irDawson==3:
        xp1 = xcolloc[0][ir]
        xp2 = xcolloc[0][ir-1]
        xp3 = xcolloc[0][ir-2]
        xp4 = xp3 - delx
        #xp1=xm(ir); xp2=xm(ir-1); xp3=xm(ir-2); xp4=xp3-delx;
    if irDawson==2:
        xp1 = xcolloc[0][ir]
        xp2 = xcolloc[0][ir-1]
        xp3 = xp2 - delx
        xp4 = xp3 - delx
        #xp1=xm(ir); xp2=xm(ir-1); xp3=xp2-delx; xp4=xp3-delx;
    if irDawson==1:
        xp1 = xcolloc[0][ir]
        xp2 = xp1 - delx
        xp3 = xp2 - delx
        xp4 = xp3 - delx
        #xp1=xm(ir); xp2=xp1-delx; xp3=xp2-delx; xp4=xp3-delx;
    #print("\n", xp1, xp2, xp3, xp4, irDawson)
    Di  =  (xp2-xp1)*(xp3-xp1)*(xp4-xp1)*(xp4-xp2)*(xp3-xp2)*(xp4-xp3)*(xp4+xp3+xp2-3*xp1);
    CDi =  ((xp2-xp1)**2)*((xp3-xp1)**2)*(xp3-xp2)*(xp3+xp2-2*xp1)/Di;
    CCi = -((xp2-xp1)**2)*((xp4-xp1)**2)*(xp4-xp2)*(xp4+xp2-2*xp1)/Di;
    CBi =  ((xp3-xp1)**2)*((xp4-xp1)**2)*(xp4-xp3)*(xp4+xp3-2*xp1)/Di;
    CAi =  -(CBi+CCi+CDi);
    p1  =  CAi; p2  =  CBi;
    p3  =  CCi; p4  =  CDi;
    #print(p1, p2, p3, p4)
    for mm in range(Nel):
        pol1=0;
        if irDawson>3:
            pol1 = p1*xalfa[ir][mm] + p2*xalfa[ir-1][mm] + p3*xalfa[ir-2][mm]+ p4*xalfa[ir-3][mm];
        Am[ir][mm] = Am[ir][mm] + pol1/kappa;

#print("Am", Am)
#print(Am[3][:])

#===============================================================================
# SOLUTION OF THE LINEAR SYSTEM
#===============================================================================
S = np.linalg.solve(Am,Bv)
#print("solution", S)

#===============================================================================
# POST-PROCESSING
#===============================================================================
Potx = -U + np.matmul(xalfa,S) #the defivative of potential in x
Poty = np.matmul(yalfa,S) #the defivative of potential in z direction
Cp = 1-(Potx**2+Poty**2)/U**2
eta = (U/gi)*(Potx[1:Nfs]+U) # free surface elevation
#print("eta", eta)

plt.plot( xcolloc[0][0:Nfs-1], eta, 'b--', label='heta')
plt.plot( xcolloc[0][:], ycolloc[0][:], 'k')
plt.plot( xcolloc[1][:], ycolloc[1][:])
plt.axis("equal")
plt.title("FIGURE.2 Free surface elevation")
plt.xlabel("x [m]")
plt.ylabel("z [m]")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
plt.grid(linestyle='-', linewidth=0.5)
plt.show()

# analytic solution for the pressure coefficient around the cylinder
#print(theta)
start = len(theta)
#print(theta[0:start-1])
#print(theta[1:start])
thi = 0.5*(theta[0:start-1] + theta[1:start])

uanal=U*(-1-np.cos(-math.pi+2*thi));
vanal=U*np.sin(-math.pi+2*thi);
Cpa=1-(uanal**2+vanal**2)/U**2;

#print(Nfs)
#print(Nbody)
#print(Cp)

plt.plot( xi[1][0:2], yi[1][0:2], 'bo')
plt.plot( xi[1][:], yi[1][:], 'k')
plt.show()

plt.plot(thi, Cp[Nfs:len(Cp)], 'bo', label="numerical")
plt.plot(thi, Cp[Nfs:len(Cp)], 'b')
plt.plot(thi, Cpa, label='analytic')
plt.title("FIGURE.1 Pressure coefficient around the cylinder")
plt.xlabel("theta (deg)")
plt.ylabel("Cp")
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
plt.grid(linestyle='-', linewidth=0.5)
plt.show()

# Calculate resistance
#?? Look up for dll
temp = np.multiply(Cp[Nfs:len(Cp)], nx[1][:])

print("1))) = ", Cp[Nfs:len(Cp)])
print("2))) = ", nx[1][:])
Resistance = dll*np.sum(temp)
print(Resistance)
#===============================================================================
