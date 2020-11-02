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
debuggTool = 0

R = 0.5 #radius of cylinder
h = 2*R #cylinder submergence
Froude = 1.2
gi = 9.81
U = Froude*np.sqrt(gi*h)
kappa = gi/U**2

Np1  = 10 #number of nodes on free-surface
Np2  = 5 #number of nodes on cylinder

# free domain creation
xa, ya, xb, yb = 3, 0, -15, 0
[x1, y1] =lineParametrization(xa, ya, xb, yb, Np1)

# cylinder surface
theta = np.linspace(2*np.pi, 0, Np2)
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
    print("\n", xp1, xp2, xp3, xp4, irDawson)
    Di  =  (xp2-xp1)*(xp3-xp1)*(xp4-xp1)*(xp4-xp2)*(xp3-xp2)*(xp4-xp3)*(xp4+xp3+xp2-3*xp1);
    CDi =  ((xp2-xp1)**2)*((xp3-xp1)**2)*(xp3-xp2)*(xp3+xp2-2*xp1)/Di;
    CCi = -((xp2-xp1)**2)*((xp4-xp1)**2)*(xp4-xp2)*(xp4+xp2-2*xp1)/Di;
    CBi =  ((xp3-xp1)**2)*((xp4-xp1)**2)*(xp4-xp3)*(xp4+xp3-2*xp1)/Di;
    CAi =  -(CBi+CCi+CDi);
    p1  =  CAi; p2  =  CBi;
    p3  =  CCi; p4  =  CDi;
    print(p1, p2, p3, p4)
    for mm in range(Nel):
        pol1=0;
        if ir>3:
            pol1 = p1*xalfa[ir][mm] + p2*xalfa[ir-1][mm] + p3*xalfa[ir-2][mm]+ p4*xalfa[ir-3][mm];
        Am[ir][mm] = Am[ir][mm] + pol1/kappa;

#print("Am", Am)
print(Am[3][:])







#===============================================================================
# SOLUTION OF THE LINEAR SYSTEM
#===============================================================================
S = np.linalg.solve(Am,Bv)
print("solution", S)
"""
kel=0;
for isec = 1:2 % for each boundary

for i = 1:N(isec) % for every panel on that boundary
kel = kel+1;
% collocation point
xp = xm(kel);
yp = ym(kel);

jel=0;
% calculate the effect from constant-source distributions
for jsec = 1:2 % from the other boundaries and from itself

    for j = 1:N(jsec) % panel-wise
        jel = jel+1;
        xi1 = xi(jsec,j);
        yi1 = yi(jsec,j);
        xi2 = xi(jsec,j+1);
        yi2 = yi(jsec,j+1);
        % from sources
            f = SL(xi1,yi1,xi2,yi2,xp,yp);
            up = f(1);
            vp = f(2);
            % save values of matrices
            Potm(kel,jel) = f(3); %from the integral expressions
            xalfa(kel,jel) = f(1); %velocities due to sinks U
            yalfa(kel,jel) = f(2); %velocities due to sinks V

        if kel == jel % for self induced velocities expression on the boundary
            % save values of matrices
            xalfa(kel,jel) = -0.5*npx(kel);%sinks
            yalfa(kel,jel) = -0.5*npy(kel);
        end

        un(kel,jel) = up*npx(kel) + vp*npy(kel);

    end;
end;

end;
end

"""





"""
if (plotBC == 1):
    plt.plot(xindex, u_bc, xindex, u_bc, 'bo')
    plt.title('Dirichlet boundary conditions')
    plt.show()

#===============================================================================
# ANALYTIC SOLUTION PREVIEW
if (plotAnalyticSolution==1):
    xanal = np.linspace(xa, xb, 100)
    yanal = np.linspace(ya, yc, 100)
    xx, yy = np.meshgrid(xanal, yanal)
    uanal = np.sinh(yy)*np.sin(xx)/math.sinh(1)
    plt.contourf(xx,yy,uanal)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(xx, yy, uanal, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    plt.show()

#===============================================================================
# BOUNDARY ELEMENT METHOD (BEM) BOYNDARY INTEGRAL EQUATION IN MATRIX FORM
sizeA = len(xcolloc)
Abij = np.zeros((sizeA, sizeA)) #G
Bbij = np.zeros((sizeA, sizeA)) #dGdn
#udirichlet = np.zeros((sizeA,1))
A = np.zeros((sizeA, sizeA))
S = np.zeros((sizeA, sizeA))

for ii in range(0, sizeA): #for each collocation point
    for jj in range(0, sizeA): #we evaluate the effect each panel
        # on each panel we have a piecewise-constat distribution of
        # sources (G) and dipoles (dG/dn)
        xp, yp = xcolloc[ii], ycolloc[ii]
        x1, y1 = xi[jj], yi[jj]
        x2, y2 = xi[jj+1], yi[jj+1]
        Abij[ii][jj] = constantStrengthSource(x1, y1, x2, y2, xp, yp)
        Bbij[ii][jj] = constantStrengthDoublet(x1, y1, x2, y2, xp, yp)
        deltaKronecker = 0
        if (ii == jj):
            deltaKronecker = 1
        A[ii][jj] = 0.5*deltaKronecker + Bbij[ii][jj]
        S[ii][jj] = Abij[ii][jj]

#===============================================================================
# SOLUTION OF THE LINEAR SYSTEM
b = np.matmul(A, u_bc)
dudn_bc = np.linalg.solve(S,b)
#===============================================================================
# COMPARISON WITH THE ANALYTIC SOLUTION
if (plotNumResultsBC==1):
    xindex = np.linspace(0, 1, len(xcolloc))
    plt.plot(xindex, u_bc)
    plt.plot(xindex, u_bc, 'bo', label='dirichlet data')
    plt.plot(xindex, dudn_bc, 'r*', label='neumann data')
    strLabel = 'Collocation points counter-clockwise indexing from (' + str(xa) + ',' + str(ya) +')'
    plt.xlabel(strLabel);
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper right', borderaxespad=0.)
    plt.grid(linestyle='-', linewidth=0.5)
    plt.show()
#===============================================================================
# DOMAIN SOLUTION u(x,y)
xtest = np.linspace(xa+0.01, xb-0.01, 20)
ytest = np.linspace(ya+0.01, yc-0.01, 20)
xnum, ynum = np.meshgrid(xtest, ytest)
uactual = np.sinh(ynum)*np.sin(xnum)/math.sinh(1)

unum = np.zeros((len(xtest), len(ytest)))

for kk in range(0, len(xtest)):
    for ll in range(0, len(ytest)):
        # for each point in the meshgrid
        utemporary = 0
        xp = xtest[ll]; yp = ytest[kk];
        for jj in range(0, len(xcolloc)):
            x1, y1 = xi[jj], yi[jj]
            x2, y2 = xi[jj+1], yi[jj+1]
            SL = constantStrengthSource(x1, y1, x2, y2, xp, yp) # (G)
            DL = constantStrengthDoublet(x1, y1, x2, y2, xp, yp) #(dG/dn)
            utemporary = utemporary + dudn_bc[jj]*SL - u_bc[jj]*DL

        unum[kk,ll]=utemporary

if (plotAnalyticSolution==1):
    plt.contourf(xtest,ytest,unum)
    plt.grid(linestyle='-', linewidth=0.5)
    plt.xlabel('x-axis')
    plt.ylabel('y-axis')
    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(xnum, ynum, unum, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    plt.show()

totalError = np.max(np.abs(unum - uactual))
message1 = 'Number of elements = ' + str(4*Np) + ', maxError = ' + str(totalError)
print(message1)
"""
#===============================================================================
