#### THE RIDDLE OF THE SIGN CONVENTION HAS BEEN SOLVED?



##### Katsikadelis - Boundary Elements and Applications (2002)

He uses Green's theorem to formulate the boundary integral equation for Dirichlet, Neumann and Robin interior boundary value problems. 

- The fundamental solution is the point source
- The normal vector points in the exterior domain and the parametrization of the boundary is counter-clockwise

$$
\upsilon = \frac{1}{2\pi}lnr
$$

In the book we have for the interior Dirichlet [see equation (3.32a), page 24]
$$
\frac{1}{2}u_{i} = - \int_{\partial \Omega} (\upsilon\frac{\partial u}{\partial n} -u\frac{\partial \upsilon}{\partial n} )ds
$$
or equivalently
$$
\frac{1}{2}u_{i} - \int_{\partial \Omega} u\frac{\partial \upsilon}{\partial n}ds =-\int_{\partial \Omega} \upsilon\frac{\partial u}{\partial n}ds
$$


##### Katz & Plotkin - Low Speed Aerodynamics (2001)

*Also the Lecture notes by Tara LaForce from Standford (2006) are in agreement with Katz & Plotkin. The Laplace tutorials for a bounded domain with various boundary conditions have been based on the following.*

The fundamental solution is the point source (but with the opposite sign)
$$
\Phi = -\frac{\sigma}{2\pi}lnr
$$
If we use the BIE from Katsikadelis and the formulas form Katz-Plotkin there is a mismatch with the signs. That is why it works for the tutorials to use
$$
\frac{1}{2}u_{i} + \int_{\partial \Omega} u\frac{\partial \Phi}{\partial n}ds =
\int_{\partial \Omega} \Phi\frac{\partial u}{\partial n}ds
$$


- The function that calculates the induced potential from a dipole is not affected by the normal vector, the velocities however are! 
- For Neumann problems typically, the boundary data are given with the normal vector pointing towards the exterior of the domain. So at hock, the sign can be altered.  