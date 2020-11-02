#### Laplace 2d Dirichlet boundary conditions

The fundamental solution for the Laplace equation considered is
$$
\Phi = -\frac{\sigma}{2\pi}lnr
$$


Collocation scheme implemented for a low-order panel method based on single-layer (source) and double-layer (dipole) potential. 

The boundary value problem with Dirichlet boundary conditions
$$
\Delta u=0
$$
has the following analytic solution
$$
u(x,y) = \sin x \sinh y - 0.25
$$
Boundary integral equation
$$
\int_{\partial D} \Phi\frac{\partial u}{\partial n} = \frac{1}{2}u_{i} + \int_{\partial D} u\frac{\partial \Phi}{\partial n}
$$
This equation gives the Neumann data with the normal vector directed towards the interior of the bounded domain.

*What is the difference between the exterior and the interior problem?*

The only difference is the normal vector. 

##### In order to validate the scheme, solve the lifting flow problem with sources and dipoles and plot the velocity potential and then validate with Hess n Smith. (OK) Thank tou theologos! 

Eάν έχω λύσει το πρόβλημα, χρησιμοποιώ το θεώρημα του green για να υπολογίσω το δυναμικό στο πεδίο. Πώς υπολογίζω όμως το πεδίο ταχυτήτων? Aριθμητικά ή μπορώ διαφορίζοντας τον τύπο του green.  (O βαγγέλης πάνω στο σύνορο υπολόγιζε τις ταχύτητες με πεπερασμένες διαφορές)

