LinAcoustics_1D_DGalerkin
=========================

LinAcoustics_1D_DGalerkin is a discontinuous Galerkin (DG) solver for one-dimensional
linear system of acoustic equations:

    dt(U) + a*dx(U) + b*dx(P) = 0,
    dt(P) + a*dx(P) + b*dx(U) = 0.

LinAcoustics_1D_DGalerkin assembles the global DG matrices for the Legendre basis, and it does for a few explicit, implicit, and semi-implicit time integration schemes. More information on the DG matrices implemented is given in the reference preprint.

The DG scheme is given in the form

    v' = a/dx * Gs * v + b/dx * Gf * v

where b >> a > 0 and v=(u^T,p^T)^T are the degrees of freedom of (u,p).

Running the code
----------------

In file dg.m, one can change the time integrators, grid size, DG polynomial order. It is possible to output the global DG matrices and print their eigenvalues to study the DG spectra. Check the reference preprint for more information on DG matrices implemented here.

Running the code is just:

    octave dg.m

Reference
---------
```bibtex
@online{BrdarKnoth2024,
  author       = {Brdar, S. and Knoth, O.},
  title        = {On Spectrum of Discontinuous Galerkin Method for Linear Acoustic System},
  date         = {2024-08-27},
  langid       = {english},
  langidopts   = {variant=american},
  eprinttype   = {ResearchGate},
  doi          = {10.13140/RG.2.2.18583.79528}
}
```
