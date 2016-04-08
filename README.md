FACTUAL: Field Abstract Calculus Type Usable and Aesthetic Library
==================================================================
This library provides abstract types to represent mathematical fields
in Fortran. These are both scalar and vector fields. It also contains
(at present) one concrete implementation of these types, using a
pseudo-spectral approach and tracking field values at Chebyshev 
collocation points.

## (Potential) Features

- [x] abstract types representing vector and scalar fields
- [ ] concrete implementations of said types
  - [ ] 1D implementation using Chebyshev pseudo-spectral techniques
        (in progress)
  - [ ] 2D implementation using Chebyshev pseudo-spectral techniques
  - [ ] 3D implementation using Chebyshev pseudo-spectral techniques
  - [ ] implementations of special "uniform" fields which provide optimal
        use of memory and CPU
- [ ] implementation with periodic boundary conditions (Fourier spectral
      approach?)
- [ ] Factory type or method to provide uniform constructor for implementations
- [ ] Front-end nonlinear solver for fields?
- [ ] Parallel implementations (for multidimensional fields)
- [ ] Derived type to represent boundary conditions, with procedure pointer to
      generate said conditions?
- [x] Equality checking of fields, within a tolerance
- [ ] Ability to set tolerance on a per-field basis

## License
FACTUAL is licensed under the GNU General Public License (GPL) v3.0 or
later. The terms are provided in the file `LICENSE`. The Lesser General
Public License (LGPL) would have been used, but the Chebyshev pseudo-spectral
implementation uses the FFTW3 library, which is licensed under the GPL.

