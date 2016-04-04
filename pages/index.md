Title: Using Fields
Author: Chris MacMackin
Date: April 2016

The field data types (descendents of [[abstract_field]]) behave much like
mathematical fields and provide many of the same behaviours. These include:

- overloaded intrinsic mathematical operations
- scalar and vector fields
- calculus functions, including vector calculus

For each implementation of fields, which will be contained in its own module,
there will be a vector and a scalar version.
There may also be uniform vector and scalar fields, special cases for 
efficiently handling fields with uniform values. In general, a field can
interact with any field of its own type, provided they have the same domain
and resolution. It is possible for a vector field to be multiplied or divided
by a scalar field, provided they have the same implementation.

"Implementation" refers to the discretization strategy used for that field.
Examples could include finite difference, finite volume, finite element, 
Chebyshev pseudo-spectral, and Fourier spectral implementations. Strategies
may each have different numbers of dimensions, or special boundary conditions.
At the present time, only the Chebyshev implementation is available.



