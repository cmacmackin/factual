project: FACTUAL
project_dir: ./src
output_dir: ./doc
author: Chris MacMackin
author_description: I am a graduate student at the University of Oxford, studying the melting and evolution of ice shelves. I enjoy programming, with my favourite languages being Fortran (for numerical work) and Python (for scripting and everything else).
website: https://cmacmackin.github.io
github: https://github.com/cmacmackin
email: cmacmackin@gmail.com
summary: Field Abstract Calculus Type Usable and Aesthetic Library: 
	 implementations of mathematical fields on which mathematical operations
	 can be performed.
project_github: https://github.com/cmacmackin/factual
display: public
         protected
	 private
graph: true
source: true
page_dir: pages

This library provides abstract types to represent mathematical fields
in Fortran. These are both scalar and vector fields. It also contains
(at present) one concrete implementation of these types, using a
pseudo-spectral approach and tracking field values at Chebyshev 
collocation points.

## License
FIAT is licensed under the GNU Lesser General Public License (LGPL) v3.0 or
later. The terms are provided in the file `LICENSE`. The LGPL make reference
to the GNU General Public License (GPL), which is provided in the file `GPL`.
In brief, the LGPL allows this library to be linked to software under any
license (with a few, minor, restrictions). However, should a modified version
of the _library itself_ be released, it must be licensed under the terms of
the LGPL or GPL.

