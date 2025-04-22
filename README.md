OpenMM XTB Plugin
=================

This project provides a connection between [OpenMM](http://openmm.org) and [XTB](https://xtb-docs.readthedocs.io/en/latest).
It lets you compute forces and energies with the semi-empirical quantum chemistry methods provided by XTB.

This plugin requires XTB and OpenMM 8.1 or later.

Installation
============

The easiest is to install the plugin into a conda environment which has already installed `openmm` and `xtb(-python)`. 
Then the provided installation script can be used, which builds, tests and installs the plugin:

```bash
# Use default EDS_OpenMM environment
./build_test_install.sh

# Use a specific conda environment
./build_test_install.sh my_env

# Skip conda and specify custom paths
./build_test_install.sh --no-conda /path/to/openmm /path/to/xtb /install/path /path/to/python
```

Using The Plugin
================

The plugin can be used in two ways: by directly creating a force and adding it to your `System`,
or by using a `ForceField`.

Creating a Force Directly
-------------------------

This is the most general and flexible approach, giving complete control over all aspects.  For example, you can use
XTB to compute forces on some parts of the system but a classical force field to compute the forces on other parts.
To use this method, simply create a `XtbForce` object and add it to your `System`.  For example to create a force with no point charges,

```Python
from openmmxtb import XtbForce
system.addForce(XtbForce(XtbForce.GFN2xTB, 0.0, 1, False, particleIndices, atomicNumbers))
```

The arguments are as follows.

- `method`: the method to use for computing forces and energy.  The currently supported options are `GFN1xTB`, `GFN2xTB`, and `GFNFF`.
- `charge`: the total charge of the system
- `multiplicity`: the spin multiplicity
- `periodic`: whether to apply periodic boundary conditions
- `particleIndices`: a list containing the indices of the particles within the System to which the force should be applied.
   This allows you to use XTB for only part of the full system.
- `atomicNumbers`: a list containing the atomic numbers of the particles to which the force should be applied.  This must
   have the same length as `particleIndices`.  Element `i` is the atomic number of the particle specified by element
  `i` of `particleIndices`.

To create a force with point charges:

```Python
from openmmxtb import XtbForce
systm.addForce(XtbForce(XtbForce.GFN2xTB, 0.0, 1, False, particleIndices, atomicNumbers, point_charges, qm_region_indices, cutoff_radius)
```

The arguments are as above with the additional
- `point_charges`: a list of lists of `XtbPointCharge`s of all the point charges to consider. The list is grouped into 
    charge groups. I.e. element `i` of `point_charges` is a list of `XtbPointCharge`s which is one charge group
- `qm_region_indices`: indices of all particles treated in the qm zone. possibly by other forces as well. this is required
    to have consistent boundary region between different forces
- `cutoff`: the cutoff whether to include point charges in the boundary region. Currently charge group based cutoff is 
    implemented, i.e. if any particle of a charge group is within the cutoff, the whole charge group is used as external charges.

`XtbPointCharge` consists of index, atomic number and charge of a particle and can be created as
```Python
from openmmxtb import XtbPointCharge
pointCharge = XtbPointCharge(particleIndex, atomicNumber, particleCharge)
```

Using a ForceField
------------------

The above method requires you to create the `System` yourself, add particles to it, and build the lists of particle indices
and atomic numbers.  In simple cases where you want to use XTB without external point charges as the only force for the entire system, it can be easier
to let a `ForceField` create the `System` for you.

```python
ff = app.ForceField('xtb/gfn2xtb.xml')
system = ff.createSystem(topology, charge=2, multiplicity=1, periodic=False)
```

The plugin provides three force field files for the available methods: `'xtb/gfn1xtb.xml'`, `'xtb/gfn2xtb.xml'` and `'xtb/gfnff.xml'`.

License
=======

Portions copyright (c) 2023-2024 Stanford University and the Authors.

Authors: Peter Eastman

Contributors:

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.
