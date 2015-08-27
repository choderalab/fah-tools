#!/bin/env python

"""
Set the PME parameters explicitly in system.xml file

Usage:

> python set-pme-parameters-in-system.py system.xml

Old file is backed up to system.xml.old
New file is created as system.xml

"""

import os
import sys
from simtk import openmm, unit
from simtk.openmm import XmlSerializer

#
# SUBROUTINES
#

def calc_pme_parameters(system):
    """Calculate PME parameters using scheme similar to OpenMM OpenCL platform.

    Parameters
    ----------
    system : simtk.openmm.System
        The system for which parameters are to be computed.

    Returns
    -------
    alpha : float
        The PME alpha parameter
    nx, ny, nz : int
        The grid numbers in each dimension

    """

    # Find nonbonded force.
    forces = { system.getForce(index).__class__.__name__ : system.getForce(index) for index in range(system.getNumForces()) }
    force = forces['NonbondedForce']
    tol = force.getEwaldErrorTolerance()
    boxVectors = system.getDefaultPeriodicBoxVectors()

    from numpy import sqrt, log, ceil
    from math import pow
    alpha = (1.0/force.getCutoffDistance())*sqrt(-log(2.0*tol))
    xsize = int(ceil(2*alpha*boxVectors[0][0]/(3*pow(tol, 0.2))))
    ysize = int(ceil(2*alpha*boxVectors[1][1]/(3*pow(tol, 0.2))))
    zsize = int(ceil(2*alpha*boxVectors[2][2]/(3*pow(tol, 0.2))))

    def findLegalDimension(minimum):
        while (True):
            # Attempt to factor the current value.
            unfactored = int(minimum)
            for factor in range(2, 8):
                while (unfactored > 1) and (unfactored%factor == 0):
                    unfactored /= factor

            if (unfactored == 1):
                return int(minimum)

            minimum += 1

    nx = findLegalDimension(xsize)
    ny = findLegalDimension(ysize)
    nz = findLegalDimension(zsize)

    return (alpha, nx, ny, nz)

def read_file(filename):
    infile = open(filename, 'r')
    contents = infile.read()
    infile.close()
    return contents

def write_file(filename, contents):
    outfile = open(filename, 'w')
    outfile.write(contents)
    outfile.close()
    return
#
# MAIN
#

def fix_system(system_xml_filename):
    """
    Set the PME parameters explicitly in a specified system XML file if they are not already set.

    The file is renamed with '.old' appended, and a corrected file written in its place.

    Parameters
    ----------
    system_xml_filename : str
        The name of the serialized system XML file to be modified

    """

    system = XmlSerializer.deserialize(read_file(system_xml_filename))

    forces = { system.getForce(force_index).__class__.__name__ : system.getForce(force_index) for force_index in range(system.getNumForces()) }
    force = forces['NonbondedForce']
    (alpha, nx, ny, nz) = force.getPMEParameters()
    if alpha == 0.0 / unit.nanometers:
        (alpha, nx, ny, nz) = calc_pme_parameters(system)
        force.setPMEParameters(alpha, nx, ny, nz)
        serialized_system = XmlSerializer.serialize(system)
        os.rename(system_xml_filename, system_xml_filename + '.old')
        write_file(system_xml_filename, serialized_system)

    return

if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise Exception('usage: python set-pme-parameters-in-system.py system.xml')
    fix_system(sys.argv[1])
