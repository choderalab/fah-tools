#!/usr/local/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Test all test systems on different platforms to ensure differences in potential energy and
forces are small among platforms.

DESCRIPTION

COPYRIGHT

@author John D. Chodera <jchodera@gmail.com>

All code in this repository is released under the GNU General Public License.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

TODO

"""

#=============================================================================================
# PYTHON 3 COMPATIBILITY CRAP
#=============================================================================================

from __future__ import print_function

#=============================================================================================
# ENABLE LOGGING
#=============================================================================================

import logging
logger = logging.getLogger(__name__)

def config_root_logger(verbose, log_file_path=None, mpicomm=None):
    """Setup the the root logger's configuration.
     The log messages are printed in the terminal and saved in the file specified
     by log_file_path (if not None) and printed. Note that logging use sys.stdout
     to print logging.INFO messages, and stderr for the others. The root logger's
     configuration is inherited by the loggers created by logging.getLogger(name).
     Different formats are used to display messages on the terminal and on the log
     file. For example, in the log file every entry has a timestamp which does not
     appear in the terminal. Moreover, the log file always shows the module that
     generate the message, while in the terminal this happens only for messages
     of level WARNING and higher.
    Parameters
    ----------
    verbose : bool
        Control the verbosity of the messages printed in the terminal. The logger
        displays messages of level logging.INFO and higher when verbose=False.
        Otherwise those of level logging.DEBUG and higher are printed.
    log_file_path : str, optional, default = None
        If not None, this is the path where all the logger's messages of level
        logging.DEBUG or higher are saved.
    mpicomm : mpi4py.MPI.COMM communicator, optional, default=None
        If specified, this communicator will be used to determine node rank.

    """

    class TerminalFormatter(logging.Formatter):
        """
        Simplified format for INFO and DEBUG level log messages.
        This allows to keep the logging.info() and debug() format separated from
        the other levels where more information may be needed. For example, for
        warning and error messages it is convenient to know also the module that
        generates them.
        """

        # This is the cleanest way I found to make the code compatible with both
        # Python 2 and Python 3
        simple_fmt = logging.Formatter('%(message)s')
        default_fmt = logging.Formatter('%(levelname)s - %(name)s - %(message)s')

        def format(self, record):
            if record.levelno <= logging.INFO:
                return self.simple_fmt.format(record)
            else:
                return self.default_fmt.format(record)

    # Check if root logger is already configured
    n_handlers = len(logging.root.handlers)
    if n_handlers > 0:
        root_logger = logging.root
        for i in xrange(n_handlers):
            root_logger.removeHandler(root_logger.handlers[0])

    # If this is a worker node, don't save any log file
    if mpicomm:
        rank = mpicomm.rank
    else:
        rank = 0

    if rank != 0:
        log_file_path = None

    # Add handler for stdout and stderr messages
    terminal_handler = logging.StreamHandler()
    terminal_handler.setFormatter(TerminalFormatter())
    if rank != 0:
        terminal_handler.setLevel(logging.WARNING)
    elif verbose:
        terminal_handler.setLevel(logging.DEBUG)
    else:
        terminal_handler.setLevel(logging.INFO)
    logging.root.addHandler(terminal_handler)

    # Add file handler to root logger
    if log_file_path is not None:
        #file_format = '%(asctime)s - %(levelname)s - %(name)s - %(message)s'
        file_format = '%(asctime)s: %(message)s'
        file_handler = logging.FileHandler(log_file_path)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(file_format))
        logging.root.addHandler(file_handler)

    # Do not handle logging.DEBUG at all if unnecessary
    if log_file_path is not None:
        logging.root.setLevel(logging.DEBUG)
    else:
        logging.root.setLevel(terminal_handler.level)

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
import os.path
import sys
import math
import glob

from simtk import unit, openmm

from simtk.openmm import XmlSerializer

#=============================================================================================
# SUBROUTINES
#=============================================================================================

# These settings control what tolerance is allowed between platforms and the Reference platform.
ENERGY_TOLERANCE = 0.06*unit.kilocalories_per_mole # energy difference tolerance
FORCE_RMSE_TOLERANCE = 0.06*unit.kilocalories_per_mole/unit.angstrom # per-particle force root-mean-square error tolerance

def assert_approximately_equal(computed_potential, expected_potential, tolerance=ENERGY_TOLERANCE):
    """
    Check whether computed potential is acceptably close to expected value, using an error tolerance.

    ARGUMENTS

    computed_potential (simtk.unit.Quantity in unit of energy) - computed potential energy
    expected_potential (simtk.unit.Quantity in unit of energy) - expected

    OPTIONAL ARGUMENTS

    tolerance (simtk.unit.Quantity in unit of energy) - acceptable tolerance

    EXAMPLES

    >>> assert_approximately_equal(0.0000 * unit.kilocalories_per_mole, 0.0001 * unit.kilocalories_per_mole, tolerance=0.06*unit.kilocalories_per_mole)

    """

    # Compute error.
    error = (computed_potential - expected_potential)

    # Raise an exception if the error is larger than the tolerance.
    if abs(error) > tolerance:
        raise Exception("Computed potential %s, expected %s.  Error %s is larger than acceptable tolerance of %s." % (computed_potential, expected_potential, error, tolerance))

    return

def compute_potential_and_force(system, positions, platform):
    """
    Compute the energy and force for the given system and positions in the designated platform.

    ARGUMENTS

    system (simtk.openmm.System) - the system for which the energy is to be computed
    positions (simtk.unit.Quantity of Nx3 numpy.array in unit of distance) - positions for which energy and force are to be computed
    platform (simtk.openmm.Platform) - platform object to be used to compute the energy and force

    RETURNS

    potential (simtk.unit.Quantity in energy/mole) - the potential
    force (simtk.unit.Quantity of Nx3 numpy.array in unit of energy/mole/distance) - the force

    """

    # Create a Context.
    kB = unit.BOLTZMANN_CONSTANT_kB
    temperature = 298.0 * unit.kelvin
    kT = kB * temperature
    beta = 1.0 / kT
    collision_rate = 90.0 / unit.picosecond
    timestep = 1.0 * unit.femtosecond
    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    context = openmm.Context(system, integrator, platform)
    # Set positions
    context.setPositions(positions)
    # Evaluate the potential energy.
    state = context.getState(getEnergy=True, getForces=True)
    potential = state.getPotentialEnergy()
    force = state.getForces(asNumpy=True)

    return [potential, force]

def compute_potential_and_force_by_force_index(system, positions, platform, force_index):
    """
    Compute the energy and force for the given system and positions in the designated platform for the given force index.

    ARGUMENTS

    system (simtk.openmm.System) - the system for which the energy is to be computed
    positions (simtk.unit.Quantity of Nx3 numpy.array in unit of distance) - positions for which energy and force are to be computed
    platform (simtk.openmm.Platform) - platform object to be used to compute the energy and force
    force_index (int) - index of force to be computed (all others ignored)

    RETURNS

    potential (simtk.unit.Quantity in energy/mole) - the potential
    force (simtk.unit.Quantity of Nx3 numpy.array in unit of energy/mole/distance) - the force

    """

    forces = [ system.getForce(index) for index in range(system.getNumForces()) ]

    # Get original force groups.
    groups = [ force.getForceGroup() for force in forces ]

    # Set force groups so only specified force_index contributes.
    for force in forces:
        force.setForceGroup(1)
    forces[force_index].setForceGroup(0) # bitmask of 1 should select only desired force

    # Create a Context.
    kB = unit.BOLTZMANN_CONSTANT_kB
    temperature = 298.0 * unit.kelvin
    kT = kB * temperature
    beta = 1.0 / kT
    collision_rate = 90.0 / unit.picosecond
    timestep = 1.0 * unit.femtosecond
    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    context = openmm.Context(system, integrator, platform)
    # Set positions
    context.setPositions(positions)
    # Evaluate the potential energy.
    state = context.getState(getEnergy=True, getForces=True, groups=1)
    potential = state.getPotentialEnergy()
    force = state.getForces(asNumpy=True)

    # Restore original force groups.
    for index in range(system.getNumForces()):
        forces[index].setForceGroup(groups[index])

    return [potential, force]

def compute_potential_and_force_by_force_group(system, positions, platform, force_group):
    """
    Compute the energy and force for the given system and positions in the designated platform for the given force group.

    ARGUMENTS

    system (simtk.openmm.System) - the system for which the energy is to be computed
    positions (simtk.unit.Quantity of Nx3 numpy.array in unit of distance) - positions for which energy and force are to be computed
    platform (simtk.openmm.Platform) - platform object to be used to compute the energy and force
    force_group (int) - index of force group to be computed (all others ignored)

    RETURNS

    potential (simtk.unit.Quantity in energy/mole) - the potential
    force (simtk.unit.Quantity of Nx3 numpy.array in unit of energy/mole/distance) - the force

    """

    forces = [ system.getForce(index) for index in range(system.getNumForces()) ]

    # Create a Context.
    kB = unit.BOLTZMANN_CONSTANT_kB
    temperature = 298.0 * unit.kelvin
    kT = kB * temperature
    beta = 1.0 / kT
    collision_rate = 90.0 / unit.picosecond
    timestep = 1.0 * unit.femtosecond
    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    context = openmm.Context(system, integrator, platform)
    # Set positions
    context.setPositions(positions)
    # Evaluate the potential energy.
    groupmask = 1 << (force_group + 1)
    state = context.getState(getEnergy=True, getForces=True, groups=groupmask)
    potential = state.getPotentialEnergy()
    force = state.getForces(asNumpy=True)

    return [potential, force]

def get_all_subclasses(cls):
    """
    Return all subclasses of a specified class.

    Parameters
    ----------
    cls : class
       The class for which all subclasses are to be returned.

    Returns
    -------
    all_subclasses : list of class
       List of all subclasses of `cls`.

    """
       
    all_subclasses = []

    for subclass in cls.__subclasses__():
        all_subclasses.append(subclass)
        all_subclasses.extend(get_all_subclasses(subclass))

    return all_subclasses
    
def get_num_runs(input_data_path):
    """Gets the number of runs from the test data and returns them"""
    runs=glob.glob(os.path.join(str(input_data_path), "RUN*"))
    n_runs= len(runs)
        
    return n_runs

def read_file(filename):
    infile = open(filename, 'r')
    contents = infile.read()
    infile.close()
    return contents


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

    print (xsize,ysize,zsize)
    def findLegalDimension(minimum):
        while (True):
            # Attempt to factor the current value.
            unfactored = minimum
            for factor in range(2, 8):
                while (unfactored > 1) and (unfactored%factor == 0):
                    unfactored /= factor
        
            if (unfactored == 1):
                return minimum

            minimum += 1
    
    nx = findLegalDimension(xsize)
    ny = findLegalDimension(ysize)
    nz = findLegalDimension(zsize)

    return (alpha, nx, ny, nz)

#=============================================================================================
# MAIN AND TESTS
#=============================================================================================

def main():
    import doctest
    import argparse

    parser = argparse.ArgumentParser(description="Check OpenMM computed energies and forces across all platforms for a suite of test systems.")
    parser.add_argument('-o', '--outfile', dest='logfile', action='store', type=str, default=None)
    parser.add_argument('-v', dest='verbose', action='store_true')
    parser.add_argument('-i', '--input', dest="input_data_path", action="store", type=str)
    parser.add_argument('-t', '--tuneplatform', dest="tune_pme_platform", action="store", type=str, default=None)
    parser.add_argument('-p', '--precision', dest="precision", action="store", type=str, default='single')
    args = parser.parse_args()

    verbose = args.verbose # Don't display extra debug information.
    config_root_logger(verbose, log_file_path=args.logfile)

    # Print version.
    logger.info("OpenMM version: %s" % openmm.version.version)
    logger.info("")

    # List all available platforms
    logger.info("Available platforms:")
    for platform_index in range(openmm.Platform.getNumPlatforms()):
        platform = openmm.Platform.getPlatform(platform_index)
        logger.info("%5d %s" % (platform_index, platform.getName()))
    logger.info("")

    # Test all systems on Reference platform.
    platform = openmm.Platform.getPlatformByName("Reference")
    print('Testing Reference platform...')
    doctest.testmod()

    # Compute energy error made on all test systems for other platforms.
    # Make a count of how often set tolerance is exceeded.
    tests_failed = 0 # number of times tolerance is exceeded
    tests_passed = 0 # number of times tolerance is not exceeded
    logger.info("%16s%16s %16s          %16s          %16s          %16s" % ("platform", "precision", "potential", "error", "force mag", "rms error"))
    reference_platform = openmm.Platform.getPlatformByName("Reference")
    n_runs=get_num_runs(args.input_data_path)
    for run in range(n_runs):
        print("Deserializing XML files for RUN%d" % run)
        state = XmlSerializer.deserialize(read_file(os.path.join(args.input_data_path,"RUN%d" % run, "state0.xml")))
        integrator = XmlSerializer.deserialize(read_file(os.path.join(args.input_data_path,"RUN%d" % run, "integrator.xml")))
        system = XmlSerializer.deserialize(read_file(os.path.join(args.input_data_path,"RUN%d" % run, "system.xml")))
        
        # Update system periodic box vectors based on state.
        system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())

        # Create test system instance.
        positions = state.getPositions()

        # Get PME parameters
        forces = [ system.getForce(force_index) for force_index in range(system.getNumForces()) ]
        force_dict = { force.__class__.__name__ : force for force in forces }
        print("PME parameters:")
        force = force_dict['NonbondedForce']
        print(force.getPMEParameters())
        (alpha, nx, ny, nz) = force.getPMEParameters()
        if alpha == 0.0 / unit.nanometers:
            # Set PME parameters explicitly.
            print("Setting PME parameters explicitly...")
            (alpha, nx, ny, nz) = calc_pme_parameters(system)
            print (alpha, nx, ny, nz)
            print(type(nx))
            force.setPMEParameters(alpha, int(nx), int(ny), int(nz))
            print(force.getPMEParameters())

        if args.tune_pme_platform:
            # Tune PME parameters for specified platform.
            from optimizepme import optimizePME
            properties = dict()
            platform = openmm.Platform.getPlatformByName(args.tune_pme_platform)
            print("Tuning PME parameters for platform '%s' precision model '%s'..." % (platform.getName(), args.precision))
            if (platform.getName() == 'OpenCL'):
                properties['OpenCLPrecision'] = args.precision
            elif (platform.getName() == 'CUDA'):
                properties['CudaPrecision'] = args.precision
            minCutoff = 0.8 * unit.nanometers
            maxCutoff = 1.2 * unit.nanometers
            optimizePME(system, integrator, positions, platform, properties, minCutoff, maxCutoff)

        class_name = 'RUN%d' % run
        logger.info("%s (%d atoms)" % (class_name, system.getNumParticles()))

        # Compute reference potential and force
        [reference_potential, reference_force] = compute_potential_and_force(system, positions, reference_platform)

        # Test all platforms.
        test_success = True
        for platform_index in range(openmm.Platform.getNumPlatforms()):
            try:
                platform = openmm.Platform.getPlatform(platform_index)
                platform_name = platform.getName()

                # Define precision models to test.
                if platform_name == 'Reference':
                    precision_models = ['double']
                else:
                    precision_models = ['single']
                    if platform.supportsDoublePrecision():
                        precision_models.append('double')

                for precision_model in precision_models:
                    # Set precision.
                    if platform_name == 'CUDA':
                        platform.setPropertyDefaultValue('CudaPrecision', precision_model)
                    if platform_name == 'OpenCL':
                        platform.setPropertyDefaultValue('OpenCLPrecision', precision_model)

                    # Compute potential and force.
                    [platform_potential, platform_force] = compute_potential_and_force(system, positions, platform)

                    # Compute error in potential.
                    potential_error = platform_potential - reference_potential

                    # Compute per-atom RMS (magnitude) and RMS error in force.
                    force_unit = unit.kilocalories_per_mole / unit.nanometers
                    natoms = system.getNumParticles()
                    force_mse = (((reference_force - platform_force) / force_unit)**2).sum() / natoms * force_unit**2
                    force_rmse = unit.sqrt(force_mse)

                    force_ms = ((platform_force / force_unit)**2).sum() / natoms * force_unit**2
                    force_rms = unit.sqrt(force_ms)

                    logger.info("%16s%16s %16.6f kcal/mol %16.6f kcal/mol %16.6f kcal/mol/nm %16.6f kcal/mol/nm" % (platform_name, precision_model, platform_potential / unit.kilocalories_per_mole, potential_error / unit.kilocalories_per_mole, force_rms / force_unit, force_rmse / force_unit))

                    # Mark whether tolerance is exceeded or not.
                    if abs(potential_error) > ENERGY_TOLERANCE:
                        test_success = False
                        logger.info("%32s WARNING: Potential energy error (%.6f kcal/mol) exceeds tolerance (%.6f kcal/mol).  Test failed." % ("", potential_error/unit.kilocalories_per_mole, ENERGY_TOLERANCE/unit.kilocalories_per_mole))
                    if abs(force_rmse) > FORCE_RMSE_TOLERANCE:
                        test_success = False
                        logger.info("%32s WARNING: Force RMS error (%.6f kcal/mol/nm) exceeds tolerance (%.6f kcal/mol/nm).  Test failed." % ("", force_rmse/force_unit, FORCE_RMSE_TOLERANCE/force_unit))
                        if verbose:
                            for atom_index in range(natoms):
                                for k in range(3):
                                    logger.info("%12.6f" % (reference_force[atom_index,k]/force_unit), end="")
                                logger.info(" : ", end="")
                                for k in range(3):
                                    logger.info("%12.6f" % (platform_force[atom_index,k]/force_unit), end="")
            except Exception as e:
                logger.info(e)

        if test_success:
            tests_passed += 1
        else:
            tests_failed += 1

        if (test_success is False):
            # Write XML files of failed tests to aid in debugging.
            
            # Place forces into different force groups.
            forces = [ system.getForce(force_index) for force_index in range(system.getNumForces()) ]
            force_group_names = dict()
            group_index = 0
            for force_index in range(system.getNumForces()):
                force_name = forces[force_index].__class__.__name__
                if force_name == 'NonbondedForce':
                    forces[force_index].setForceGroup(group_index+1)
                    force_group_names[group_index] = 'NonbondedForce (direct)'
                    group_index += 1
                    forces[force_index].setReciprocalSpaceForceGroup(group_index+1)
                    force_group_names[group_index] = 'NonbondedForce (reciprocal)'
                    group_index += 1
                else:
                    forces[force_index].setForceGroup(group_index+1)
                    force_group_names[group_index] = force_name
                    group_index += 1
            ngroups = len(force_group_names)

            # Test by force group.
            logger.info("Breakdown of discrepancies by Force component:")
            nforces = system.getNumForces()
            for force_group in range(ngroups):
                force_name = force_group_names[force_group]
                logger.info(force_name)
                [reference_potential, reference_force] = compute_potential_and_force_by_force_group(system, positions, reference_platform, force_group)
                logger.info("%16s%16s %16s          %16s          %16s          %16s" % ("platform", "precision", "potential", "error", "force mag", "rms error"))

                for platform_index in range(openmm.Platform.getNumPlatforms()):
                    try:
                        platform = openmm.Platform.getPlatform(platform_index)
                        platform_name = platform.getName()
                        
                        # Define precision models to test.
                        if platform_name == 'Reference':
                            precision_models = ['double']
                        else:
                            precision_models = ['single']
                            if platform.supportsDoublePrecision():
                                precision_models.append('double')

                        for precision_model in precision_models:
                            # Set precision.
                            if platform_name == 'CUDA':
                                platform.setPropertyDefaultValue('CudaPrecision', precision_model)
                            if platform_name == 'OpenCL':
                                platform.setPropertyDefaultValue('OpenCLPrecision', precision_model)
                                
                            # Compute potential and force.
                            [platform_potential, platform_force] = compute_potential_and_force_by_force_group(system, positions, platform, force_group)

                            # Compute error in potential.
                            potential_error = platform_potential - reference_potential

                            # Compute per-atom RMS (magnitude) and RMS error in force.
                            force_unit = unit.kilocalories_per_mole / unit.nanometers
                            natoms = system.getNumParticles()
                            force_mse = (((reference_force - platform_force) / force_unit)**2).sum() / natoms * force_unit**2
                            force_rmse = unit.sqrt(force_mse)

                            force_ms = ((platform_force / force_unit)**2).sum() / natoms * force_unit**2
                            force_rms = unit.sqrt(force_ms)

                            logger.info("%16s%16s %16.6f kcal/mol %16.6f kcal/mol %16.6f kcal/mol/nm %16.6f kcal/mol/nm" % (platform_name, precision_model, platform_potential / unit.kilocalories_per_mole, potential_error / unit.kilocalories_per_mole, force_rms / force_unit, force_rmse / force_unit))

                    except Exception as e:
                        logger.info(e)
                        pass
        logger.info("")

    logger.info("%d tests failed" % tests_failed)
    logger.info("%d tests passed" % tests_passed)

    if (tests_failed > 0):
        # Signal failure of test.
        sys.exit(1)
    else:
        sys.exit(0)

if __name__ == "__main__":
    main()
