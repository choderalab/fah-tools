#!/usr/bin/env python

# Perform quality control checking of OpenMM FAH data.

import simtk.unit as units

def extract_trajectory(results_filename, reference_pdb_filename, pdb_trajectory_filename, atomSubset=None):
    """
    Extract PDB trajectory from specified compressed (bz2) payload.

    ARGUMENTS

    results_filename (string) - name of compressed results file
    reference_pdb_filename (string) - name of PDB reference file
    pdb_trajectory_filename (string) - name of PDB filename to write trajectory to

    TODO

    * Add support for writing only a subset of atoms
    * Add support for concatenating multiple gens into a single trajectory

    """
    import numpy    
    import mdtraj # MDTraj: https://github.com/rmcgibbo/mdtraj
    
    # Read reference PDB file.
    import simtk.openmm 
    import simtk.unit as units
    import simtk.openmm.app as app
    pdb = app.PDBFile(reference_pdb_filename)

    # Copy results to temporary directory.
    import os, os.path, tempfile, shutil
    cwd = os.getcwd()
    tmpdir = tempfile.mkdtemp()
    print tmpdir
    shutil.copyfile(results_filename, os.path.join(tmpdir, 'results.tar.bz2'))
    os.chdir(tmpdir)

    # Extract XTC file from payload.
    import commands
    command = 'bzcat results.tar.bz2 | tar x positions.xtc' 
    commands.getoutput(command)
    
    # Read XTC file.
    from mdtraj.xtc import XTCReader
    xtc = XTCReader('positions.xtc')

    # Clean up temporary directory.
    os.chdir(cwd)
    for filename in os.listdir(tmpdir):
        os.unlink(os.path.join(tmpdir, filename))
    os.removedirs(tmpdir)

    # Write multi-model PDB file of trajectory.
    outfile = open(pdb_trajectory_filename, 'w')

    #app.PDBFile.writeModel(pdb.topology, pdb.positions, outfile, 0, atomSubset=atomSubset)
    
    nframes = 0
    for (frame_index, frame) in enumerate(xtc):
        [xyz, time, step, box, prec] = frame
        print (time[0], step[0])
        positions = units.Quantity(xyz[0,:,:], units.nanometers)
        #app.PDBFile.writeModel(pdb.topology, positions, outfile, frame_index+1, atomSubset=atomSubset)
        nframes += 1
    print "trajectory has %5d frames" % nframes
    outfile.close()
    return

def read_file(filename):
    infile = open(filename, 'r')
    contents = infile.read()
    infile.close()
    return contents


def check_xtc(filename):
    from mdtraj.xtc import XTCReader
    xtc = XTCReader(filename)

    # Check XTC.
    import numpy
    nframes = 0
    for (frame_index, frame) in enumerate(xtc):
        [xyz, time, step, box, prec] = frame    
        
        if numpy.any( numpy.isnan(xyz) ):
            raise Exception("positions are NaN at frame %d, time %f" % frame, time)
        
        nframes += 1
        
    print "XTC contains %d frames" % nframes        
    return

def check_state(system_xml_filename, integrator_xml_filename, state_xml_filename, energy_error_tolerance=0.6*units.kilocalories_per_mole):
    import simtk.openmm as mm

    # Read System.
    system_xml = read_file(system_xml_filename)
    system = mm.XmlSerializer.deserialize(system_xml)

    # Read Integrator.
    integrator_xml = read_file(integrator_xml_filename)
    integrator = mm.XmlSerializer.deserialize(integrator_xml)

    # Read state.
    state_xml = read_file(state_xml_filename)
    import string
    state_xml = string.replace(state_xml, '<checkpoint', '<State')
    serialized_state = mm.XmlSerializer.deserialize(state_xml)
        
    # Create Context.
    context = mm.Context(system, integrator)
    context.setPositions(serialized_state.getPositions())
    context.setVelocities(serialized_state.getVelocities())
    box_vectors = serialized_state.getPeriodicBoxVectors()
    context.setPeriodicBoxVectors(*box_vectors)

    # Check initial state.
    import numpy
    state = context.getState(getForces=True, getEnergy=True)
    serialized_potential = serialized_state.getPotentialEnergy()
    computed_potential = state.getPotentialEnergy()
    delta_potential = serialized_potential - computed_potential
    print "potential: serialized | computed : %12.3f %12.3f kcal/mol : delta = %.3f kcal/mol" % (serialized_potential / units.kilocalories_per_mole, computed_potential / units.kilocalories_per_mole, delta_potential / units.kilocalories_per_mole)
    if abs(delta_potential) > energy_error_tolerance:
        raise Exception("Energy discrepancy exceeds tolerance of %s" % str(energy_error_tolerance))
    
    # Check initial energy from serialized coordinates.
    if numpy.isnan(state.getPotentialEnergy() / units.kilocalories_per_mole):
        raise Exception("Initial energy is NaN.")
    
    # Integrate.
    nsteps = 10
    integrator.step(nsteps)

    # Check final state.
    import numpy
    state = context.getState(getPositions=True, getVelocities=True, getForces=True, getEnergy=True)
    if numpy.isnan(state.getPotentialEnergy() / units.kilocalories_per_mole):
        raise Exception("Final energy is NaN")
    
    # Clean up.
    del context, state, integrator

def test_resume_from_xtc(system_xml_filename, integrator_xml_filename, state_xml_filename, xtc_filename, energy_error_tolerance=0.6*units.kilocalories_per_mole):
    """
    Check if we can resume from XTC frames.

    """

    import simtk.openmm as mm

    # Read System.
    system_xml = read_file(system_xml_filename)
    system = mm.XmlSerializer.deserialize(system_xml)

    # Read Integrator.
    integrator_xml = read_file(integrator_xml_filename)
    integrator = mm.XmlSerializer.deserialize(integrator_xml)

    # Read state.
    state_xml = read_file(state_xml_filename)
    import string
    state_xml = string.replace(state_xml, '<checkpoint', '<State')
    serialized_state = mm.XmlSerializer.deserialize(state_xml)
        
    # Create Context.
    context = mm.Context(system, integrator)
    context.setPositions(serialized_state.getPositions())
    context.setVelocities(serialized_state.getVelocities())
    box_vectors = serialized_state.getPeriodicBoxVectors()
    context.setPeriodicBoxVectors(*box_vectors)

    # Check initial state.
    import numpy
    state = context.getState(getForces=True, getEnergy=True)
    serialized_potential = serialized_state.getPotentialEnergy()
    computed_potential = state.getPotentialEnergy()
    delta_potential = serialized_potential - computed_potential
    print "potential: serialized | computed : %12.3f %12.3f kcal/mol : delta = %.3f kcal/mol" % (serialized_potential / units.kilocalories_per_mole, computed_potential / units.kilocalories_per_mole, delta_potential / units.kilocalories_per_mole)
    if abs(delta_potential) > energy_error_tolerance:
        raise Exception("Energy discrepancy exceeds tolerance of %s" % str(energy_error_tolerance))

    # Check initial energy from serialized coordinates.
    if numpy.isnan(state.getPotentialEnergy() / units.kilocalories_per_mole):
        raise Exception("Initial energy is NaN.")    

    # Read XTC file.
    from mdtraj.xtc import XTCReader
    xtc = XTCReader(xtc_filename)
    import numpy
    nframes = 0
    for (frame_index, frame) in enumerate(xtc):
        print "Attempting to resume from frame %d" % frame_index

        [xyz, time, step, box, prec] = frame    
        
        positions = units.Quantity(xyz[0,:,:], units.nanometers)
        box_vectors = units.Quantity((mm.Vec3(box[0,0,0], 0, 0), mm.Vec3(0, box[0,1,1], 0), mm.Vec3(0, 0, box[0,2,2])), units.nanometers)

        # Set positions and box vectors.
        context.setPositions(serialized_state.getPositions())
        context.setPeriodicBoxVectors(*box_vectors)

        # Randomize velocities.
        temperature = integrator.getTemperature()
        context.setVelocitiesToTemperature(temperature)

        # Integrate.
        nsteps = 50
        integrator.step(nsteps)

        # Check final state.
        import numpy
        state = context.getState(getEnergy=True)
        print "final energy %12.3f kcal/mol" % (state.getPotentialEnergy() / units.kilocalories_per_mole)
#        if numpy.isnan(state.getPotentialEnergy() / units.kilocalories_per_mole):
#            raise Exception("Final energy is NaN")
        
        nframes += 1
        
    print "XTC contains %d frames" % nframes        

    # Clean up.
    del context, state, integrator

def check_results(results_filename):
    """
    Test FAH OpenMM results packet.

    ARGUMENTS

    results_filename (string) - name of compressed results file to test

    TODO

    * Add support for writing only a subset of atoms
    * Add support for concatenating multiple gens into a single trajectory

    """
    import numpy    
    import mdtraj # MDTraj: https://github.com/rmcgibbo/mdtraj
    
    # Read reference PDB file.
    import simtk.openmm 
    import simtk.unit as units
    from simtk.openmm import app

    # Create temporary directory.
    import os, os.path, tempfile, shutil
    cwd = os.getcwd()
    tmpdir = tempfile.mkdtemp()

    # Extract source directory.
    [directory, filename] = os.path.split(results_filename)

    # Copy results to temporary directory.
    shutil.copyfile(results_filename, os.path.join(tmpdir, 'results.tar.bz2'))

    # Copy payload to the temporary directory.
    gen = 0
    payload_filename = os.path.join(directory, 'payload-%03d.tar.bz2' % gen)
    while not os.path.exists(payload_filename):
        gen += 1
        payload_filename = os.path.join(directory, 'payload-%03d.tar.bz2' % gen)
    shutil.copyfile(payload_filename, os.path.join(tmpdir, 'payload.tar.bz2'))

    # Change to temporary directory.
    os.chdir(tmpdir)

    # Extract payload and results.
    import commands
    command = 'bzcat payload.tar.bz2 | tar x' 
    commands.getoutput(command)
    command = 'bzcat results.tar.bz2 | tar x' 
    commands.getoutput(command)
    
    # Check XTC file.
    check_xtc('positions.xtc')

    # Check state
    check_state('system.xml', 'integrator.xml', 'checkpointState.xml')

    # Test if we can resume from XTC frames
    #test_resume_from_xtc('system.xml', 'integrator.xml', 'checkpointState.xml', 'positions.xtc')

    # Clean up temporary directory.
    os.chdir(cwd)
    for filename in os.listdir(tmpdir):
        os.unlink(os.path.join(tmpdir, filename))
    os.removedirs(tmpdir)

    return

if __name__ == '__main__':
    # Name of project directory to check.
    project_directory = 'PROJ8900'
    write_trajectories = False
    natoms = 4091 # number of atoms to write
    
    # Trap SIGXCPU (for CUDA 5.0 bug / SGE reasons).
    ignore_signals = ['SIGXCPU']
    if len(ignore_signals) > 0:
        import signal
        # Create a dummy signal handler.
        def signal_handler(signal, frame):
            print 'Signal %s received and ignored.' % str(signal)
        # Register the dummy signal handler.
        for signal_name in ignore_signals:
            print "Will ignore signal %s" % signal_name
            signal.signal(getattr(signal, signal_name), signal_handler)

    # Build list of directories.
    import os, os.path
    runs = os.listdir(project_directory)
    for run in runs:
        run_directory = os.path.join(project_directory, run)
        clones = os.listdir(run_directory)
        for clone in clones:
            clone_directory = os.path.join(run_directory, clone)
            
    # Get all results packets.
    import glob
    results_list = glob.glob('%s/RUN*/CLONE*/results-*.tar.bz2' % project_directory)
    print results_list
    for results_filename in results_list:
        # Extract run, clone, and gen number.
        import re
        match = re.search('RUN(\d+)/CLONE(\d+)/results-(\d+).tar.bz2', results_filename)
        [run,clone,gen] = (int(match.groups(1)[0]), int(match.groups(1)[1]), int(match.groups(1)[2]))
        print [run,clone,gen]

        print "Checking %s..." % results_filename
        check_results(results_filename)

        if write_trajectories:
            # Write trajectory of protein only.
            import os.path
            reference_pdb_filename = '%s/RUN0/system.pdb' % project_directory
            [results_directory, filename] = os.path.split(results_filename)
            pdb_trajectory_filename = 'trajectory-RUN%03d-CLONE%03d-GEN%03d.pdb' % (run, clone, gen)
            atomSubset = range(natoms)
            extract_trajectory(results_filename, reference_pdb_filename, pdb_trajectory_filename, atomSubset)
        
