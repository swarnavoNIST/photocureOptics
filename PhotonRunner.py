import sys, os
from optparse import OptionParser
from MC_light_absorption import MC_light_absorption

def processOptions():
    parser = OptionParser()
    parser.add_option("-i", dest="inputfile", help="Name of the file containing the inputs for a photon transport simulation. Use option -f to find the format.", default="")
    parser.add_option("-f", dest="inputFormat", help="Use 'show' to see the input format. Or specify the name to create an empty inputfile with the placeholders.", default="")
    
    [options, args] = parser.parse_args()

    if options.inputFormat == "show":
        print >> sys.stdout , '\n\nThe inputfile must contain the following quantities. Use only SI units.\n'
        print >> sys.stdout , '\n# Material Properties \n'
        print >> sys.stdout , 'Absorption coefficient - '
        print >> sys.stdout , 'Chemical absorption - '
        print >> sys.stdout , 'Scattering coefficient - '
        print >> sys.stdout , 'Scattering anisotropy - '

        print >> sys.stdout , '\n# Simulation details \n'
        print >> sys.stdout , 'Initial photon weight - '
        print >> sys.stdout , 'No. of photons - '
        print >> sys.stdout , 'Display intervals - '
        print >> sys.stdout , 'Weight threshold - '
        print >> sys.stdout , 'Roulette Constant - '

        print >> sys.stdout , '\n# Domain and boundary specifications\n'
        print >> sys.stdout , 'Sample height - '
        print >> sys.stdout , 'Sample radius - '
        print >> sys.stdout , 'Refractive index (at Z=0) - '
        print >> sys.stdout , 'Refractive index (at Z=height) - '
        print >> sys.stdout , 'Refractive index (at X^2 + Y^2 = R^2) - '

        print >> sys.stdout , '\n# Multiprocessing option'
        print >> sys.stdout , 'Total threads - '
        
        print >> sys.stdout , '\n# Output details'
        print >> sys.stdout , 'Output filename - '

        print >> sys.stdout , '\n# Steady-State or Transient Simulation'
        print >> sys.stdout , 'Simulation type - '
        print >> sys.stdout , 'Tolerance - '
        print >> sys.stdout , 'Total number of steps - '
        
    elif len(options.inputFormat) != 0:
        emptyInputFile = open(options.inputFormat,'w')
        originalStdOut = sys.stdout
        sys.stdout = emptyInputFile
        print >> sys.stdout , '# Input for the photon transport code. Use only SI units.'
        print >> sys.stdout , '\n# Material Properties'
        print >> sys.stdout , 'Adsorption coefficient - '
        print >> sys.stdout , 'Scattering coefficient - '
        print >> sys.stdout , 'Scattering anisotropy - '
        print >> sys.stdout , 'Refractive index - '

        print >> sys.stdout , '\n# Simulation details'
        print >> sys.stdout , 'Initial photon weight - '
        print >> sys.stdout , 'No. of photons - '
        print >> sys.stdout , 'Display intervals - '
        print >> sys.stdout , 'Weight threshold - '
        print >> sys.stdout , 'Roulette constant - '
        print >> sys.stdout , 'Grid resolution - '

        print >> sys.stdout , '\n# Domain and boundary specifications'
        print >> sys.stdout , 'Sample height - '
        print >> sys.stdout , 'Sample radius - '
        print >> sys.stdout , 'Refractive index (at Z=0) - '
        print >> sys.stdout , 'Refractive index (at Z=height) - '
        print >> sys.stdout , 'Refractive index (at X^2 + Y^2 = R^2) - '

        print >> sys.stdout , '\n# Multiprocessing option'
        print >> sys.stdout , 'Total threads - '
        
        print >> sys.stdout , '\n# Output details'
        print >> sys.stdout , 'Output filename - '

        print >> sys.stdout , '\n# Steady-State or Transient Simulation'
        print >> sys.stdout , 'Simulation type - '
        print >> sys.stdout , 'Tolerance - '
        print >> sys.stdout , 'Total number of steps - '

        emptyInputFile.close()
        sys.stdout = originalStdOut
    elif len(options.inputfile) != 0:
        try:
            InputFile = open(options.inputfile, 'r')
        except IOError:
            print >> sys.stderr , "ERROR : Cannot open inputfile. Check inputfile name."
            
        InputLines = InputFile.readlines()
        
        removeSet = []
        for l in InputLines:
            if l[0] == '#' or l[0] == '\n':
                removeSet.append(l)
        
        for rem in removeSet:
            InputLines.remove(rem)

        return InputLines

    else:
        print >> sys.stderr , "ERROR : Please try 'python PhotoRunner.py -h' to learn about correct usage"

if __name__ == '__main__':
    InputLines = processOptions()
        
    if InputLines:
        photonSimulator = MC_light_absorption(InputLines)
        photonSimulator.initializeStructs()
        photonSimulator.timeLoop()
        photonSimulator.write_outputs()
        #print "Total Absorption : ", total_absorption
        #photonSimulator.writeVTKfile()
        #photonSimulator.writePVD()