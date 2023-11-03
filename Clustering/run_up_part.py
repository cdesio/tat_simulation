import argparse

# currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
# sys.path.insert(0, parentdir)

import runClustering_up_part
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--filename", type=str,
                        help="Root file to run clustering on")
    parser.add_argument("--filenamePhoton", type=str,
                        help="Root file to run clustering on")
    parser.add_argument("--output", type=str, help="outfile name")
    parser.add_argument("--min", type=float,
                        help="Emin direct damage, default 5eV")
    parser.add_argument("--max", type=float,
                        help="Emax direct damage, default 37.5eV")
    parser.add_argument("--indirect", type=float,
                        help="Indirect damage percentage, default 0.405")
    parser.add_argument("--sugar", type=str, help="sugar geometry file")
    # parser.add_argument("--primary", type=str, help="select particle to save to output file", required=False)
    parser.add_argument("--sepR", type=bool, help="save to separate files for R", default=True)
    parser.add_argument("--ndiv_R", type=int, help="no. of divisions in radius", default=80)
    parser.add_argument("--continuous", type=bool, help="flag for continuous geometry", default = False)
    # parser.add_argument("--ndiv_Z", type=int, help="no. of divisions in Z", default=40)
    # parser.add_argument("--startR", type=float,
    #                     help="starting radius for boxes (um)", default=10.5)
    # parser.add_argument("--spacing", type=float,
    #                     help="boxes spacing in R (um)", default=1)
    
    args = parser.parse_args()

    if args.filename:
        filename = args.filename
    else:
        print("Missing filename")

    if args.filenamePhoton:
        filenamePhoton = args.filenamePhoton
    else:
        filenamePhoton = False

    if args.output:
        outputFilename = args.output
    else:
        print("Missing output filename")

    if args.min:
        fEMinDamage = args.min
    else:
        fEMinDamage = 5

    if args.max:
        fEMaxDamage = args.max
    else:
        fEMaxDamage = 37.5

    if args.indirect:
        probIndirect = args.indirect
    else:
        probIndirect = 0.405

    if args.sugar:
        sugarPosFilename = args.sugar
    else:
        print("Missing sugar geometry file")

        
    

    runClustering_up_part.runClustering(filename, outputFilename, fEMinDamage, fEMaxDamage,
                                probIndirect, sugarPosFilename, filenamePhoton=filenamePhoton, separate_r=args.sepR, 
                                ndiv_R=args.ndiv_R, continuous=args.continuous) #, ndiv_Z=ndiv_Z, start_R=10.5, spacing=1, primaryParticle=args.primary)
