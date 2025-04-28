import sys
import getopt

def parse_param():
    long_opts_list = ['ref_dir=', 'bim_prefix=', 'sst_file=', 'a=', 'b=', 'phi=', 'n_gwas=', 'pop=',
                      'n_iter=', 'n_burnin=', 'thin=', 'out_dir=', 'out_name=', 'chrom=', 'meta=', 'write_pst=', 'seed=', 'help']

    param_dict = {'ref_dir': None, 'bim_prefix': None, 'sst_file': None, 'a': 1, 'b': 0.5, 'phi': None, 'n_gwas': None, 'pop': None,
                  'n_iter': None, 'n_burnin': None, 'thin': 5, 'out_dir': None, 'out_name': None, 'chrom': range(1,23),
                  'meta': 'FALSE', 'write_pst': 'FALSE', 'seed': None}

    print('\n')

    if len(sys.argv) > 1:
        try:
            opts, args = getopt.getopt(sys.argv[1:], "h", long_opts_list)          
        except getopt.GetoptError:
            print('* Option not recognized.')
            print('* Use --help for usage information.\n')
            sys.exit(2)

        for opt, arg in opts:
            if opt == "-h" or opt == "--help":
                print(__doc__)
                sys.exit(0)
            elif opt == "--ref_dir": 
                param_dict['ref_dir'] = arg
            elif opt == "--bim_prefix": 
                param_dict['bim_prefix'] = arg
            elif opt == "--sst_file": 
                param_dict['sst_file'] = arg.split(',')
            elif opt == "--a": 
                param_dict['a'] = float(arg)
            elif opt == "--b": 
                param_dict['b'] = float(arg)
            elif opt == "--phi": 
                param_dict['phi'] = float(arg)
            elif opt == "--n_gwas": 
                param_dict['n_gwas'] = list(map(int,arg.split(',')))
            elif opt == "--pop": 
                param_dict['pop'] = arg.split(',')
            elif opt == "--n_iter": 
                param_dict['n_iter'] = int(arg)
            elif opt == "--n_burnin": 
                param_dict['n_burnin'] = int(arg)
            elif opt == "--thin": 
                param_dict['thin'] = int(arg)
            elif opt == "--out_dir": 
                param_dict['out_dir'] = arg
            elif opt == "--out_name": 
                param_dict['out_name'] = arg
            elif opt == "--chrom": 
                param_dict['chrom'] = arg.split(',')
            elif opt == "--meta": 
                param_dict['meta'] = arg.upper()
            elif opt == "--write_pst": 
                param_dict['write_pst'] = arg.upper()
            elif opt == "--seed": 
                param_dict['seed'] = int(arg)
    else:
        print(__doc__)
        sys.exit(0)


    if param_dict['ref_dir'] is None:
        print('* Please specify the directory to the reference panel using --ref_dir\n')
        sys.exit(2)
    elif param_dict['bim_prefix'] is None:
        print('* Please specify the directory and prefix of the bim file for the target (validation/testing) dataset using --bim_prefix\n')
        sys.exit(2)
    elif param_dict['sst_file'] is None:
        print('* Please provide at least one summary statistics file using --sst_file\n')
        sys.exit(2)
    elif param_dict['n_gwas'] is None:
        print('* Please provide the sample size of the GWAS using --n_gwas\n')
        sys.exit(2)
    elif param_dict['pop'] is None:
        print('* Please specify the population of the GWAS sample using --pop\n')
        sys.exit(2)
    elif param_dict['out_dir'] is None:
        print('* Please specify the output directory using --out_dir\n')
        sys.exit(2)
    elif param_dict['out_name'] is None:
        print('* Please specify the prefix of the output file using --out_name\n')
        sys.exit(2)
    elif (len(param_dict['sst_file']) != len(param_dict['n_gwas']) or 
          len(param_dict['sst_file']) != len(param_dict['pop'])):
        print('* Length of sst_file, n_gwas and pop does not match\n')
        sys.exit(2)

    n_pop = len(param_dict['pop'])
    if param_dict['n_iter'] is None or param_dict['n_burnin'] is None:
        param_dict['n_iter'] = n_pop*1000
        param_dict['n_burnin'] = n_pop*500

    for key in param_dict:
        print('--%s=%s' % (key, param_dict[key]))

    print('\n')
    return param_dict
