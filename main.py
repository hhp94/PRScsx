import os

import parse_genet
import mcmc_gtb
import parse_param

def main():
    param_dict = parse_param.parse_param()
    n_pop = len(param_dict['pop'])

    if n_pop == 1:
        print('*** only %d discovery population detected ***\n' % n_pop)
    else:
        print('*** %d discovery populations detected ***\n' % n_pop)


    for chrom in param_dict['chrom']:
        print('##### process chromosome %d #####' % int(chrom))

        if os.path.isfile(param_dict['ref_dir'] + '/snpinfo_mult_1kg_hm3'):
            ref = '1kg'
            ref_dict = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_mult_1kg_hm3', int(chrom), ref)
        elif os.path.isfile(param_dict['ref_dir'] + '/snpinfo_mult_ukbb_hm3'):
            ref = 'ukbb'
            ref_dict = parse_genet.parse_ref(param_dict['ref_dir'] + '/snpinfo_mult_ukbb_hm3', int(chrom), ref)

        vld_dict = parse_genet.parse_bim(param_dict['bim_prefix'], int(chrom))

        sst_dict = {}
        for pp in range(n_pop):
            sst_dict[pp] = parse_genet.parse_sumstats(ref_dict, vld_dict, param_dict['sst_file'][pp], param_dict['pop'][pp], param_dict['n_gwas'][pp])

        ld_blk = {}
        blk_size = {}
        for pp in range(n_pop):
            ld_blk[pp], blk_size[pp] = parse_genet.parse_ldblk(param_dict['ref_dir'], sst_dict[pp], param_dict['pop'][pp], int(chrom), ref)

        snp_dict, beta_dict, frq_dict, idx_dict = parse_genet.align_ldblk(ref_dict, vld_dict, sst_dict, n_pop, int(chrom))

        mcmc_gtb.mcmc(param_dict['a'], param_dict['b'], param_dict['phi'], snp_dict, beta_dict, frq_dict, idx_dict, param_dict['n_gwas'], ld_blk, blk_size,
            param_dict['n_iter'], param_dict['n_burnin'], param_dict['thin'], param_dict['pop'], int(chrom),
            param_dict['out_dir'], param_dict['out_name'], param_dict['meta'], param_dict['write_pst'], param_dict['seed'])

        print('\n')


if __name__ == '__main__':
    main()
