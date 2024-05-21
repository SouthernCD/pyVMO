import argparse
from pyvmo import VMO
import pandas as pd
from pyvmo.distance import get_IBS_matrix
import shutil

class Job(object):
    def __init__(self):
        pass

    def run_arg_parser(self):
        # argument parse
        parser = argparse.ArgumentParser(
            prog='pyVMO',
        )

        subparsers = parser.add_subparsers(
            title='subcommands', dest="subcommand_name")

        # argparse for converter
        parser_a = subparsers.add_parser('converter',
                                         description='Convert the vcf file to vmo\n')

        parser_a.add_argument('input_path', type=str,
                              help='input file or directory')
        parser_a.add_argument('output_path', type=str,
                              help='output file or directory')
        parser_a.add_argument('-vmo2bimbam', action='store_true',
                                help='convert vmo to bimbam', default=False)
        parser_a.add_argument('-vcf2vmo', action='store_true',
                                help='convert vcf to vmo', default=True)
        parser_a.add_argument('-vmo2vcf', action='store_true',
                                help='convert vmo to vcf', default=False)
        parser_a.set_defaults(func=converter_main)

        # argparse for extractor
        parser_a = subparsers.add_parser('extractor',
                                         description='Extraction of useful submatrices by sample listing, variant site quality control\n')

        parser_a.add_argument('sample_id_list_file', type=str,
                              help='a list file which store samples id')
        parser_a.add_argument('input_vmo_path', type=str,
                              help='input vmo directory')
        parser_a.add_argument('output_vmo_path', type=str,
                              help='output vmo directory')
        parser_a.add_argument('-m', '--maf_thr', type=float,
                              help='min minor allele frequency', default=0.05)
        parser_a.add_argument('-s', '--mis_thr', type=float,
                              help='max proportion of missing data', default=0.5)
        parser_a.add_argument('-c', '--chunk_size', type=int,
                              help='chunk size', default=2000)
        parser_a.add_argument('-n', '--n_jobs', type=int,
                              help='number of parallel jobs', default=20)
        parser_a.set_defaults(func=extractor_main)

        # argparse for distance
        parser_a = subparsers.add_parser('distance',
                                         description='Calculate the genetic distance between samples\n')

        parser_a.add_argument('input_vmo_path', type=str,
                              help='input vmo directory')
        parser_a.add_argument('output_path', type=str,
                                help='output file')
        parser_a.add_argument('-t', '--distance_type', type=str,
                                help='distance type', default='IBS')
        parser_a.add_argument('-c', '--chunk_size', type=int,
                              help='chunk size', default=200)
        parser_a.add_argument('-n', '--n_jobs', type=int,
                              help='number of parallel jobs', default=20)
        parser_a.set_defaults(func=distance_main)

        self.arg_parser = parser

        self.args = parser.parse_args()
        
        # parser.set_defaults(func=parser.print_help())

    def run(self):
        self.run_arg_parser()
        # self.args = self.arg_parser.parse_args()
        # if hasattr(self, "subcommand_name"):
        #     self.args.func(self.args)
        # else:
        #     self.arg_parser.print_help()

        args_dict = vars(self.args)

        if args_dict["subcommand_name"] == "converter":
            converter_main(self.args)
        elif args_dict["subcommand_name"] == "extractor":
            extractor_main(self.args)
        elif args_dict["subcommand_name"] == "distance":
            distance_main(self.args)
        else:
            self.arg_parser.print_help()



def converter_main(args):
    if args.vmo2bimbam:
        vmo = VMO(args.input_path)
        vmo.to_bimbam(args.output_path)
    elif args.vcf2vmo:
        vmo = VMO(args.output_path)
        vmo.store_vcf_to_vmo(args.input_path)
    elif args.vmo2vcf:
        vmo = VMO(args.input_path)
        vmo.to_vcf(args.output_path)

def extractor_main(args):
    vmo = VMO(args.input_vmo_path)

    # load sample list
    sample_id_list = [line.strip() for line in open(args.sample_id_list_file)]
    
    # extract by sample list
    spl_idx_list = vmo.get_samples_index(sample_id_list)
    spl_vmo_path = args.output_vmo_path + "_spl"
    spl_vmo = vmo.extract_sub_vmo(spl_vmo_path, spl_idx_list=spl_idx_list)

    # extract by variant site quality control
    var_idx_list = spl_vmo.site_filter(maf_thr=args.maf_thr, mis_thr=args.mis_thr, chunk_size=args.chunk_size, n_jobs=args.n_jobs)
    var_vmo_path = args.output_vmo_path
    var_vmo = spl_vmo.extract_sub_vmo(var_vmo_path, var_idx_list=var_idx_list)

    # remove spl_vmo
    shutil.rmtree(spl_vmo_path)

def distance_main(args):
    vmo = VMO(args.input_vmo_path)
    m = vmo.get_matrix()
    if args.distance_type == "IBS":
        results = get_IBS_matrix(m, chunk_size=args.chunk_size, n_jobs=args.n_jobs)
    else:
        raise ValueError("Unknown distance type %s" % args.distance_type)
    
    with open(args.output_path, 'w') as f:
        for row in results:
            f.write("%s\n" % ', '.join(["%.5f" % x for x in row]))

def main():
    job = Job()
    job.run()


if __name__ == '__main__':
    main()