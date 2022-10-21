# CLI for MetaQuad

import os
import sys
import time
from optparse import OptionParser, OptionGroup
from metaquad.version import __version__
from vireoSNP.utils.io_utils import read_cellSNP
from metaquad.metaquad import Metaquad

START_TIME = time.time()


def main():
    # With warnings
    parser = OptionParser()
    parser.add_option("--Data", "-i", dest="import_data", default=None,
                      help=("The cellSNP folder with AD and DP matirx."))
    parser.add_option("--outDir", "-o", dest="out_dir", default=None,
                      help=("Directory for output files [default: $Data/metaquad_out]"))

    group1 = OptionGroup(parser, "Different parameters")

    group1.add_option("--cutoff", "-t", type="int", dest="cutoff", default=10,
                      help=("User-defined deltaBIC cutoff [default:10] "))
    group1.add_option("--export_heatmap", "-d", dest="export_heatmap", default=False,
                      help=("Drew heatmap based on informative mutations [default:False] "))
    group1.add_option("--randSeed", type="int", dest="rand_seed", default=None,
                      help="Seed for random initialization [default: %default]")
    group1.add_option("--nproc", "-p", type="int", dest="nproc", default=10,
                      help="Number of subprocesses [default: 10]")
    group1.add_option("--minDP", type="int", dest="minDP", default=10,
                      help="Minimum DP to include for modelling [default: 10]")
    group1.add_option("--minSample", type='int', dest="minSample", default=5,
                      help=("Minimum number of samples in minor component [default: 5]"))
    group1.add_option("--batchSize", type='int', dest="batch_size", default=128,
                      help=("Number of variants in one batch, cooperate with --nproc for speeding up [default: 128]"))


    parser.add_option_group(group1)
    (options, args) = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        print("Welcome to MetaQuad v%s!\n" % (__version__))
        print("use metaquad -h or metaquad --help for help on parameters.")
        sys.exit(1)

    ## output directory
    if options.out_dir is None:
        print("Warning: no output directory provided, $Data/metaquad_out will be used")
        out_dir = os.path.dirname(os.path.abspath(options.cell_file)) + "/metaquad_out"
    elif os.path.dirname(options.out_dir) == "":
        out_dir = "./" + options.out_dir
    else:
        out_dir = options.out_dir
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    ## input data (cellsnp folder)
    if (options.import_data is None):
        print("Error: need data in cellSNP output folder")
        sys.exit(1)

    else:
        print("[MetaQuad] Loading data ...")
        input_data = read_cellSNP(options.import_data)

    ## other parameters
    nproc = options.nproc
    minDP = options.minDP
    batch_size = options.batch_size
    cutoff = options.cutoff
    minSample = options.minSample
    export_heatmap = options.export_heatmap

    ## Main functions

    mdphd = Metaquad(
                AD=input_data['AD'],
                DP=input_data['DP'],
                variant_names=input_data['variants']
            )
    df = mdphd.fit_deltaBIC(
                out_dir=out_dir,
                minDP=minDP,
                nproc=nproc,
                batch_size=batch_size,
            )
    best_ad, best_dp = mdphd.selectInformativeVariants(min_samples=minSample, out_dir=out_dir,
                                                       cutoff=cutoff, export_heatmap=export_heatmap)

    run_time = time.time() - START_TIME
    print("[MetaQuad] Time usage: %d min %.1f sec" % (int(run_time / 60),
                                                 run_time % 60))
    print()


if __name__ == "__main__":
    main()
