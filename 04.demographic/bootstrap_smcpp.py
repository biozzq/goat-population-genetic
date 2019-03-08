#!/stor9000/apps/users/NWSUAF/2012010954/Software/Anaconda4.4_py3.6/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 15:22:36 2017

"""

import click
import argparse
import os
import sys
import random
import gzip
import time


@click.command()
@click.option('--nr_bootstraps',  type=int, help="nr of bootstraps [20]", default=20)
@click.option("--chunk_size", type=int, help="size of bootstrap chunks [5000000]", default=5000000)
@click.option("--chunks_per_chromosome", type=int,
              help="nr of chunks to put on one chromosome in the bootstrap [20]", default=20)
@click.option("--nr_chromosomes", type=int, help="nr of chromosomes to write [30]", default=30)
@click.option("--seed", type=int, help="initialize the random number generator", default=None)
@click.argument("out_dir_prefix")
@click.argument("files", nargs=-1)
def main(nr_bootstraps, chunk_size, chunks_per_chromosome, nr_chromosomes, seed, out_dir_prefix, files):
    chunks = []
    offset = 0
    chunks_in_chrom = []
    if not seed:
        seed = int(time.time())
    random.seed(seed)
    print('seed: %s' % seed)
    for file in files:
        with gzip.open(file, 'rb') as f:
            header = f.readline().decode()
            for line in f:
                line = line.decode()
                line = [int(x) for x in line.strip().split()]
                pos = line[0] + offset
                offset += line[0]
                chunk_index = (pos - 1) // chunk_size
                if chunk_index > len(chunks_in_chrom)-1:
                    chunks_in_chrom.append([])
                chunks_in_chrom[chunk_index].append(line)
        chunks.extend(chunks_in_chrom)

    for bootstrap_id in range(1, nr_bootstraps +1):
        for chr_ in range(1, nr_chromosomes + 1):
            chr_dir = "{}_{}".format(out_dir_prefix, bootstrap_id)
            if not os.path.exists(chr_dir):
                os.makedirs(chr_dir)
            chr_file = "{}/bootstrap_chr{}.gz".format(chr_dir, chr_)
            print("writing", chr_file, file=sys.stderr)
            with gzip.open(chr_file, 'wb') as f:
                f.write(header.encode())
                for i in range(chunks_per_chromosome):
                    chunk_id = random.randrange(len(chunks))
                    for line in chunks[chunk_id]:
                        line = ' '.join([str(x) for x in line]) + '\n'
                        f.write(line.encode())



if __name__ == "__main__":
    main()
