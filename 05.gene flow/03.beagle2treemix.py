# -*- coding: utf-8 -*-
"""
Created on Sat Nov 24 20:56:36 2018

@author: YudongCai

@Email: yudongcai216@gmail.com
"""

import gzip
import click
import random


def random_select(probs):
    """
    probs: list
    beagle文件中的某个个体的概率, 3个一组
    如[0.941178, 0.058822, 0]
    Genotype likelihood for the major/major, major/minor, minor/minor
    """
    cutoff1 = probs[0]
    cutoff2 = cutoff1 + probs[1]
    random_val = random.random()
    if random_val < cutoff1:
        return '2,0'
    elif random_val < cutoff2:
        return '1,1'
    else:
        return '0,2'


def load_sites(keep):
    markers = []
    with gzip.open(keep, 'rb') as f:
        for line in f:
            marker = line.decode().strip()
            markers.append(marker)
    return set(markers)



@click.command()
@click.option('--beagle', help='beagle.gz file')
@click.option('--outfile', help='output file(.gz)')
@click.option('--markerfile', help='keep only markers listed in this gzipped file', default=None)
def main(beagle, outfile, markerfile):
    if markerfile:
        markers = load_sites(markerfile)
    with gzip.open(beagle, 'rb') as f1, gzip.open(outfile, 'wb') as f2:
        header = f1.readline().decode().strip().split()
        samples = header[3::3]
        outheader = ' '.join(samples) + '\n'
        f2.write(outheader.encode())
        if markerfile: # 筛点
            for line in f1:
                tline = line.decode().strip().split()
                marker = tline[0]
                if marker in markers:
                    counts = []
                    tline = [float(i) for i in tline[3:]]
                    for probs in zip(tline[::3], tline[1::3], tline[2::3]):
                        count = random_select(probs)
                        counts.append(count)
                    outline = ' '.join(counts) + '\n'
                    f2.write(outline.encode())
        else: # 不筛
            for line in f1:
                counts = []
                tline = [float(i) for i in line.decode().strip().split()[3:]]
                for probs in zip(tline[::3], tline[1::3], tline[2::3]):
                    count = random_select(probs)
                    counts.append(count)
                outline = ' '.join(counts) + '\n'
                f2.write(outline.encode())


if __name__ == '__main__':
    main()







