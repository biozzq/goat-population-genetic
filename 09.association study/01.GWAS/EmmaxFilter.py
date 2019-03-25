# -*- coding: utf-8 -*-
#!/usr/bin/env python3
'''
Created on Wed Sep 26 14:01:01 CST 2018
@Mail: minnglee@163.com
@Author: Ming Li
'''

import sys,os,logging,click

logging.basicConfig(filename='{0}.log'.format(os.path.basename(__file__).replace('.py','')),
                    format='%(asctime)s: %(name)s: %(levelname)s: %(message)s',level=logging.DEBUG,filemode='w')
logging.info(f"The command line is:\n\tpython3 {' '.join(sys.argv)}")

@click.command()
@click.option('-i','--input',type=click.File('r'),help='The input file',required=True)
@click.option('-o','--output',type=click.File('w'),help='The output file',required=True)
def main(input,output):
    line = input.readline()
    line = line.strip().split()
    output.write(f'{line[0]}\t{line[2]}\t{line[3]}\n')
    MaxChr = 0
    for line in input:
        line = line.strip().split()
        line[0] = line[0].replace('CHR','')
        line[0] = line[0].replace('Chr','')
        line[0] = line[0].replace('chr','')
        try:
            line[0] = int(line[0])
            if MaxChr < line[0] : MaxChr = line[0]
            output.write(f'{line[0]}\t{line[2]}\t{line[3]}\n')
        except:
            if line[0] == 'X': output.write(f'{MaxChr+1}\t{line[2]}\t{line[3]}\n')
            elif line[0] == 'X1': output.write(f'{MaxChr+1}\t{line[2]}\t{line[3]}\n')
            elif line[0] == 'X2': output.write(f'{MaxChr+2}\t{line[2]}\t{line[3]}\n')
            elif line[0] == 'Y': output.write(f'{MaxChr+2}\t{line[2]}\t{line[3]}\n')
            elif line[0] == 'W': output.write(f'{MaxChr+1}\t{line[2]}\t{line[3]}\n')
            elif line[0] == 'Z': output.write(f'{MaxChr+2}\t{line[2]}\t{line[3]}\n')
if __name__ == '__main__':
    main()
