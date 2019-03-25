# -*- coding: utf-8 -*-
"""
Created on Fri May  4 21:43:23 2018
@Mail: minnglee@163.com
@Author: Ming Li
"""

import sys,os,logging,click

logging.basicConfig(filename='{0}.log'.format(os.path.basename(__file__).replace('.py','')),
                    format='%(asctime)s: %(name)s: %(levelname)s: %(message)s',level=logging.DEBUG,filemode='w')
logging.info('The command line is:\n\t{0}'.format('python3 {0}'.format(' '.join(sys.argv))))

@click.command()
@click.option('-i','--input',help='Tne divergence site file',type=click.File('r'),required=True)
@click.option('-o','--output',help='Output file',type=click.File('w'),required=True)
def main(input,output):
    '''
    Input Format:
    #Chr    Start   End     DivergenceSiteNum       35      36      38      42      44      47
    28      0       1000000 1853                    0|0     0|1     0|1     0|1     0|0     0|0
    28      1000000 2000000 2099                    0|0     0|1     0|1     0|1     0|1     0|0
    '''
    for line in input:
        if line.startswith('#'):
            line = line.strip().split()
            PreviousLine = line
            continue
        line = line.strip().split()
        if line[0] != PreviousLine[0]:
            output.write('{0}\n'.format('\t'.join(PreviousLine)))
            PreviousLine = line
            continue
        NotSame = 0
        for i in range(4,len(line)):
            if line[i] != PreviousLine[i]: NotSame += 1
        if NotSame != 0:
            output.write('{0}\n'.format('\t'.join(PreviousLine)))
            PreviousLine = line
        else:
            PreviousLine[2] = line[2]
    output.write('{0}\n'.format('\t'.join(PreviousLine)))
if __name__=='__main__':
    main()
