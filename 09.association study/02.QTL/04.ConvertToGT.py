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
@click.option('-i','--Input',help='Tne divergence site file',type=click.File('r'),required=True)
@click.option('-o','--Output',help='Output file',type=click.File('w'),required=True)
@click.option('-t','--Threshold',help='Threshold',type=float,default=0.25)
def main(input,output,threshold):
    '''
    #Chr    Start   End     DivergenceSiteNum       35      36      38      42      44
    28      0       1000000 1853                    1560|0  685|592 977|810 736|643 1811|0
    28      1000000 2000000 2099                    1779|3  838|686 989|987 782|644 1422|622
    '''
    for line in input:
        if line.startswith('#'):
            output.write(line)
            continue
        line = line.strip().split()
        output.write('\t'.join(line[:4]))
        for i in range(4,len(line)):
            ReadsCount = line[i].split('|')
            Ratio = int(ReadsCount[0])/(int(ReadsCount[0])+int(ReadsCount[1]))
            if Ratio > 1 - threshold:
                output.write('\tA')  #0|0
            elif threshold <= Ratio <= 1 - threshold:
                output.write('\tH')  #0|1
            elif Ratio < threshold:
                output.write('\tB')  #1|1
        output.write('\n')
if __name__ == '__main__':
    main()
