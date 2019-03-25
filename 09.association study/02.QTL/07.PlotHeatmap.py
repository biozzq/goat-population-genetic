# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Created on Tue May 15 22:18:08 2018
@Mail: minnglee@163.com
@Author: Ming Li
"""

import sys,os,logging,click
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns; sns.set()
import pandas as pd

logging.basicConfig(filename='{0}.log'.format(os.path.basename(__file__).replace('.py','')),
                    format='%(asctime)s: %(name)s: %(levelname)s: %(message)s',level=logging.DEBUG,filemode='w')
logging.info('The command line is:\n\t{0}'.format('python3 {0}'.format(' '.join(sys.argv))))

def LoadAllele(Input):
    '''
    #Chr	Start	End	DivergenceSiteNum	35	36	38	42	44	47	48	50	52	53	54	55	57
    1		0	100000	323			H	H	H	A	A	H	H	A	H	A	A	A	A
    1		100000	200000	256			H	H	H	A	A	H	H	A	H	A	A	A	A
    '''
    AlleleList = []
    Name = []
    PosList = []
    Pos = 0
    ChrStartList = []
    PreviousChr = ''
    for line in Input:
        if line.startswith('#'):
            line = line.strip().split()
            Name = line[4:]
            continue
        line = line.strip().split()
        Allele = []
        Pos += 1
        PosList.append(Pos)
        if line[0] != PreviousChr:
            ChrStartList.append(Pos)
            PreviousChr = line[0]
        for i in range(4,len(line)):
            if line[i] == 'A': Allele.append(0)
            elif line[i] == 'H': Allele.append(1)
            elif line[i] == 'B': Allele.append(2)
            else: print(line[i])
        AlleleList.append(Allele)
    return AlleleList,Name,PosList,ChrStartList
@click.command()
@click.option('-o','--output',help='the output file',type=str,required=True)
@click.option('-i','--input',help='the input file',type=click.File('r'),required=True)
def main(input,output):
    plt.switch_backend('Agg')
    lightyellow = (247/255, 248/255, 213/255)
    darkred = (191/255, 56/255, 38/255)
    black = (0/255, 0/255, 0/255)
    Allele,Name,Pos,ChrStartList = LoadAllele(input)
    LenList = []
    for i in Allele:
        LenList.append(len(i))
#    print(LenList)
#    print(Allele[0])
    Allele = np.array(Allele,dtype=np.int32).T
    Individual = np.repeat(Name,len(Pos))
    Position = Pos * len(Name)
    value = []
    for i in Allele:
        for j in i:
            value.append(j)
    df=pd.DataFrame({'Chromosome': Position, 'Sample ID': Individual, 'value': value })
    
    figsize = 20,6
#    FontSize = 80
    figure, ax = plt.subplots(figsize=figsize)
    df_wide=df.pivot_table( index='Sample ID', columns='Chromosome', values='value' )
    sns.heatmap(df_wide,cmap=[lightyellow,darkred,black],ax=ax,xticklabels=False)
    for ChrStart in ChrStartList[1:]:
        plt.axvline(ChrStart,linewidth=1, color='black',linestyle = "solid") 
    plt.tight_layout() #防止x,y轴超出画布或重叠
    ####保存图片
    format = 'pdf'
    plt.savefig(f'{output}.{format}',format = format)
    plt.show()
if __name__ == '__main__':
    main()
