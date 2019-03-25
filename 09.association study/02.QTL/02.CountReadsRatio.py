# -*- coding: utf-8 -*-
"""
Created on Fri May  4 21:43:23 2018
@Mail: minnglee@163.com
@Author: Ming Li
"""

import sys,os,logging,click

logging.basicConfig(filename=os.path.basename(__file__).replace('.py','.log'),
                    format='%(asctime)s: %(name)s: %(levelname)s: %(message)s',level=logging.DEBUG,filemode='w')
logging.info(f"The command line is:\n\tpython3 {' '.join(sys.argv)}")

def LoadVCF(File,Output):
    Dict = {}
    Chr = set()
    VCF = os.popen(f'less {File}')
    for line in VCF:
        line = line.strip().split()
        if line[0].startswith('##'): continue
        if line[0].startswith('#'):
            Output.write('#Chr\tStart\tEnd\tDivergenceSiteNum\t{0}\n'.format('\t'.join(line[9:])))
            IndividualNum =len(line[9:])
            continue
        Chr.add(line[0])
        line[1] = int(line[1])
        GTList = [line[1]]
        for GT in line[9:]:
            GT = GT.split(':')
            RefAltCount = GT[1].split(',')
            if GT[0] == './.': GTCount = [0,0]
            elif GT[0] == '0/0': GTCount = [int(RefAltCount[0]),0]
            elif GT[0] == '1/1': GTCount = [0,int(RefAltCount[-1])]
            elif GT[0] == '0/1': GTCount = [int(RefAltCount[0]),int(RefAltCount[-1])]
            GTList.append(GTCount)
        Dict.setdefault(line[0],[]).append(GTList)
    return Dict,Chr,IndividualNum
def SlideWindow(start,end,step,Index,InList):
    SubList = []
    for i in range(Index,len(InList)):
        if InList[i][0] <= start + step:
            SubList.append(InList[i])
            Index = i+1
        elif InList[i][0] <= end: SubList.append(InList[i])
        else: return SubList,Index,True
    return SubList,Index,False
def Calculate_Windowed_Reads(CHR,WinSnpList,IndividualNum):
    TotalDivergenceSite = len(WinSnpList)
    CountDict = {}
    if TotalDivergenceSite == 0:
        return 0,0
    for snp in WinSnpList:
        for i in range(IndividualNum):
            if i not in CountDict: CountDict[i] = [0,0]
            CountDict[i][0] += snp[i+1][0]
            CountDict[i][1] += snp[i+1][1]
    GTCount = []
    for i in range(IndividualNum):
        GTCount.append('{0}|{1}'.format(CountDict[i][0],CountDict[i][1]))
    return TotalDivergenceSite,GTCount
@click.command()
@click.option('-w','--window',help='Tne window size',type=int,default=100000)
@click.option('-s','--step',help='Step Size',type=int,default=100000)
@click.option('-v','--vcf',help='Tne VCF file;need .gz file',type=str,required=True)
@click.option('-o','--output',help='Output file',type=click.File('w'),required=True)
def main(window,step,vcf,output):
    GTDict,ChrSet,IndividualNum = LoadVCF(vcf,output)
    for CHR in ChrSet:
        start,end = 0,window
        Index = 0
        NotEnd = True
        while NotEnd :
            WinSnpList,Index,NotEnd = SlideWindow(start,end,step,Index,GTDict[CHR])
            TotalDivergenceSite,GTCount = Calculate_Windowed_Reads(CHR,WinSnpList,IndividualNum)
            if TotalDivergenceSite > 0:
                output.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(CHR,start,end,TotalDivergenceSite,'\t'.join(GTCount)))
            start += step
            end += step
if __name__=='__main__':
    main()
