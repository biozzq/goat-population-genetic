# -*- coding: utf-8 -*-
"""
Created on Fri May  4 21:43:23 2018
@Mail: minnglee@163.com
@Author: Ming Li
"""

import sys,os,logging,click
import gzip
import msprime

logging.basicConfig(filename='{0}.log'.format(os.path.basename(__file__).replace('.py','')),
                    format='%(asctime)s: %(name)s: %(levelname)s: %(message)s',level=logging.DEBUG,filemode='w')
logging.info('The command line is:\n\t{0}'.format('python3 {0}'.format(' '.join(sys.argv))))

@click.command()
@click.option('-s','--seed',help='The number of seed (we used seeds 1-100)',type=int,default=1)
#@click.option('-o','--output',help='The prefix of output file',type=str,required=True)
@click.option('-r','--recombination',help='The recombination ratio',type=float,default=1e-8)
@click.option('-i','--introgression',help='The introgression ratio',type=float,default=0.001)
def main(seed,recombination,introgression):
    if not os.path.exists('vcf'): os.system('mkdir vcf')
    OutputPrefix = f'vcf/Simulation.s{seed}.r{recombination}.i{introgression}'

    nbez = 48 # num. sampled wild goat haplotypes 
    ndom = 320 # num. sampled domestic goat haplotypes 
    nibex = 20 # num. sampled C.ibex haplotypes 
    nbp = 10000000 # simulated 10 Mb 
    rho = recombination # recombination rate 
    mu = 4.32e-9 # mutation rate of base simulation (can downsample to mu = 1.2e-8 for analysis) 
 
    Ne_ibex = 20000 # C.ibex effective size 
    Ne_bez = 100000 # wild goat eff. size 
    Ne_dom = 15000 # domestic goat eff. size 
    r = 0 # growth rate since domestication 
 
    # times are in generations before present 
    TIBEX = 750000 # C.ibex-goats split time 
    TDom = 5000 # wild - dom split time
    Tintrostart = 4000 # start of introgression 
    Tintroend = Tintrostart-10 # end introgression 
#    Tgrowth = 100 # time of recent growth  
 
    # migration rates are proportion of population made of new immigrants each generation 
#    mWIDO = 5.0e-8 # mig. rate betw. wild. and dom. 
    mWIDO = 0
    mintro = introgression # introgression from C.ibex into domestic goat 
 
    pop_config = [ 
            msprime.PopulationConfiguration(sample_size = None, initial_size = Ne_dom,growth_rate = r), 
            msprime.PopulationConfiguration(sample_size = None, initial_size = Ne_bez, growth_rate = r), 
            msprime.PopulationConfiguration(sample_size = None, initial_size = Ne_ibex)]
 
    samples = [msprime.Sample(0,0) for i in range(ndom)] + [msprime.Sample(1,0) for i in range(nbez)] + [msprime.Sample(2,0) for i in range(nibex)]
 
    # migration between wild and dom 
    mig_mat = [[0,mWIDO,0],[mWIDO,0,0],[0,0,0]]
    
 
    # onset of growth with agriculture 
#    growth_event = [msprime.PopulationParametersChange(time = Tgrowth, growth_rate = 0.0)] 
 
    # neanderthal introgression 
    intro_event = [msprime.MigrationRateChange(time = Tintroend, rate = mintro, matrix_index = (0,2)), 
                   msprime.MigrationRateChange(time = Tintrostart, rate = 0, matrix_index = (0,2))]
 
    # domestication event
    dom_event = [msprime.MigrationRateChange(time=TDom-.01, rate=0), msprime.MassMigration(time=TDom, 
                 source=0, destination=1)] 
 
    # the split of C.ibex and bezoer
    ibex_event = [msprime.MassMigration(time=TIBEX, source=1, destination=2)] 
 
    events =  intro_event + dom_event + ibex_event
    
    '''
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    dd = msprime.DemographyDebugger(
        population_configurations=pop_config,
        migration_matrix=None,
        demographic_events=events)
    dd.print_history()
    '''
    # make the simulation 
    treeseq = msprime.simulate(population_configurations = pop_config, samples = samples, migration_matrix = mig_mat,
                               demographic_events = events, length = nbp, recombination_rate = rho,
                               mutation_rate = mu, random_seed = seed)
    
    # Get genetic map
    '''
    mapfile = open(f'{OutputPrefix}.map',"w")
    breakpoints = np.array(list(treeseq.breakpoints()))
    num_bins = 500
    y, bin_edges = np.histogram(breakpoints, num_bins, density=True)
    TotalMap = 0
    for i in range(1,len(bin_edges)):
        Map = y[i-1] * (bin_edges[i] - bin_edges[i-1])
        TotalMap += Map
        mapfile.write('1\t.\t{0}\t{1}\n'.format(round(TotalMap,6),int(bin_edges[i])))
    '''
    
    # output the vcf 
    vcffile = open(f'{OutputPrefix}.vcf',"w")
    treeseq.write_vcf(vcffile,2)
    os.system(f'gzip {OutputPrefix}.vcf')
 
    # determine intro. status for dom. hapl. 
    def node_get_pop(tree,node,admix_time, split_time): 
        # if the time of a node is < admix_time, its population is found by tracing up the tree 
        while tree.get_time(node) <= admix_time: 
            node = tree.get_parent(node) 
        # if the time of the parent node is > split_time, the returned pop. is -1 (unknown)  
        if tree.get_time(node) > split_time: 
            return -1 
        else: return tree.get_population(node)
     
    # output intro. status for dom. hapl.
    introfile = open(f'{OutputPrefix}.intro',"w") 
    eurohaps = range(ndom) 
    for tree in treeseq.trees(): 
        introgressed = [x for x in eurohaps if node_get_pop(tree,x,Tintrostart,TIBEX) == 2]
        if len(introgressed) > 0: 
            start,stop=tree.get_interval() 
            introfile.write('{0}\t{1}\t'.format(round(start,1),round(stop,1)))
            for x in introgressed:
                introfile.write('{0}\t'.format(x))
            introfile.write('\n')
if __name__ == '__main__':
    main()
