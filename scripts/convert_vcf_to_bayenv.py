import argparse
import numpy as np
import vcf
import sys


def get_args():
    parser = argparse.ArgumentParser(description='convert vcf.gz to a bayenv SNPSFILE', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--vcf', type=str, required=True, help='vcf.gz file')
    parser.add_argument('--samples',type=str, required=True, help='text file of all samples in vcf.gz file', nargs='*')
    parser.add_argument('--o', type=str, required=True, help='output filename')

    return parser.parse_args()


def write_info(ids, out):

    with open(out+'_info.txt', 'w') as f:
        for snp_id in ids:
            f.write(snp_id)


def write_snps(snps, out):
    
    with open(out+'.SNPS', 'w') as f:
        for i in range(snps.shape[0]):
            line = ''
            for j in range(snps.shape[1]):
                if j < snps.shape[1] - 1:
                    line += f'{snps[i,j]}\t'
                else:
                    line += f'{snps[i,j]}\n'
            f.write(line)


def get_snps(v, pop_samples):
    
    snps = {i:[] for i in range(len(pop_samples))}
    snps_info = {i:[] for i in range(len(pop_samples))}
    print(len(pop_samples)) 
    
    for record in v:
            
        #we only want polymorphic sites that are also SNPs per reqs of bayenv
        if record.is_monomorphic or record.is_snp is False:
            continue
        
        fixed_refs=0
        fixed_alts=0
        #The counts of allele 1 and allele 2 are assumed to sum to the 
        #sample size typed at this SNP in this population 
        #(i.e. the total sample size excluding missing data)
        for i, samples in enumerate(pop_samples):
            refs = 0
            alts = 0

            for sample in samples:
                try:
                    res = record.genotype(sample)['GT']
                except:
                    #there's going to be some missing data which is expected
                    continue

                ref_count = res.count('0')
                alt_count = res.count('1')
                
                if ref_count is not None:
                    refs += ref_count
                if alt_count is not None:
                    alts += alt_count
            
            if refs == 0:
                fixed_refs += 1
            if alts == 0:
                fixed_alts += 1
            
            snps_info[i].append(f'{record.CHROM}\t{record.POS}\t{record.REF},\t{record.ALT}\n')
            #since we're going to stack them as columns later
            snps[i].append(refs)
            snps[i].append(alts)
            
            #ideally we fix this earlier in the pipeline
            #check that it isn't fixed
            #also we know that we're at the last population for this to be true
            if fixed_refs == len(pop_samples) or fixed_alts == len(pop_samples):
                for i in snps:
                    #remove refs and alts that were appended from all populations
                    #maybe we change this whole things to use .extend instead of popping
                    snps[i] = snps[i][:-2]
                    snps_info[i].pop()

    #    count+=1
    #    if count == cut:
    #        break

    return snps, snps_info


def read_samples(samples):
    
    ids = []
    with open(samples, 'r') as f:
        for s in f.readlines():
            ids.append(s.strip())
    
    return ids


def main():
    args = get_args()
    
    #probably need to use tbi and not use 'rb'
    
    pop_samples = [read_samples(samples) for samples in args.samples]
    #print(samples[:5])

    v = vcf.Reader(open(args.vcf, 'rb'))
    
    pop_snps, pop_snps_info = get_snps(v, pop_samples)
        
    for n in pop_snps_info:
        write_info(pop_snps_info[n], args.o+f'_{n}')

    #stack snps arrays
    #in order of population
    snps_fmt = np.column_stack([pop_snps[i] for i in pop_snps])
    print(snps_fmt.shape)
    
    write_snps(snps_fmt, args.o)


if __name__ == '__main__':
    main()

