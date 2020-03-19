""" Adapted from WASP/find_interesecting_snps.py """

import argparse
import glob
import gzip
from itertools import product, groupby
import numpy as np
import os
import pysam
import sys
import tables
import pandas as pd
#import importlib

import snptable
import util

#importlib.reload(snptable)


#BASEQ_DEFAULT = 10

class DataFiles(object):
    """Object to hold names and filehandles for all input / output 
    datafiles"""
    
    def __init__(self, bam_filename, is_sorted, is_paired,
                 output_dir=None, snp_dir=None):
        # flag indicating whether reads are paired-end
        self.is_paired = is_paired
        
        # prefix for output files
        self.prefix = None

        # name of input BAM filename
        self.bam_filename = bam_filename        
        # name of sorted input bam_filename
        # (new file is created if input file is not
        #  already sorted)
        self.bam_sort_filename = None
        # pysam file handle for input BAM
        self.input_bam = None

        # name of output  txt file for imbalance
        self.post_remapping_AI_filename = None

        # pysam file handles for output txt filename
        self.post_remapping_AI_txt = None

        # name of directory to read SNPs from
        self.snp_dir = snp_dir
            
        # separate input directory and bam filename
        tokens = self.bam_filename.split("/")
        bam_dir = "/".join(tokens[:-1])
        filename = tokens[-1]

        if output_dir is None:
            # if no output dir specified, use same directory as input
            # bam file
            output_dir = bam_dir
        else:
            if output_dir.endswith("/"):
                # strip trailing '/' from output dir name
                output_dir = output_dir[:-1]
                
        name_split = filename.split(".")
        if len(name_split) > 1:
           self.prefix = output_dir + "/" + ".".join(name_split[:-1])
        else:
            self.prefix = output_dir + "/" + name_split[0]
        
        # create output dir if does not exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        if not is_sorted:
            util.sort_bam(self.bam_filename, self.prefix)
            self.bam_sort_filename = self.prefix + ".sort.bam"
        else:
            self.bam_sort_filename = self.bam_filename

        self.post_remapping_AI_filename = self.prefix + ".post_remapping_AI.txt"

        sys.stderr.write("reading reads from:\n  %s\n" %
                         self.bam_sort_filename)
        sys.stderr.write("writing output files to:\n")

        self.input_bam = pysam.AlignmentFile(self.bam_sort_filename, "rb")
        self.post_remapping_AI_txt = open(self.post_remapping_AI_filename, "w+")
        sys.stderr.write("  %s\n  " % self.post_remapping_AI_filename)


    
        
    def close(self):
        """close open filehandles"""
        filehandles = [self.post_remapping_AI_txt]

        for fh in filehandles:
            if fh:
                fh.close()

class ReadStats(object):
    """Track information about reads and SNPs that they overlap"""

    def __init__(self):
        # number of read matches to reference allele
        self.ref_count = 0
        # number of read matches to alternative allele
        self.alt_count = 0
        # number of reads that overlap SNP but match neither allele
        self.other_count = 0
        
        # number of reads discarded becaused not mapped
        self.discard_unmapped = 0
        
        # number of reads discarded because not proper pair
        self.discard_improper_pair = 0

        # number of reads discarded because mate unmapped
        self.discard_mate_unmapped = 0

        # paired reads map to different chromosomes
        self.discard_different_chromosome = 0

        # number of faked reads that dont match to original chrom/position
        self.discard_different_position = 0

        # number of reads discarded because overlap an indel
        self.discard_indel = 0

        # number of reads discarded because secondary match
        self.discard_secondary = 0

        # number of chimeric reads discarded
        self.discard_supplementary = 0

        # number of reads discarded because of too many overlapping SNPs
        self.discard_excess_snps = 0
        
        # number of reads discarded because too many allelic combinations
        self.discard_excess_reads = 0

        # when read pairs share SNP locations but have different alleles there
        self.discard_discordant_shared_snp = 0
        
        # reads where we expected to see other pair, but it was missing
        # possibly due to read-pairs with different names
        self.discard_missing_pair = 0

        # number of single reads that need remapping
        self.remap_single = 0
        # number of read pairs kept
        self.remap_pair = 0

    def write(self, file_handle):
        sys.stderr.write("DISCARDED reads:\n"
                         "  unmapped: %d\n"
                         "  mate unmapped: %d\n"
                         "  improper pair: %d\n"
                         "  different chromosome: %d\n"
                         "  diffrent chromosome/position as orinal read: %d\n"
                         "  indel: %d\n"
                         "  secondary alignment: %d\n"
                         "  supplementary alignment: %d\n"
                         "  excess overlapping snps: %d\n"
                         "  excess allelic combinations: %d\n"
                         "  read pairs with discordant shared SNPs: %d\n"
                         "  missing pairs (e.g. mismatched read names): %d\n"
                         "REMAP reads:\n"
                         "  single-end: %d\n"
                         "  pairs: %d\n" %
                         (self.discard_unmapped,
                          self.discard_mate_unmapped,
                          self.discard_improper_pair,
                          self.discard_different_chromosome,
                          self.discard_different_position,
                          self.discard_indel,
                          self.discard_secondary,
                          self.discard_supplementary,
                          self.discard_excess_snps,
                          self.discard_excess_reads,
                          self.discard_discordant_shared_snp,
                          self.discard_missing_pair,
                          self.remap_single,
                          self.remap_pair))

        file_handle.write("read SNP ref matches: %d\n" % self.ref_count)
        file_handle.write("read SNP alt matches: %d\n" % self.alt_count)
        file_handle.write("read SNP mismatches: %d\n" % self.other_count)
        
        total = self.ref_count + self.alt_count + self.other_count
        if total > 0:
            mismatch_pct = 100.0 * float(self.other_count) / total
            if mismatch_pct > 10.0:
                sys.stderr.write("WARNING: many read SNP overlaps do not match "
                                 "either allele (%.1f%%). SNP coordinates "
                                 "in input file may be incorrect.\n" %
                                 mismatch_pct)
    



def parse_options():
    
    parser = argparse.ArgumentParser(description="Looks for SNPs and indels "
                                     "overlapping reads and counts reads mapping to ref or "
                                     "alt allele for each SNP to determine "
                                     "allelic imbalance. "
                                     "Reads that overlap indels are not "
                                     "counted.")
                                   

    parser.add_argument("--is_paired_end", "-p", action='store_true',
                        dest='is_paired_end', 
                        default=False,
                        help=("Indicates that reads are paired-end "
                              "(default is single)."))
    
    parser.add_argument("--is_sorted", "-s", action='store_true',
                        dest='is_sorted', 
                        default=False,
                        help=('Indicates that the input BAM file'
                              ' is coordinate-sorted (default '
                              'is False).'))
    
    parser.add_argument("--output_dir", default=None,
                        help="Directory to write output files to. If not "
                        "specified, output files are written to the "
                        "same directory as the input BAM file.")

    parser.add_argument("--snp_dir", action='store', default=None,
                        help="Directory containing SNP text files "
                        "This directory should contain one file per "
                        "chromosome named like <chr>.snps.txt.gz. "
                        "Each file should contain 3 columns: position "
                        "RefAllele AltAllele. ")
                        
    parser.add_argument("bam_filename", action='store',
                        help="Coordinate-sorted input BAM file "
                        "containing mapped reads.")
    
        
    options = parser.parse_args()
    
    return options


def filter_count_ref_alt_matches(read, read_stats, snp_tab, snp_idx, read_pos):
    
    #base_idx = [i for i in range(len(base_q)) if base_q[i] >= base_qual]

    # update based on qual
    snp_idx = [snp_idx[x] for x in base_idx]
    snp_read_pos = [read_pos[x] for x in base_idx]
    #base_q = [base_q[x] for x in base_idx]

    ref_alleles = snp_tab.snp_allele1[snp_idx]
    alt_alleles = snp_tab.snp_allele2[snp_idx]

    to_rem = []
    
    for i in range(len(snp_idx)):
        ref = ref_alleles[i].decode("utf-8")
        alt = alt_alleles[i].decode("utf-8")
        
        if ref == read.query_sequence[snp_read_pos[i]-1]:
            # read matches reference allele
            read_stats.ref_count += 1
        elif alt == read.query_sequence[snp_read_pos[i]-1]:
            # read matches non-reference allele
            read_stats.alt_count += 1
        else:
            # read matches neither ref nor other
            read_stats.other_count += 1
            to_rem.append(i)

    if len(to_rem):
        for i in sorted(to_rem, reverse=True): # remove in descending order so indices are not 
            del snp_idx[i]
            del snp_read_pos[i]
            #del base_q[i]
            

    return snp_idx, snp_read_pos #, base_q



def count_reads(chrom, reads,  pair_snp_read_pos, pair_snp_idx, ref_alleles, alt_alleles, genome_pos, df, read_stats):
    """ Count reads mapping to ref or alt allele for those reads that map to the same chromosome and position as original read.
    Handles the possibility of shared SNPs among the pairs (ie doesn't treat them as independent, has to be the same base).
    Returns None when the original read pair has discordant alleles at shared SNPs. Also keeps records for reads overlapping ref and alt alleles by updating df """

    snp_read_pos = pair_snp_read_pos[:]

    # paired reads only
    if len(snp_read_pos) == 2:
        for i in range(len(snp_read_pos)):
            # first check if the same snp is in both reads and if so if it has the same allele, otherwise discard read
            # get the indices of the SNP indices that are in both reads
            idx_idxs = np.nonzero(np.in1d(pair_snp_idx[i], pair_snp_idx[(i+1) % 2]))[0]
            # now, use the indices in idx_idxs to get the relevant snp positions
            # and convert positions to indices
            snp_read_pos[i] = np.array(snp_read_pos[i], dtype=int)[idx_idxs] - 1

        # check: are there discordant alleles at the shared SNPs?
        # if so, discard these reads
        if  slice_read(reads[0], snp_read_pos[0]) != slice_read(reads[1], snp_read_pos[1]) :
            read_stats.discard_discordant_shared_snp += 1
            return 
        
  # shared SNPs among reads, avoid double-counting
           
    for i in range(len(pair_snp_read_pos)):    
        if len(pair_snp_read_pos[i]) != 0:
            for x in range(len(pair_snp_read_pos[i])):
                idx = pair_snp_read_pos[i][x]-1
                if reads[i][idx] == ref_alleles[i][0][x].decode("utf-8"):  # index 0 is to unlist
                    # add counts to df, both for ref and alt one for each read, except when read2 is matching the same snp in read 1
                    if not (len(genome_pos[0]) and len(genome_pos[1])) :
                        df.loc[(chrom, genome_pos[i][0][x] ), ['NREF'] ] +=1
                    else:
                        if not (i == 1 and genome_pos[i][0][x] in genome_pos[0][0]):
                            df.loc[(chrom, genome_pos[i][0][x] ), ['NREF'] ] +=1
                    
                else:               
                    if not (len(genome_pos[0]) and len(genome_pos[1])) :
                        df.loc[(chrom, genome_pos[i][0][x] ), ['NALT'] ] +=1
                    else:
                        if not (i == 1 and genome_pos[i][0][x] in genome_pos[0][0]):
                            df.loc[(chrom, genome_pos[i][0][x] ), ['NALT'] ] +=1
    return 


def filter_reads(files):
    cur_chrom = None
    cur_tid = None
    seen_chrom = set([])

    snp_tab = snptable.SNPTable()
    read_stats = ReadStats()
    read_pair_cache = {}
    cache_size = 0
    read_count = 0

    # data frame to collect AI across all chromosomes
    df = pd.DataFrame(columns=["POS", "CHROM", "REF", "ALT", "NREF", "NALT"])

    for read in files.input_bam:
        read_count += 1
        # if (read_count % 100000) == 0:
        #     sys.stderr.write("\nread_count: %d\n" % read_count)
        #     sys.stderr.write("cache_size: %d\n" % cache_size)

        # TODO: need to change this to use new pysam API calls
        # but need to check pysam version for backward compatibility
        if read.tid == -1:
            # unmapped read
            read_stats.discard_unmapped += 1
            continue

        if (cur_tid is None) or (read.tid != cur_tid):
            # when a read is "fake"  (read name: orginal_name.chrom.pos-pos, name can have "." or for single reads orginal_name.chrom.pos-) check whether it matches to the same chromosome and positions as original, otherwise discard
            l = str.split(read.qname, sep=".")           
            if len(l) >= 3 and "-" in l[-1]:
                coord = str.split(l[-1], sep="-")             
                if '' not in coord: # paired read
                    if not read.is_paired:
                        continue
                    if not read.is_proper_pair:
                        continue
                    lpos = [int(x) for x in coord]
                
                    if l[-2] != read.reference_name or read.pos+1 < lpos[0] or read.pos+1 > lpos[1]:
                        read_stats.discard_different_position += 1
                        continue

                else: # single end read
                
                    if l[-2] != read.reference_name or read.pos+1 != int(coord[0]):
                        read_stats.discard_different_chromosome += 1
                        continue     
                
            # this is a new chromosome
            cur_chrom = files.input_bam.getrname(read.tid)

            if len(read_pair_cache) != 0:
                sys.stderr.write("WARNING: failed to find pairs for %d "
                                 "reads on this chromosome\n" %
                                 len(read_pair_cache))
                read_stats.discard_missing_pair += len(read_pair_cache)
            read_pair_cache = {}
            cache_size = 0
            read_count = 0

            if cur_chrom in seen_chrom:
                # sanity check that input bam file is sorted
                raise ValueError("expected input BAM file to be sorted "
                                 "but chromosome %s is repeated\n" % cur_chrom)
            seen_chrom.add(cur_chrom)
            cur_tid = read.tid
            sys.stderr.write("starting chromosome %s\n" % cur_chrom)

            # use text files from SNP dir
            #snp_filename = "%s/%s.snps.txt.gz" % (files.snp_dir, cur_chrom)
            snp_filename = glob.glob(files.snp_dir + "/*" +  cur_chrom + "*")[0]
            sys.stderr.write("reading SNPs from file '%s'\n" % snp_filename)
            snp_tab.read_file(snp_filename)

            # add new fSNPs to df
            tempDF = pd.DataFrame(np.hstack((snp_tab.snp_pos[:,None])), columns=["POS"]) 
            tempDF['CHROM'] =  cur_chrom
            tempDF["REF"] = snp_tab.snp_allele1.astype('U13')
            tempDF["ALT"] = snp_tab.snp_allele2.astype('U13')
            tempDF['NREF'] = np.zeros(tempDF.shape[0], dtype=int)
            tempDF['NALT'] = np.zeros(tempDF.shape[0], dtype=int)
           
            df = df.append(tempDF)
            df.set_index(["CHROM", "POS"], inplace = True)
            df.sort_index(inplace=True)
            
            sys.stderr.write("processing reads\n")

        if read.is_secondary:
            # this is a secondary alignment (i.e. read was aligned more than
            # once and this has align score that <= best score)
            read_stats.discard_secondary += 1
            continue

        if read.is_supplementary:
            # this is a supplementary alignment (ie chimeric and not the representative alignment)
            read_stats.discard_supplementary += 1
            continue

        if read.is_paired:
            if read.mate_is_unmapped:
                # other side of pair not mapped
                # we could process as single... but these not likely
                # useful so discard
                # process_single_read(read, read_stats, files,
                #                     snp_tab, max_seqs, max_snps)
                read_stats.discard_mate_unmapped += 1
            elif(read.next_reference_name == cur_chrom or
                 read.next_reference_name == "="):
                # other pair mapped to same chrom

                # sys.stderr.write("flag: %s" % read.flag)
                if not read.is_proper_pair:
                    # sys.stderr.write(' => improper\n')
                    read_stats.discard_improper_pair += 1
                    continue
                # sys.stderr.write(' => proper\n')
                if read.qname in read_pair_cache:
                    # we already saw prev pair, retrieve from cache keeping track of original fastq files
                    if read.is_read1:
                        read1 = read
                        read2 = read_pair_cache[read.qname]
                    else:
                        read1 = read_pair_cache[read.qname]
                        read2 = read
                    del read_pair_cache[read.qname]
                    cache_size -= 1

                    if read2.next_reference_start != read1.reference_start:
                        sys.stderr.write("WARNING: read pair positions "
                                         "do not match for pair %s\n" %
                                         read.qname)
                    else:
                        # if 23115001 in read1.get_reference_positions() or 23115001 in read2.get_reference_positions():
                        #     print(read1.qname)
                        #     break
                        # if read.qname == 'ERR188315.9027537':
                        #     break
                        process_paired_read(read1, read2, read_stats, files, snp_tab, df)
                        

                else:
                    # we need to wait for next pair
                    read_pair_cache[read.qname] = read

                    cache_size += 1


            else:
                # other side of pair mapped to different
                # chromosome, discard this read
                read_stats.discard_different_chromosome += 1

        else:
            process_single_read(read, read_stats, files, snp_tab) ## make sure it works

    if len(read_pair_cache) != 0:
        sys.stderr.write("WARNING: failed to find pairs for %d "
                         "reads on this chromosome\n" %
                         len(read_pair_cache))
        read_stats.discard_missing_pair += len(read_pair_cache)

    read_stats.write(sys.stderr)

    # save df
    df.to_csv(files.post_remapping_AI_txt, sep=" ", index=True)




def slice_read(read, indices):
    """slice a read by an array of indices"""
    return "".join(np.array(list(read))[indices])

def process_paired_read(read1, read2, read_stats, files,
                        snp_tab,df):
    """Checks if either end of read pair overlaps SNPs or indels
    and writes read pair (or generated read pairs) to appropriate
    output files. Also add counts ref and alt alleles mapping reads"""

    pair_snp_idx = []
    pair_snp_read_pos = []
    ref_alleles = []
    alt_alleles = []
    genome_pos = []  # keep record to count mapping reads to snps
    

    for read in (read1, read2):
        # check if either read overlaps SNPs or indels
        # check if read overlaps SNPs or indels
        snp_idx, snp_read_pos, \
            indel_idx, indel_read_pos, base_q = snp_tab.get_overlapping_snps(read)

        if len(indel_idx) > 0:
            # for now discard this read pair, we want to improve this to handle
            # the indel reads appropriately
            read_stats.discard_indel += 2
            # TODO: add option to handle indels instead of throwing out reads
            return
            
        if len(snp_idx) > 0:
            # Generate a read pair with the opposite allele as the one in read but discarding reads that
            # dont overlap ref or alt allele
                ref_alleles.append([snp_tab.snp_allele1[snp_idx]])
                alt_alleles.append([snp_tab.snp_allele2[snp_idx]])
                genome_pos.append([snp_tab.snp_pos[snp_idx]])
 
                #new_reads.append(read_seqs)
                pair_snp_idx.append(snp_idx)
                pair_snp_read_pos.append(snp_read_pos)

        else:
            # no SNPs overlap this read by quality
            pair_snp_idx.append([])
            pair_snp_read_pos.append([])
            ref_alleles.append([])
            alt_alleles.append([])
            genome_pos.append([])
            
    if len(pair_snp_idx[0]) !=0 or len(pair_snp_idx[1]) != 0:
        reads=[read1.seq, read2.seq][:]
        count_reads(read1.reference_name, reads, pair_snp_read_pos, pair_snp_idx, ref_alleles, alt_alleles, genome_pos, df, read_stats)



def process_single_read(read, read_stats, files, snp_tab):
    """Check if a single read overlaps SNPs or indels, and writes
    this read (or generated read pairs) to appropriate output files"""

    ref_alleles = []
    alt_alleles = []
    genome_pos = []  # keep record to count mapping reads to snps
                
    # check if read overlaps SNPs or indels
    snp_idx, snp_read_pos, \
        indel_idx, indel_read_pos, base_q = snp_tab.get_overlapping_snps(read)

    
    if len(indel_idx) > 0:
        # for now discard this read, we want to improve this to handle
        # the indel reads appropriately
        read_stats.discard_indel += 1
        # TODO: add option to handle indels instead of throwing out reads
        return

    if len(snp_idx) > 0:
        # Generate a read pair with the opposite allele as the one in read but discarding reads that
        # dont overlap ref or alt allele and base quality below BASEQ
        
        ref_alleles.append([snp_tab.snp_allele1[snp_idx]])
        alt_alleles.append([snp_tab.snp_allele2[snp_idx]])
        genome_pos.append([snp_tab.snp_pos[snp_idx]])
        
        count_reads(read1.reference_name, [read.query_sequence],  snp_read_pos, snp_idx,
                                        ref_alleles, alt_alleles, genome_pos, df)

         



            
def main(bam_filenames, is_paired_end=False,
         is_sorted=False,
         output_dir=None,
         snp_dir=None):
    
    files = DataFiles(bam_filenames,  is_sorted, is_paired_end,
                      output_dir=output_dir,
                      snp_dir=snp_dir)
    
    filter_reads(files)

    files.close()
    
    

if __name__ == '__main__':
    sys.stderr.write("command line: %s\n" % " ".join(sys.argv))
    sys.stderr.write("python version: %s\n" % sys.version)
    sys.stderr.write("pysam version: %s\n" % pysam.__version__)
    sys.stderr.write("pytables version: %s\n" % tables.__version__)

    util.check_pysam_version()
    util.check_pytables_version()
    util.check_python_version()
        
    options = parse_options()
    
    main(options.bam_filename,
         is_paired_end=options.is_paired_end,
         is_sorted=options.is_sorted,
         output_dir=options.output_dir,
         snp_dir=options.snp_dir)


# bam_filenames = "/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/STAR2/HG00120/Aligned.sortedByCoord.out.bam"
# bam_filenames = "/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/inputs/bam/HG00120.chrom22.MapUnique.bam"

# is_sorted=True
# is_paired_end =True
# output_dir="/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias"
# snp_dir="/mrc-bsu/scratch/ev250/EGEUV1/quant/refbias/inputs/fSNP"
# base_qual = 10

#/mrc-bsu/scratch/ev250/EGEUV1/RNA_seq_fastq
# 'ERR188091.16512776'
# zcat HG00256.chrom22.MapUnique.remap.fq1.gz | awk 'NR==418282'
# CACGGACCTGCCAACATGGGGCGCGTTCGGACCAAAACCGTGAAGAAGGCGGCCCGGGTCATCATAGAAAAGGACC
# zcat HG00256.chrom22.MapUnique.remap.fq2.gz | awk 'NR==418282'
# CGGAGCTTTTTGCTGGGGATAATGGCGATCTCCTCGCACACGCGCTTGTTCGTGTGGGAGTCGTTGCCGAGGCGCG

# zcat ERR188091_1.fastq.gz | awk 'NR==66051102'
# CGGAGCTTTTTGCTGGGGATAATGGCGATCTCCTCGCACACGCGCTTGTTCGTGTGGGAGTCGTTGCCGAGGCGCG

# zcat ERR188091_2.fastq.gz | awk 'NR==66051102' 
# CACGGACCTGCCAACATGGGGCGCGTTCGGACCAAAACCGTGAAGAAGGCGGCCCGGGTCATCATAGAAAAGGACC

# # bam file STAR2
# CGGAGCTTTTTGCTGGGGATAATGGCGATCTCCTCGCACACGCGCTTGTTCGTGTGGGAGTCGTTGCCGAGGCGCG #read1
# GGTCCTTTTCTATGATGACCCGGGCCGCCTTCTTCACGGTTTTGGTCCGAACGCGCCCCATGTTGGCAGGTCCGTG #read2 rev-comp

# # bam file input
# CACGGACCTGCCAACATGGGGCGCGTTCGGACCAAAACCGTGAAGAAGGCGGCCCGGGTCATCATAGAAAAGGACC #read2
# CGCGCCTCGGCAACGACTCCCACACGAACAAGCGCGTGTGCGAGGAGATCGCCATTATCCCCAGCAAAAAGCTCCG #read1 rev-comp
