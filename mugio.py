# -*- coding: utf-8 -*-
"""
Created on Mon May  6 21:02:55 2019
# mugio - low as a verb

Version 1.0 - 11.20.19 - Public Release - (Mind Television)
    
@author: Pieter Spealman ps163@nyu.edu
"""
import os
import numpy as np
import argparse
import scipy.stats as stats
import random
import pickle
import math
import re
import subprocess
import pandas as pd
import json

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fastqfile_name')
parser.add_argument('-s', '--samfile_name')
parser.add_argument('-o', '--outfile_name')
parser.add_argument('-bam', '--bamfile_name')

parser.add_argument('-bkgd', '--background', action='store_true')
parser.add_argument('-lt', '--lower_threshold')
parser.add_argument('-mt', '--middle_threshold')
parser.add_argument('-ut', '--upper_threshold')

parser.add_argument('-burnout', '--burnout', action='store_true')
parser.add_argument('-min', '--minimum_reads')
parser.add_argument('-un', '--unchanged')
parser.add_argument('-step', '--steps')

parser.add_argument('-bprd', '--breakpoint_retrieval_and_definition', action='store_true')

parser.add_argument('-blast', '--reflection_point_blast', action='store_true')
parser.add_argument('-rp', '--reflection_point_seq')

parser.add_argument('-pp', '--plot_phred', action='store_true')
parser.add_argument('-nt_min', '--nt_min')
parser.add_argument('-nt_max', '--nt_max')
parser.add_argument('-mm', '--make_median', action='store_true')

parser.add_argument('-e', '--evaluate', action='store_true')
parser.add_argument('-n', '--name')
parser.add_argument('-b', '--breakpoint')
parser.add_argument('-bpf', '--breakpoint_file')
parser.add_argument('-snf', '--sniffle_file')
parser.add_argument('-x', '--filter_list_file')

parser.add_argument('-demo', '--run_demo', action='store_true')

args = parser.parse_args()

def os_mkdir(in_name):    
    if '/' in in_name:
        directory_name = in_name.rsplit('/',1)[0]
        
        if not os.path.exists(directory_name):
            os.makedirs(directory_name)

def run_demo():
    '''runs demo given a known set of files'''
    outline = ('\nMugio reguires the following programs:\n'+
               '\tgzip\t1.5+\n'+
               '\tpython\t2.7+ or 3.6+\n'+
               '\tsamtools\t1.6+\n'+
               '\tbedtools\t2.26.0+\n'+
               '\t(optional) The --blast command requires blast+ 2.9.0+\n'+
               '\n'+
               'Mugio reguires the following python packages:\n'+
               '\tos\n'+
               '\tnumpy\n'+
               '\targparse\n'+
               '\tscipy.stats\n'+
               '\trandom\n'+
               '\tpickle\n'+
               '\tmath\n'+
               '\tre\n'+
               '\tsubprocess\n'+
               '\tpandas\n'+
               '\tjson\n')
    print(outline)
    
    if not os.path.isfile('demo/demo.fastq'):
        outline = ('GUnzipping the demo.fastq file\n')
        print(outline)
        bash_command = ('\tgunzip -f demo/demo.fastq.gz')
        print(bash_command)
        output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
        if output:
            print(output)
        
    outline = ('Converting bam to sam\n')
    print(outline)
    bash_command = ('\tsamtools view -h -o demo/demo.sam demo/demo.bam')
    print(bash_command)
    output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
    if output:
        print(output)
        
    outline = ('Command --background calculates the frequency of low_phred scoring reads.\n'+
               '\tHere we use the --burnout option to halt the process if no change greater than 1% is observed over 100 reads')
    print(outline)
    bash_command = ('\tpython mugio.py --background -f demo/demo.fastq -s demo/demo.sam -o demo_output/demo_bkgrd')
    print(bash_command)
    output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
    if output:
        print(output)
        
    1/0
        
    outline = ('Command --bprd (breakpoint retrieval and definition) identifies loci likely to be reflection points in inverted duplications.')
    print(outline)
    bash_command = ('\tpython mugio.py -bprd -f demo/demo.fastq -s demo/demo.sam -bam demo/demo.bam -o demo_output/demo_bprd')
    print(bash_command)
    output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
    if output:
        print(output)

    outline = ("Command --evaluate calculates the correlation (Spearman's rho) between pre-breakpoint seqeunce length and post-breakpoint low scoring region length.")
    print(outline)
    bash_command = ('\tpython mugio.py --evaluate -bpf demo_output/demo_bprd_bprd.bed -f demo/demo.fastq -s demo/demo.sam -o demo/demo_bprd_lengths')
    print(bash_command)
    output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
    if output:
        print(output)
        
    outline = ("Command --evaluate calculates the correlation (Spearman's rho) between pre-breakpoint seqeunce length and post-breakpoint low scoring region length.")
    print(outline)
    bash_command = ('\tpython mugio.py --evaluate -snf demo/demo_sniffles.vcf -f demo/demo.fastq -s demo/demo.sam -o demo/demo_vcf_lengths')
    print(bash_command)
    output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
    if output:
        print(output)
          
def get_char(gc_cigar, runmode):
    #This identifies the hard or soft clip number preceeding 
    # the first and last matches in a sam file cigar line
    prechar = 0
    if runmode == 'first':
        for each_char in gc_cigar:
            prechar += 1
            if each_char.isalpha():
                if each_char == 'M':
                    return(False, 0)
                else:
                    return(True, prechar)
                
    if runmode == 'last':
        if gc_cigar[-1] == 'M':
            return(False, 0)
            
        else:
            gc_cigar = gc_cigar[:-1]
            for each_char in gc_cigar[::-1]:
                prechar += 1
                if each_char.isalpha():
                    return(True, prechar)
                            
def distance_in(strand, cigar, runmode):
    match = 0 
        
    temp_list = re.split('M|S|D|I|H', cigar)
            
    char_per_step = []
    for each_char in cigar:
        if each_char.isalpha():
            char_per_step.append(each_char)
                        
    index = -1       
    for each_step in range(len(char_per_step)):
        index +=1
        if runmode == 'reference':
            if (char_per_step[each_step] != 'I') and (char_per_step[each_step] != 'H') and (char_per_step[each_step] != 'S'):
                match += int(temp_list[each_step])
                
        if runmode == 'cigar':
            if (char_per_step[each_step] != 'D') and (char_per_step[each_step] != 'H') and (char_per_step[each_step] != 'S'):
                match += int(temp_list[each_step])
                
    return(match)
        
def distance_cigar_to_phred(strand, cigar, threshold, prechar, postchar):
    match = 0
    projection = 0
    index = -1
        
    temp_list = re.split('M|S|D|I|H', cigar)
            
    char_per_step = []
    for each_char in cigar:
        if each_char.isalpha():
            char_per_step.append(each_char)
    
    if strand == 0:
        for each_step in range(len(char_per_step)):
            index +=1
            
            if (char_per_step[each_step] != 'H') and (char_per_step[each_step] != 'S'):
                match += int(temp_list[each_step])
                
            if (char_per_step[each_step] != 'D'):
                projection += int(temp_list[each_step])
                if match >= threshold:
                    correction = (match-threshold)
                    return(projection-correction)
                        
                    
        correction = (postchar)
        return(projection-correction)
        
    if strand == 1:
        temp_list = temp_list[::-1]
        temp_list = temp_list[1:]
        char_per_step = char_per_step[::-1]

        for each_step in range(len(char_per_step)):
            index +=1
            
            if (char_per_step[each_step] != 'H') and (char_per_step[each_step] != 'S'):
                match += int(temp_list[each_step])
            
            if (char_per_step[each_step] != 'D'):
                projection += int(temp_list[each_step])
                if match >= threshold:
                    correction = (match-threshold)
                    return(projection-correction)
                    
        correction = (prechar)
        return(projection-correction)
        
def parse_starpoints(start, stop, chromo, cigar, starpoints_by_chromo):    
    starred = False
    starpower = 0
    char_per_step = []
    for each_char in cigar:
        if each_char.isalpha():
            char_per_step.append(each_char)
    
    if char_per_step[0] != 'M':
        starpower =  int(cigar.split(char_per_step[0])[0])      
        if starpower > 100:
            starred = True
            if chromo in starpoints_by_chromo:
                starpoints_by_chromo[chromo].add(start)
            else:
                starpoints_set = set()
                starpoints_set.add(start)
                starpoints_by_chromo[chromo] = starpoints_set
    
    if char_per_step[-1] != 'M':
        starpower = int(cigar.rsplit(char_per_step[-2], 1)[1].split(char_per_step[-1])[0])
        if starpower > 100:
            starred = True
            if chromo in starpoints_by_chromo:
                starpoints_by_chromo[chromo].add(stop)
            else:
                starpoints_set = set()
                starpoints_set.add(stop)
                starpoints_by_chromo[chromo] = starpoints_set
                
    return(starred, starpower, starpoints_by_chromo)
        
def convert_sort(outfile_name):
    case = True
    bash_command = ('samtools view -Sb {output_file}.sam > {output_file}_u.bam').format(output_file=outfile_name)
    output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
    if output:
        print(output)
        case = False

    bash_command = ('samtools sort -T tmp_sort -o {output_file}.bam {output_file}_u.bam').format(output_file=outfile_name)
    output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
    if output:
        print(output)
        case = False
        
    bash_command = ('samtools index {output_file}.bam').format(output_file=outfile_name)
    output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
    if output:
        print(output)
        case = False
        
    bash_command = ('rm {output_file}.sam {output_file}_u.bam').format(output_file=outfile_name)
    output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
    if output:
        print(output)
        case = False

    bash_command = ('bedtools genomecov -ibam {output_file}.bam -d > {output_file}.depth').format(output_file=outfile_name)
    output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
    if output:
        print(output)
        case = False
        
    return(case)

def distance_ref_to_cigar(strand, cigar, threshold, prechar, postchar):
    match = 0
    index = -1
        
    temp_list = re.split('M|S|D|I|H', cigar)
            
    char_per_step = []
    for each_char in cigar:
        if each_char.isalpha():
            char_per_step.append(each_char)
    
    if strand == 0:
        projection = -1*prechar
        for each_step in range(len(char_per_step)):
            index +=1
            projection += int(temp_list[each_step])
            if (char_per_step[each_step] != 'I') and (char_per_step[each_step] != 'H') and (char_per_step[each_step] != 'S'):
                match += int(temp_list[each_step])
                if match >= threshold:
                    correction = (match-threshold)
                    return(projection-correction)
                        
        correction = (postchar)
        return(projection-correction)
        
    if strand == 1:
        projection = -1*prechar
        for each_step in range(len(char_per_step)):
            index +=1
            
            projection += int(temp_list[each_step])
            if (char_per_step[each_step] != 'I') and (char_per_step[each_step] != 'H') and (char_per_step[each_step] != 'S'):
                match += int(temp_list[each_step])
                if match >= threshold:
                    correction = (match-threshold)
                    return(projection-correction)
    
        correction = (postchar)
        return(projection-correction)
                                
def parallel_construction(strand, cigar, len_phred_line, sam_start, sam_stop, brk_start, brk_stop):
    '''
    Use the coordinates of the breakpoint in the reference genome to lookup where the breakpoint is in the phred scores. 
                            |----B------|
    samfile coordinates:    s----B------e 
    cigar                  SMDIMIDIMIMDIMH
    phred             ?(*&^)!@#$^&*%*==--!@#$     
    
    The we know where the reference start (s_ref) is in the cigar, it is the first M. The distance between s_ref and B_ref in 
    the reference space is related to their distance in cigar space by counting the matches (M) and deletions (D) cigar scores. 
    The location of B_phred is the sum of the preceding soft (S) or hard (H) clip before the s_cigar and all (M) and inserts (I) 
    between the s_cigar and B_cigar.
    
    # Is this strand specific? 
    
                            |----B------|
    samfile coordinates:    s----B------e 
    cigar                  SMDIMIDIMIMDIMH
    phred             ?(*&^)!@#$^&*%*==--!@#$
    #                      
    cigar_rev              HMIDMIMIDIMIDMS
    phred_rev         $#@!--==*%*&^$#@!)^&*(?
    
    In the opposite direction the samfile retains + strand orientation, so too does the cigar, but the relevant phred scroe in the 
    fastq file is reversed relative to this positive mapping. By reversing the direction of the phred (now phred_rev) we can use the 
    same incremental steps as described previously.   
    
    '''
    if strand == 0:
        preceede, first_cut = get_char(cigar, 'first')
        if not preceede:
            prechar = 0
        else:
            prechar = int(cigar[:first_cut-1])
                        
        preceede, last_cut = get_char(cigar, 'last')
        if preceede:
            postchar = int(cigar[-1*last_cut:-1])
        else:
            postchar = 0
            
        ref_distance_between_samstart_and_breakpoint = abs(sam_start-brk_start)
        cigar_distance_between_samstart_and_breakpoint = distance_ref_to_cigar(strand, cigar, ref_distance_between_samstart_and_breakpoint, prechar, postchar)
            
        phred_distance_between_samstart_and_breakpoint = distance_cigar_to_phred(strand, cigar, cigar_distance_between_samstart_and_breakpoint, prechar, postchar)

        ref_distance_between_samstop_and_breakpoint = abs(sam_start-brk_stop)
        cigar_distance_between_samstop_and_breakpoint = distance_ref_to_cigar(strand, cigar, ref_distance_between_samstop_and_breakpoint, prechar, postchar)

        phred_distance_between_samstop_and_breakpoint = distance_cigar_to_phred(strand, cigar, cigar_distance_between_samstop_and_breakpoint, prechar, postchar)
        
        return(phred_distance_between_samstart_and_breakpoint, phred_distance_between_samstop_and_breakpoint)

    if strand == 1:
        preceede, first_cut = get_char(cigar, 'first')
        if not preceede:
            prechar = 0
        else:
            prechar = int(cigar[:first_cut-1])

        preceede, last_cut = get_char(cigar, 'last')
        if preceede:
            postchar = int(cigar[-1*last_cut:-1])
        else:
            postchar = 0
            
        ref_distance_between_samstart_and_breakpoint = abs(sam_stop-brk_start)
        cigar_distance_between_samstart_and_breakpoint = distance_ref_to_cigar(strand, cigar, ref_distance_between_samstart_and_breakpoint, prechar, postchar)
        
        phred_distance_between_samstart_and_breakpoint = distance_cigar_to_phred(strand, cigar, cigar_distance_between_samstart_and_breakpoint, prechar, postchar)
        ref_distance_between_samstop_and_breakpoint = abs(sam_stop-brk_stop)
        cigar_distance_between_samstop_and_breakpoint = distance_ref_to_cigar(strand, cigar, ref_distance_between_samstop_and_breakpoint, prechar, postchar)
            
        phred_distance_between_samstop_and_breakpoint = distance_cigar_to_phred(strand, cigar, cigar_distance_between_samstop_and_breakpoint, prechar, postchar)
        
        return(phred_distance_between_samstop_and_breakpoint,  phred_distance_between_samstart_and_breakpoint)
 
def ord_phred(phred_line):
    phred_list = []
            
    for each_nt in phred_line:
        phred_list.append((ord(each_nt)-33))
        
    return(phred_list)

def unpackbits(x,num_bits=12):
    xshape = list(x.shape)
    x = x.reshape([-1,1])
    to_and = 2**np.arange(num_bits).reshape([1,num_bits])
    upb = (x & to_and).astype(bool).astype(int).reshape(xshape + [num_bits])

    #0  (rp)    read_paired
    #1  (rmp)    read_mapped_in_proper_pair
    #2  (ru)    read_unmapped
    #3  (mu)    mate_unmapped
    #4  (rrs)    read_reverse_strand
    #5  (mrs)    mate_reverse_strand
    #6  (fip)    first_in_pair
    #7  (sip)    second_in_pair
    #8  (npa)    not_primary_alignment
    #9  (rfp)    read_fails_platform
    #10 (pcr)    read_is_PCR_or_optical_duplicate
    #11 (sa)    supplementary_alignment
    
    """ DISCORDANT definition (from samblaster)
        Both side of the read pair are mapped (neither FLAG 0x4 or 0x8 is set).
        The properly paired FLAG (0x2) is not set.
        Note: We implemented an additional criteria to distinguish between strand re-orientations and distance issues
        Strand Discordant reads must be both on the same strand.
    """
        
    """ SPLIT READS
        Identify reads that have between two and --maxSplitCount [2] primary and supplemental alignments.
        Sort these alignments by their strand-normalized position along the read.
        Two alignments are output as splitters if they are adjacent on the read, and meet these criteria:
            each covers at least --minNonOverlap [20] base pairs of the read that the other does not.
            the two alignments map to different reference sequences and/or strands. 
            the two alignments map to the same sequence and strand, and represent a SV that is at least --minIndelSize [50] in length, 
            and have at most --maxUnmappedBases [50] of un-aligned base pairs between them.
        Split read alignments that are part of a duplicate read will be output unless the -e option is used.
    """
    return(upb) 
    
def rolling_median(x, y, runmode = 'graphical_curve'):
    if runmode == 'graphical_curve':
        new_x = []
        new_y = []
        
        for nt in range(len(y)-100):
            if len(y[nt:nt+100]) > 0:
                new_y.append(np.median(y[nt:nt+100]))
                new_x.append(x[nt])
                    
        return(new_x, new_y)    
        
    if runmode == 'plot_phred':
        new_x = []
        new_y = []
                
        for nt in range(len(y)-1000):
            if len(y[nt:nt+1000]) > 0:
                new_y.append(np.median(y[nt:nt+1000]))
                new_x.append(x[nt])
                    
        return(new_x, new_y) 

def graphical_curve(curve_dict, uid_ct, strand, brkpt, is_outfile_name, output_all=True):   
    
    if not output_all:
        check = curve_dict[(len(curve_dict)-1)][0]
        if check != 'g':
            for index in [0, len(curve_dict)-1]:
                color = curve_dict[index][0]
                plt.plot(curve_dict[index][1], curve_dict[index][2], c=color, alpha=0.1)
                c_x, c_y = rolling_median(curve_dict[index][1], curve_dict[index][2])
                plt.plot(c_x, c_y, c=color)
                
                if color == 'b':
                    rel_brk_start = curve_dict[index][3]
                    rel_brk_stop = curve_dict[index][4]
                    simple_break = round((rel_brk_start+rel_brk_stop)/2)
                    outline = ('M.pt:{}\n{}').format(simple_break, brkpt)
                    plt.text(simple_break, 0, outline)
                    plt.axvline(x=simple_break, linewidth=3, color='red', alpha=0.2)
                    plt.scatter(simple_break, 0, c='red', s=23)
                
                if color != 'b':
                    outline = ('Pre {}\nPost {}').format(curve_dict[index][5], max(curve_dict[index][1])-min(curve_dict[index][1]))
                    plt.text(0, 25, outline)
                        
            plt.savefig(is_outfile_name + '_' + uid_ct + '_' + str(strand) +'.png')
            plt.show()
            plt.close()
    else:
        for index in [0, len(curve_dict)-1]:
            color = curve_dict[index][0]
            plt.plot(curve_dict[index][1], curve_dict[index][2], c=color, alpha=0.1)
            c_x, c_y = rolling_median(curve_dict[index][1], curve_dict[index][2])
            plt.plot(c_x, c_y, c=color)
            
            if color == 'b':
                rel_brk_start = curve_dict[index][3]
                rel_brk_stop = curve_dict[index][4]
                simple_break = round((rel_brk_start+rel_brk_stop)/2)
                outline = ('M.pt:{}\n{}').format(simple_break, brkpt)
                plt.text(simple_break, 0, outline)
                plt.axvline(x=simple_break, linewidth=3, color='red', alpha=0.2)
                plt.scatter(simple_break, 0, c='red', s=23)
            
            if color != 'b':
                outline = ('Pre {}\nPost {}').format(curve_dict[index][5], max(curve_dict[index][1])-min(curve_dict[index][1]))
                plt.text(0, 25, outline)
                    
        plt.savefig(is_outfile_name + '_' + uid_ct + '_' + str(strand) +'.png')
        plt.show()
        plt.close()
            
def phred_at_random_simulation(obs_median, length, phred_list, correction_factor, phred_rates, simnum=1000):
    expected = 0
    if length > 10:
        for each in range(simnum):
            random_list = random.sample(phred_rates, length)
            if random_list:
                if np.median(random_list) <= obs_median:
                    expected+=1

        #Bonferroni correction: (correction_factor*expected)/float(simnum)
        return((correction_factor*expected)/float(simnum))
    else:
        return(1)
        
def phred_line_to_nt_list(phred_line, read_threshold, step_size=25):    
    phred_list = ord_phred(phred_line)
        
    window_median = []
    nt_list = []
    zones = []
    
    for nt in range(int(len(phred_list)-step_size)):
        temp_list = phred_list[nt:nt+step_size]
        win_med = read_threshold
        
        if temp_list:
            try:
                win_med = np.median(temp_list)
                if win_med < read_threshold:
                    zones.append(nt)
                    
            except Exception as ex:
                print(ex)
            
        window_median.append(win_med)
        nt_list.append(nt)
        
    return(nt_list, window_median, zones)
    
def purity_check(region, zones):
    processing = False
    temp_ct = 0
    region_length = 0
        
    for each in range(min(region), max(region)+1):
        region_length += 1
        if each in zones:
            temp_ct += 1
            
    if (region_length > 10):
        if temp_ct/float(region_length) > 0.9:
            processing = True
            
    return(processing)
    
def low_zones(phred_line, read_threshold, read_upper_threshold, phred_rates):
    _nt_list, _window_median, zones = phred_line_to_nt_list(phred_line, read_threshold)
    phred_list = ord_phred(phred_line)
    
    cold_zones = {}
    region = set()
    region_ct = 0
    stitch = {}
    stitch_pass = {}
    case = False
    resolved_zone_dict = {}

    #fine grain 10 nt region building
    for z_nt in zones:
        if z_nt not in region:
            
            # previous region is sufficient size to add to cold zones and next step is sufficiently far away
            if (len(region) >= 10) and (z_nt > (max(region)+10)):
                if purity_check(region, zones):
                    cold_zones[region_ct] = region
                    region_ct+=1
                region = set()

            # scan downstream for next hit, add to region
            for next_nt in range(z_nt, z_nt + 10):
                if next_nt in zones:
                    region.add(next_nt)
            
    if (len(region) >= 10):  
        if purity_check(region, zones):
            cold_zones[region_ct] = region
            
    if cold_zones:    
        for each_zone, regions in cold_zones.items():
            region_length = max((max(regions)-min(regions))/float(10), 50)
            group_set = set()
            group_set.add(each_zone)
                    
            if each_zone-1 in cold_zones:
                pre_scores = phred_list[max(cold_zones[each_zone-1]):min(cold_zones[each_zone])]
                            
                if len(pre_scores) > 20:
                    if np.median(pre_scores)+np.std(pre_scores) < read_upper_threshold:
                        group_set.add(each_zone-1)
                
                if (min(regions) - region_length <= max(cold_zones[each_zone-1])):
                    group_set.add(each_zone-1)
                    
            if each_zone+1 in cold_zones:
                pre_scores = phred_list[max(cold_zones[each_zone]):min(cold_zones[each_zone+1])]
                
                if len(pre_scores) > 20:
                    if np.median(pre_scores)+np.std(pre_scores) < read_upper_threshold:
                        group_set.add(each_zone+1)
                
                if (max(regions) + region_length >= min(cold_zones[each_zone+1])):
                    group_set.add(each_zone+1)
                    
            if len(group_set) >= 2:
                add_new = True
                for each_gs in group_set:
                    for stitch_num, stitch_group in stitch.items():
                        if each_gs in stitch_group:
                            for each_member in group_set:
                                stitch[stitch_num].add(each_member)
                                add_new = False
                      
                if add_new:
                    stitch[len(stitch)]=group_set
                    
            if len(group_set) == 1:
                stitch[len(stitch)]=group_set
    
    if stitch:
        resolved_zone_dict = {}
        case = False 
    
        already_processed = []
    
        for stitch_num, stitch_group in stitch.items():
            if stitch_num not in already_processed:
                zone_start = min(cold_zones[min(stitch_group)])
                zone_stop = max(cold_zones[max(stitch_group)])
                
                mod_len = max((zone_stop - zone_start)/float(10),500)
                mod_zone_start = zone_start - mod_len
                mod_zone_stop = zone_stop + mod_len
                
                if (zone_stop-zone_start) >= 50:
                    stitch_pass[stitch_num] = stitch_group
                    already_processed.append(stitch_num)
                                        
                    if (stitch_num+1 in stitch) and (stitch_num+1 not in already_processed):
                        post_stitch_group = stitch[stitch_num+1]
                        post_zone_start = min(cold_zones[min(post_stitch_group)])
                        post_zone_stop = max(cold_zones[max(post_stitch_group)])
                        post_mod_len = (post_zone_stop - post_zone_start)/float(2)
                        post_mod_zone_start = int(post_zone_start - post_mod_len)
                        post_mod_zone_stop = int(post_zone_stop + post_mod_len)
                        post_zone_phred = phred_list[post_mod_zone_start:post_mod_zone_stop]                            
                                               
                        if (len(phred_list[post_zone_stop:]) > 100) and post_zone_phred: 
                            exceed_list = [score for score in post_zone_phred if score >= read_stop_threshold]
                            if len(exceed_list) < (len(post_zone_phred)/10):
                                if (np.median(phred_list[post_zone_stop:]) >= read_mid_threshold) and (np.median(post_zone_phred) <= read_threshold):        
                                    if mod_zone_stop >= post_mod_zone_start:
                                        for each_group in post_stitch_group:
                                            stitch_group.add(each_group)
                                               
                                        stitch_pass[stitch_num] = stitch_group
                                        already_processed.append(stitch_num+1)
                    
                    if (stitch_num-1 in stitch) and (stitch_num-1 not in already_processed):
                        pre_stitch_group = stitch[stitch_num-1]
                        pre_zone_start = min(cold_zones[min(pre_stitch_group)])
                        pre_zone_stop = max(cold_zones[max(pre_stitch_group)])
                        pre_mod_len = (pre_zone_stop - pre_zone_start)/float(2)
                        pre_mod_zone_start = int(pre_zone_start - pre_mod_len)
                        pre_mod_zone_stop = int(pre_zone_stop + pre_mod_len)
                        
                        pre_zone_phred = phred_list[pre_mod_zone_start:pre_mod_zone_stop]
                                               
                        if (len(phred_list[pre_zone_stop:]) > 100) and pre_zone_phred:
                            exceed_list = [score for score in pre_zone_phred if score >= read_stop_threshold]
                            if len(exceed_list) < (len(pre_zone_phred)/10):
                                if (np.median(phred_list[:pre_zone_start]) >= read_mid_threshold) and (np.median(pre_zone_phred) <= read_threshold):        
                                    if mod_zone_start <= pre_mod_zone_stop:                                        
                                        for each_group in pre_stitch_group:
                                            stitch_group.add(each_group)
                                    
                                        stitch_pass[stitch_num] = stitch_group
                                        already_processed.append(stitch_num-1)
                                                                                                            
    if stitch_pass:           
        case = False
        for stitch_num, stitch_group in stitch_pass.items():
            zone_start = min(cold_zones[min(stitch_group)])
            zone_stop = max(cold_zones[max(stitch_group)])
            
            if len(phred_list[zone_start:zone_stop]) > 50:
                zone_median = np.median(phred_list[zone_start:zone_stop])
                
                pval = phred_at_random_simulation(zone_median, (zone_stop-zone_start), phred_list, len(stitch_pass), phred_rates, 1000)
                if pval <= 0.05:
                    case = True
                    resolved_zone_dict[len(resolved_zone_dict)]=[zone_median, abs(zone_start-zone_stop)]
     
    return(case, resolved_zone_dict)

def drop_anchor(start, stop, anchor_seed, anchor_set, anchor_point_dict, anchor_weight):
    
    if (start,stop) not in anchor_seed:
                
        pA = len(anchor_point_dict)
        anchor_point_dict[pA] = [start, start]
        
        pB = len(anchor_point_dict)
        anchor_point_dict[pB] = [stop, stop]
        
        anchor_set.add((pA, pB))
        
        anchor_seed[(start,stop)] = (pA, pB)
        
        if (pA, pB) not in anchor_weight:
            anchor_weight[(pA, pB)] = 1
        else:
            anchor_weight[(pA, pB)] += 1
                            
    else:
        (pA, pB) = anchor_seed[(start,stop)]
        
        if (pA, pB) not in anchor_weight:
            anchor_weight[(pA, pB)] = 1
        else:
            anchor_weight[(pA, pB)] += 1
                
    return(anchor_seed, anchor_set, anchor_point_dict, anchor_weight)
    
def weigh_anchor(starpoints_list, anchor_point_dict, anchor_weight, thumb=2):
    for star in starpoints_list:
        for (pA, pB) in anchor_weight:
            rA = anchor_point_dict[pA]
            rB = anchor_point_dict[pB]
        
            if star in range(min(rA), max(rA)):
                anchor_weight[(pA, pB)] += thumb
                
            if star in range(min(rB), max(rB)):
                anchor_weight[(pA, pB)] += thumb

    return(anchor_weight)
        
def expand_anchor(anchor_seed, anchor_set, anchor_point_dict):
    # similar to reduce but adds range
    for (start, stop), (pC, pD) in anchor_seed.items():
        for pA, pB in anchor_set:
            if pA != pC:
                rA = anchor_point_dict[pA]
                rB = anchor_point_dict[pB]
                
                if start in range(min(rA)-50, max(rA)+50):
                    if stop in range(min(rB)-50, max(rB)+50):
                        #first update points, if need be
                        if start not in range(min(rA),max(rA)):
                            rA = [min(min(rA),start), max(max(rA),start)]
                            anchor_point_dict[pA] = rA
                            
                        if stop not in range(min(rB),max(rB)):
                            rB = [min(min(rB),stop), max(max(rB),stop)]
                            anchor_point_dict[pB] = rB
                            
    return(anchor_set, anchor_point_dict)
            
def reduce_anchor(anchor_set, anchor_point_dict, anchor_weight, step=50):
    for pA, pB in anchor_set:
        for pC, pD in anchor_set:
            if pA != pC:
                
                #if pA in anchor_point_dict:
                rA = anchor_point_dict[pA]
                rB = anchor_point_dict[pB]
                rC = anchor_point_dict[pC]
                rD = anchor_point_dict[pD]
                        
                if min(rA) in range(min(rC)-step, max(rC)+step) or max(rA) in range(min(rC)-step, max(rC)+step):
                    if min(rB) in range(min(rD)-step, max(rD)+step) or max(rB) in range(min(rD)-step, max(rD)+step):

                        #first update points, if need be
                        rC = [min(rC[0],rA[0]), max(rC[1],rA[1])]
                        anchor_point_dict[pC] = rC
                            
                        rD = [min(rD[0],rB[0]), max(rD[1],rB[1])]
                        anchor_point_dict[pD] = rD
                        
                        anchor_weight[(pC,pD)] += anchor_weight.pop((pA,pB))
                        
                        anchor_set.remove((pA, pB))
                        
                        return(True, anchor_set, anchor_point_dict, anchor_weight)
    else:            
        return(False, anchor_set, anchor_point_dict, anchor_weight)

def build_weights_in_range(weight, each_range, weight_nt_list, nt_weight):    
    if min(each_range) == max(each_range):
        each_nt = min(each_range)
        weight_nt_list.append(weight)
            
        if each_nt in nt_weight:
            nt_weight[each_nt] += weight
        else:
            nt_weight[each_nt] = weight
        
        return(weight_nt_list, nt_weight)
        
    else:
        for each_nt in range(min(each_range),max(each_range)):
            weight_nt_list.append(weight)
            
            if each_nt in nt_weight:
                nt_weight[each_nt] += weight
            else:
                nt_weight[each_nt] = weight
                
        return(weight_nt_list, nt_weight)
               
def derive_nt_weights(each_range, average, nt_weight):  
    contested = True
    max_votes = ['', 0]
    weight_list = []
    
    if min(each_range) == max(each_range):
        each_nt = min(each_range)
        if each_nt in nt_weight:
            weight = nt_weight[each_nt]
            
        return(each_nt, weight, 0)
        
    else:
        for each_nt in range(min(each_range),max(each_range)):
            if each_nt in nt_weight:
                weight = nt_weight[each_nt]
                weight_list.append(weight)
                
                if weight > max_votes[1]:
                    max_votes = [each_nt, weight]       
        
        while contested:
            contested = False
            for each_nt in range(min(each_range),max(each_range)):
                if each_nt in nt_weight:
                    weight = nt_weight[each_nt]
                    
                    if weight == max_votes[1] and each_nt != max_votes[0]:
                        if min((abs(each_nt)-average),(abs(max_votes[0])-average)) == (abs(each_nt)-average):
                            max_votes = [each_nt, weight]
                            contested = True
        
        if weight_list:
            if np.std(weight_list) == 0:
                std_weight = round((max_votes[1]-np.median(weight_list)),2)
            else:
                std_weight = round((max_votes[1]-np.median(weight_list))/np.std(weight_list),2)
        
            return(max_votes[0], max_votes[1], std_weight)
        else:
            return(max_votes[0], max_votes[1], average)

def find_anchorpoint(starpoints_set, anchor_set, anchor_point_dict, anchor_weight):
    anchorpoints = {}
    weight_nt_list = []
    nt_weight = {}
    
    for (pA, pB) in anchor_set:
        weight = anchor_weight[(pA, pB)]
        rA = anchor_point_dict[pA]
        rB = anchor_point_dict[pB]
        
        weight_nt_list, nt_weight = build_weights_in_range(weight, rA, weight_nt_list, nt_weight)
        weight_nt_list, nt_weight = build_weights_in_range(weight, rB, weight_nt_list, nt_weight)
        
        star = False
        for a_nt in range(min(rA),max(rA)):
            if a_nt in starpoints_set:
                star = True

        if star:
            star = False
            for b_nt in range(min(rB),max(rB)):
                if b_nt in starpoints_set:
                    star = True
                
        if star:
            anchorpoints[(pA, pB)] = [0,0]
            
    mode_nt = max(set(weight_nt_list), key=weight_nt_list.count)
            
    for (pA, pB) in anchorpoints:
        rA = anchor_point_dict[pA]
        rB = anchor_point_dict[pB]
        
        anchor_A, aA_median, aA_zscore = derive_nt_weights(rA, mode_nt, nt_weight)
        anchor_B, aB_median, aB_zscore = derive_nt_weights(rB, mode_nt, nt_weight)
        anchorpoints[(pA, pB)] = [anchor_A, aA_median, aA_zscore, anchor_B, aB_median, aB_zscore]
        
    return(anchorpoints)
            
def peak_size_test(phred_list, score_threshold, length_threshold, purity):
    temp_list = []

    for step in range(len(phred_list)):
        above_threshold = 0
        temp_list.append(phred_list[step])
        
        if len(temp_list) > length_threshold:
            temp_list = temp_list[1:]
        
        for isval in temp_list:
            if isval >= score_threshold:
                above_threshold += 1
        
        if len(temp_list) >= length_threshold:
            if above_threshold/float(len(temp_list)) >= purity:
                return(False)
                
    return(True)   
    
def scan_flanking(stitch_pass, already_processed, stitch, stitch_num, stitch_group, cold_zones, zone_stop, zone_start, mod_zone_start, mod_zone_stop, phred_list, read_mid_threshold, read_stop_threshold):
    if (zone_stop-zone_start) >= 50:
        if np.median(phred_list[zone_start:zone_stop]) <= read_mid_threshold :
            # length test needs the low_phred region to be bounded on 
            # both sides. Here we make that test by requiring the boundaries to be 
            # at least 50nt long and exceed the read_threshold
            head_median = 0
            tail_median = 0
            
            if len(phred_list[:zone_start]) > 50:
                head_median = np.median(phred_list[:zone_start])
                
            if len(phred_list[zone_stop:]) > 50:    
                tail_median = np.median(phred_list[zone_stop:])
                
            if (head_median >= read_mid_threshold) and (tail_median >= read_mid_threshold):
                stitch_pass[stitch_num] = stitch_group
                
                # what about flanking regions? We also want to include these if they are 
                # nearby, low_phred, and still abides by the boundary condition
                if (stitch_num+1 in stitch) and (stitch_num+1 not in already_processed):
                    post_stitch_group = stitch[stitch_num+1]
                    post_zone_start = zone_stop
                    post_zone_stop = max(cold_zones[max(post_stitch_group)])
                    post_mod_len = (post_zone_stop - post_zone_start)/float(2)
                    post_mod_zone_start = int(post_zone_start - post_mod_len)
                    post_mod_zone_stop = int(post_zone_stop + post_mod_len)
                    post_zone_phred = phred_list[post_mod_zone_start:post_mod_zone_stop]                            
                                           
                    if (len(phred_list[post_zone_stop:]) > 100):
                        pass_excess = True
                        exceed_list = [score for score in post_zone_phred if score >= read_stop_threshold]
                        
                        if len(exceed_list) > (len(post_zone_phred)/10):
                            pass_excess = False
                               
                            if peak_size_test(post_zone_phred, read_stop_threshold, 25, 0.9):
                                if peak_size_test(post_zone_phred, read_threshold, 100, 0.7):
                                    pass_excess = True
                                
                        if pass_excess:  
                            if (np.median(phred_list[post_zone_stop:]) >= read_mid_threshold) and (np.median(post_zone_phred) <= read_threshold):        
                                if mod_zone_stop >= post_mod_zone_start:
                                    for each_group in post_stitch_group:
                                        stitch_group.add(each_group)
    
                                    stitch_pass[stitch_num] = stitch_group
                                    already_processed.append(stitch_num+1)
            
                if (stitch_num-1 in stitch) and (stitch_num-1 not in already_processed):
                    pre_stitch_group = stitch[stitch_num-1]
                    pre_zone_start = min(cold_zones[min(pre_stitch_group)])
                    pre_zone_stop = zone_start
                    pre_mod_len = (pre_zone_stop - pre_zone_start)/float(2)
                    pre_mod_zone_start = int(pre_zone_start - pre_mod_len)
                    pre_mod_zone_stop = int(pre_zone_stop + pre_mod_len)
                    
                    pre_zone_phred = phred_list[pre_mod_zone_start:pre_mod_zone_stop]
                                           
                    if (len(phred_list[pre_zone_stop:]) > 100):
                        pass_excess = True
                        exceed_list = [score for score in pre_zone_phred if score >= read_stop_threshold]
                        
                        if len(exceed_list) > (len(pre_zone_phred)/10):
                            pass_excess = False
                             
                            if peak_size_test(pre_zone_phred, read_stop_threshold, 25, 0.9):
                                if peak_size_test(pre_zone_phred, read_threshold, 100, 0.7):
                                    pass_excess = True
                    
                        if pass_excess:
                            if (np.median(phred_list[:pre_zone_stop]) >= read_mid_threshold) and (np.median(pre_zone_phred) <= read_threshold):        
                                if mod_zone_start <= pre_mod_zone_stop:
                                    for each_group in pre_stitch_group:
                                        print(each_group)
                                        stitch_group.add(each_group)
                                
                                    stitch_pass[stitch_num] = stitch_group
                                    already_processed.append(stitch_num-1)
                                    
        
    return(stitch_pass, already_processed)       
    
def resolve_zones(phred_line, read_threshold, read_mid_threshold, read_upper_threshold, brk_start, brk_stop, read_stop, phred_rates, runmode='resolve'):
    _nt_list, _window_median, zones = phred_line_to_nt_list(phred_line, read_threshold)
    phred_list = ord_phred(phred_line)
    
    cold_zones = {}
    region = set()
    region_ct = 0
    stitch = {}
    stitch_pass = {}
    case = False
    resolved_zone_dict = {}

    #fine grain 10 nt region building
    for z_nt in zones:
        if z_nt not in region:
            # previous region is sufficient size to add to cold zones and next step is sufficiently far away
            if (len(region) >= 10) and (z_nt > (max(region)+10)):
                if purity_check(region, zones):
                    cold_zones[region_ct] = region
                    region_ct+=1
                region = set()

            # scan downstream for next hit, add to region
            for next_nt in range(z_nt, z_nt + 10):
                if next_nt in zones:
                    region.add(next_nt)
            
    if (len(region) >= 10):  
        if purity_check(region, zones):
            cold_zones[region_ct] = region
            
    if cold_zones:    
        for each_zone, regions in cold_zones.items():
            region_length = max((max(regions)-min(regions))/float(10), 50)
            group_set = set()
            group_set.add(each_zone)
                                
            if each_zone-1 in cold_zones:
                pre_scores = phred_list[max(cold_zones[each_zone-1]):min(cold_zones[each_zone])]
                            
                if len(pre_scores) > 20:
                    if np.median(pre_scores)+np.std(pre_scores) < read_upper_threshold:
                        group_set.add(each_zone-1)
                        
                    else:
                        if peak_size_test(pre_scores, read_stop_threshold, 10, 0.9):
                            group_set.add(each_zone-1)
                
                if (min(regions) - region_length <= max(cold_zones[each_zone-1])):
                    group_set.add(each_zone-1)
                    
            if each_zone+1 in cold_zones:
                pre_scores = phred_list[max(cold_zones[each_zone]):min(cold_zones[each_zone+1])]
                
                if len(pre_scores) > 20:
                    if np.median(pre_scores)+np.std(pre_scores) < read_upper_threshold:
                        group_set.add(each_zone+1)
                    else:
                        if peak_size_test(pre_scores, read_stop_threshold, 10, 0.9):
                            group_set.add(each_zone+1)
                    
                if (max(regions) + region_length >= min(cold_zones[each_zone+1])):
                    group_set.add(each_zone+1)
                    
            if len(group_set) >= 2:
                add_new = True
                for each_gs in group_set:
                    for stitch_num, stitch_group in stitch.items():
                        if each_gs in stitch_group:
                            for each_member in group_set:
                                stitch[stitch_num].add(each_member)
                                add_new = False
                      
                if add_new:
                    stitch[len(stitch)]=group_set
                    
            if len(group_set) == 1:
                stitch[len(stitch)]=group_set
        
    if stitch:
        resolved_zone_dict = {}
        case = False 
    
        already_processed = []
    
        for stitch_num, stitch_group in stitch.items():
            if stitch_num not in already_processed:
                zone_start = min(cold_zones[min(stitch_group)])
                zone_stop = max(cold_zones[max(stitch_group)])
                
                mod_len = max((zone_stop - zone_start)/float(10),500)
                mod_zone_start = zone_start - mod_len
                mod_zone_stop = zone_stop + mod_len
                
                if runmode == 'resolve':
                    # Roadmap
                    #   start===brstart===stop-----brstop                               brstart----start===brststop===stop                            brstart----start===stop---brstop                               start===brstart===brstop===stop
                    if (mod_zone_start <= brk_start and mod_zone_stop >= brk_start) or (mod_zone_start <= brk_stop and mod_zone_stop >= brk_stop) or (mod_zone_start >= brk_start and mod_zone_stop <= brk_stop) or (mod_zone_start <= brk_start and mod_zone_stop >= brk_stop):
                        stitch_pass, already_processed = scan_flanking(stitch_pass, already_processed, stitch, stitch_num, stitch_group, cold_zones, zone_stop, zone_start, mod_zone_start, mod_zone_stop, phred_list, read_mid_threshold, read_stop_threshold)
                else:
                    stitch_pass, already_processed = scan_flanking(stitch_pass, already_processed, stitch, stitch_num, stitch_group, cold_zones, zone_stop, zone_start, mod_zone_start, mod_zone_stop, phred_list, read_mid_threshold, read_stop_threshold)

    if stitch_pass:           
        case = False
        for stitch_num, stitch_group in stitch_pass.items():
            zone_start = min(cold_zones[min(stitch_group)])
            zone_stop = max(cold_zones[max(stitch_group)])
            
            if len(phred_list[zone_start:zone_stop]) > 0:
                zone_median = np.median(phred_list[zone_start:zone_stop])
                
                pval = phred_at_random_simulation(zone_median, (zone_stop-zone_start), phred_list, len(stitch_pass), phred_rates, 1000)
                if pval <= 0.05:
                    case = True
                
                resolved_zone_dict[len(resolved_zone_dict)]=[zone_median, zone_start, zone_stop, pval, phred_list[zone_start:zone_stop]]
     
    return(case, resolved_zone_dict)
    
def short_run(low_list):
    low_list.sort()
    
    run = 0
    step = 0
    while run < 1000:
        step+=1
        if step+1 >= len(low_list):
            return(False)
        else:
            if low_list[step+1] == low_list[step]:
                run+=1
                
    if run >= 100:
        return(True)
        
def length_function(brkpt, threshold_deets, first_fastq_dict, supplemental_align_set, read_ct, read_dict_prefilter, is_outfile_name):
    starting_distances = []
    zone_lengths = []
    cor_lengths = []
    uid_brkpt_sites = {}
    cross_brkpt_set = set()
    
    brk_chromo, nts = brkpt.split(':')
    brk_start = int(nts.split('-')[0])
    brk_stop = int(nts.split('-')[1])
    
    outstats_name = ('{}_stats.log').format(is_outfile_name)
    outstats = open(outstats_name,'w')
    
    outline = ('For sample: {}\nFor breakpoint {}\n').format(is_outfile_name, brkpt)
    outstats.write(outline)
    
    phred_rates, read_threshold, read_mid_threshold, read_upper_threshold, read_stop_threshold = threshold_deets
                                
    outline = ('Median global phred score: {}\nLower read threshold: {}\nMiddle read threshold: {}\nUpper read threshold: {}\n').format(read_upper_threshold, read_threshold, read_mid_threshold, read_stop_threshold)
    outstats.write(outline)
    
    read_dict = {}
    for uid, max_ct in read_ct.items():
        temp_chromo = {}
        for ct in range(max_ct+1):
            uid_ct = ('{}.{}').format(uid, ct)
            strand = read_dict_prefilter[uid_ct][0]
            chromo = read_dict_prefilter[uid_ct][5]

            if chromo == brk_chromo:
                if chromo in temp_chromo:
                    temp_chromo[chromo].add(strand)
                else:
                    strand_set = set()
                    strand_set.add(strand)
                    temp_chromo[chromo] = strand_set
        
        process = False
        if brk_chromo in temp_chromo:
            strand_set = temp_chromo[brk_chromo]
            if len(strand_set) > 1:
                process = True
                
        if process:
            for ct in range(max_ct+1):
                uid_ct = ('{}.{}').format(uid, ct)
                read_dict[uid_ct] = read_dict_prefilter[uid_ct]
            
    for uid_ct, deets in read_dict.items():
        uid = uid_ct.split('.')[0]
        if uid in supplemental_align_set:
            outline = ('For uid_ct: {}\n').format(uid_ct)
            outstats.write(outline)
            strand, start, stop, local_phred_line, cigar, chromo = deets
                        
            if chromo == brk_chromo:
                if (start <= brk_start and stop >= brk_start) or (start <= brk_stop and stop >= brk_stop) or (start >= brk_start and stop <= brk_stop):
                    cross_brkpt_set.add(uid)
                    outline = ('uid added to cross_brkpt_set: {}\n').format(uid_ct)
                    outstats.write(outline)
                    
                    uid_brkpt_dict[uid_ct] = parallel_construction(strand, cigar, len(local_phred_line), start, stop, brk_start, brk_stop)
                    
                    if uid not in uid_brkpt_sites:
                        uid_brkpt_sites[uid] = [uid_brkpt_dict[uid_ct]]
                    else:
                        uid_brkpt_sites[uid] += [uid_brkpt_dict[uid_ct]]     
    samfile.close()

    outline = ('{} aligned regions. {} alignments overlap breakpoint. {} supplemental alignments to multiple loci.\n').format(len(read_dict), len(cross_brkpt_set), len(supplemental_align_set)) 
    print(outline)
    outstats.write(outline)
    
    fastq_dict = {}

    for uid, deets in first_fastq_dict.items():            
        if (uid in read_ct) and (uid in supplemental_align_set) and (uid in cross_brkpt_set):
            fastq_dict[uid] = deets
            
    outline = ('{} reads with supplemental alignments with at least one overlapping the breakpoint.\n').format(len(fastq_dict)) 
    print(outline)
    outstats.write(outline)
         
    for uid, phred_line in fastq_dict.items():
        outline = ('Processing uid {}...').format(uid) 
        print(outline)
        outstats.write(outline)
        
        if (uid in read_ct) and (uid in supplemental_align_set) and (uid in cross_brkpt_set):
            map_dict = {}
            phred_line = fastq_dict[uid]
            read_stop = len(phred_line)
                        
            uid_ct_to_section = {}
            uid_ct_to_coordinates = {}
            
            for ct in range(read_ct[uid]+1):
                uid_ct = ('{}.{}').format(uid, ct)
                outline = ('Processing read section {}\n').format(uid_ct)
                print(outline)
                outstats.write(outline)
                
                if uid_ct in read_dict:
                    section_phred = read_dict[uid_ct][3]
                    
                    if len(section_phred) == len(phred_line):
                        uid_ct_to_section[uid_ct] = len(phred_line)
                        uid_ct_to_coordinates[uid_ct] = [len(phred_line.split(section_phred)[0]), len(phred_line.split(section_phred)[0])+len(section_phred)]
                    
                    else:
                        if section_phred in phred_line:
                            uid_ct_to_section[uid_ct] = len(phred_line.split(section_phred)[0])
                            uid_ct_to_coordinates[uid_ct] = [len(phred_line.split(section_phred)[0]), len(phred_line.split(section_phred)[0])+len(section_phred)]

                        else:
                            if section_phred[::-1] in phred_line:
                                uid_ct_to_section[uid_ct] = len(phred_line.split(section_phred[::-1])[0])
                                uid_ct_to_coordinates[uid_ct] = [len(phred_line.split(section_phred[::-1])[0]), len(phred_line.split(section_phred)[0])+len(section_phred)]

                    for each_nt in range(uid_ct_to_coordinates[uid_ct][0], uid_ct_to_coordinates[uid_ct][1]+1):
                        map_dict[each_nt]=uid_ct
                
            sorted_uid_ct_to_section = sorted(uid_ct_to_section.items(), key=lambda kv: kv[1])
            
            start_sites = []
            stop_sites = []
            for sites in uid_brkpt_sites[uid]:
                start_sites.append(sites[0])
                stop_sites.append(sites[1])
                
            curve_dict = {}
                
            nt_list, window_median, _Z = phred_line_to_nt_list(phred_line, read_threshold)
            
            if start_sites:
                median_starts = np.median(start_sites)
            else:
                median_starts = window_median
                
            
            if stop_sites:
                median_stops = np.median(stop_sites)
            else:
                median_stops = window_median

            curve_dict[0]=['b', nt_list, window_median, median_starts, median_stops]
            
            for index in range(int(math.floor(len(sorted_uid_ct_to_section)/(2)))):                
                uid_ct_f = sorted_uid_ct_to_section[2*index][0]
                uid_ct_r = sorted_uid_ct_to_section[(2*index)+1][0]
                outline = ('Sorted_uid_ct_to_section {}\t{}\n').format(uid_ct_f,uid_ct_r)
                print(outline)
                outstats.write(outline)
                #assign to the reverse
                for each_nt in range(len(phred_line[uid_ct_to_coordinates[uid_ct_f][1]:uid_ct_to_coordinates[uid_ct_r][0]])):
                    if (each_nt + uid_ct_to_coordinates[uid_ct_f][1]) not in map_dict:
                        map_dict[(each_nt + uid_ct_to_coordinates[uid_ct_f][1])]=uid_ct_r
                        
            mapped_uid_ct = {}
            for map_nt, uid_ct in map_dict.items():
                if uid_ct not in mapped_uid_ct:
                    map_set = set()
                    map_set.add(map_nt)
                    mapped_uid_ct[uid_ct]=map_set

                else:
                    mapped_uid_ct[uid_ct].add(map_nt)
            
            outline = ('Reads with significant low zone:\n')
            outstats.write(outline)
            
            for s_uid in sorted_uid_ct_to_section:
                uid_ct = s_uid[0]
                strand = read_dict[uid_ct][0]
                
                outline = ('{}\t{}\n').format(s_uid,uid_ct)
                print(outline)
                outstats.write(outline)
                
                if (uid_ct in uid_brkpt_dict):
                    rel_brk_start, rel_brk_stop = uid_brkpt_dict[uid_ct]
                    print('Running resolved zones')
                    case = False
                    outline = ('In Resolved Zones: {}\t{}\t{}\n').format(s_uid, rel_brk_start, rel_brk_stop)
                    print(outline)
                    outstats.write(outline)
                
                    case, resolved_zone_dict = resolve_zones(phred_line, read_threshold, read_mid_threshold, read_upper_threshold, rel_brk_start, rel_brk_stop, read_stop, phred_rates)

                    if len(resolved_zone_dict) > 0:
                        for each_zone, res_stats in resolved_zone_dict.items():
                            pval = res_stats[3]
                            
                            if pval <= 0.05 and case:
                                zone_type = 'r'
                                outline = ('{}\n').format(uid)
                                outstats.write(outline)
                                
                            else:
                                zone_type = 'g'
                                                                
                            low_zone = res_stats[4]
                            lz_start = res_stats[1]
                            
                            x_is = []
                            for each_lz in range(len(low_zone)):
                                x_is.append(each_lz+lz_start)
                            
                            curve_dict[len(curve_dict)]=[zone_type, x_is, low_zone, rel_brk_start, rel_brk_stop, lz_start]

                            zone_lengths.append(len(low_zone))
                            starting_dist = min(lz_start-rel_brk_start, lz_start-rel_brk_stop)
                            starting_distances.append(starting_dist)

                            outline = ('Distance from breakpoint: {}\tLow zone length: {}\n').format(starting_dist, len(low_zone))
                            outstats.write(outline)

                        graphical_curve(curve_dict, uid, strand, brkpt, is_outfile_name, False)
                                            
                    if case:
                        for each_zone, res_stats in resolved_zone_dict.items(): 
                            zone_median, zone_start, zone_stop, pval, p_list = res_stats
                            lengths = [zone_start, (zone_stop-zone_start)]
                            if lengths not in cor_lengths:
                                cor_lengths.append(lengths)
    
    rho, pval = stats.spearmanr(cor_lengths)
            
    outline = ('Sample: {} Rho: {} pval: {} from {} sequences.\n').format(is_outfile_name, round(float(rho),4), round(float(pval),6), len(cor_lengths))
    print(outline)
    outstats.write(outline)
            
    if starting_distances:
        outline = ('Distances from breakpoint:\nMedian: {}\tLeast: {}\tGreatest: {}\n').format(np.median(starting_distances), min(starting_distances), max(starting_distances))
        outstats.write(outline)
    
    if zone_lengths:
        outline = ('Low Phred Score zone lengths:\nMedian: {}\tLeast: {}\tGreatest: {}\n').format(np.median(zone_lengths), min(zone_lengths), max(zone_lengths))
        outstats.write(outline)
    
    outstats.close()
    
    resource_pickle_name = ('{}_lengths.p').format(is_outfile_name)         
    with open(resource_pickle_name, 'wb') as file:
        pickle.dump(cor_lengths, file)
           
def test_for_close_inversion(uid, read_ct, read_dict, uid_starpower):
    
    for uid_is_1 in range(read_ct[uid]+1):
        for uid_is_2 in range(read_ct[uid]+1):
            if uid_is_1 != uid_is_2:
                uid_ct_1 = ('{}.{}').format(uid, uid_is_1)
                uid_ct_2 = ('{}.{}').format(uid, uid_is_2)

                deets_1 = read_dict[uid_ct_1]
                deets_2 = read_dict[uid_ct_2]

                #Does it have a significant inversion? It should be mapped in two directions.
                if deets_1[0] != deets_2[0] and deets_1[1] == deets_2[1]:
                    starpower = max(uid_starpower[uid_ct_1], uid_starpower[uid_ct_2])

                    #Test if there is are two aligned segments that are within the starpower of each other. 
                    if (abs(deets_1[2]-deets_2[2]) <= starpower) or (abs(deets_1[3]-deets_2[3]) <= starpower):
                        return(True)
    return(False)
    
def test_for_long_phred(phred_list):
    low_list = [score for score in phred_list if score < read_threshold]
    if len(low_list) >= 1000:
        if short_run(low_list):
            return(True)
    
    return(False)
    
def filter_bprd(start, stop, coor_list):
    for fstart, fstop in coor_list: 
        for ent in [start, stop]:
            if (ent >= fstart) and (ent <= fstop):
                return(True)
    else:
        return(False)
        
if args.run_demo:
    run_demo()
            
if args.plot_phred:
    os_mkdir(args.outfile_name)
    fastqfile_name = (args.fastqfile_name)
    phred_file = open(fastqfile_name)
    
    line_ct = 4
    for line in phred_file:
        line_ct += 1
        
        if line[0] == '@':
            uid = line.split('@')[1].split(' ')[0].strip()
            line_ct = 0
                                        
        if line_ct == 3:
            phred_line = line.strip()
    
    phred_file.close()
    
    if args.nt_min and args.nt_max:
        phred_list = ord_phred(phred_line[int(args.nt_min):int(args.nt_max)])

    if args.nt_min and not args.nt_max:
        phred_list = ord_phred(phred_line[int(args.nt_min):])    

    if not args.nt_min and args.nt_max:
        phred_list = ord_phred(phred_line[:int(args.nt_max)])
        
    if not args.nt_min and not args.nt_max:
        phred_list = ord_phred(phred_line)
    
    x = []
    y = []
    
    index=0
    for score in phred_list:
        x.append(index)
        y.append(score)
        index+=1
        
    plt.plot(x, y)
    
    nx, ny = rolling_median(x, y, 'plot_phred')
    plt.plot(nx, ny)
    
    plt.savefig(args.outfile_name)
    plt.close()

if args.breakpoint_retrieval_and_definition: 
    os_mkdir(args.outfile_name)
        
    fastqfile_name = args.fastqfile_name
    samfile_name = args.samfile_name
    outfile_name = args.outfile_name
    
    print('Processing fastq file to build complete fastq phred set ...')
    infile = open(fastqfile_name)
    
    filter_dict = {}

    if args.filter_list_file:     
        filter_list = open(args.filter_list_file)
        
        for line in filter_list:
            if line[0] != '#':
                line = line.strip()
                chromo = line.split('\t')[0]
                start = int(line.split('\t')[1])
                stop = int(line.split('\t')[2])
                
                if chromo not in filter_dict:
                    filter_dict[chromo] = [[start,stop]]
                else:
                    filter_dict[chromo] += [start,stop]

    phred_dict = {}
    anchor_point_dict = {}
    anchor_set = set()
    
    all_phred_list = []
    line_ct = 4
    total_read_ct = 0

    for line in infile:
        line_ct += 1
        
        if line[0] == '@':
            uid = line.split('@')[1].split(' ')[0].strip()
            line_ct = 0
            total_read_ct += 1
                                        
        if line_ct == 3:
            phred_line = line.strip()
            if phred_line:
                phred_list = ord_phred(phred_line)
                phred_dict[uid] = phred_list

                if phred_list:
                    all_phred_list.append(np.median(phred_list))
            
    infile.close()
    
    read_threshold = np.median(all_phred_list)-(5*np.std(all_phred_list))
    
    print('Processing sam file for breakpoint retrieval and definition (bprd)...')
    samfile = open(samfile_name)
    multimapper_align_set = set()
    supplemental_align_set = set()
    read_ct = {}
    read_dict = {}
    chromo_size_dict = {}
    starpoints_by_chromo = {}
    uid_starpower = {}

    for line in samfile:
        if line[0] == '@':
            if 'LN:' in line:
                line = line.strip()
                chromo = line.split('SN:')[1].split('\t')[0]
                size = int(line.split('LN:')[1])
                            
                if chromo not in chromo_size_dict:
                    chromo_size_dict[chromo] = size
                else:
                    print('chromosome name duplication in sam file')
                
        if line[0] != '@':
            multimapper = unpackbits(np.array([int(line.split('\t')[1])]))[0][8]
            if multimapper == 0:    
                if unpackbits(np.array([int(line.split('\t')[1])]))[0][11] == 1:
                    supplemental_align_set.add(uid)
                    
                uid = line.split('\t')[0]
                chromo = line.split('\t')[2]
                
                start = int(line.split('\t')[3])
                strand = unpackbits(np.array([int(line.split('\t')[1])]))[0][4]
                cigar = line.split('\t')[5]
                
                if cigar != '*':
                    stop = start + distance_in(strand, cigar, 'reference') - 1                        
                    filter_line = False
                    
                    if filter_dict:
                        if chromo in filter_dict:
                            coor_list = filter_dict[chromo]
                            filter_line = filter_bprd(start, stop, coor_list)
                    
                    if not filter_line:
                        starred, starpower, starpoints_by_chromo = parse_starpoints(start, stop, chromo, cigar, starpoints_by_chromo)
                                                    
                        if starred:
                            if uid not in read_ct:
                                read_ct[uid] = 0
                            else:
                                read_ct[uid] += 1
                                
                            uid_ct = ('{}.{}').format(uid, read_ct[uid])
                            read_dict[uid_ct] = [strand, chromo, start, stop]
                            uid_starpower[uid_ct] = starpower
                                        
    samfile.close()
        
    if not chromo_size_dict:
        print('Failed to load chromosome sizes from sam file. Please ensure sam file contains chromosome names.')
        chromo_size_dict = {'NC_001133.9':230218,'NC_001134.8':813184,'NC_001135.5':316620,'NC_001136.10':1531933,'NC_001137.3':576874,'NC_001138.5':270161,'NC_001139.9':1090940,'NC_001140.6':562643,'NC_001141.2':439888,'NC_001142.9':745751,'NC_001143.9':670191,'NC_001144.5':1078177,'NC_001145.3':924431,'NC_001146.8':784333,'NC_001147.6':1091291,'NC_001148.4':948066,'NC_001224.1':85779}
        
    same_chromo_list = []
    for uid, ct in read_ct.items():
        uid_base = ('{}.0').format(uid)
        base_chromo = read_dict[uid_base][1]
        
        for index in range(1, ct+1):
            uid_base = ('{}.{}').format(uid, index)
            other_chromo = read_dict[uid_base][1]
            if base_chromo == other_chromo:
                same_chromo_list.append(uid)
    
    hit_dict = {}
    pair_dict = {}
    cleared_uid_list = []
    
    breakpoints_by_chromo = {}
    
    #Does it have supplemental reads? It should be in supplemental_align_set
    for uid in supplemental_align_set:
        process = False        
        if uid in same_chromo_list:                
            if test_for_close_inversion(uid, read_ct, read_dict, uid_starpower):
                #Does it have a low phred score region at least 1000 nt long?
                
                if test_for_long_phred(phred_dict[uid]):                                                
                    for index in range(read_ct[uid]+1):
                        uid_ct = ('{}.{}').format(uid, index)
                        chromo = (read_dict[uid_ct][1])
                        rstart = int(read_dict[uid_ct][2])
                        rstop = int(read_dict[uid_ct][3])
                        
                        if chromo in chromo_size_dict:                            
                            if (rstart >= 500) and (rstop < chromo_size_dict[chromo]-500):
                                
                                if chromo not in breakpoints_by_chromo:
                                    breakpoints_by_chromo[chromo] = [[rstart, rstop, uid_ct]]
                                else:
                                    breakpoints_by_chromo[chromo] += [[rstart, rstop, uid_ct]]
                                    
                                cleared_uid_list.append(uid)
       
                                    
    if breakpoints_by_chromo:
        phase_one_outfile_name = outfile_name + '_phaseone.bed'
        phase_one_outfile = open(phase_one_outfile_name, 'w')
        anchor_to_brkpt = {}
        
        for chromo, breakpoint_list in breakpoints_by_chromo.items():
            starpoints_set = starpoints_by_chromo[chromo]
            anchor_seed = {}
            anchor_set = set()
            anchor_point_dict = {}
            anchor_weight = {}
            
            for rstart, rstop, uid_ct in breakpoint_list:
                anchor_seed, anchor_set, anchor_point_dict, anchor_weight = drop_anchor(rstart, rstop, anchor_seed, anchor_set, anchor_point_dict, anchor_weight)
    
            anchor_weight = weigh_anchor(starpoints_set, anchor_point_dict, anchor_weight, 2)
    
            reduced = True
            while reduced:
                reduced, anchor_set, anchor_point_dict, anchor_weight = reduce_anchor(anchor_set, anchor_point_dict, anchor_weight, 50)
            
            anchor_set, anchor_point_dict = expand_anchor(anchor_seed, anchor_set, anchor_point_dict)
                
            reduced = True
            while reduced:
                reduced, anchor_set, anchor_point_dict, anchor_weight = reduce_anchor(anchor_set, anchor_point_dict, anchor_weight, 1000)
                
            anchorpoints = find_anchorpoint(starpoints_set, anchor_set, anchor_point_dict, anchor_weight)
                                
            anchor_index = 0
            for (pA, pB), (anchor_start, start_median, start_zscore, anchor_stop, stop_median, stop_zscore) in anchorpoints.items():
                anchor_uid = str(hash(('{},{}').format(chromo, anchor_index)))
                anchor_index += 1  
                
                if (start_median >= 5):
                    outline = ('{}\t{}\t{}\t{}.0\t{}\t.\n').format(chromo, anchor_start-10, anchor_start+10, anchor_uid, start_median)
                    anchor_to_brkpt[anchor_uid+'.0']=outline
                    phase_one_outfile.write(outline)
                    
                if (stop_median >= 5):
                    outline = ('{}\t{}\t{}\t{}.1\t{}\t.\n').format(chromo, anchor_stop-10, anchor_stop+10, anchor_uid, stop_median)
                    anchor_to_brkpt[anchor_uid+'.1']=outline
                    phase_one_outfile.write(outline)

        phase_one_outfile.close()
        
        if anchor_to_brkpt:
            print('Running Anchorpoints')
            bamfile_name = args.bamfile_name         
            
            temp_genome_file_name = outfile_name + '_temp_genome.tab'
            temp_genome_file = open(temp_genome_file_name, 'w')
            
            for chromo, stop in chromo_size_dict.items():
                outline = ('{}\t{}\n').format(chromo, stop)
                temp_genome_file.write(outline)
            temp_genome_file.close()
                        
            temp_flank_file_name = outfile_name + '_temp_flank.tab'
            bash_command = ('bedtools flank -b 1000 -i {} -g {} > {}').format(phase_one_outfile_name, temp_genome_file_name, temp_flank_file_name)
            output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
            if output:
                print(output)
            
            temp_count_file_name = outfile_name + '_temp_count.tab'
            bash_command = ('bedtools coverage -counts -a {} -b {} > {}').format(temp_flank_file_name, bamfile_name, temp_count_file_name)
            output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
            if output:
                print(output)
                
            temp_count_file = open(temp_count_file_name)
            
            eval_brkpts = {}
            passed_brkpts = {}
            quant_score = {}
            chromo_regions = {}
                        
            #TODO fixx file name
            outfile = open(outfile_name + '_bprd.bed', 'w')
            
            for line in temp_count_file:
                anchor_uid, score, reads = line.split('\t')[3], int(line.split('\t')[4]), int(line.split('\t')[6])
                
                #psuedocount for log safety
                if reads == 0:
                    reads = 1

                if anchor_uid not in eval_brkpts:
                    eval_brkpts[anchor_uid] = reads
                else:
                    prev_reads = eval_brkpts[anchor_uid]

                    # greater than 1.25-fold or 0.22314355131420976 = np.log(125/100)
                    check_1 = stats.binom_test([prev_reads, reads])
                    check_2 = abs(np.log(prev_reads/float(reads)))
                    
                    if ((prev_reads + reads) >= 20) and (check_1 <= 0.05) and (check_2 >= 0.22314355131420976):
                        if anchor_uid in anchor_to_brkpt:
                            
                            deets = anchor_to_brkpt[anchor_uid]
                            
                            chromo, start, stop, anchor_uid, start_median, _dot = deets.split('\t')
                            start, stop, start_median = int(start), int(stop), int(start_median)
                            
                            anchor_id = anchor_uid.split('.')[0]
                            if anchor_id not in quant_score:
                                quant_score[anchor_id] = [start_median, check_1, check_2]
                            else:
                                prev_start_median, prev_check_1, prev_check_2 = quant_score[anchor_id]
                                quant_score[anchor_id] += [start_median+prev_start_median, prev_check_1+check_1, prev_check_2+check_2]
                                
                            if chromo not in chromo_regions:
                                chromo_regions[chromo] = []
                                
                            if chromo not in passed_brkpts:
                                passed_brkpts[chromo] = {}
                                
                            juggle_dict = passed_brkpts[chromo]
                            for nt in range(start, stop+1):
                                if nt not in juggle_dict:
                                    juggle_dict[nt]=[anchor_id]
                                else:
                                    juggle_dict[nt].append(anchor_id)
                                
                                region_range = range(start-50, stop+51) 
                                if region_range not in chromo_regions[chromo]:
                                    chromo_regions[chromo].append(region_range)
                            passed_brkpts[chromo] = juggle_dict
                            
            temp_count_file.close()
            
            anchor_index = 0
            for chromo, region_set in chromo_regions.items():
                juggle_dict = passed_brkpts[chromo]
                
                for region_range in region_set:
                    hit_dict = {}
                    max_anchor = [0,'']
                    for nt in region_range:
                        if nt in juggle_dict:
                            anchor_list = juggle_dict[nt]
                            
                            if len(anchor_list) > 1:
                                for anchor_id in anchor_list:
                                    if anchor_id in hit_dict:
                                        hit_dict[anchor_id] += quant_score[anchor_id][0]
                                    else:
                                        hit_dict[anchor_id] = quant_score[anchor_id][0]
                            else:
                                if anchor_id in hit_dict:
                                    hit_dict[anchor_id] += quant_score[anchor_id][0]
                                else:
                                    hit_dict[anchor_id] = quant_score[anchor_id][0]
                    
                    for anchor_id, score in hit_dict.items():
                        if score > max_anchor[0]:
                            max_anchor = [score, anchor_id]
                        
                        else: 
                            if (score == max_anchor[0]) and (anchor_id != max_anchor[1]):
                                if quant_score[anchor_id][2] != quant_score[max_anchor[1]][2]:
                                    if max(quant_score[anchor_id][2],quant_score[max_anchor[1]][2]) == quant_score[anchor_id][2]:
                                        max_anchor = [score, anchor_id]
                                else:
                                    if quant_score[anchor_id][1] != quant_score[max_anchor[1]][1]:
                                        if min(quant_score[anchor_id][1],quant_score[max_anchor[1]][1]) == quant_score[anchor_id][1]:
                                            max_anchor = [score, anchor_id]
                    
                    if max_anchor[0] != 0:
                        start = min(region_range)-50
                        stop = max(region_range)+50
                        
                        outline = ('{}\t{}\t{}\t{}\t{}\t.\n').format(chromo, start, stop, anchor_index, quant_score[max_anchor[1]][0])
                        anchor_index += 1
                        print(outline)
                        outfile.write(outline)
                                                    
            outfile.close()
        
    else:
        print('No significant breaks identified.')
        
    print('Processing sam file for samfile recreation...')
    samfile = open(samfile_name)
    outsamfile_name = outfile_name + '_background'
    outsamfile = open(outsamfile_name + '.sam', 'w')
    make_out = False

    for line in samfile:
          
        if line[0] == '@':
            outsamfile.write(line)
        
        if line[0] != '@':
            uid = line.split('\t')[0]
            
            if uid in cleared_uid_list:
                make_out = True
                outsamfile.write(line)
                                        
    samfile.close()
    outsamfile.close()
    
    if make_out:
        _ = convert_sort(outsamfile_name)
            
def convert_to_fasta(fastq_filename, spliton, uid_list=[], filter_uid=False):
    fasta_file_name = fastq_filename.split(spliton)[0] + '_temp.fa'
    fasta_file = open(fasta_file_name,'w')
    
    fastq_file = open(fastq_filename)
    line_ct = 4
    
    for line in fastq_file:
        line_ct += 1
        
        if line[0] == '@':
            line = line.strip()
            
            uid = line.split('@')[1].split(' ')[0]
            if filter_uid:
                if uid in uid_list:
                    line_ct = 0
            else:
                line_ct = 0
            
        if line_ct == 1:
            outline = ('>{}\n{}\n').format(uid, line.strip())
            fasta_file.write(outline)
            
    fasta_file.close()
    fastq_file.close()
    
    return(fasta_file_name)
    
def gff_from_json(json_file_name, gff_file_name, rp_seq, max_eval, min_coverage):
    active = False
    query_deets_dict = {}
    chromo_dict={}
        
    try:
        data = json.load(open(json_file_name))
    except:
        print('json file ' + json_file_name + ' not found')
        
    for report_index in range(len(data["BlastOutput2"])):
        data_dict = (data["BlastOutput2"][report_index])
        for each_report in data_dict.items():
            for a_key, a_value in enumerate(each_report):
                if type(a_value)==dict:                   
                    for b_key, b_value in a_value.items():
                        if type(b_value)==dict:
                            for c_key, c_value in b_value.items():
                                if ('bl2seq' in c_key) and (type(c_value)==list):
                                    hit_dict = c_value[0]
                                    
                                    for d_key, d_value in hit_dict.items():                                        
                                        q_title = str(hit_dict['query_title'])
                                          
                                        if ('hits' in d_key) and (type(d_value)==list) and (len(d_value)>0):
                                            for each_hits in d_value:
                                                
                                                for e_key, e_value in each_hits.items():
                                                
                                                    base = q_title + '.'+str(each_hits['num']) 
                                                    chromo = each_hits['description']
                                                    chromo = chromo[0]
                                                    chromo_dict[base] = str(chromo['id'])
                                                        
                                                    if (e_key == 'hsps') and (type(e_value)==list):
                                                        for e_index in range(len(e_value)):
                                                            each_hsps = e_value[e_index]
    
                                                            numb = str(base)+'.'+str(each_hsps['num'])
                                                            
                                                            if len(numb)>1:
                                                                active = True
                                                                
                                                                hit_from = int(each_hsps["hit_from"])                                                                
                                                                hit_to = int(each_hsps["hit_to"])                                                                
                                                                query_from = str(each_hsps["query_from"])                                                                
                                                                query_to = str(each_hsps["query_to"])                                                            
                                                                bit_score = float(each_hsps["bit_score"])
                                                                evalue_score = float(each_hsps["evalue"])
                                                                query_strand = str(each_hsps["query_strand"])                                                                
                                                                hit_strand = str(each_hsps["hit_strand"])                                                                
                                                                qseq = str(each_hsps["qseq"])                                                            
                                                                hseq = str(each_hsps["hseq"])
                                                            
                                                            if (hit_to - hit_from)/float(len(rp_seq)) < float(min_coverage):
                                                                active = False
                                                            
                                                            if evalue_score > max_eval:
                                                                active = False
                                                            
                                                            if active:
                                                                active = False
                                                                query_deets_dict[numb] = ['q_id','hit_from','hit_to','query_from','query_to','bit_score','query_strand','hit_strand','qseq', 'hseq','q_title']
                                                                query_deets_dict[numb][0] = base
                                                                query_deets_dict[numb][1] = hit_from
                                                                query_deets_dict[numb][2] = hit_to
                                                                query_deets_dict[numb][3] = query_from
                                                                query_deets_dict[numb][4] = query_to
                                                                query_deets_dict[numb][5] = bit_score
                                                                query_deets_dict[numb][6] = query_strand
                                                                query_deets_dict[numb][7] = hit_strand
                                                                query_deets_dict[numb][8] = qseq
                                                                query_deets_dict[numb][9] = hseq
                                                                query_deets_dict[numb][10] = q_title
                                                                numb = 0

    if query_deets_dict:
        gff_file = open(gff_file_name, 'w')
        
        ct = 0
        outline = ('Identified {} reads containing a near match to the reflection point sequence.').format(len(query_deets_dict))
        print(outline)
        
        for numb, deets in query_deets_dict.items():
            q_title = query_deets_dict[numb][10]
                    
            base = query_deets_dict[numb][0]
            reduced_base = query_deets_dict[numb][0]
            hypo_id = str(reduced_base + '_' + str(ct))
                    
            chromo = chromo_dict[base]
            
            hit_from = query_deets_dict[numb][1]
            hit_to = query_deets_dict[numb][2]
                        
            if hit_from < hit_to:
                start = hit_from
                stop = hit_to
            else:
                start = hit_to
                stop = hit_from
                        
            if query_deets_dict[numb][7] == 'Plus':
                sign = '+'
            else:
                sign = '-'
                
            bit_score = float(query_deets_dict[numb][5])
            
            if query_deets_dict[numb][6] != query_deets_dict[numb][7]:
                orient = 'reverse'
            else:
                orient = 'forward'
                
            mod_seq = ''
            qseq = query_deets_dict[numb][8]
            hseq = query_deets_dict[numb][9]
            
            for index in range(len(qseq)):
                if qseq[index] == hseq[index]:
                    mod_seq += qseq[index].upper()
                if qseq[index] != hseq[index]:
                    mod_seq += qseq[index].lower()
            
            node_loci = ('{},{},{},{}').format(q_title, query_deets_dict[numb][3], query_deets_dict[numb][4], float(query_deets_dict[numb][5]))
            gff_line = ('{}\tmugio\t{}_reflecton_point\t{}\t{}\t.\t{}\t{}\tnode_uid={}; orient={}; alt_seq={}; ref_seq={}; contig={}\n').format(chromo, hypo_id, start, stop, sign, int(round(bit_score)), node_loci, orient, mod_seq, hseq, rp_seq)
            gff_file.write(gff_line)
                        
            ct+=1
            
        gff_file.close()
        
    return(len(query_deets_dict))
    
if args.reflection_point_blast:
    #TODO blast
    #The use of this module is to identify reads that contain close matches to the reflection point identified by short-read sequencing
    os_mkdir(args.outfile_name)
    
    outfile_name = args.outfile_name
    rp_seq = args.reflection_point_seq
    
    if ('.fastq' in args.fastqfile_name) or ('.fq' in args.fastqfile_name):
        if ('.fastq' in args.fastqfile_name):
            fasta_file_name = convert_to_fasta(args.fastqfile_name, '.fastq')
        else:
            fasta_file_name = convert_to_fasta(args.fastqfile_name, '.fq')
    else:
        fasta_file_name = args.fastqfile_name
        
    query_file_name = ('{}_temp_query.fa').format(outfile_name)
    query_file = open(query_file_name,'w')
    
    outline = ('>temp_query\n{}').format(rp_seq)
    query_file.write(outline)
    query_file.close()
    
    json_file_name = ('{}_blast.json').format(outfile_name)
    gff_file_name = ('{}_reflection_point.gff').format(outfile_name)
    
    bash_command = ('blastn -task "blastn-short" -query {} -subject {} -outfmt 15 -out {}').format(query_file_name, fasta_file_name, json_file_name)
    output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
    if output:
        print(output)
    
    _hit_ct = gff_from_json(json_file_name, gff_file_name, rp_seq, 0.01, 0.51)

if args.background:
    os_mkdir(args.outfile_name)
    #TODO background
                    
    fastqfile_name = (args.fastqfile_name)
    samfile_name = (args.samfile_name)
    phred_stats_file_name = ('{}_phred_stats.log').format(args.outfile_name)
    phred_stats_file = open(phred_stats_file_name,'w')
    
    if args.burnout:
        burnout = True
        minimum_reads = 1000
        unchanged = 10
        steps = 100
        
        if args.minimum_reads:
            minimum_reads = int(args.minimum_reads)
            
        if args.unchanged:
            unchanged = int(args.unchanged)
            
        if args.steps:
            steps = int(args.steps)
            
    else:
        burnout = False
        steps = 100
        
    print('Processing fastq file to build complete fastq phred set ...')
    infile = open(fastqfile_name)
    first_fastq_dict = {}
    line_ct = 4
    
    all_phred_list = []
    each_phred_list = []
    longest_line = 0

    character_count = {}

    for line in infile:
        line_ct += 1
        
        if line[0] == '@':
            uid = line.split('@')[1].split(' ')[0].strip()
            line_ct = 0
                                        
        if line_ct == 3:
            phred_line = line.strip()
            if phred_line:
                        
                if len(phred_line) >= longest_line:
                    longest_line = len(phred_line)
                
                first_fastq_dict[uid] = phred_line
                
                phred_list = ord_phred(phred_line)
                
                if phred_list:
                    all_phred_list.append(np.median(phred_list))    
            
                for each in phred_list:
                    each_phred_list.append(each)
                    if each not in character_count:
                        character_count[each] = 1
                    else:
                        character_count[each] += 1
            
    infile.close()

    phred_rates = random.sample(each_phred_list, 10*longest_line)
    
    if args.lower_threshold:
        read_threshold= float(args.lower_threshold)
    else:
        read_threshold = np.median(all_phred_list)-(5*np.std(all_phred_list))
        
    if args.middle_threshold:
        read_mid_threshold = float(args.middle_threshold)
    else:
        read_mid_threshold = np.median(all_phred_list)
        
    if args.upper_threshold:
        read_stop_threshold = float(args.upper_threshold)
    else:
        read_stop_threshold = np.median(all_phred_list)+(2*np.std(all_phred_list))
                    
    read_deets = [phred_rates, read_threshold, read_mid_threshold, read_stop_threshold]
    
    outline = ('{} median phred score of {} bases. {} minimum threshold for discovery.').format(np.median(all_phred_list), len(each_phred_list), read_threshold) 
    print(outline)
    phred_stats_file.write(outline)
    
    #short stats:
    char_list = []
    for char, value in character_count.items():
        char_list.append(char)
    char_list.sort()
    
    val_list = []
    for each in char_list:
        val_list.append(character_count[each])
        
    y_pos = np.arange(len(char_list))
    
    plt.bar(y_pos, val_list, align='center', alpha=0.5)
    plt.xticks(y_pos, char_list)
    plt.ylabel('Count')
    plt.title('Phred score')
    plt.savefig('phred_score_count.png')
            
    low_zone_uids = {}
    isct = 0
    pct_history = []
    stopping = False
    
    for uid, phred_line in first_fastq_dict.items():
        if not stopping:
            isct += 1
            pct_so_far = round(len(low_zone_uids)/float(isct),2)
            
            if isct % steps == 0:
                outline = ('{} reads processed out of {}. {} ({}%) reads have siginificant low zones').format(isct, len(first_fastq_dict), len(low_zone_uids), pct_so_far*100)
                print(outline)
                phred_stats_file.write(outline)
            
            if burnout:
                pct_history.append(pct_so_far)

                if isct >= minimum_reads:                    
                    if len(pct_history) >= unchanged and pct_history:
                        if round(np.median(pct_history),2) == pct_so_far:
                            outline = ("Stopping collection, median of reads with low zones unchanged over {} cycles of {}. \n{}% of reads have low phred zones.").format(unchanged, steps, pct_so_far*100)
                            print(outline)
                            phred_stats_file.write(outline)
                            stopping = True
                    
            case, resolved_zone_dict = low_zones(phred_line, read_threshold, read_mid_threshold, phred_rates)
            
            if case:
                low_zone_uids[uid] = resolved_zone_dict
            
    print('Processing sam file visualization...')
    samfile = open(samfile_name)
    genome_file = open('temp.genome','w')
    
    outfile_name = (args.outfile_name + '_low')
    outfile = open(outfile_name +'.sam', 'w')
        
    for line in samfile:
        if line[0] == '@':
            outfile.write(line)
                                   
        if line[0] != '@':                    
            uid = line.split('\t')[0]
            if uid in low_zone_uids:
                outfile.write(line)                         
    samfile.close()
    outfile.close()
                         
    region_number_per_read = []
    region_length = []
    median_phred = []
    for uid, resolved_zones in low_zone_uids.items():
        region_number_per_read.append(len(resolved_zones))
        
        for unid, deets in resolved_zones.items():
            region_length.append(deets[1])
            median_phred.append(deets[0])
    
    print('Processing low phred regions file visualization...')
    if region_number_per_read:
        outline = ('Total reads: {}\nReads with low phred score regions: {}\nMedian regions per read: {}\n').format(isct,len(low_zone_uids),np.median(region_number_per_read))
        print(outline)
        phred_stats_file.write(outline)

    if region_length:
        outline = ('Median region length: {}\n Minimum Length: {}\n Maximum Length: {}\n').format(np.median(region_length), min(region_length), max(region_length))
        print(outline)
        phred_stats_file.write(outline)
    
    if median_phred:
        outline = ('Median median phred score: {}\n Minimum median phred score: {}\n Maximum median phred score: {}\n').format(np.median(median_phred), min(median_phred), max(median_phred))
        print(outline)
        phred_stats_file.write(outline)
    
    phred_stats_file.close()
    
    print('Converting sam and generating depth file.')
    if convert_sort(outfile_name):
        depth_filename = ('{output_file}.depth').format(output_file=outfile_name)
        depth_df = pd.read_csv(depth_filename, sep='\t')
            

if args.evaluate:
    os_mkdir(args.outfile_name)
    
    brkpt_list = []
    
    if args.breakpoint_file:
        breakpoint_file = open(args.breakpoint_file)
        for line in breakpoint_file:
            if line[0]!='#':
                line = line.strip()
                chromo = line.split('\t')[0]
                start = line.split('\t')[1]
                stop = line.split('\t')[2]
                #TODO remove filter
                if (chromo != 'NC_001224.1'):
                    outline = ('{}:{}-{}').format(chromo, start, stop)
                    brkpt_list.append(outline)
            
    if args.sniffle_file:
        '''because of the way structural variants are defined by sniffles some 
        of SV sites are 100-10Ks long. Because of the read criteria these 
        overlong regions don't seem to induce artifacts but it could be a potential problem.   
        '''
        sniffle_file = open(args.sniffle_file)
        for line in sniffle_file:
            if line[0]!='#':
                if 'PASS' in line.split('\t')[6]:
                    line = line.strip()
                    chromo1 = line.split('\t')[0]
                    chromo2 = line.split('\t')[0]
                    start = int(line.split('\t')[1])
                    stop = int(line.split('\t')[1])

                    deets = line.split('\t')[7]                    
                    if 'CHR2=' in deets:
                        chromo2 = deets.split('CHR2=')[1].split(';')[0]
                    if 'END=' in deets:
                        stop = int(deets.split('END=')[1].split(';')[0])
        
                    if (chromo1 == chromo2) and (start != stop):
                        #TODO remove filter
                        if (chromo1 != 'NC_001224.1'):
                            outline = ('{}:{}-{}').format(chromo1, start, stop)
                            brkpt_list.append(outline)
        
    if not args.breakpoint_file and not args.sniffle_file:
        brkpt = args.breakpoint
        brkpt_list.append(brkpt)
            
    outline = ('The following breakpoints will be evaluated: {}').format(brkpt_list)
    print(outline)
        
    fastqfile_name = (args.fastqfile_name)
    samfile_name = (args.samfile_name)
                    
    #build sets of reads that span the breakpoint
    uid_brkpt_dict = {}
    uid_brkpt_sites = {}
    
    print('Processing fastq file to build complete fastq phred set ...')
    infile = open(fastqfile_name)
    first_fastq_dict = {}
    line_ct = 4
    
    all_phred_list = []
    each_phred_list = []
    longest_line = 0

    for line in infile:
        line_ct += 1
        
        if line[0] == '@':
            uid = line.split('@')[1].split(' ')[0].strip()
            line_ct = 0
                                        
        if line_ct == 3:
            phred_line = line.strip()
            if phred_line:
                if len(phred_line) >= longest_line:
                    longest_line = len(phred_line)
                
                first_fastq_dict[uid] = phred_line
                
                phred_list = ord_phred(phred_line)
                
                if phred_list:
                    all_phred_list.append(np.median(phred_list))    
            
                for each in phred_list:
                    each_phred_list.append(each)                    
            
    infile.close()
        
    phred_rates = random.sample(each_phred_list, 10*longest_line)
    read_threshold = np.median(all_phred_list)-(5*np.std(all_phred_list))
    read_mid_threshold = np.median(all_phred_list)-(3*np.std(all_phred_list))
    read_upper_threshold = np.median(all_phred_list)
    read_stop_threshold = np.median(all_phred_list)+(2*np.std(all_phred_list))
                
    threshold_deets = [phred_rates, read_threshold, read_mid_threshold, read_upper_threshold, read_stop_threshold]
    
    outline = ('{} median phred score of {} bases. {} minimum threshold for discovery.').format(np.median(all_phred_list), len(each_phred_list), read_threshold) 
    print(outline)
    
    print('Processing sam file for breakpoint retrieval and definition (bprd)...')
    samfile = open(samfile_name)
    
    supplemental_align_set = set()
    read_ct = {}
    read_dict = {}
    chromo_size_dict = {}

    for line in samfile:
        if line[0] == '@':
            if 'LN:' in line:
                line = line.strip()
                chromo = line.split('SN:')[1].split('\t')[0]
                size = int(line.split('LN:')[1])
                if chromo not in chromo_size_dict:
                    chromo_size_dict[chromo] = size
                else:
                    print('chromosome name duplication in sam file')
                
        if line[0] != '@':
            if unpackbits(np.array([int(line.split('\t')[1])]))[0][8] != 1:
                uid = line.split('\t')[0]
                chromo = line.split('\t')[2]
                
                if unpackbits(np.array([int(line.split('\t')[1])]))[0][11] == 1:
                    supplemental_align_set.add(uid)                        
                    
                start = int(line.split('\t')[3])
                strand = unpackbits(np.array([int(line.split('\t')[1])]))[0][4]
                cigar = line.split('\t')[5]
                
                stop = start + distance_in(strand, cigar, 'reference') - 1
                
                if uid not in read_ct:
                    read_ct[uid] = 0
                else:
                    read_ct[uid] += 1
                    
                uid_ct = ('{}.{}').format(uid, read_ct[uid])
                read_dict[uid_ct] = [strand, start, stop, line.split('\t')[10], cigar, chromo]
                                                        
    samfile.close()           
            
    if brkpt_list:
        outdir = args.outfile_name.rsplit('/', 1)[0]
        bash_command = ('mkdir -p {}/traces').format(outdir)
        output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
        if output:
            print(output)
        
    for brkpt in brkpt_list:
        brkpt_clean = brkpt.replace(':','_')

        bash_command = ('mkdir -p {}/traces/{}').format(outdir, brkpt_clean)
        output = subprocess.check_output([bash_command],stderr=subprocess.STDOUT,shell=True)
        if output:
            print(output)
        
        outline = ('{}/traces/{}/{}').format(outdir, brkpt_clean, brkpt_clean)
        is_outfile_name = outline
        
        length_function(brkpt, threshold_deets, first_fastq_dict, supplemental_align_set, read_ct, read_dict, is_outfile_name)
