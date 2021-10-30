# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 12:15:23 2021

@author: liuh
"""

import pysam   
import os,re
import argparse, sys
parser = argparse.ArgumentParser(description='Removing PCR duplicates with mapping coordinates and cell barcode.')
parser.add_argument('infile', action = 'store', nargs = '?',
                    type = str, default = sys.stdin, help = 'Input mapping paired bam')
parser.add_argument('outfile',action = 'store', nargs = '?',
                    type = argparse.FileType('w'), default = sys.stdout,help = 'Output file name')
args = parser.parse_args()

in_file = pysam.AlignmentFile(args.infile, "rb")
out_file = pysam.AlignmentFile(args.outfile, "wb", template=in_file)

in_file = pysam.AlignmentFile("BoneMarrow_62016.bam", "rb")
out_file = pysam.AlignmentFile("Test.bam", "wb", template=in_file)

header_chr = in_file.references
read_dict = {}

#sample_name =os.path.basename(args.infile)
#print(sample_name + ' is removing duplicates with cell barcodes, chr, chr_start, and fragment_length.')

for read in in_file.fetch('chr1', 1, 5000000):
    is_valid = (read.mapq >= 20 and read.reference_name == read.next_reference_name and read.flag in [83, 163, 99,147])
    if is_valid:
        if read.query_name not in read_dict:
            read_dict[read.query_name] = [read]
        else:
            read_dict[read.query_name].append(read)
                       
uniq_alignments = {}

for k in read_dict.keys():
    if len(read_dict[k]) == 2:
        soft_clip_5 = re.compile(r'^(\d+)S')
        soft_clip_3 = re.compile(r'(\d+)S$')
        frag_end_pos = 0
        
        ## create a pair of new reads
        global left_r
        global right_r
        left_r =  pysam.AlignedSegment(header = in_file.header)
        right_r = pysam.AlignedSegment(header = in_file.header)
        
        for r in read_dict[k]:
            ## accounting for soft-clipping at both ends
            ## shifting reads mapped + strand by +4, and - strand by -5
            
            cigar = r.cigarstring
            if r.flag in [99, 163]:  # + strand alignment
                left_r.query_name = r.query_name
                left_r.flag = r.flag
                left_r.reference_id = r.reference_id
                left_r.mapping_quality = r.mapping_quality
                left_r.next_reference_id = r.next_reference_id
                left_r.next_reference_start = r.next_reference_start
                match_5 = re.search(soft_clip_5, cigar)
                if match_5:
                    left_r.reference_start = r.reference_start - match_5.group(1) + 4
                else:
                    left_r.reference_start = r.reference_start + 4
                left_r.query_sequence = r.query_sequence[4:]
                left_r.query_qualities = r.query_qualities[4:]
                left_r.cigar =((0, r.reference_end - left_r.reference_start),)  # Returned cigar is only matches
            else:
                right_r.query_name = r.query_name
                right_r.flag = r.flag
                right_r.reference_id = r.reference_id
                right_r.mapping_quality = r.mapping_quality
                right_r.next_reference_id = r.next_reference_id
                right_r.reference_start = r.reference_start
                
                match_3 = re.search(soft_clip_3, cigar)
                if match_3:
                    frag_end_pos = r.reference_end + match_3.group(1) - 5
                else:
                    frag_end_pos = r.reference_end - 5
                right_r.query_sequence = r.query_sequence[:-5]
                right_r.query_qualities = r.query_qualities[:-5]
                right_r.next_reference_start = left_r.reference_start
                right_r.cigar = ((0, frag_end_pos - r.reference_start),)
            left_r.template_length = frag_end_pos - left_r.reference_start
            right_r.template_length = - left_r.template_length
            left_r.set_tag(tag = "CB", value = r.query_name.split(":")[0], value_type = "Z")
            right_r.set_tag(tag = "CB", value = r.query_name.split(":")[0], value_type = "Z")
                
            
    # output alignments from unique fragments        
    frag_id = ":".join([left_r.get_tag(tag = "CB"), left_r.reference_name, 
                        str(left_r.reference_start), str(abs(left_r.template_length))])
    print frag_id
    
    if frag_id not in uniq_alignments:
        uniq_alignments[frag_id] = 1
        for r in [left_r, right_r]:
            out_file.write(r)
            
    else:
        uniq_alignments[frag_id] +=  1
        
out_file.close()
        
print("Done!")
