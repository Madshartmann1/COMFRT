######## 
# Processing of NUMTs
# Written by Kirstine Tersbøl Melsen, s215096
# Modified for pipeline integration
########

import numpy as np
import pandas as pd 
import os
import re
import sys

# Check for command line arguments
if len(sys.argv) != 3:
    print("Usage: python FragmentsMappingTwicePY.py <input_file> <mt_scaffold_name>")
    sys.exit(1)

input_file = sys.argv[1]
mt_scaffold = sys.argv[2]

### Opening the file
with open(input_file, 'r', encoding='utf-8') as infile:
    content = [line.split() for line in infile.readlines()]

### Parse the string of the 20th column into a tuple of 2 lists: 1) 'scaffold', 2) 'mismatches'

def parse_xa_field(xa_field: str):
    
    if not xa_field.startswith('XA:Z:'):
        return []
    payload = xa_field[5:]  # strip 'XA:Z:'
    entries = [e for e in payload.split(';') if e] # only include if the string is non-empty
    out = []
    for e in entries:
        parts = e.split(',')
        if len(parts) != 4:
            return []
        scaffold = parts[0]
        try:
            mismatches = int(parts[3])
        except ValueError:
            return []
        out.append({'scaffold': scaffold, 'mismatches': mismatches})
    return out

def findingNUMTs(content, mt_scaffold):
    mapping_mt = []
    mapping_mt_nucl = []

    for line in content:

        # The scaffold of the primary mapping is defined
        primary_scaf = line[2]

        # The scaffold(s) of the secondary mappings are extracted into a tuple
        sec_list = []
        if len(line) == 20:
            sec_list = parse_xa_field(line[19])

        # Sequences mapping only to the MT genome
        if (primary_scaf == mt_scaffold) and (len(sec_list) == 0 or all((s['scaffold'] == mt_scaffold) for s in sec_list)):
            mapping_mt.append(line)

        ## Sequences mapping to the MT and the nuclear genome
        # Occurences of mapp11ing to the MT genome
        # If mt_count = 0, then there are no MT mappings
        # If mt_count = 1, then there is exactly one MT mapping - which is what we want.
        # If mt_count >= 2, then there is a mix of MT and nucl mappings, but more than one MT mappings, which is not good. 
        mt_count = (1 if primary_scaf == mt_scaffold else 0) + sum(1 for s in sec_list if (s['scaffold'] == mt_scaffold))
        if len(line) == 20 and mt_count == 1:
            mapping_mt_nucl.append(line)
        else:
            pass

    records = [{'Reads that map only to the mitogenome': len(mapping_mt), 'Reads that map once to the mitogenome and at least once to the nuclear genome': len(mapping_mt_nucl)}]
    
    numts = []
    betterMT_thanNucl = []

    for line in mapping_mt_nucl:
        # Again parse the 20th column 
        sec_list = []
        if len(line) == 20:
            sec_list = parse_xa_field(line[19])
        
        # Collect the scaffolds into one list and all mismatch parameters into another list
        scaffolds = [line[2]]
        mismatch_para = [line[12].split(':')[-1]]
        for s in sec_list:
            scaffolds.append(s['scaffold'])
            mismatch_para.append(s['mismatches'])
        
        if len(scaffolds) != len(mismatch_para):
            pass

        idx = 0
        for k in scaffolds:
            if mt_scaffold in k:
                mt_idx = idx
            else: 
                pass
            idx += 1
        
        mismatch_para_array = np.array(mismatch_para)


        if (mismatch_para_array[mt_idx]) == 0 or (mismatch_para_array[mt_idx] == min(mismatch_para_array)):
            betterMT_thanNucl.append(line[0])
        else:
            numts.append(line[0])

    records.append({'Reads that map better to the nuclear genome (numts)': len(numts)})
    records.append({'Reads that map better or equal to the mitogenome than to the nuclear genome': len(betterMT_thanNucl)})

    return records, numts, betterMT_thanNucl


# Run the analysis
records, numts, betterMT_thanNucl = findingNUMTs(content, mt_scaffold)

# Output results to stdout
# First output the records (will be captured as first 4 lines)
for record in records:
    for key, value in record.items():
        print(f"{key}: {value}")

# Then output the NUMTS read IDs (one per line, starting from line 5)
for read_id in numts:
    print(read_id)