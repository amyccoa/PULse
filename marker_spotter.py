#!/usr/bin/env python3

import sys
from itertools import islice
import time
from collections import deque

#~#~#~#~#~#~#~#~#
# Searches for susD and tonB genes in the organised protein file, writes to output with the start+end of both genes together to create a list of PUL markers
#~#~#~#~#~#~#~#~#

infi = sys.argv[1]  
outfi = sys.argv[2]  
start = time.perf_counter()

#~#~#~#~#~#~#~#~#

with open(infi, 'r') as fi, open(outfi, 'w') as ou:
    lines = iter(fi)
    #stores prev lines for checking later
    buffer = deque(maxlen=10)

    for ln in lines:
        Parts = ln.strip().split()
        if len(Parts) < 3:
            print(f"Skipping malformed line: {ln.strip()}")
            continue  

        susdID = Parts[0]
        try:
            susdStart = int(Parts[2])  
            susdEnd = int(Parts[3])    
        except ValueError:
            print(f"Skipping line with invalid start/end positions: {ln.strip()}")
            continue  #skip lines with invalid start/end positions

        #check for susD
        if "SusD" in ln:
            print(f"Found SusD: {ln.strip()}")  #debugging output

            #check above for tonb
            found_above = None
            for line in buffer:
                if "TonB" in line:
                    TonB_Parts = line.strip().split()
                    tonbID = TonB_Parts[0]
                    try:
                        tonbStart = int(TonB_Parts[2])
                        tonbEnd = int(TonB_Parts[3])
                        if abs(susdStart - tonbEnd) <= 500:  #check distance is between 500bp
                            found_above = line
                            break
                    except ValueError:
                        print(f"Skipping malformed TonB line: {line.strip()}")
                        continue

            # Check below for tonb
            check_lines = list(islice(lines, 10))
            found_below = None
            for line in check_lines:
                if "TonB" in line:
                    TonB_Parts = line.strip().split()
                    tonbID = TonB_Parts[0]
                    try:
                        tonbStart = int(TonB_Parts[2])
                        tonbEnd = int(TonB_Parts[3])
                        if abs(tonbStart - susdEnd) <= 500:  #check distance is between 500bp
                            found_below = line
                            break
                    except ValueError:
                        print(f"Skipping malformed TonB line: {line.strip()}")
                        continue

            #write to file, if found and within 500bp
            if found_above:
                ou.write(f"SusD+TonB pair found\t{susdID}\t{tonbID}\t{tonbStart}\t{susdEnd}\n")
                #print(f"SusD+TonB pair found: {susdID}\t{tonbID}")  #debugging output
            if found_below:
                ou.write(f"SusD+TonB pair found\t{susdID}\t{tonbID}\t{susdStart}\t{tonbEnd}\n")
                #print(f"SusD+TonB pair found: {susdID}\t{tonbID}")  #debugging output

            #add peeked lines back to the buffer
            for line in reversed(check_lines):
                buffer.append(line)

        #add current line to the buffer
        buffer.append(ln)

#~#~#~#~#~#~#~#~#

end = time.perf_counter()
elapsed = end - start
print(f'PULs written, time taken: {elapsed:.6f} seconds')