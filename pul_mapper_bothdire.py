#!/usr/bin/env python3

import sys
import time

#~#~#~#~#~#~#~#~#
# This uses the marker file created by marker_spotter.py to search for the locations of the markers and then check up and downstream for cazymes, stopping if it cannot find any within 500bp or within 2 extra genes, in attempt to catch regulatory genes or transporters
#~#~#~#~#~#~#~#~#

in_faa = sys.argv[1]
in_marker = sys.argv[2]
in_cazyme = sys.argv[3]
in_pfam = sys.argv[4]
outfi = sys.argv[5]
start = time.perf_counter()

markerhits = []
faahits = []

#~#~#~#~#~#~#~#~#

#process marker genes and find nearby genes in faa
with open(in_marker, 'r') as markerlist, open(in_faa, 'r') as faa, open(in_cazyme, 'r') as cazyme, open(in_pfam, 'r') as pfam, open(outfi, 'w') as ou:
    marker_lines = list(markerlist)  #read all marker lines
    faa_lines = list(faa) #read all faa lines
    cazyme_lines = list(cazyme)  #read all cazyme lines
    pfam_lines = list(pfam)  #read all pfam lines

    for marker_line in marker_lines:
        Parts = marker_line.strip().split('\t')
        if len(Parts) < 5:
            #print(f"Skipping malformed marker line: {marker_line.strip()}")  #used for debugging, remove later
            continue

        try:
            M_susD_ID = Parts[1]
            M_tonB_ID = Parts[2]
            MarkerStart = int(Parts[3])
            MarkerEnd = int(Parts[4])
        except ValueError:
            #print(f"Skipping marker line with invalid start/end positions: {marker_line.strip()}") #used for debugging, remove later
            continue

        #print(f"Processing marker: {M_susD_ID}, Start: {MarkerStart}, End: {MarkerEnd}")

        found_marker = False
        for i, faa_line in enumerate(faa_lines):
            if faa_line.startswith(">"):  #check if it's a protein line
                faaParts = faa_line.strip().split()
                faaID = faaParts[0].lstrip(">")  #remove ">" from the FAA ID
                if faaID == M_susD_ID or faaID == M_tonB_ID:  #match the exact MarkerID
                    found_marker = True

                    #lists for storing results
                    upstream_results = []
                    downstream_results = []

                    #Downstream search
                    downstream_counter = 0  
                    current_marker_end = MarkerEnd  
                    for j in range(i + 1, len(faa_lines)):
                        next_faa_line = faa_lines[j]
                        if next_faa_line.startswith(">"):  #check for protein description lines
                            next_faaParts = next_faa_line.strip().split()
                            try:
                                next_faaID = next_faaParts[0].lstrip(">")
                                next_faaStart = int(next_faaParts[2])
                                next_faaEnd = int(next_faaParts[4])
                            except (IndexError, ValueError):
                                continue

                            #checks if the next protein starts within 500bp of the current marker end
                            if abs(next_faaStart - current_marker_end) <= 500:
                                #search for the faaID in the cazyme file
                                cazyme_matches = []
                                for cazyme_line in cazyme_lines:
                                    cazymeParts = cazyme_line.strip().split('\t')
                                    if len(cazymeParts) < 4:
                                        continue
                                    cazymeID = cazymeParts[0]
                                    if cazymeID == next_faaID:
                                        cazyme_matches.append(cazymeParts[3])
                                cazyme_string = ", ".join(cazyme_matches) if cazyme_matches else "None"

                                #search for PFAM matches
                                pfam_matches = []
                                for pfam_line in pfam_lines:
                                    pfamParts = pfam_line.strip().split('\t')
                                    if len(pfamParts) < 4:
                                        continue
                                    pfamID = pfamParts[0]
                                    pfamDescription = pfamParts[4]  
                                    if pfamID == next_faaID:
                                        for domain in pfamDescription.split(','):
                                            pfam_matches.append(domain.strip())

                                pfam_string = ", ".join(pfam_matches) if pfam_matches else "None"

                                downstream_results.append(f"Gene: {next_faaID}, Start: {next_faaStart}, End: {next_faaEnd}, CAZymes: {cazyme_string}, pfam: {pfam_string}")

                                if cazyme_matches or any(
                                    any(x in pfam for x in ["SusD", "TonB", "CarboxypepD_reg", "TPR", "Glyco_hydro", "HTH_18", "SASA", "HisKA", "AraC", "DUF"])
                                    for pfam in pfam_matches
                                ):
                                    downstream_counter = 0  # Reset counter if we find a qualifying gene
                                else:
                                    downstream_counter += 1  # Increment only for non-qualifying genes
                                    if downstream_counter >= 4:  # Break only after 4 consecutive non-qualifying genes
                                        break

                                #update current_marker_end to the end of the current protein
                                current_marker_end = next_faaEnd
                            else:
                                #stop checking if the next protein is not within 500bp
                                break

                    # Upstream search
                    upstream_counter = 0  
                    current_marker_start = MarkerStart  
                    for j in range(i - 1, -1, -1):  #iterates backward = upstream
                        prev_faa_line = faa_lines[j]
                        if prev_faa_line.startswith(">"):  
                            prev_faaParts = prev_faa_line.strip().split()
                            try:
                                prev_faaID = prev_faaParts[0].lstrip(">")
                                prev_faaStart = int(prev_faaParts[2])
                                prev_faaEnd = int(prev_faaParts[4])
                            except (IndexError, ValueError):
                                continue

                            #checks if the previous protein ends within 500bp of the current marker's start
                            if abs(current_marker_start - prev_faaEnd) <= 500:
                                #search for the faaID in the cazyme file
                                cazyme_matches = []
                                for cazyme_line in cazyme_lines:
                                    cazymeParts = cazyme_line.strip().split('\t')
                                    if len(cazymeParts) < 4:
                                        continue
                                    cazymeID = cazymeParts[0]
                                    if cazymeID == prev_faaID:
                                        cazyme_matches.append(cazymeParts[3])
                                cazyme_string = ", ".join(cazyme_matches) if cazyme_matches else "None"

                                #searches for PFAM matches
                                pfam_matches = []
                                for pfam_line in pfam_lines:
                                    pfamParts = pfam_line.strip().split('\t')
                                    if len(pfamParts) < 4:
                                        continue
                                    pfamID = pfamParts[0]
                                    pfamDescription = pfamParts[4]  
                                    if pfamID == prev_faaID:
                                        pfam_matches.append(pfamDescription)
                                pfam_string = ", ".join(pfam_matches) if pfam_matches else "None"

                                upstream_results.append(f"Gene: {prev_faaID}, Start: {prev_faaStart}, End: {prev_faaEnd}, CAZymes: {cazyme_string}, pfam: {pfam_string}")

                                if cazyme_matches or any(
                                    any(x in pfam for x in ["SusD", "TonB", "CarboxypepD_reg", "TPR", "Glyco_hydro", "HTH_18", "SASA", "HisKA", "AraC", "DUF"])
                                    for pfam in pfam_matches
                                ):
                                    upstream_counter = 0  # Reset counter if we find a qualifying gene
                                else:
                                    #increment the counter if no CAZyme matches are found
                                    upstream_counter += 1  # Increment only for non-qualifying genes
                                    if upstream_counter >= 4:  # Break only after 4 consecutive non-qualifying genes
                                        break

                                # update current_marker_start to the start of the current protein
                                current_marker_start = prev_faaStart
                            else:
                                #stops checking if the previous protein is not within 500bp
                                break

                    # Write results for BOTH upstream and downstream
                    ou.write(f"Marker: {M_susD_ID}, Start: {MarkerStart}, End: {MarkerEnd}\n")
                    
                    # Downstream results
                    ou.write(f"  Downstream:\n")
                    for result in downstream_results:
                        ou.write(f"\t{result}\n")
                    

                    # Upstream results  
                    ou.write(f"  Upstream:\n")
                    for result in upstream_results:
                        ou.write(f"\t{result}\n")

                    
                    ou.write(f"\n")  #spacing between markers

                    print(f"Processed {M_susD_ID}: {len(upstream_results)} upstream, {len(downstream_results)} downstream results.")
                    

        if not found_marker:
            print(f"Marker {M_susD_ID} not found in FAA file.")

end = time.perf_counter()
elapsed = end - start
print(f'Annotated all proteins, time taken: {elapsed:.6f} seconds')