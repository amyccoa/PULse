#!/usr/bin/env python3

import sys
import os
import time
from collections import deque, defaultdict  

def log(msg):
    now = time.strftime('%H:%M:%S')
    print(f"[{now}] {msg}")


### everything is now indexed, this should speed up the code but also make it more safe for big searches
def parse_faa_headers(faa_lines):
    parsed = []
    for line in faa_lines:
        if line.startswith(">"):
            parts = line.strip().split()
            try:
                gene_id = parts[0].lstrip(">")
                start = int(parts[2])
                end = int(parts[4])
                parsed.append((gene_id, start, end, line))
            except (IndexError, ValueError) as e:
                log(f"Failed to parse line: {line.strip()} — {e}")
                continue
    return parsed

def index_annotation_lines(lines, min_fields=4, field_index=3, sep=','):
    data = defaultdict(list)
    for line in lines:
        parts = line.strip().split('\t')
        if len(parts) >= min_fields:
            key = parts[0]
            values = parts[field_index].split(sep)
            data[key].extend([v.strip() for v in values if v.strip()])
    return data

def find_nearby_genes(parsed_faa, index, direction, limit_bp, max_hits, cazy_dict, pfam_dict, relevant_pfams):
    results = []
    recent_hits = deque(maxlen=max_hits)

    if direction == "down":
        rng = range(index + 1, len(parsed_faa))
        get_distance = lambda cur_end, next_start: abs(next_start - cur_end)
        update_anchor = lambda next_end: next_end
        anchor = parsed_faa[index][2]  # end
    else:
        rng = range(index - 1, -1, -1)
        get_distance = lambda cur_start, prev_end: abs(cur_start - prev_end)
        update_anchor = lambda prev_start: prev_start
        anchor = parsed_faa[index][1]  # start

    for j in rng:
        gene_id, start, end, _ = parsed_faa[j]
        distance = get_distance(anchor, start if direction == "down" else end)
        if distance > limit_bp:
            break

        cazy_hits = cazy_dict.get(gene_id, [])
        pfam_hits = pfam_dict.get(gene_id, [])
        cazy_str = ", ".join(cazy_hits) if cazy_hits else "None"
        pfam_str = ", ".join(pfam_hits) if pfam_hits else "None"

        results.append(f"{gene_id}, {start}, {end}, CAZymes: {cazy_str}, pfam: {pfam_str}")

        hit = bool(cazy_hits) or any(any(k in pf for k in relevant_pfams) for pf in pfam_hits)
        recent_hits.append(hit)
        if len(recent_hits) == max_hits and not any(recent_hits):
            log(f"Stopped early at {gene_id}, no relevant hits in last {max_hits}")
            break

        anchor = update_anchor(end if direction == "down" else start)

    return results

###upstream and downstream are now combined compared to prev version
if __name__ == "__main__":
    in_faa, in_marker, in_cazyme, in_pfam, outfi = sys.argv[1:6]
    start = time.perf_counter()

    genome = os.path.basename(in_faa).split(".")[0]

    relevant_pfams = ["SusD", "TonB", "CarboxypepD_reg", "TPR", "Glyco_hydro", "HTH_18", "SASA", "HisKA", "AraC", "DUF"] #this checks for extra susC/susD pairs, DUF (domain of unknown function) and regulatory genes and extends the search as these are often found inside PULs

    with open(in_marker) as f: marker_lines = f.readlines()
    with open(in_faa) as f: faa_lines = f.readlines()
    with open(in_cazyme) as f: cazyme_lines = f.readlines()
    with open(in_pfam) as f: pfam_lines = f.readlines()

    parsed_faa = parse_faa_headers(faa_lines)
    cazy_dict = index_annotation_lines(cazyme_lines, min_fields=4, field_index=3)
    pfam_dict = index_annotation_lines(pfam_lines, min_fields=5, field_index=4)

    with open(outfi, 'w') as ou:
        pul_counter = 1
        for marker_line in marker_lines:
            parts = marker_line.strip().split('\t')
            if len(parts) < 5:
                continue

            try:
                M_susD_ID = parts[1]
                M_tonB_ID = parts[2]
                marker_start = int(parts[3])
                marker_end = int(parts[4])
            except ValueError:
                continue

            found = False
            for i, (gene_id, start, end, header) in enumerate(parsed_faa):
                if gene_id == M_susD_ID or gene_id == M_tonB_ID:
                    found = True
                    down = find_nearby_genes(parsed_faa, i, "down", 500, 4, cazy_dict, pfam_dict, relevant_pfams)
                    up = find_nearby_genes(parsed_faa, i, "up", 500, 4, cazy_dict, pfam_dict, relevant_pfams)

                    ou.write(f"Genome: {genome}_PUL_{pul_counter}, Marker: {gene_id}, {marker_start}, {marker_end}\n")
                    ou.write("  Downstream:\n")
                    for d in down:
                        ou.write(f"\t{d}\n")
                    ou.write("  Upstream:\n")
                    for u in up:
                        ou.write(f"\t{u}\n")
                    ou.write("\n")

                    pul_counter += 1
                    break

            if not found:
                log(f"Marker {M_susD_ID} not found in FAA file.")

    end = time.perf_counter()
    log(f"Completed pul mapping annotation in {end - start:.2f} seconds.")
