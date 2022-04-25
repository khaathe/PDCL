#!/usr/bin/python3 

# Command line to retrieve pathway informations from the Reactome Database
# curl -X GET --header 'Accept: application/json' https://reactome.org/ContentService/data/query/enhanced/R-HSA-70171

import os
import sys
import requests
import json

# Path to reactome2go file to write
reactome2go_path = "/home/spinicck/PhD/tmp/reactome2go.txt"

# reactome db url
base_reactome_url = "https://reactome.org/ContentService/data"

# Path to gmt file to collect Reactome IDs
gmt_path = "/home/spinicck/PhD/PDCL/Data/gene-set/ReactomePathways.gmt"

# Reactome maximum number of IDs per POST request
reactome_max_nb_ids = 20

# Wrapper for get request to the reactome db
def reactome_get(resources):
    url = "{0}/{1}".format(base_reactome_url, '/'.join(resources) )
    response = requests.get(url)
    response.raise_for_status()
    return(response)

# Wrapper for POST request to the reactome db
def reactome_post(resources, data, headers = dict() ):
    url = "{0}/{1}".format(base_reactome_url, "/".join(resources))
    custom_headers = { "Content-Type" : "text/plain"}
    custom_headers.update(headers)
    response = requests.post(url = url, data = data, headers=custom_headers)
    response.raise_for_status()
    return(response)
    
# Sent a request and parse the response
# Get the list of all the human pathways in db
# Return a dict with pathway id as key and pathway desc
# as value
def get_list_human_pathways(id_col = 0):
    human_pathways = list()
    with open(gmt_path, "r") as fin:
        [human_pathways.append(line.strip().split("\t")[id_col]) for line in fin]
    return(human_pathways)
    
# Write the reactome2go file at the specified path
def write_reactome2go(pathway_ids, source_db = "Reactome", sep = "\t"):
    no_go_id_pathways = 0 # Number of pathways without mapped go terms
    resources = ["query", "ids"]
    i = 0 # Incrementer to keep track of progress
    nb_pathways = len(pathway_ids) # Number of pathway to compute progress
    with open(reactome2go_path, "w") as fout:
        header = sep.join(["source_db", "source_id", "source_desc", "go_id", "go_desc"]).strip() + os.linesep
        fout.write(header)
        # Loop while the list of pathways IDs isn't empty
        while(len(pathway_ids) > 0):
            # Retrieve information for the N-th first pathways.
            # Reactome cannot process more than a certain number of pathways
            # Hence we need to process the list by parts.
            response = reactome_post(resources, ", ".join(pathway_ids[:reactome_max_nb_ids]))
            pathways = response.json()
            for pathway in pathways:
                # Check if there is a GO BP associated
                if "goBiologicalProcess" in pathway:
                    go_object = pathway["goBiologicalProcess"]
                    path_id = pathway["stId"]
                    path_desc = pathway["displayName"]
                    go_id = f'{go_object["databaseName"]}:{go_object["accession"]}'
                    go_desc = go_object["displayName"]
                    # Write the line
                    line = sep.join([source_db, path_id, path_desc, go_id, go_desc]) + os.linesep
                    fout.write(line)
                else:
                    no_go_id_pathways = no_go_id_pathways + 1
                i = i + 1
                print_progress(i/nb_pathways)
            # Remove the pathways we just mapped from the list
            del pathway_ids[:nb_pathways]
    print_progress(1.0)
    print()
    print(f'Pathway with no GO Term mapped : {no_go_id_pathways}')

# Generator that loop over an iterable and
# print progress as we loop through
# Show a progress bar like : 
# Progress: [##############################] - 100.00%
def progress_bar(it):
    count = len(it)
    out = sys.stdout
    print_progress(0.0)
    for i, item in enumerate(it):
        yield item
        print_progress(i/count)
    print_progress(1.0)
    out.write(os.linesep)
    out.flush()

# Print the progress of a task to
# keep track of a time consumming task
def print_progress(progress, length = 30, progress_char = '#', empty_char = ' ', out = sys.stdout):
    progress_str = "Progress: [{:" + empty_char + "<" + f'{length}' + "s}] - {:.2%}\r"
    out.write(progress_str.format(progress_char * int(progress * length), progress))
    out.flush()

def main():
    human_pathways = get_list_human_pathways(1)
    write_reactome2go(human_pathways)

main()
