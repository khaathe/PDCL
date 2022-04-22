#!/usr/bin/python3 

# Command line to retrieve pathway informations from the Kegg Database
# curl http://rest.kegg.jp/get/hsa00010

import os
import sys
import requests
import re

# Path to kegg2go file to write
kegg2go_path = "/home/spinicck/PhD/tmp/kegg2go.txt"

# Kegg db url
base_kegg_url = "http://rest.kegg.jp"

# Wrapper for get request to the kegg db
def kegg_get(operation, args):
    url = "{0}/{1}/{2}".format(base_kegg_url, operation, '/'.join(args) )
    r = requests.get(url)
    return(r)
    
# Sent a request and parse the response
# Get the list of all the human pathways in db
# Return a dict with pathway id as key and pathway desc
# as value
def get_list_human_pathways():
    human_pathways = dict()
    r = kegg_get("list", ["pathway", "hsa"])
    for line in progress_bar(r.text.strip().splitlines()):
        (pathway_id, pathway_desc) = line.split("\t")
        human_pathways[pathway_id] = pathway_desc
    return(human_pathways)
    
# Return GO ID corresponding to Kegg ID
def get_kegg_2_go(path_id):
    go_ids = list()
    currentState = "START"
    response = kegg_get("get", [path_id])
    for line in response.text.strip().splitlines():
        m = re.search("^(\w+)", line, flags=re.IGNORECASE)
        if (m != None) : currentState = m.group(1)
        if (currentState.upper() == "DBLINKS") :
            m = re.search("GO:\s*([\d\s]*)", line)
            if (m != None):
                go_ids.extend(m.group(1).split())
    return(go_ids)

# Write the kegg2go file at the specified path
def write_kegg2go(pathways, source_db = "KEGG", sep="\t"):
    with open(kegg2go_path, "w") as out:
        header = sep.join(["source_db", "source_id", "source_desc", "go_id", "go_desc"]).strip() + os.linesep
        out.write(header)
        for (path_id, path_desc) in progress_bar( pathways.items() ):
            for go_id in get_kegg_2_go(path_id):
                str = sep.join([source_db, path_id, path_desc, f'GO:{go_id}', "NA"]).strip() + os.linesep
                out.write(str)

# Generator that loop over an iterable and
# print progress as we loop through
# Show a progress bar like : 
# Progress: [##############################] - 100.00%
def progress_bar(it):
    count = len(it)
    out = sys.stdout
    # Print the progress of a task to
    # keep track of a time consumming task
    def print_progress(progress, length = 30, progress_char = '#', empty_char = ' '):
        progress_str = "Progress: [{:" + empty_char + "<" + f'{length}' + "s}] - {:.2%}\r"
        out.write(progress_str.format(progress_char * int(progress * length), progress))
        out.flush()
    print_progress(0.0)
    for i, item in enumerate(it):
        yield item
        print_progress(i/count)
    print_progress(1.0)
    out.write(os.linesep)
    out.flush()



def main():
    print("##### Download Kegg2Go file #####")
    print("Get Kegg Human Pathways ...")
    human_pathways = get_list_human_pathways()
    print("Done !")
    print(f'Writting {kegg2go_path} ...')
    write_kegg2go(human_pathways)
    print("Done !")
    
main()
