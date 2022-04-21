#!/usr/bin/python3 

# Command line to retrieve pathway informations from the Kegg Database
# curl http://rest.kegg.jp/get/hsa00010

import requests
import warnings
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
    for line in r.text.strip().split("\n"):
        (pathway_id, pathway_desc) = line.split("\t")
        human_pathways[pathway_id] = pathway_desc
    return(human_pathways)
    
# Map each kegg pathways to a list
# of GO Term ID
def link_kegg_2_go(pathways):
    # TODO: Test
    kegg2go = dict()
    for (path_id) in pathways.keys():
        currentState = "START"
        response = kegg_get("get", [path_id])
        for line in response.text.strip().splitlines():
            m = re.search("^(\w+)", line, flags=re.IGNORECASE)
            if (m != None) : currentState = m.group(1)
            if (currentState.upper() == "DBLINKS") :
                m = re.search("GO:\s(\.*)", line)
                if (m != None):
                    go_ids = m.group(1).split()
                    kegg2go[path_id] = go_ids
    return(kegg2go)

# Write the kegg2go file at the specified path
def write_kegg2go():
    # TODO: implement
    # with open(gmt_path, "w") as out:
    #     for k, v in pathway_2_genes.items():
    #         out.write("{}\t{}\t{}\n".format(k, pathways[k], '\t'.join(v) ))
    pass

def main():
    human_pathways = get_list_human_pathways()
    write_kegg2go()
    
main()
