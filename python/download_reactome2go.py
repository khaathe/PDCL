#!/usr/bin/python3 

# Command line to retrieve pathway informations from the Reactome Database
# curl -X GET --header 'Accept: application/json' https://reactome.org/ContentService/data/query/enhanced/R-HSA-70171

import requests
import warnings
import re

# Path to reactome2go file to write
reactome2go_path = "/home/spinicck/PhD/tmp/reactome2go.txt"

# reactome db url
base_reactome_url = "https://reactome.org/ContentService/data"

# Wrapper for get request to the reactome db
def reactome_get(operations, args):
    url = "{0}/{1}/{2}".format(base_reactome_url, '/'.join(operations), '/'.join(args) )
    r = requests.get(url)
    return(r)
    
# Sent a request and parse the response
# Get the list of all the human pathways in db
# Return a dict with pathway id as key and pathway desc
# as value
def get_list_human_pathways():
    human_pathways = dict()
    # TODO : Implement
    
# Write the reactome2go file at the specified path
def write_reactome2go():
    # TODO: implement
    # with open(gmt_path, "w") as out:
    #     for k, v in pathway_2_genes.items():
    #         out.write("{}\t{}\t{}\n".format(k, pathways[k], '\t'.join(v) ))
    pass

def main():
    human_pathways = get_list_human_pathways()
    write_reactome2go()
    
main()
