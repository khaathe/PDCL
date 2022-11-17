#!/usr/bin/python3 

# Use this script to download the gene from the Kegg database.
# To filter out small and large pathways (nb of genes <15 or >1000)
# Use the awk program with a command such as :
# awk -F "\t" '{if(NF-2 > 15 || NF-2<1000) print}' kegg.gmt > kegg_no_less_15_no_more_1000.gmt 
import requests
import warnings
import re

# Path to gmt file to write
kegg_gmt_path = "/home/spinicck/PhD/tmp/kegg.gmt"

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
    
# Sent a request and parse the response
# Get the list of all the human genes in db
# Return a dict with gene id as key and as 
# value a dict](symbols=list(), desc=string)
def get_list_human_genes():
    human_genes = dict()
    r = kegg_get("list", ["hsa"])
    no_symbol_genes = 0
    for line in r.text.splitlines():
        (gene_id, gene_info) = line.split("\t")
        match = re.search(";", gene_info)
        # Check that the gene has symbols
        if (match != None):
            (gene_symbols_str, gene_desc) = gene_info.strip().split(";")
            gene_symbols = re.split(",\s", gene_symbols_str)
            human_genes[gene_id] = dict(symbols = gene_symbols, desc = gene_desc.strip() )
        else:
            human_genes[gene_id] = dict(symbols = None, desc = gene_desc.strip() )
            no_symbol_genes = no_symbol_genes + 1
    if (no_symbol_genes>0):
        warnings.warn("Number of genes without symbols : {}".format(no_symbol_genes) )
    return(human_genes)

# Sent a request and parse the response
# Get the list of all the human genes in db
# Return a dict[pathway_id] = list(gene_id)
def get_link_pathway_2_genes():
    link_pathway_2_genes = dict()
    r = kegg_get("link", ["hsa", "pathway"])
    for line in r.text.strip().split("\n"):
        (pathway_id, gene_id) = line.split("\t")
        gene_list = link_pathway_2_genes.get(pathway_id, list() )
        gene_list.append(gene_id)
        link_pathway_2_genes.update({pathway_id : gene_list})
    return(link_pathway_2_genes)

# Convert a list of gene id ot their corresponding
# Symbols.
# Return a list of symbols
def convert_gene_id_to_symbol(gene_info, gene_ids):
    gene_symbols = [ gene_info[id]["symbols"][0] if (gene_info[id]["symbols"] != None) else id  for id in gene_ids ]
    return(gene_symbols)

# Write the gmt file at the specified path
def write_gmt(gmt_path, pathways, pathway_2_genes):
    with open(gmt_path, "w") as out:
        for k, v in pathway_2_genes.items():
            out.write("{}\t{}\t{}\n".format(k, pathways[k], '\t'.join(v) ))

def main():
    human_pathways = get_list_human_pathways()
    human_genes = get_list_human_genes()
    link_pathway_2_genes = get_link_pathway_2_genes()
    for (path_id, gene_list) in link_pathway_2_genes.items():
        gene_symbols = convert_gene_id_to_symbol(human_genes, gene_list)
        link_pathway_2_genes[path_id] = gene_symbols
    write_gmt(kegg_gmt_path, human_pathways, link_pathway_2_genes)
    
    

main()
