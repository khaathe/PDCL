# Thanks to Jorge Bretonnes-Santamarina who created the algorithm to assign categories, defined
# the regular expression in word_list and underesired_words.

######## Load Library needed
library(ggplot2)
library(ggh4x)
library(tidyr)
library(purrr)

# Made by Jorge Bretonnes-Santamarina
#Define undesired words
undesired_words <- c("developmental","valproic acid","disease","\\breproduction\\b","^(?=.*stem)(?!.*recept).*$","hemostasis","congenital","mouse","\\srat(\\b|s\\b)","ataxia","^(?=.*(\\b|ch)l1(\\b|cam))(?!.*(adhesion|signal)).*$",
                     "\\bcrcp\\b","\\bvpr\\b","fibrin clot","\\bclotting\\b","bacteri","\\bhost","diabete","defective","parkinson","microb","alzheimer","salmonella",
                     "\\bhiv","influenza","infect(ed|ion)","vir(al|us|ion)","pathogen","invasion","^(?=.*tumor)(?!.*(suppress|factor)).*$","shigell","nemia","toxoplasm","legionel","\\bmusc(ular|le)",
                     "myo(phagy|card)","ltr circle","injury","ischemia","cancer","leukemia","keratinization","\\bmhc\\b","prion","\\blupus\\b","dementia","rodent","inflammat(ory|ion)","\\btruncat",
                     "neurofibrillary","vasospasm","\\bepilepti","follicle","^(?=.*oma\\b)(?!.*retinoblastom).*$","neoplasm","antigen","(\\b(mast|t|b|beta)(-cell| cell|-help|reg|-reg))","endometriosis","(angio|myo)genesis",
                     "\\bv(\\(d\\)|d)j\\b","polarity","gastrulation","coagul(um|ation)","parasit(e|ism)","\\bexrna\\b","crcp cell","placenta","^(?=.*immun)(?!.*(apopto|metabo)).*$","axon","^(?=.*toxi(c|n))(?!.*metabol).*$","^(?=.*blast)(?!.*retinoblastom).*$",
                     "hypothe(sis|tical)","retinal","dendritic","\\bocular\\b","melan(osome|in)","\\beph-ephrin\\b","papilla","neurotransmit","\\bhmc-[1-9]","cystic fibrosis","mucopolysaccharidoses","\\bneur(o|al|e)","sperm",
                     "^(?=.*syndrome)(?!.*molecule).*$","^(?=.*card)(?!.*cardiolipin).*$","botulinum","hypertroph","fertil","leishmania","brain","measles","alcoholism",
                     "pertussis","amoeb","sclero","killer","\\bnk\\b","^(?=.*drug)(?!.*metabol).*$","malaria","thyroidism","macroph","^(?=.*hepato)(?!.*receptor).*$","granulo","ophil(\\b|s\\b)","pathy","disorder","\\bovari","pecia","pathie","trophy",
                     "asthma","nephro","neutro","^(?=.*synap(tic|s))(?!.*(\\bmeio|\\bdna\\b|strand)).*$","\\bgout\\b","itis\\b","ectasia","obesity","fever","plasia\\b","ilia\\b","rhea\\b","\\bhpv","vascular","vitiligo","\\bacne\\b",
                     "(keratin|seb|my|chondr|adip|leuk)ocyte","cataract","\\bear\\b","intest","algia\\b","^(?=.*\\brev)(?=.*prote).*$","\\btherap","theli(um|a)\\b","spinal","dorsal","^(?=.*potassium)(?=.*channel).*$",
                     "ventral","chemotaxis","penile","wound","\\bheal","heart","\\bskin\\b","apical","^(?=.*oderm)(?!.*(\\begf(r\\b| receptor)|growth factor)).*$","myel(in|oyd)","(\\boo|lympho)cyte(\\b|s)",
                     "\\bcd4\\b","cd8","saliva","mut(ant|tion|ated)","physiolog","\\blearn","epilepsy","\\bgland\\b","proximal","periodism","\\bcomplement\\b","amyloid","hem(opo|ato)","mesenchymal","leads to",
                     "fibrino","\\bampa\\b","density lipoprotein","nep/ns2","\\bdectin","\\btube\\b","isotype","hodgkin and reed-sternberg","insulin resistance","hypertension",
                     "dysfunction","potentiation","\\bduct\\b","platelet","\\bph\\b","\\blectin\\b","\\bspike\\b","ossification","abacavir","phenylketonuria")

# Made by Jorge Bretonnes-Santamarina
#Define list of regexp that determine how pathways will be assigned into Reactome categories
word_list <- list(
  #Autophagy
  Autophagy = c("phag"),
  #Cell cycle
  `Cell Cycle` = c("cell (cycle|division)","^(?=.*\\bgrowth\\b)(?!.*(factor|signal)).*$","(?=.*\\bcdk[1,2,4,6]\\b)","\\bcyclin\\b","^(?=.*spindle)(?=.*microtubul).*$","checkpoint","\\bflap\\b","chromatid","meio","mito(tic|sis)","^(?=.*phase)(?!.*(apopto|attenuation|strand|homologous)).*$","\\bg1\\b","\\bg2\\b","\\bt-circle","retinoblastom","\\brb1\\b","^(?=.*\\be2f(\\b|[1-9]))(?!.*transcript).*$","\\bapc/c\\b",
                   "centriole","^(?=.*chromosome)(?!.*(protein|localiz|organiz|\\sx\\s)).*$","cleavage furrow","\\bnek2","^(?=.*(telomer|centrosome|kinetochore))(?!.*(protein|localiz)).*$","telomerase","cytokinesis","nuclear division","^(?=.*nuclear)(?=.*(lamina|envelope|pore))(?=.*(sealing|formati|polymeri|breakdown|assembly)).*$","^(?=.*tumor)(?=.*suppressor)(?!.*(signal|activat)).*$"),
  #Cell-cell communication
  `Cell-Cell communication` = c("^(?=.*junction)(?!.*(exon|signal)).*$","cell-cell","^(?=.*(cell|focal))(?=.*(communication|spreading|adhesion))(?!.*(signal|apopto|matrix)).*$","^(?=.*(nephrin|\\bnetrin|\\bnectin))(?!.*signal).*$"),
  #Cellular response to external stimuli
  `Cellular responses to external stimuli` = c("senescence","\\bhsf1\\b","deprivation","longevity","^(?=.*(reactive oxygen|\\bros\\b|zink|attenuation|cadmium|superoxide|response(|s) (to|of)|\\buv\\b))(?!.*signal).*$","detoxificat","(hyp|norm)oxia","stress","^(?=.*heat)(?!.*(\\batp\\b|electron|respiratory)).*$","^(?=.*damage)(?=.*response)(?!.*(protein|ylation)).*$","stimuli"),
  #Chromatin organization
  `Chromatin organization` = c("^(?=.*histone)(?!.*rna(\\b|s)).*$","^(?=.*chromatin)(?!.*localiz)(?!.*protein).*$","polycomb","nucleosome","^(?=.*epigenetic)(?=.*(alteration|instability)).*$"),
  #Circadian
  `Circadian Clock` = c("circadian","clock"),
  #DNA repair
  `DNA Repair`= c("^(?=.*repair)(?!.*(transcript|protein)).*$","ap site","\\sap\\s","\\bhdr\\b","translesion","path replacement","excision","incision","^(?=.*joining)(?=.*end).*$","abasic","mismatch","^(?=.*(dna\\b|strand))(?=.*(ligation|break|annealing|pairing|exchange))(?!.*(transcript|ylation|protein|replicat)).*$",
                  "^(?=.*(dna\\b|strand|purine))(?=.*damage)(?!.*(transcript|protein)).*$","depurination"),
  #"strand (break|annealing|exchange|damage)"
  #DNA replication
  `DNA Replication` = c("^(?=.*(\\bdna))(?!.*(\\brna|transcript|metabol|ylation|repression|stress))(?=.*(re(|du)plicat|elongation|synthesis|polymerase|winding)).*$","(pre|post)-replicati","recombination","\\bfork\\b","^(?=.*strand)(?!.*\\brna).*$","d-loop",
                        "licensing factor","\\bhdr\\b","^(?=.*(polymerase|origin))(?=.*switching).*$","\\borc\\b"),
  #Extracellular matrix organisation
  `Extracellular matrix organization` = c("^(?=.*(integrin|elastin|laminin|cadherin|ecan\\b|fibronectin|collagen))(?!.*signal).*$","fibril","^(?=.*(extracellular|\\bcell))(?=.*matrix).*$","\\becm\\b"),
  #Transcription
  `Gene expression (Transcription)` = c("^(?=.*(expression|silencing|transcription))(?!.*(signal|kinase|translat|factor)).*$","^(?=.*transcription (co-|co|)factor)(?!.*(ylation\\b|activat|kinase)).*$","^(?=.*elongation)(?!.*(translat|protein|process|peptide|chain|ribosom|fatty|golgi|glycan)).*$","promoter","rna pol(\\b|ymerase)","rna biosynthe",
                                        "\\b(pi|mi|si|small )rna(s\\b|\\b)","^(?=.*methylat)(?=.*(dna|cytosine)).*$","interference","enhancer","^(?=.*\\brna\\b)(?=.*regulat)(?!.*(metaboli|\\bbind|(ex|im|trans)port|splicing)).*$","^(?=.*\\bgene(\\b|s))(?=.*regulat)(?!.*splicing).*$",
                                        "^(?=.*(\\bgene(\\b|s)|\\sx\\s))(?=.*activat).*$"),
  #Metabolism
  `Metabolism` = c("^(?=.*(metabol|anabol|catabol|synthe|modificat|\\bchange|cataly|formation|transfer|degradation|isomer(ase|i)))(?!.*(rna(\\s|$)|crosslink|actin|destruction complex|filament|proteasome|ubiquitin|translat|\\ber\\b|endoplasmic|\\berad\\b|golgi|\\bprotein|peptid(e|yl)|subunit|\\b[0-9]{2}s\\b)).*$","\\btca\\b","tricarboxylic","electron",
                   "^(?=.*((fatty |amino )acid|glucose|lip(id|ase)|glycogen|mevalonate|carbohy|pyruvate|methionine|p(uri|yri)midine|\\burea\\b|(pent|hex)ose|itol\\b|sterol|nicotin|sulfur|\\bgaba\\b|glutathione|folate|\\batp|\\badp))(?!.*(signal|recept|(trans|im|ex)port|\\bprotein)).*$"
                   ,"(gluconeo|lipo|keto|thermo|adipo)genesis","glycoly(sis|tic)","oxidat","salvag(e|ing)","^(?=.*regulat)(?=.*activity)(?=.*(\\base|ylation))(?!.*(prot|kinase|signal)).*$","nutrient",
                   "breakdown","^(?=.*hydro(lys|lase|genase))(?!.*(ribosom|rna\\b)).*$","conversion","^(?=.*((acetyl|methyl|phosphoryl|alkyl|glucuronid|sulfon)ation))(?!.*(rna\\b|signal|cascade|kinase|activat|pathway|protein|peptid(e|yl)|activity|target|\\bbad\\b)).*$",
                   "^(?=.*cycle)(?!.*(calnexin|calreticulin)).*$","cyclase","xenobiotic","shunt","\\benos\\b","^(?=.*substrate)(?!.*(foldin|\\bcct\\b|\\btric\\b)).*$","imprinting",
                   "^(?=.*insulin)(?=.*(secret|metabol|regulat|process))(?!.*(signal|recept)).*$","respiration","carboxylase"),
  #Metabolism of proteins
  `Metabolism of proteins` = c("foldin","\\bcct\\b","\\btric\\b","\\berad\\b","unfold","n-glycan","calnexin","calreticulin","chaperone","\\btranslation\\b","^(?=.*translation)(?=.*(initiation|elongation|termination|regulat|recycling|silencing)).*$","peptide chain","ribosom","(?=.*(\\b[0-9]{2}s\\b))",
                               "ubiquitin","proteas(e|ome)","proteinase","^(?=.*proteoly)(?!.*signal).*$","^(?=.*(degradation|ylation))(?!.*(rna\\b|phosphoryl|destruction complex|receptor)).*$","sumoylat","\\beif[1-6]",
                               "^(?=.*(protein|\\bpepti|chain))(?=.*(process|metaboli|linkage|binding|modeling|synthesis|modif|stabil|merization|phosphoryl))(?!.*(rna\\b|nucleotide|signal|transport)).*$","citrullin","cytochrome","gpi anchor","^(?=.*removal)(?=.*amino).*$",
                               "insulinotropic","quality control","cardiolipin","\\bcl\\b","^(?=.*(golgi|\\ber\\b|endoplasmic))(?=.*(modific|elongation|trimming)).*$"),
  #"^(?=.*ribosom)(?!.*((ex|im|trans)port|nucleus)).*$"
  #Metabolism of RNA
  `Metabolism of RNA` = c("exon","intron","editing","nonsense","snrnp","non-coding","splic(eosome\\b|ing\\b)","microprocessor","^(?=.*trimming)(?!.*rna\\b).*$",
                          "\\b(nc|sn|sno|r|t|m|pre-m)rna(s\\b|\\b)","^(?=.*(rna\\b|transcr))(?=.*(processing|modif|hydrol|(meta|ana|cata)boli|binding|cap(ping\\b|\\b)|degradat|decay|cleavage|unwind|secondary struct|(ex|im|trans)port|stabilizat|loading|ylation)).*$","(peptide|protein) cross","^(?=.*\\brna\\b)(?=.*surveillance).*$"),
  #Organelle biogenesis and maintenance
  `Organelle biogenesis and maintenance` = c("^(?=.*(\\bactin\\b|microtubule|\\bactomyosin\\b|\\bmyosin|\\b(bio|)genesis))(?!.*(vesicle|ligand|transport|protein|osome\\b)).*$","^(peroxisome)$","^(lysosome)$","basal body","^(?=.*organelle)(?!.*protein).*$",
                                             "^(?=.*(mito|cytoskelet|vacuole|nucle|golgi|projection|ruffle|cilium|chromosome|endoplasmic|ribbon|membrane|lysosome|peroxi|component|fiber|lamellipod|multivesic|microvillus))(?=.*(assembl|organi|locali|format|migrat|fission|fusion|division|anchor|inherit|invagination|maintenance|tubulation))(?!.*(vesicle|protein|destruction complex)).*$",
                                             "furrow","spindle","^(?=.*(cili(um|ary)|flagel))(?!.*transport).*$","^(?=.*assembl)(?!.*destruction complex).*$","geometric change","cell size"),
  #Apoptosis
  `Programmed Cell Death`= c("apopto","caspase","^(?=.*cell)(?=.*death).*$","necro(s|pto)","inflammasome","anoikis","ferroptosis","\\b(c-|)flip(|-)","\\bbad\\b","^(?=.*lamina)(?=.*cleavage).*$"),
  #Protein localization
  `Protein localization` = c("^(?=.*(protein|\\bran\\b))(?=.*(\\ber\\b|mitoch|peroxis|cili(um|ary)|nucle|cytoskelet|reticulum|nucleus|lysosome|vacuole|membrane)).*$","dynein",
                             "^(?=.*(protein|peptide|subunit))(?=.*(import|export|exit|traffic|entry|locali|sequest|target|transport))(?!.*vesic).*$","^(?=.*receptor)(?=.*cluster).*$"),
  #vesicle mediated transport)
  `Vesicle-mediated transport` = c("\\brabgaps\\b","(endo|exo)som","vesic(ul|le)","\\bcopi","golgi","clathrin","^(?=.*(membrane|rab))(?=.*traffic).*$","(endo|exo|pino|trans)cyt","^(?=.*scaveng)(?=.*receptor).*$","kinesin","\\bcargo\\b","adp-ribosyl","\\bglut[1-9]"),
  #Transport of small molecules
  `Transport of small molecules` = c("clearance","\\babc","^(?=.*((im|ex|trans)port|exit|entry|locali|sequest|exchange|release|translocat|\\bbind|(disso|asso)ciat))(?!.*(signal|activit)).*$","\\bslc","aquapo","heme",
                                     "^(?=.*plasma)(?=.*lipoprotein).*$","\\biron\\b","(?=.*(\\bfe[0-9]\\b))","^(?=.*(stimul|\\btrp))(?=.*channel).*$","shuttle","\\bion\\b","\\bsmdt1\\b"),
  #"^(?=.*target)(?!.*stream).*$"
  #Signal transduction
  `Signal Transduction` = c("signal","phosphoryl","destruction complex","target","kinase","\\bmtor","\\bp53\\b","\\beif2\\b","\\begf","\\bfgf","\\ntrk","\\bpdgf","\\bvegf","\\bpten\\b","\\bras\\b","\\berbb","\\bgpcr","\\bnotch",
                            "\\bwnt","\\bhippo","\\bhedge","\\bmapk","\\bpip3","rho","\\bmtor","cascade","trigger","\\bon\\b","\\boff\\b","transduction","\\smet\\s",
                            "acti","inhib","recycl","pathway","effect","signal","events","network","\\bnicd\\b","\\bnotch","\\bikk\\b","nod-like","\\bnlr\\b",
                            "^(?=.*(nuclear|secretin|estrogen))(?=.*receptor).*$","downstream","g-protein","regulate","motility",
                            "^(?=.*insulin)(?=.*(recept|stimulus|process)).*$","insulin-like","\\bp38","oncostatin")
)

# Generic function turn a list into a data.frame.
# Used by the collapse.gem.list function.
map.list.to.df <- function(x, var_name){
  map2_dfr(x, names(x), function(el, el_name){
    el[,var_name] = el_name
    el
  })
}

method_vs_method <- function(collapsed_by_method, m1_name, m2_name){
  m1_df <- collapsed_by_method %>%
    filter(method == m1_name)
  m2_df <- collapsed_by_method %>%
    filter(method == m2_name)
  m_vs_m <- m1_df %>% 
    full_join(m2_df, by = c("pathway_id", "pdcl", "db", "category")) %>% 
    dplyr::rename(p_value_m1 = p_value.x, p_value_m2 = p_value.y, fdr_m1 = fdr.x, fdr_m2 = fdr.y) %>%
    select(pathway_id, pdcl, db, p_value_m1, fdr_m1, p_value_m2, fdr_m2, category) %>%
    mutate(method_1 = m1_name, method_2 = m2_name)
  m_vs_m
}

# Collapse a list of GEM results into one data.frame with :
# - all the column from a GEM file (https://enrichmentmap.readthedocs.io/en/latest/FileFormats.html#generic-results-files)
# - db : Source of pathway informations (ex: Reactome)
# - pdcl : ID of the pdcl used for the enrichment analysis
# - method : the enrichment method that has generated this result (ex: GSEA, G:Profiler) 
collapse.enrichment <- function(all_enrichment){
  all_enrichment_collapsed <- all_enrichment
  
  all_enrichment_collapsed <-  map_depth(.x = all_enrichment_collapsed, .depth = 2, .f = function(x){
    map.list.to.df(x, "sample")
  })
  
  all_enrichment_collapsed <-  map_depth(.x = all_enrichment_collapsed, .depth = 1, .f = function(x){
    map.list.to.df(x, "db")
  })
  
  all_enrichment_collapsed <-  map_depth(.x = all_enrichment_collapsed, .depth = 0, .f = function(x){
    imap(all_enrichment_collapsed, function(x,y){
      if( all(isEmpty(x))){
        stop(y, " result is empty.")
      }
    })
    map.list.to.df(x, "method")
  })
  
  all_enrichment_collapsed$sample <- as.factor(all_enrichment_collapsed$sample)
  all_enrichment_collapsed$db <- as.factor(all_enrichment_collapsed$db)
  all_enrichment_collapsed$method <- as.factor(all_enrichment_collapsed$method)
  
  all_enrichment_collapsed
}

# Create a data.frame where one line contain the P-Value and FDR result for
# 2 different enrichment method (ex: G:Profiler vs GSEA). Usefull to compare
# the result of two methods. The collapsed_by_method is a data.frame given by the function
# collapse.gem.list .
collapse.method.vs.method <- function(collapsed_by_method){
  method_names <- unique(collapsed_by_method$method)
  collapsed_by_pathways <- data.frame()
  for (m1 in head(method_names, n=-1)){
    index <- which(method_names == m1)
    for (m2 in tail(method_names, n = -index)){
      m_vs_m <- method_vs_method(collapsed_by_method, m1, m2)
      collapsed_by_pathways <- rbind(collapsed_by_pathways, m_vs_m)
    }
  }
  collapsed_by_pathways
}

# Find the pathway common between 2 method in a data.frame for each threshold specified
# in the threshold_vect arg. The collapsed_method_vs_method is a data.frame given by the function
# collapse.method.vs.method.
find.method.common.pathways <- function(collapsed_method_vs_method, threshold_vect){
  common <- collapsed_method_vs_method
  for ( threshold in threshold_vect){
    var_name <- paste0("alpha_", threshold)
    is_common <- ( common$fdr_m1<threshold & common$fdr_m2<threshold )
    is_common[is.na(is_common)] <- F
    is_common <- factor(
      is_common, 
      levels = c("TRUE", "FALSE"), 
      labels = c("Common Pathways", "Specific Pathways")
    )
    common[,var_name] <- is_common
  }
  col_to_keep <- c("pdcl", "pathway_id", "method_1", "method_2", paste0("alpha_", threshold_vect) )
  common %>% select(all_of(col_to_keep))
}

# Assign categories to pathways by performing a regular expression test using the regular 
# expression defined in word_list.
# This method is based on Jorge Bretonnes-Santamarina works.
assign.categories <- function(x){
  data <- x
  undesired_words <- paste0(undesired_words,collapse = "|")
  data$undesired <- grepl(undesired_words, data$pathway_id, perl = TRUE, ignore.case = T)
  data$category <- NA
  for (i in 1:length(word_list)){
    pattern <- paste0(word_list[[i]], collapse = "|")
    match_index <- grepl(pattern, data$pathway_id, perl = T, ignore.case = T)
    if( any(match_index) ){
      data$category[match_index & is.na(data$category)] <- names(word_list)[i]
    }
  }
  data
}

# Plot the FDR value of all pathways in the enrichment analysis
# results of each PDCL.
# Usefull to compare correlation between 2 method.
plot_method_vs_method <- function(collapsed_method_vs_method, m1, m2, threshold){
  data <- collapsed_method_vs_method %>% filter(method_1 == m1, method_2 == m2)
  alpha_col <- paste0("alpha_", threshold)
  data$is.common <- data[,alpha_col]
  ggplot(
    data = data, 
    aes(
      x = -log10(fdr_m1), 
      y = -log10(fdr_m2),
      colour = is.common
    )
  ) +
    geom_point(
      size = 0.5
    ) +
    scale_colour_manual(values = c("red", "grey30")) +
    facet_wrap(vars(pdcl)) +
    labs(
      title = paste0( m1, " against ", m2),
      subtitle = paste0("Threshold alpha : ", threshold),
      x = paste0( "-log10(fdr ",  m1, ")"),
      y = paste0( "-log10(fdr ",  m2, ")"),
      colour = ""
    ) +
    theme(
      plot.title = element_text(hjust = 0.5,face = "bold",color = "red"),
      plot.subtitle = element_text(hjust = 0.5,face = "italic",color = "blue"),
      axis.text.x = element_text(angle = 90)
    ) 
}

# Plot the number of enriched pathways (pathways with an associated FDR below
# a threshold) of all PDCL for each method present in the data.frame collapsed.by.method.
plot.count.pathways <- function(collapsed.by.method, threshold){
  data <- collapsed.by.method %>% filter(fdr<threshold)
  ggplot(
    data = data, 
    aes(x = pdcl)
  ) +
    geom_bar(fill = "steelblue") +
    labs(
      title = "Number of enriched pathways per PDCL for each enrichment method",
      subtitle = paste0("Threshold alpha : ", threshold),
      x = "PDCL",
      y = "Counts",
      fill = ""
    ) +
    facet_wrap(vars(method)) +
    theme(
      plot.title = element_text(hjust = 0.5,face = "bold",color = "red"),
      plot.subtitle = element_text(hjust = 0.5,face = "italic",color = "blue"),
      axis.text.x = element_text(angle = 90)
    )  +
    #We need to specify this. If not, unpopulated categories do not appear
    scale_x_discrete(drop = F)
}

# Plot the number of enriched pathways (pathways with an associated FDR below
# a threshold) of all PDCL for two methods specified by m1 and m2 args present in the 
# data.frame collapsed.by.method. Highlight the pathways commons between the 2 methods.
plot.count.pathways.m.vs.m <- function(collapsed.by.method, common.pathways, threshold, m1, m2){
  common.pathways <- common.pathways %>%
    filter(method_1 == m1, method_2 == m2)
  collapsed.by.method <- collapsed.by.method %>% 
    filter( method %in% c(m1, m2) & fdr < threshold ) %>%
    inner_join(common.pathways, by=c("pdcl", "pathway_id"))
  collapsed.by.method$is.common <- collapsed.by.method[,paste0("alpha_", threshold)]
  
  ggplot(
    data = collapsed.by.method, 
    aes(x = pdcl, fill = is.common )
  ) +
    geom_bar() +
    scale_fill_manual(values = c("red", "steelblue")) + 
    labs(
      title = "Number of enriched pathways per PDCL for each enrichment method",
      subtitle = paste0("Threshold alpha : ", threshold),
      x = "PDCL",
      y = "Counts",
      fill = ""
    ) +
    facet_wrap(vars(method)) +
    theme(
      plot.title = element_text(hjust = 0.5,face = "bold",color = "red"),
      plot.subtitle = element_text(hjust = 0.5,face = "italic",color = "blue"),
      axis.text.x = element_text(angle = 90)
    )  +
    #We need to specify this. If not, unpopulated categories do not appear
    scale_x_discrete(drop = F)
}


# Plot the number of enriched pathways (pathways with an associated FDR below
# a threshold) of all categories for each method present in the data.frame collapsed.by.method.
plot.count.categories <- function(collapsed.by.method, threshold){
  data <- collapsed.by.method %>% filter(fdr<threshold)
  ggplot(
    data,
    aes(x = category)
  ) +
    geom_bar(fill = "steelblue") +
    facet_wrap(~method) +
    labs(
      x = "Category name",
      y = "Counts",
      title = "Number of enriched pathways in a category for each enrichment method",
      subtitle = paste0("Threshold alpha : ", threshold)
    ) +
    theme(
      plot.title = element_text(hjust = 0.5,face = "bold",color = "red"),
      plot.subtitle = element_text(hjust = 0.5,face = "italic",color = "blue"),
      axis.text.x = element_text(angle = 90)
    ) +
    #We need to specify this. If not, unpopulated categories do not appear
    scale_x_discrete(drop = F)
}

# Plot the number of enriched pathways (pathways with an associated FDR below
# a threshold) of all categories for two methods specified by m1 and m2 args present in the 
# data.frame collapsed.by.method. Highlight the pathways commons between the 2 methods.
plot.count.categories.m.vs.m <- function(collapsed.by.method, common.pathways, threshold, m1, m2){
  common.pathways <- common.pathways %>%
    filter(method_1 == m1, method_2 == m2)
  collapsed.by.method <- collapsed.by.method %>% 
    filter( method %in% c(m1, m2) & fdr < threshold ) %>%
    inner_join(common.pathways, by=c("pdcl", "pathway_id"))
  collapsed.by.method$is.common <- collapsed.by.method[,paste0("alpha_", threshold)]
  ggplot(
    collapsed.by.method,
    aes(x = category, fill = is.common)
  ) +
    geom_bar() +
    scale_fill_manual(values = c("red", "steelblue")) +
    facet_wrap(~method) +
    labs(x = "Category name",
         y = "Counts",
         title = "Number of enriched pathways in a category for each enrichment method",
         subtitle = paste0("Threshold alpha : ", threshold),
         fill = ""
    ) +
    theme(
      plot.title = element_text(hjust = 0.5,face = "bold",color = "red"),
      plot.subtitle = element_text(hjust = 0.5,face = "italic",color = "blue"),
      axis.text.x = element_text(angle = 90)
    ) +
    #We need to specify this. If not, unpopulated categories do not appear
    scale_x_discrete(drop = F)
}

# Generate a vector of breaks for ggplot graph with n values
# inferior and superior a midpoint.
generate.breaks <- function(min, max, midpoint, n){
  left <- seq(from = min, to = midpoint, length.out = n ) 
  right <- seq(from = midpoint, to = max, length.out = n )[-1]
  c(left, right)
}

# Turn the numeric FDR col into a factor with correct labels. 
fdr.to.factor <- function(x, n){
  break.vect <- generate.breaks(0, 1, threshold, n)
  cut.vect <- cut(x$fdr, breaks = break.vect, labels = break.vect[-1])
  cut.vect <- factor(cut.vect, levels = rev(levels(cut.vect)))
  x$fdr <- cut.vect
  x
}

# Heatmap of the FDR values for pathways common in at least one PDCL between the 2
# method specified by m1 and m2 args. Common pathways have an FDR value below the
# threshold specified.
heatmap.pathways <- function(collapsed.by.method, common.pathways, threshold, m1, m2){
  common_col_name <- paste0("alpha_", threshold)
  common.pathways <- common.pathways %>%
    filter(method_1 == m1, method_2 == m2) %>%
    dplyr::rename(is.common = all_of(common_col_name), method = method_1) %>%
    select(pdcl, pathway_id, method, is.common )
  
  collapsed.by.method <- collapsed.by.method %>% 
    left_join(common.pathways, by=c("pdcl", "pathway_id", "method"))
  
  common.pathways <- common.pathways %>%  
    filter(is.common == "Common Pathways") %>%
    select(!method)
  pathway.ids <- unique(common.pathways$pathway_id)
  data <- collapsed.by.method %>% filter(pathway_id %in% pathway.ids)
  data <- fdr.to.factor(data, 5)
  
  category_df <- data %>% distinct(category, pathway_id) %>% arrange(pathway_id)
  
  plot_title <- "FDR value each pathways by PDCL and method "
  plot_subtitle <- paste0("Pathways are common between : ", m1, " and ", m2)
  heatmap_plot <- ggplot(
    data = data, 
    aes(x = pdcl, y = pathway_id, fill =  fdr)
  ) +
    geom_tile() +
    geom_tile(
      data = common.pathways,
      mapping = aes(x = pdcl, y = pathway_id),
      inherit.aes = F,
      colour = 'black',
      fill = NA,
      size = 0.5
    ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "PDCL",
      y = "Pathway",
      fill = "FDR"
    ) +
    facet_wrap(vars(method)) +
    theme(axis.text.x = element_text(angle = 90) ) +
    scale_fill_brewer(
      palette = "RdBu",
      na.value = "gray30",
      direction = -1
    ) +
    #We need to specify this. If not, unpopulated categories do not appear
    scale_x_discrete(drop = F) +
    # By default it is sorted in desc order, we reverse it so it is in asc order
    scale_y_discrete(limits=rev) + 
    guides(
      y.sec =  guide_axis_manual(
        labels = rev(category_df$category),
        title = "Category"
      )
    ) +
    theme(panel.background = element_rect(fill = "white"))
  heatmap_plot
}

# Heatmap of the number of pathways found enriched in all categories.
# Enriched pathways have a FDR value below the threshold specified  
heatmap.categories <- function(collapsed.by.method, threshold){
  count <- collapsed.by.method %>% 
    filter(fdr < threshold) %>%
    dplyr::count(method, pdcl, category, .drop = F)
  ggplot(
    data = count,
    aes(x = pdcl, y = category, fill = n)
  ) + 
    geom_raster() +
    labs(
      title = "Number of pathways enriched in each categories for each PDCL",
      x = "PDCL",
      y = "Categories",
      fill = "Count"
    ) +
    facet_wrap(~method) +
    scale_fill_viridis_c(direction = 1) +
    theme(axis.text.x = element_text(angle = 90) ) +
    # By default it is sorted in desc order, we reverse it so it is in asc order
    scale_y_discrete(limits=rev)
}

# Turn a list of rank for GSEA analysis into a data.frame object.
collapse.rnk <- function(rnk){
  rnk <- imap_dfr(rnk, .f = function(rank,pdcl){
    data.frame(gene = names(rank), rank = rank, pdcl = pdcl, row.names = NULL)
  })
  rnk$pdcl <- factor(rnk$pdcl)
  rnk
}

# Get the genes in the leadingEdge column from GEM result.
# Return a data.frame where one line correspond to one gene
# in a pathway for one PDCL. All informations from the original
# data.frame are kept.
get.genes.in.pathway <- function(collapsed.by.method){
  genes.in.common.pathways <- collapsed.by.method %>%
    mutate(genes = strsplit(genes, ",")) %>%
    unnest(genes) %>%
    dplyr::rename(gene = genes)
}

# Plot the genes present in the leadingEdge of enriched pathways
# for one PDCL specified by the pdcl.name arg. Genes are ranked using
# their abs(rank) value in decreasing order. Only the nb.max.gene first
# genes are displayed in the plot.
# The higher the abs(rank) is, the bigger the point in the plot.
# Colour of the point indicate if the gene is up/down-regulated.
# A warning is printed if there is no pathway below the threshold specified and
# and a NULL value is returned by the function.
plot.leading.edge <- function(rnk.collapsed, collapsed.by.method, pdcl.name, threshold, gsea.method.name = "GSEA", nb.max.genes=100, pathways = NULL){
  data <- get.genes.in.pathway(collapsed.by.method) %>% 
    filter(method == gsea.method.name, fdr < threshold, pdcl == pdcl.name)
  
  if ( !is.null(pathways) ){
    data <- filter(data, pathway_id %in% pathways)
  }
  
  if( nrow(data) == 0){
    warning("No pathway found enriched for : ", pdcl.name, " at alpha = ", threshold)
    return(NULL)
  }
  
  rnk <- rnk.collapsed %>% filter(pdcl == pdcl.name)
  
  data <- data %>% 
    inner_join(rnk, by=c("pdcl", "gene")) %>%
    mutate(abs_rank = abs(rank), dereg = sign(rank) )
  data$dereg <- factor(data$dereg, levels = c(1, -1), labels = c("upregulated", "downregulated"))
  
  nb_gene_pdcl <- n_distinct(data$gene)
  if ( nb_gene_pdcl > nb.max.genes ){
    keep <- data %>% 
      distinct(gene, rank) %>% 
      arrange(desc(abs(rank))) %>% 
      slice(1:nb.max.genes) %>% 
      pull(gene)
    data <- filter(data, gene %in% keep)
  }
  
  plot_title <- paste0("Top ", nb.max.genes , " genes in enriched pathways for ", pdcl.name)
  plot_subtitle <- paste0("Alpha : ", threshold, ", Total number of different genes in pathways leading edge : ", nb_gene_pdcl)
  plot_leading_edge <- ggplot(data, aes(gene, pathway_id, size = abs_rank, colour = dereg) ) +
    geom_point() + 
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Gene Symbol",
      y = "Pathways",
      size = "abs(rank)",
      colour = ""
    ) +
    scale_colour_manual( values = c("upregulated" = "red", "downregulated"="blue") ) +
    theme(axis.text.x = element_text(angle = 90) ) +
    scale_y_discrete(limits = rev)
  plot_leading_edge
}

# Same as the plot.leading.edge function, but only plot
# genes found for pathways common between the method specified by m1 and m2.
plot.common.pathway.leading.edge <- function(rnk.collapsed, collapsed.by.method, pdcl.name, threshold, all.common.pathways, m1="G:Profiler", m2="GSEA", gsea.method.name = m2, nb.max.genes=100){
  alpha_col <- paste0("alpha_", threshold)
  common_pathways <- all_common_pathways %>%
    filter( .data[[alpha_col]] == "Common Pathways", method_1 == m1, method_2 == m2, pdcl == pdcl.name) %>%
    distinct(pathway_id) %>%
    pull(pathway_id)
  if( length(common_pathways) == 0){
    warning("No common pathway found between : ", m1, " and ", m2, " for ", pdcl.name)
  }
  plot.leading.edge(rnk.collapsed, collapsed.by.method, pdcl.name, threshold, gsea.method.name, nb.max.genes, common_pathways)
}

# Volcano plot to compare the deregulation type (up or down)
# between DESeq2 and PENDA 
plot.volcano <- function(data){
  ggplot(
    data, 
    aes( log2FoldChange, -log10(padj), colour = penda_predict, alpha = penda_predict )
  ) +
    geom_point() +
    scale_colour_manual( values = c("1" = "red", "-1"="blue", "0" = "grey50") )  +
    scale_alpha_manual(values = c("1" = 0.8, "-1" = 0.8, "0" = 0.4) ) +
    labs(
      title = "Comparison DESeq2 fold change vs PENDA Result",
      colour = "Penda Prediction", 
      alpha = "Penda Prediction"
    )
}

# Create a data.frame to compare the deregulation type (up or down)
# between DESeq2 and PENDA for the PDCL specified in args.
# Call the plot.volcano method to create the plot.
volcano.deseq2.vs.penda <- function(deseq2_res_list, penda_res, pdcl){
  data <- merge(as.data.frame(deseq2_res_list[[pdcl]]), penda_res, by="row.names")
  data <- data %>% select(log2FoldChange, padj, all_of(pdcl)) %>%
    dplyr::rename( penda_predict = all_of(pdcl) )
  data$penda_predict <- factor(data$penda_predict)
  
  plot_title <- paste0("Comparison DESeq2 fold change vs PENDA Result for ", pdcl)
  plot.volcano(data) + labs(title = plot_title )
}

# Create a data.frame to compare the deregulation type (up or down)
# between DESeq2 and PENDA for ALL PDCL.
# Call the plot.volcano method to create the plot.
volcano.deseq2.vs.penda.all <- function(deseq2_res_list, penda_res, pdcl_names){
  df <- map_dfr(pdcl_names, function(pdcl){
    data <- merge(as.data.frame(deseq2_res_list[[pdcl]]), penda_res, by="row.names")
    data <- data %>% 
      select(log2FoldChange, padj, all_of(pdcl) ) %>%
      dplyr::rename( penda_predict = all_of(pdcl) ) %>%
      mutate(pdcl = pdcl)
  })
  df$penda_predict <- as.factor(df$penda_predict)
  plot.volcano(df) + facet_wrap(~pdcl)
}
