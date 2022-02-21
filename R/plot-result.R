######## Load Library needed
library(ggplot2)
library(ggh4x)

# Made by Jorge
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

# Made by Jorge
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

collapse.data.by.pathways <- function(gost.penda.gem, gsea.deseq2.gem, pdcl.samples.name, db.source){
  collapsed <- data.frame()
  for (pdcl in pdcl.samples.name){
    for (db in db.source) {
      gost <- gost.penda.gem[[pdcl]][[db]]
      gsea <- gsea.deseq2.gem[[pdcl]][[db]]
      join <- gost %>% 
        full_join(gsea, by = "pathway_id") %>%
        dplyr::rename(p_value_gprofiler = p_value.x, p_value_gsea = p_value.y, fdr_gprofiler = fdr.x, fdr_gsea = fdr.y) %>%
        select(pathway_id, p_value_gprofiler, fdr_gprofiler, p_value_gsea, fdr_gsea) %>%
        mutate(pdcl = pdcl, db = db) 
      collapsed <- rbind(collapsed, join)
    }
  }
  collapsed
}

collapse.data.by.method <- function(collapsed.by.pathways){
  gprofiler <- collapsed.by.pathways %>%
    select(pathway_id, fdr_gprofiler, pdcl, db, category) %>%
    dplyr::rename(fdr = fdr_gprofiler) %>%
    mutate(method = "gprofiler")
  
  gsea <- collapsed.by.pathways %>%
    select(pathway_id, fdr_gsea, pdcl, db, category) %>%
    dplyr::rename(fdr = fdr_gsea) %>%
    mutate(method = "gsea")
  
  rbind(gprofiler, gsea)
}

collapse.gem.list <- function(gem_list){
  collapsed_by_method <- gem_list
  collapsed_by_method <-  map_depth(.x = collapsed_by_method, .depth = 2, .f = function(x){
    map.list.to.df(x, "db")
  })
  
  collapsed_by_method <-  map_depth(.x = collapsed_by_method, .depth = 1, .f = function(x){
    map.list.to.df(x, "pdcl")
  })
  
  collapsed_by_method <-  map_depth(.x = collapsed_by_method, .depth = 0, .f = function(x){
    map.list.to.df(x, "method")
  })
  
  collapsed_by_method <- assign.categories(collapsed_by_method)
  
  collapsed_by_method <- collapsed_by_method %>%
    filter(!(undesired | is.na(category)))
  
  collapsed_by_method$category <- factor(
    collapsed_by_method$category,
    levels = names(word_list)
  )
  
  collapsed_by_method$pdcl <- factor(
    collapsed_by_method$pdcl,
    levels = pdcl_names
  )
  
  collapsed_by_method$method <- factor(
    collapsed_by_method$method,
    levels = unique(collapsed_by_method$method)
  )
  
  collapsed_by_method
}

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

find.common.pathways <- function(collapsed.by.pathways, threshold){
  common <- collapsed.by.pathways 
  common <- common %>%
    mutate(is.common = ( common$fdr_gprofiler<threshold & common$fdr_gsea<threshold ) ) %>%
    select(pdcl, pathway_id, is.common)
  common$is.common[is.na(common$is.common)] <- F
  
  common$is.common <- factor(
    common$is.common, 
    levels = c("TRUE", "FALSE"), 
    labels = c("Common Pathways", "Specific Pathways")
  )
  
  common
}

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

plot.gprofiler.vs.gsea <- function(x){
  ggplot(
    data = x, 
    aes(
      x = -log10(fdr_gprofiler), 
      y = -log10(fdr_gsea),
      colour = is.common
    )
  ) +
    geom_point(
      size = 0.5
    ) +
    scale_colour_manual(values = c("red", "grey30")) +
    facet_wrap(vars(pdcl)) +
    labs(
      title = "G:Profiler against GSEA",
      subtitle = paste0("Threshold alpha : ", threshold),
      x = "-log10(pvalue G:Profiler)",
      y = "-log10(pvalue GSEA)",
      colour = ""
    ) +
    theme(
      plot.title = element_text(hjust = 0.5,face = "bold",color = "red"),
      plot.subtitle = element_text(hjust = 0.5,face = "italic",color = "blue"),
      axis.text.x = element_text(angle = 90)
    ) 
}

plot_method_vs_method <- function(x, m1, m2, threshold){
  data <- x %>% filter(method_1 == m1, method_2 == m2)
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
      x = paste0( "-log10(fdr ",  m2, ")"),
      y = paste0( "-log10(fdr ",  m2, ")"),
      colour = ""
    ) +
    theme(
      plot.title = element_text(hjust = 0.5,face = "bold",color = "red"),
      plot.subtitle = element_text(hjust = 0.5,face = "italic",color = "blue"),
      axis.text.x = element_text(angle = 90)
    ) 
}

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

generate.breaks <- function(min, max, midpoint, n){
  left <- seq(from = min, to = midpoint, length.out = n ) 
  right <- seq(from = midpoint, to = max, length.out = n )[-1]
  c(left, right)
}

fdr.to.factor <- function(x, n){
  break.vect <- generate.breaks(0, 1, threshold, n)
  cut.vect <- cut(x$fdr, breaks = break.vect, labels = break.vect[-1])
  cut.vect <- factor(cut.vect, levels = rev(levels(cut.vect)))
  x$fdr <- cut.vect
  x
}

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
    )
}

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

collapse.rnk <- function(rnk){
  rnk <- imap_dfr(rnk, .f = function(rank,pdcl){
    data.frame(gene = names(rank), rank = rank, pdcl = pdcl, row.names = NULL)
  })
  rnk$pdcl <- factor(rnk$pdcl)
  rnk
}

summarize.rank <- function(rnk.collaspe){
  rnk.collaspe %>% 
    group_by(gene) %>% 
    summarise(mean = mean(rank), sd = sd(rank), abs.mean = mean(abs(rank))) %>%
    arrange(desc(abs.mean))
}

get.genes.in.pathway <- function(collapsed.by.method){
  genes.in.common.pathways <- collapsed.by.method %>%
    mutate(genes = strsplit(genes, ",")) %>%
    unnest(genes) %>%
    dplyr::rename(gene = genes)
}

get.genes.in.common.pathways <- function(collapsed.by.method, all.common.pathways, threshold, m1, m2, common.method.col){
  alpha_col <- paste0("alpha_", threshold)
  common_pathways <- all_common_pathways %>%
    filter( .data[[alpha_col]] == "Common Pathways", method_1 == m1, method_2 == m2) %>%
    dplyr::rename(method = all_of(common.method.col))
  
  genes.in.common.pathways <- get.genes.in.pathway(collapsed.by.method) %>%
    inner_join(common_pathways, by = c("pdcl", "pathway_id", "method"))
  genes.in.common.pathways
}

plot.common.pathway.genes <- function(rnk.collapsed, collapsed.by.method, all.common.pathways, threshold, m1, m2, gsea.col, nb.gene){
  genes.in.common.pathways <- get.genes.in.common.pathways(collapsed.by.method, all.common.pathways, threshold, m1, m2, gsea.col) %>%
    add_count(pathway_id, gene)
  
  rnk <- rnk.collapsed %>% filter(gene %in% unique(genes.in.common.pathways$gene)) 
  rnk_summary <- summarize.rank(rnk)
  rnk_summary <- rnk_summary %>% dplyr::slice(1:nb.gene)
  
  top.genes <- rnk_summary$gene
  
  rnk <- filter(rnk, gene %in% top.genes)
  genes.in.common.pathways <- filter(genes.in.common.pathways, gene %in% top.genes)
  
  plot_title <- paste0("Top ", nb.gene, " genes in common pathways")
  plot_subtitle <- paste0("Pathways are common between : ", m1, " and ", m2)
  plot <- ggplot(genes.in.common.pathways, aes(gene, pathway_id, size = n) ) +
    geom_point() +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Gene Symbol",
      y = "Pathways",
      size = "Number of PDCL"
    ) +
    theme(axis.text.x = element_text(angle = 90) ) +
    scale_y_discrete(limits = rev)
  plot
}

heatmap.common.gene.ranks <- function(rnk.collapsed, collapsed.by.method, all.common.pathways, threshold, m1, m2, gsea.col, nb.gene){
  genes.in.common.pathways <- get.genes.in.common.pathways(collapsed.by.method, all.common.pathways, threshold, m1, m2, gsea.col) %>%
    add_count(pathway_id, gene)
  
  rnk <- rnk.collapsed %>% filter(gene %in% unique(genes.in.common.pathways$gene)) 
  rnk_summary <- summarize.rank(rnk)
  rnk_summary <- rnk_summary %>% dplyr::slice(1:nb.gene)
  
  top.genes <- rnk_summary$gene
  
  rnk <- filter(rnk, gene %in% top.genes)
  genes.in.common.pathways <- filter(genes.in.common.pathways, gene %in% top.genes)
  
  plot_title <- paste0("Top ", nb.gene, " genes rank in common pathways")
  plot_subtitle <- paste0("Pathways are common between : ", m1, " and ", m2)
  heatmap <- ggplot(rnk, aes(pdcl, gene, fill = rank)) +
    geom_tile() +
    scale_fill_gradient2(
      high = "red",
      low = "blue",
      na.value = "gray30"
    ) +
    geom_tile(
      data = genes.in.common.pathways,
      mapping = aes(x = pdcl, y = gene),
      inherit.aes = F,
      colour = 'black',
      fill = NA,
      size = 0.5
    ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "PDCL",
      y = "Gene Symbol",
      fill = "Rank"
    ) +
    theme(axis.text.x = element_text(angle = 90) ) +
    scale_y_discrete(limits = rev)
  heatmap
}


plot.gsea.pathway.genes <- function(rnk.collapsed, collapsed.by.method, gsea.method.name, threshold, nb.gene){
  genes_to_pathway <- get.genes.in.pathway(collapsed.by.method) %>%
    filter(fdr<threshold, method == gsea.method.name) %>%
    add_count(pathway_id, gene)
  
  rnk <- rnk.collapsed %>% filter(gene %in% unique(genes_to_pathway$gene)) 
  rnk_summary <- summarize.rank(rnk)
  rnk_summary <- rnk_summary %>% dplyr::slice(1:nb.gene)
  
  top.genes <- rnk_summary$gene
  
  rnk <- filter(rnk, gene %in% top.genes)
  genes_to_pathway <- filter(genes_to_pathway, gene %in% top.genes)
  
  plot_title <- paste0("Top ", nb.gene, " genes in enriched pathways for ", gsea.method.name)
  plot_subtitle <- paste0("Alpha : ", threshold)
  plot <- ggplot(genes_to_pathway, aes(gene, pathway_id, size = n) ) +
    geom_point() +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "Gene Symbol",
      y = "Pathways",
      size = "Number of PDCL"
    ) +
    theme(axis.text.x = element_text(angle = 90) ) +
    scale_y_discrete(limits = rev)
  plot
}

heatmap.gsea.gene.ranks <- function(rnk.collapsed, collapsed.by.method, gsea.method.name, threshold, nb.gene){
  genes_to_pathway <- get.genes.in.pathway(collapsed.by.method) %>%
    filter(fdr<threshold, method == gsea.method.name)
  
  rnk <- rnk.collapsed %>% filter(gene %in% unique(genes_to_pathway$gene)) 
  rnk_summary <- summarize.rank(rnk)
  rnk_summary <- rnk_summary %>% dplyr::slice(1:nb.gene)
  
  top.genes <- rnk_summary$gene
  
  rnk <- filter(rnk, gene %in% top.genes)
  genes_to_pathway <- filter(genes_to_pathway, gene %in% top.genes)
  
  plot_title <- paste0("Top ", nb.gene, " gene ranks in enriched pathways for ", gsea.method.name, " for each PDCL")
  plot_subtitle <- paste0("Alpha : ", threshold)
  heatmap <- ggplot(rnk, aes(pdcl, gene, fill = rank)) +
    geom_tile() +
    scale_fill_gradient2(
      high = "red",
      low = "blue",
      na.value = "gray30"
    ) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "PDCL",
      y = "Gene Symbol",
      fill = "Rank"
    ) +
    theme(axis.text.x = element_text(angle = 90) ) +
    scale_y_discrete(limits = rev)
  heatmap
}