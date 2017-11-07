library(data.table)

reads = "/data/processed_data/RB7/nanopore/asmset_corrected/canu_corr150x.min20kb.renamed.correctedReads.fasta.gz"
valreads = "/data/processed_data/RB7/nanopore/valset/validation.fastq.gz"
threads = 55

run_minimap_ava = function(k, f, threads, reads, output){
    cmd = paste("minimap2 -x ava-ont -k", k, " -f", f, " -t", threads, " -I200G", reads reads, " | gzip -1  > ", output)
    system(cmd)    
}

run_miniasm = function(params, reads, threads, output){
    "miniasm -f $reads -m$M -i$i -s$S -c$C -h$H -I$I -g$G -d$D $paffile > $dir"assembly.gfa""
}

run_minimap_val = function(){
    
}

run_quals = function(){
    
}

propose_new_params = function(params, param, change, lower, upper){
    # take a list of params
    # return up to 6 new lists
    
    r = list()
    if(params[[param]] + change < upper){
        p1 = params
        p1[[param]] = p1[[param]] + change
        r$p1 = p1
    }
    
    if(params[[param]] + 2*change < upper){
        p2 = params
        p2[[param]] = p2[[param]] + 2*change
        r$p2 = p2
    }
    
    if(params[[param]] + 10*change < upper){
        p3 = params
        p3[[param]] = p3[[param]] + 10*change
        r$p3 = p3
    }

    if(params[[param]] - change > lower){
        p4 = params
        p4[[param]] = p4[[param]] - change
        r$p4 = p4
    }
    
    if(params[[param]] - 2*change > lower){
        p5 = params
        p5[[param]] = p5[[param]] - 2*change
        r$p5 = p5
    }
    
    if(params[[param]] - 10*change > lower){
        p6 = params
        p6[[param]] = p6[[param]] - 10*change
        r$p6 = p6
    }
    
    r = as.data.frame(rbindlist(r))
    
    return(r)
}

get_new_assembly_params = function(params){
    
    # make a list of lists, containing all new sets of assemblies to analyse
    M = propose_new_params(params, "M", 10, 0, 100000) 
    i = propose_new_params(params, "i", 0.005, 0, 1) 
    S = propose_new_params(params, "S", 10, 0, 100000) 
    C = propose_new_params(params, "C", 1, 0, 150)
    H = propose_new_params(params, "H", 10, 0, 100000) 
    I = propose_new_params(params, "I", 0.005, 0, 1) 
    G = propose_new_params(params, "G", 10, 0, 100000) 
    #D = propose_new_params(params, "D", 500, 0, 1000000000) # max is 2x genome size
    
    r = rbind(M, i, S, C, H, I, G)
    
    return(r)
}


ava_params = list(k=21,     # kmer length
                  f=0.012  # ignore top fraction of most occurring minimizers in minimap, typically these are repeats we want to ignore [default 0.001]
)

assembly_params = list(   M=1000,   # min match length [default 100]
                          i=0.10,  # min identity [0.05]
                          S=2000,  # min span [2000]
                          C=3,     # min coverage [3]
                          H=1000,  # max overhang length [ 1000]
                          I=0.8,   # max end-to-end match ratio [0.8]
                          G=1000,  # max gap differences between reads for trans-reduction [1000]
                          D=500000000 # max distance for bubble popping [50000]
)

best_mqs = 391960095458

for(i in 1:100){
    
    new_assembly_params = get_new_assembly_params(assembly_params)
    results = run_assemblies(new_assembly_params)
    new_best_mqs = get_best_mqs(results)
    
    if(new_best_mqs>best_mqs){
        params = get_best_params(results)
    }else{
        # try some new mapping 
        get_new_mapping_params(get_mapping_params)
        results = run_mappings(mapping_params)
        new_best_mqs = get_best_mqs(results)
        if(new_best_mqs>best_mqs){
            params = get_best_params(results)
        }else{
            # done
            break
        }
    }
}
