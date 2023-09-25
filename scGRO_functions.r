
# inputs = count matrix and the number of
# iterations (permutations) to run
permuted_correlatePairs = function( countMatrix, Niters ) {

    # sampling probability = read count / total reads in each cell
    pmatrix = t(t(countMatrix) / colSums(countMatrix));
    # average across cells
    pvector = rowMeans(pmatrix);

    Ncells = ncol(countMatrix);
    Ngenes = nrow(countMatrix);
    Nreads = colSums(countMatrix);

    # Binarize counts
    obsx = countMatrix;
    obsx[ which(obsx>0) ] = 1;

    # Multiply each Ngenes x Ncells matrix by its
    # transpose within each iteration.
    # Because data is binary, this counts co-occurence
    # of 1's among all genes pairs.
    obsx = obsx %*% t(obsx) / Ncells;
    dim(obsx)

    # repeat each cellID by its read count
    cell_index = rep(1:Ncells, Nreads);
    samplesize = sum(Nreads);
    emp_p = matrix(0, nrow=Ngenes, ncol=Ngenes);

    for( n in 1:Niters ) {
        if( n %% 100 == 0 ) {
            message(n);
        }

        # randomly sample from genes with replacement
        simx = sample.int( Ngenes, size=samplesize, replace=T, prob=pvector);
        # assign sampled genes to cells based on read count
        simx = cbind(simx, cell_index);
        # binarize
        simx = unique(simx);
        simx = sparseMatrix(i=simx[,1], j=simx[,2], x=1, dims=c(Ngenes, Ncells));

        # compute coexpression
        simx = simx %*% t(simx) / Ncells;
        # compare to observed coexpression
        emp_p = emp_p + (simx >= obsx);
        gc();
    }
    
    results = data.frame(
        gene1 = rep(1:Ngenes, each=Ngenes),
        gene2 = rep(1:Ngenes, Ngenes),
        Pboth = as.vector(obsx),
        emp_p = as.vector(emp_p) / Niters
    );
    # remove redundant pair
    results = results[ results$gene1 < results$gene2, ];
    results$gene1 = rownames(countMatrix)[results$gene1];
    results$gene2 = rownames(countMatrix)[results$gene2];
    
    return(chr_corr);
}

######################################################################
######################################################################


# distance between molecules within same gene & cell
distance_to_neighbor = function(reads) {
	if(length(reads) < 2) {
		return(data.frame(gene=NA, type=NA, distance=NA, prev=NA));
	}
    #reads = reads[ order(reads$cellID, reads$gene, start(reads)) ];
    prevpos = reads[1:length(reads)-1];
    nextpos = reads[2:length(reads)  ];
	distance  = abs(start(nextpos) - start(prevpos));
    strpair = paste0(strand(prevpos), strand(nextpos));
    results = data.frame(
        gene = prevpos$gene,
        type = strpair,
        distance,
        prev = ifelse( strand(prevpos) == "+", start(prevpos), start(nextpos) )
    ) %>%
        filter(prevpos$gene == nextpos$gene & prevpos$cellID == nextpos$cellID);
	if(nrow(results) > 0) {
		return( results );
	}
	return(data.frame(gene=NA, type=NA, distance=NA, prev=NA));
}


# distance between all molecules within same gene & cell
distance_to_all = function(reads) {
	if(length(reads) < 2) {
		return(data.frame(gene=NA, type=NA, distance=NA));
	}
    reads = reads[ order(reads$cellID, reads$gene, start(reads)) ];
    # find first polymerase within each gene
    first  = reads[!duplicated(reads$gene)];
    others = reads[ duplicated(reads$gene)];
    gcount = table(reads$gene) %>% as.integer;
    # repeat first polymerase N times, where
    # N is the # of other polymerases in this gene
    first  = rep(first, times=gcount-1);
	distance = abs(start(others) - start(first));
    strpair = paste0(strand(first), strand(others));
    results = data.frame(gene = first$gene, type = strpair, distance);
	if(nrow(results) > 0) {
		return( results );
	}
	return(data.frame(gene=NA, type=NA, distance=NA));
}


# compute time between events from positions
time_between_events = function(reads, features, pol2speed=2500*60) {
    reads = reads[ order(reads$gene, reads$cellID, start(reads)) ] %>%
        mutate(count=1);
    allgenes = unique(reads$gene);
    allcells = unique(reads$cellID);
    
    prevpos = reads[1:length(reads)-1];
    nextpos = reads[2:length(reads)  ];
	distance  = abs(start(nextpos) - start(prevpos));
    results = data.frame(
        gene = prevpos$gene,
        time = distance/pol2speed
    ) %>%
        filter(prevpos$gene == nextpos$gene & prevpos$cellID == nextpos$cellID);
    #return est_times;
}

######################################################################
######################################################################


compute_distances = function( reads, cellIDs ) {
	# iterate through each cellID,
	# and combine results using cbind (column bind)
	distances = foreach(
		cell = unique(cellIDs),
		.combine="rbind"
	) %dopar% {
		dists = reads[cellIDs == cell] %>%
			distance_to_neighbor();
		return(dists);
	}
	return(distances);
}

compute_distances2 = function( reads, cellIDs ) {
	# iterate through each cellID,
	# and combine results using cbind (column bind)
	distances = foreach(
		cell = unique(cellIDs),
		.combine="rbind"
	) %dopar% {
		dists = reads[cellIDs == cell] %>%
			distance_to_neighbor_bygene();
		return(dists);
	}
	return(distances);
}

######################################################################
######################################################################

plot.cell.counts = function( reads, gene ) {
    Ncells = length(levels(reads$cellID));
	hits = subsetByOverlaps(reads, gene);

	real  = data.frame(ID=hits$cellID, source="observed");
    simd0 = data.frame(ID=hits$id0, source="permuted");

	output = rbind(real, simd0) %>%
        count(source, ID, .drop=F) %>%
        count(source, n, name="count") %>%
        mutate(fraction = round(count/Ncells, 3)) %>%
		na.omit;

    output %>%
        ggplot(aes(x=n, y=fraction)) +
        geom_col(fill="#dddddd", color="black") +
        geom_text(aes(y=0.6, label=count, angle=90)) +
        scale_x_continuous(labels=0:10, breaks=0:10) +
        ylim(0, 1) +
        ggtitle(paste( sub("GN-", "", names(gene)) ) ) +
        xlab("Molecules per cell") +
        ylab(paste0("Fraction of cells (N=", Ncells, ")")) +
		facet_wrap(~source)
}


######################################################################
######################################################################

plot.kinetics = function( reads, gene, genelen = NA ) {
    if( is.na(genelen) )
        genelen = 0.9*width(gene);
    genetime = genelen / 2500;

	allcells = unique(reads$cellID);
	reads = subsetByOverlaps(reads, gene);
	reads = reads[ order(start(reads)) ];

	dists = compute_distances( reads, reads$cellID );
	dists$source = "observed";
	
	simd0 = compute_distances( reads, reads$id0 );
	simd1 = compute_distances( reads, reads$id1 );
	simd2 = compute_distances( reads, reads$id2 );
	simd3 = compute_distances( reads, reads$id3 );
	simd4 = compute_distances( reads, reads$id4 );
	simd5 = compute_distances( reads, reads$id5 );
	simd6 = compute_distances( reads, reads$id6 );
	simd7 = compute_distances( reads, reads$id7 );
	simd8 = compute_distances( reads, reads$id8 );
	simd9 = compute_distances( reads, reads$id9 );
    simd0$source = "permuted0";
    simd1$source = "permuted1";
    simd2$source = "permuted2";
    simd3$source = "permuted3";
    simd4$source = "permuted4";
    simd5$source = "permuted5";
    simd6$source = "permuted6";
    simd7$source = "permuted7";
    simd8$source = "permuted8";
    simd9$source = "permuted9";

	distances = rbind(dists, simd0, simd1, simd2, simd3, simd4, simd5, simd6, simd7, simd8, simd9) %>%
		na.omit;
	
    # estimate bursts per minute from cell counts
	counts = table(reads$cellID);
    bpm = mean(counts[counts>1]) / genetime;

    real = distances %>%
		filter(source == "observed") %>%
        mutate(time=distance/2500);
	
	sim = distances %>%
		filter(source != "observed") %>%
        mutate(time=distance/2500);

    # test whether real and simulated distances are from same distribution
    kspval = ks.test( real$distance, sim$distance, alternative="greater", exact=F )$p.value;

    freq = table( round(real$time,0) );
    freq = data.frame( bin=names(freq), count=as.integer(freq)/sum(freq) );
    rownames(freq)=NULL;
    freq$bin = as.integer(freq$bin);
	freq$gen = "1. scGRO";

	sim = table( round(sim$time,0) );
    sim = data.frame( bin=names(sim), count=as.integer(sim)/sum(sim) );
    rownames(sim)=NULL;
    sim$bin = as.integer(sim$bin);
	sim$gen = "2. Permuted cells";

    rbind(freq, sim) %>%
        ggplot(aes(x=bin, y=count)) +
        geom_col() +
        geom_function( fun=dexp, args=list(rate=bpm), aes(linetype="Exponential") ) +
        scale_linetype_discrete(name = paste0(round(bpm, 3), " bursts/min")) +
        theme(legend.position = c(0.35,0.8)) +
        ggtitle(paste( sub("GN-", "", names(gene)), "p =", round(kspval, 3) ) ) +
        xlab("Time between events (min)") +
        ylab(paste0("Fraction of events (N=", nrow(real), ")")) +
        xlim(-0.5, genetime+0.5) +
		facet_wrap(~gen)
}

# to limit permutation to cells with more than 1 Pol II
# the cellID permutation is done fresh
plot.kinetics2 = function( reads, gene, genelen = NA ) {
    if( is.na(genelen) )
        genelen = 0.9*width(gene);
    genetime = genelen / 2500;

	reads = subsetByOverlaps(reads, gene); 
	reads = reads[ order(start(reads)) ];
    reads$cellID = droplevels(reads$cellID);
	allcells = levels(reads$cellID);
    
    reads$id0 = sample(reads$cellID);
    reads$id1 = sample(reads$cellID);
    reads$id2 = sample(reads$cellID);
    reads$id3 = sample(reads$cellID);
    reads$id4 = sample(reads$cellID);
    reads$id5 = sample(reads$cellID);
    reads$id6 = sample(reads$cellID);
    reads$id7 = sample(reads$cellID);
    reads$id8 = sample(reads$cellID);
    reads$id9 = sample(reads$cellID);

	dists = compute_distances( reads, reads$cellID );
	dists$source = "observed";
    
	simd0 = compute_distances( reads, reads$id0 );
	simd1 = compute_distances( reads, reads$id1 );
	simd2 = compute_distances( reads, reads$id2 );
	simd3 = compute_distances( reads, reads$id3 );
	simd4 = compute_distances( reads, reads$id4 );
	simd5 = compute_distances( reads, reads$id5 );
	simd6 = compute_distances( reads, reads$id6 );
	simd7 = compute_distances( reads, reads$id7 );
	simd8 = compute_distances( reads, reads$id8 );
	simd9 = compute_distances( reads, reads$id9 );
    simd0$source = "permuted0";
    simd1$source = "permuted1";
    simd2$source = "permuted2";
    simd3$source = "permuted3";
    simd4$source = "permuted4";
    simd5$source = "permuted5";
    simd6$source = "permuted6";
    simd7$source = "permuted7";
    simd8$source = "permuted8";
    simd9$source = "permuted9";

	distances = rbind(dists, simd0, simd1, simd2, simd3, simd4, simd5, simd6, simd7, simd8, simd9) %>%
		na.omit;
	
    # estimate bursts per minute from cell counts
	counts = table(reads$cellID);
    bpm = mean(counts[counts>1]) / genetime;

    real = distances %>%
		filter(source == "observed") %>%
        mutate(time=distance/2500);
	
	sim = distances %>%
		filter(source != "observed") %>%
        mutate(time=distance/2500);

    # test whether real and simulated distances are from same distribution
    kspval = ks.test( real$distance, sim$distance, alternative="greater", exact=F )$p.value;

    freq = table( round(real$time,0) );
    freq = data.frame( bin=names(freq), count=as.integer(freq)/sum(freq) );
    rownames(freq)=NULL;
    freq$bin = as.integer(freq$bin);
	freq$gen = "1. scGRO";

	sim = table( round(sim$time,0) );
    sim = data.frame( bin=names(sim), count=as.integer(sim)/sum(sim) );
    rownames(sim)=NULL;
    sim$bin = as.integer(sim$bin);
	sim$gen = "2. Permuted cells";

    rbind(freq, sim) %>%
        ggplot(aes(x=bin, y=count)) +
        geom_col() +
        geom_function( fun=dexp, args=list(rate=bpm), aes(linetype="Exponential") ) +
        scale_linetype_discrete(name = paste0(round(bpm, 3), " bursts/min")) +
        theme(legend.position = c(0.35,0.8)) +
        ggtitle(paste( sub("GN-", "", names(gene)), "p =", round(kspval, 3) ) ) +
        xlab("Time between events (min)") +
        ylab(paste0("Fraction of events (N=", nrow(real), ")")) +
        xlim(-0.5, genetime+0.5) +
		facet_wrap(~gen)
}

######################################################################
######################################################################

plot.ecdf = function( reads, gene ) {
    if( is.na(genelen) )
        genelen = 0.95*width(gene);
    genetime = genelen / 2500;

	allcells = unique(reads$cellID);
	reads = subsetByOverlaps(reads, gene);
	reads = reads[ order(start(reads)) ];

	dists = compute_distances( reads, reads$cellID );
	dists$gen = "scGRO";
	
	simd1 = compute_distances( reads, reads$id1 );
	simd2 = compute_distances( reads, reads$id2 );
	simd3 = compute_distances( reads, reads$id3 );
	simd4 = compute_distances( reads, reads$id4 );
	simd5 = compute_distances( reads, reads$id5 );
	simd6 = compute_distances( reads, reads$id6 );
	simd7 = compute_distances( reads, reads$id7 );
	simd8 = compute_distances( reads, reads$id8 );

	distances = rbind(dists, simd1, simd2, simd3, simd4, simd5, simd6, simd7, simd8) %>%
		na.omit;
	
    # estimate bursts per minute from cell counts
	counts = table(reads$cellID);
    bpm = mean(counts[counts>1]) / genetime;

    real = distances %>%
		filter(gen == "scGRO") %>%
        mutate(time=distance/2500);
	
	sim = distances %>%
		filter(gen == "sim") %>%
        mutate(time=distance/2500);

    # test whether real and simulated distances are from same distribution
    kspval = ks.test( real$distance, sim$distance, exact=F, alter )$p.value;

    freq = table( round(real$time,0) );
    freq = data.frame( bin=names(freq), count=as.integer(freq)/sum(freq) );
    rownames(freq)=NULL;
    freq$bin = as.integer(freq$bin);
	freq$gen = "1. scGRO";

	sim = table( round(sim$time,0) );
    sim = data.frame( bin=names(sim), count=as.integer(sim)/sum(sim) );
    rownames(sim)=NULL;
    sim$bin = as.integer(sim$bin);
	sim$gen = "2. Permuted cells";

	exp_norm = function(x, rate) {
		out = dexp(x, rate);
		return( out / sum(out) );
	}
    
    distances %>%
        filter(distance < 5000) %>%
        filter(type %in% c("++", "--")) %>%
        ggplot(aes(x=distance)) +
        stat_ecdf(geom = "step") +
        # geom_histogram(binwidth=100, boundary=0) +
        ggtitle("Consecutive Pol2") +
        xlab("Distance between molecules (bp)");
}


######################################################################
# signal = scGRO-seq reads object (GRanges)
# grange = Genomic range to be plotted (GRange)
# max.cells = max number of cells to plot (default 100)
# min.rpc = minimum reads per cell required to be shown (default 1)
# sortcells = sort cells based on polymerase position (default: true)
######################################################################

plot_polymerase_view = function( signal, grange, max.cells=100, min.rpc=1, sortcells=T ) {
    signal = subsetByOverlaps(signal, grange, ignore.strand=T);
    strsig = subsetByOverlaps(signal, grange, ignore.strand=F);
    #anchor= resize(grange, width=0, fix="start");

    #counts = table(signal$cellID);
    counts = table(strsig$cellID);
    counts = counts[counts >= min.rpc];
    signal = signal %>%
        filter( cellID %in% names(counts) );
    signal$cellID = droplevels(signal$cellID);

    ncells = n_distinct(signal$cellID);
    if( ncells > max.cells ) {
        signal = signal[ signal$cellID %in% sample(unique(signal$cellID))[1:max.cells] ];
    }
    if(!sortcells) {
        signal = signal[order(signal$cellID)];
        levels(signal$cellID) = 1:ncells;
    }
    signal$cellID = as.integer(signal$cellID);
    ncells = n_distinct(signal$cellID);
    all_cell_str = c(
        paste("-", unique(signal$cellID)),
        paste("+", unique(signal$cellID))
    );
    signal = signal %>%
        mutate(cellID = paste(strand, cellID)) %>%
        mutate(cellID=factor(cellID, levels=all_cell_str));
    
    #print(head(as.data.frame(signal)))
    
    as.data.frame(signal) %>%
        ggplot( aes(x=start, y=cellID, color=strand) ) +
        ggtitle(as.character(grange)) +
        scale_color_manual(values = c("#B10000", "#0099FF")) +
        geom_line(color="#e0bf00", lwd = 0.5) +
        geom_point(size=0.5, shape=15) +
        scale_y_discrete(drop=F) +
        xlim(start(grange), end(grange)) +
        xlab("Position (bp)") +
        ylab(paste0("Cells (n=", ncells, ")")) +
        theme(
            legend.position="none",
            panel.spacing = unit(0, "lines"),
            axis.line=element_blank()
        )
}


######################################################################
# reads = scGRO-seq reads object (GRanges)
# gene = Genomic range to be plotted (GRanges)
# enh = Genomic range to be plotted (GRanges)
# dreg = dreg peaks to be plotted (GRanges)
######################################################################

plot_position_pairs = function( reads, gene, enh, dreg=GRanges() ) {
    query = c(gene, enh);
    qstart= promoters(query, upstream=0, downstream=1);
    qreads = subsetByOverlaps(reads, query) %>%
        mutate( pos = distanceToNearest(., qstart) );
    qreads$pos.queryHits = NULL;
    qreads$str = as.character(strand(qreads));
    #mn_reads = which(qreads$str == "-");
    #if(length(mn_reads)>0)
    #    qreads$pos[mn_reads] = -1*qreads$pos[mn_reads];
    
    # only plot co-occuring polymerases
    hasGene  = qreads$cellID[ qreads$pos.subjectHits == 1 ];
    hasEnh   = qreads$cellID[ qreads$pos.subjectHits == 2 ];
    labels   = c( gene$name, enh$name );

    cotrans = as.data.frame(mcols(qreads)) %>%
        filter( cellID %in% hasGene ) %>%
        filter( cellID %in% hasEnh  ) %>%
        mutate( element = labels[ pos.subjectHits ] ) %>%
        arrange_at("cellID");
    cotrans$pos.subjectHits = NULL;
    
    models = data.frame(
        pos.distance = c(0, 0, width(query)),
        element = c(0.8, 2.2, 0.8, 2.2),
        cellID=c(1:2, 1:2)
    );
    
    dreg = subsetByOverlaps(dreg, query) %>%
        mutate( pos = distanceToNearest(., qstart[2]) ) %>%
        mutate( element=2.3, cellID=1 ) %>%
        as.data.frame;
    
    cotrans %>%
        arrange_at("pos.distance") %>%
        ggplot( aes(x=pos.distance, y=element, group=cellID) ) +
        # line connecting the Pol IIs
        geom_line(alpha=0.5, size=1, color="#e0bf00") +
        geom_line(data=models, size=2) +
        # Pol II dots:
        geom_beeswarm(size=2, color="#0099FF", groupOnX=F) + # alpha=0.5,
        # to disable the dREG, either make size 0 instead of 4 or color white:
        geom_point(data=dreg, size=0, pch='|', color="green") +
        xlab("Distance transcribed (bp)") +
        ylab("") +
        xlim(0, max(c(width(gene), width(enh))))
}

#####################

plot_2strand_position_pairs = function( reads, gene, enh, dreg ) {
    strand(enh) = "+";
    pl = plot_position_pairs( reads, gene, enh, dreg );
    strand(enh) = "-";
    mn = plot_position_pairs( scGRO, gene, enh, dreg );
    
    return(pl);
    return(mn);
}

#####################

# get matrix of Pol II positions. Nate wrote this but never used it:
get_position_matrix = function(reads, window, binsize=100) {
    tss = promoters(window, upstream=0, downstream=1);
    nbins = ceiling(width(window)/binsize);
    reads = subsetByOverlaps(reads, window) %>%
        mutate( pos = distanceToNearest(., tss) ) %>%
        mutate( bin = factor(ceiling(pos.distance/binsize), levels=1:nbins) );
    cells = levels(reads$cellID);
    #out = Matrix( nrow=nbins, ncol=length(cells), sparse=T );
    out = stats::xtabs( ~ bin + cellID, data=reads, sparse=T );
    colnames(out) = cells;
    
    return(out);
}

#####################

# get the table with position differnce between two features
get_position_table = function(reads, gene, enh) {
    gtss = promoters(gene, upstream=0, downstream=1);
    etss = promoters(enh,  upstream=0, downstream=1);
    
    greads = subsetByOverlaps(reads, gene) %>%
        mutate( pos = distanceToNearest(., gtss) ) %>%
        as.data.frame %>%
        select(cellID, pos=pos.distance);
    
    ereads = subsetByOverlaps(reads, enh ) %>%
        mutate( pos = distanceToNearest(., etss) ) %>%
        as.data.frame %>%
        select(cellID, pos=pos.distance);
    
    out = left_join(greads, ereads,
                    by=c("cellID"), suffix = c(".gene", ".enh"),
                    multiple="all") %>%
        na.omit %>%
        mutate(diff = pos.enh - pos.gene);
    #out = c(greads, ereads);
    #out = out %>%
    #    group_by(cellID) %>%
    #    mutate(pos.distance = pos.distance-pos.distance[1]) %>%
    #    ungroup;
    
    return(out);
}

#####################

# plot the table obtained from above function
plot_position_table = function(reads, gene, enh) {
    strand(enh) = "+"
    pl = get_position_table(reads, gene, enh);
    strand(enh) = "-"
    mn = get_position_table(reads, gene, enh);
    table = rbind(pl, mn);
    out = table %>%
    ggplot(aes(x=diff/1000)) +
    geom_histogram(binwidth=1, fill = "#426872") + 
    # scale_fill_manual(values=c()) +
    xlab("Difference in Polymerase position") +
    ylab(paste0("Number of cells", " - (", nrow(table), ")")) +
    # xlim(-30, 30) +
    ggtitle(paste(names(gene), "-", names(enh)));
    
    return(out);
}

#############################
#############################

getGOgenes = function( GOresult, GOterm ){
    selGenes = data.frame(GOresult) %>%
    filter( ID == GOterm );
    geneIDs = lapply(selGenes$geneID, function(x) {
        Ids = unlist(strsplit(gsub("/", ',', x), ","));
        # Ids = mapIds(org.Mm.eg.db, Ids, 'SYMBOL', 'ENTREZID');
        Ids = c(unname(Ids));
        return(Ids);
    })
    geneIDs = unlist(geneIDs)
        
    selGeneIDs = features %>%
        mutate( name = names) %>%
        filter( sub("GN-", "", names) %in% geneIDs) %>%
        data.frame();

    markerGenes = data.frame(
    Type = GOterm,
    Shape = "box",
    Chr = selGeneIDs$seqnames,
    Start = selGeneIDs$start,
    End = selGeneIDs$end,
    color = "6a3d9a");
    rownames(markerGenes) = sub("GN-", "", selGeneIDs$name);
    
    return(markerGenes);
}

#############################

# to get bed files for all genes in a GOterm, unlike below function that only gets if the genes in GOterm are present as a pair in corrF corrleation matrix:

# bedGOgenes = function( GOresult, clusterNumber ){
bedGOgenes = function( GOresult, GOterm ){
    selGenes = data.frame(GOresult) %>%
    # filter( Cluster == clusterNumber );
    filter( ID == GOterm );
    geneIDs = lapply(selGenes$geneID, function(x) {
        Ids = unlist(strsplit(gsub("/", ',', x), ","));
        # Ids = mapIds(org.Mm.eg.db, Ids, 'SYMBOL', 'ENTREZID');
        Ids = c(unname(Ids));
        return(Ids);
    })
    geneIDs = unlist(geneIDs)
        
    selGeneIDs = features %>%
        mutate( name = names) %>%
        filter( sub("GN-", "", names) %in% geneIDs) %>%
        # anchor_5p() %>%
        # mutate( start = start - 750 ) %>%
        # mutate( end = start + 1000 ) %>%
        promoters( upstream = 750, downstream = 250 ) %>%
        mutate( names = sub("GN-", "", names) ) %>% 
        data.frame();
    
    return(selGeneIDs);
}

#############################

# get bed files for genes in GOterm ONLY for genes if the co-expressed partner is also in the GOterm:
bedGOgenes = function( corrMatrix, GOresult, GOterm ){
    genesList = getGOgenes(GOresult, GOterm);
    net = corrMatrix %>%
        # Even if the corr instead of corrF matrix is used, it ensures that we consider relatively strongly co-Ex genes
        filter( corr >= 0.1 & pAdj <= 0.05 ) %>%
        filter( sub("GN-", "", geneA) %in% rownames(genesList) & sub("GN-", "", geneB) %in% rownames(genesList) );
        
    selGeneIDs = features %>%
        mutate( name = sub("GN-", "", names)) %>%
        filter( names %in% c(net$geneA, net$geneB)) %>%
        promoters( upstream = 1000, downstream = 500 ) %>%
        data.frame();
    
    return(selGeneIDs);
}

#############################


# gets both genes and enhancers:
getGOfeatures = function( corrMatrix, GOresult, GOterm ){
    selGenes = data.frame(GOresult) %>%
    filter( ID %in% GOterm );
    geneIDs = lapply(selGenes$geneID, function(x) {
        Ids = unlist(strsplit(gsub("/", ',', x), ","));
        # Ids = mapIds(org.Mm.eg.db, Ids, 'SYMBOL', 'ENTREZID');
        Ids = c(unname(Ids));
        return(Ids);
    })
    geneIDs = unlist(geneIDs)
        
    # select corr pairs of geneIDs:
    net = corrMatrix %>%
        filter(  pAdj < 0.05 ) %>% #  corr > 0.075 &
        filter( sub("GN-", "", Gene) %in% geneIDs );
    
    selGeneIDs = features %>%
        mutate( name = names) %>%
        filter( names %in% c(net$Gene, net$Enhancer)) %>%
        data.frame();

    markerGenes = data.frame(
    Type = ifelse(substr(selGeneIDs$name, 0, 3) == "GN-", "Gene", ifelse(substr(selGeneIDs$name, 0, 3) == "chr", "Enhancer", "SE")),
    Shape = ifelse(substr(selGeneIDs$name, 0, 3) == "GN-", "circle", ifelse(substr(selGeneIDs$name, 0, 3) == "chr", "triangle", "box")),
    # Shape = ifelse(substr(selGeneIDs$name, 0, 3) == "GN-", "circle", "triangle"),
    Chr = selGeneIDs$seqnames,
    Start = selGeneIDs$start,
    End = selGeneIDs$end,
    # color = ifelse(substr(selGeneIDs$name, 0, 3) == "GN-", "3a4664", "cc8921")
    color = ifelse(substr(selGeneIDs$name, 0, 3) == "GN-", "3a4664", ifelse(substr(selGeneIDs$name, 0, 3) == "chr", "cc8921", "cf4a49")))
    
    rownames(markerGenes) = sub("GN-", "", selGeneIDs$name);
    
    return(markerGenes);
}

#############################

# get bed files for genes in GOterm ONLY for genes if the co-expressed partner is also in the GOterm:
bedGOfeatures = function( corrMatrix, GOresult, GOterm ){
    genesList = getGOfeatures(corrMatrix, GOresult, GOterm);
    net = corrMatrix %>%
        # Even if the corr instead of corrF matrix is used, it ensures that we consider relatively strongly co-Ex genes
        filter( corr > 0.075 & pAdj < 0.05 ) %>%
        filter( sub("GN-", "", Gene) %in% rownames(genesList) & Enhancer %in% rownames(genesList) );
        
    selGeneIDs = features %>%
        # mutate( name = sub("GN-", "", names)) %>%
        filter( names %in% c(net$Gene, net$Enhancer)) %>%
        promoters( upstream = 750, downstream = 750 ) %>%
        data.frame();
    
    return(selGeneIDs);
}


#############################

# to make network plot using the gens from GOterm
networkPlot = function(corrMatrix, GOresult, GOterm ){
    genesList = getGOgenes(GOresult, GOterm);
    net = corrMatrix %>%
        filter( corr > 0.075 & pAdj <= 0.05 ) %>%
        mutate( Gene = sub("GN-", "", Gene) ) %>%
        filter( Gene %in% rownames(genesList) );
    
    net = graph_from_data_frame(net, directed = F) 

    # conditional color based on pAdj:
    # E(net)$color = ifelse(E(net)$pAdj <= 0.05 & E(net)$corr >= 0.1, "#ce968b", "gray")
    # conditional on corr value:
    E(net)$color = ifelse(E(net)$corr > 0.1, "tomato", "gray")
    # V(net)$label.color = ifelse(substr(V(net), 0, 3) == "GN-", "black", "white")

    netPlot = plot(net, 
                   vlayout = layout_on_sphere(net),
                   vertex.label.family = "Helvetica",
                   vertex.label.font = 1,
                   edge.arrow.size = .1, 
                   # edge.color = "gray", 
                   edge.width = 3,
                   vertex.color = "#cee2f4", 
                   vertex.label.color = "black",
                   vertex.frame.color = "white")
    
    return(netPlot);
    
    # ggsave(filename=sprintf("../plots/scGROv2p8_GxGmodules_max10kbp_binary_corr_%s%_network.pdf", GOterm), width=12, height=12, units="in")
}

#############################

# get bed files for genes in GOterm ONLY for genes if the co-expressed partner is also in the GOterm:
bedModules = function( corr, expressedFeatures, submodule ){
    # get features in the submodule:
    proms =  expressedFeatures[ unlist(submodule) ] %>%
        # filter features such that they have to be present in the correlated matrix
        filter( names %in% c(corr$Gene, corr$Enhancer)) %>%
        # 750 of gene is -500 to +250, as the gene starts are trimmed 250 nt
        # It is redundant for enhancers and the two 750 nt blocks overlap 500 nt
        promoters( upstream = 750, downstream = 0 ) %>%
        # removing the excess 250 nt
        # the resulting 500 nt is just promoter of genes and center of enhancers
        anchor_5p() %>%
        resize( width = 500 ) %>%
        # promoters( ifelse(substr(names, 0, 3) == "GN", 
        #                   (upstream = 750, downstream = 0), 
        #                   (upstream = 500, downstream = 0)) ) %>%
        data.frame();
        return(proms); 
}


#############################