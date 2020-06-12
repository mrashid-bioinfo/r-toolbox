vector_to_color_map = function( vec = NULL )
{
    vec_uniq = unique( vec )
    return(sample(length(colors()), length(vec_uniq), replace=F))
}



## ----------------------------------------------------------------- ##
##      Overrides lolliplot function from track-viewer package   
##      Only does it better 
##  
## ----------------------------------------------------------------- ##


tSNE_pos_neg_ratio_plot = function( matrix = NULL, zero_conversion_thrs = 1e-4, round_digit = 5)
{

    # matrix = tSNE_annotated_enhancer_pos_neg$ratio_smoothed
    # zero_conversion_thrs = 1e-4
    # round_digit = 5
    # image_file = "tSNE.pdf"
    # # legend_file = "legend_file.pdf"
    # require(ggplot2)

    mat_dim = dim(matrix)

    ## ------- ##
    ## Extract ratio matrix and transform to vector 
        matrix_log2 = log2(matrix)
        matrix_log2_vec = as.vector( matrix_log2 )
        ## Need some manual fiddling. Data points too close to zero are assigned 
        matrix_log2_vec_1 = ifelse( matrix_log2_vec < zero_conversion_thrs & matrix_log2_vec > -1*(zero_conversion_thrs), 0 , matrix_log2_vec)
        matrix_log2_vec_2 = round( matrix_log2_vec_1, digits = round_digit )

    ## ------- ##
    ## Bin the data points to nearest quantile values
    ##
    ##  call function : quantile_transform
    ##
        pos_ix = matrix_log2_vec_2 > 0
        pos_values = matrix_log2_vec_2[ pos_ix ]
        pos_quants = quantile( pos_values, (0:20)*0.05 )
        pos_values_q = sapply( pos_values, quantile_transform, pos_quants, "greater" )
        matrix_log2_vec_2[ pos_ix ] = pos_values_q

        neg_ix = matrix_log2_vec_2 < 0 
        neg_values = matrix_log2_vec_2[ neg_ix ]
        neg_quants = quantile( neg_values, (0:20)*0.05 )
        neg_values_q = sapply( neg_values, quantile_transform, neg_quants, "smaller" )
        matrix_log2_vec_2[ neg_ix ] = neg_values_q

        dim(matrix_log2_vec_2) = mat_dim

    ## ------- ##
    ## Very small values are converted to zero
        matrix_log2_vec_2_m = melt( matrix_log2_vec_2 )
        matrix_log2_vec_2_m$value = ifelse( matrix_log2_vec_2_m$value == 0, NA, matrix_log2_vec_2_m$value )

    ## ------- ##
    ## Create color pallete
        pos_color_values = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(length(unique(pos_values_q)))
        names(pos_color_values) = sort(unique(pos_values_q))
        neg_color_values = colorRampPalette(rev(brewer.pal(n = 9, name ="YlGnBu")))(length(unique(neg_values_q)))
        names(neg_color_values) = sort(unique(neg_values_q))
        color_code = c( neg_color_values, c("0"="white"), pos_color_values  )
        names(color_code) = round(as.numeric(names(color_code)), digits = round_digit )
        
    ## ------- ##
    ## Plot ratio image and legend
        
        point = ggplot( matrix_log2_vec_2_m, aes( Var1, Var2) ) + geom_point(aes(colour = factor(value))) + scale_color_manual( values = color_code ) + theme_bw() + theme( legend.position = "none" ) + xlab("Dimension 1") + ylab("Dimension 2")
        
        sub_ix = c(1,3,5,7,10,15,22,25,27,33,38,40,43)
        legend_df = data.frame( color = color_code , id = 1:length(color_code), value = round(as.numeric(names(color_code)),digits=5) )
        bar = ggplot( legend_df[sub_ix,], aes( x = id, y = value, fill = factor(value) ) ) + geom_bar(stat = "identity") + scale_fill_manual(values = color_code[sub_ix]) + theme_bw() + coord_flip() + theme(legend.box ="vertical")

        return( list( point, bar ) )

}



## ----------------- ##
##
##  Function :: Plot two vector with error range 
##  
##   Example :: Two Classifier Fscore / AUC score
##
##  Input : 
##      1. Data frame with various run / negative sets
##      2. Column name holding classifier name [ Still hard coded ]
##      3. Classifier one name 
##      4. Classifier two name 
##

clf_pair_metric = function( df, clf_column="classifier", clf1 = "RFC", clf2 = "RFC_Smote" )
{
    temp1 = as.data.frame(df %>% filter(classifier == clf1))
    temp2 = as.data.frame(df %>% filter(classifier == clf2))
    
    temp_df = data.frame( classifier = temp1$classifier , negative = temp1$negative, 
    temp1.med = temp1$med, temp1.min = temp1$min, temp1.max = temp1$max, 
    temp2.med = temp2$med, temp2.min = temp2$min, temp2.max = temp2$max )
    temp_df$negative = factor( temp_df$negative, levels = unique( temp_df$negative ), ordered = T  )
    
    p = ggplot( temp_df, aes(x=temp1.med, y = temp2.med,  color = factor(negative)) ) + geom_errorbar( aes(ymin=temp2.min, ymax=temp2.max) ) + geom_errorbarh( aes(xmin=temp1.min, xmax=temp1.max) ) + theme_bw() + xlab(clf1)+ ylab(clf2) + geom_smooth(aes(group=factor(classifier)), colour=c("green4"), size = 0.5 , alpha=0.2 ) + theme( axis.text.x = element_text(angle=90,vjust = 0.5) )  + geom_abline(linetype="dotted",intercept = 0, slope = 1)
        # + scale_x_continuous( limits=c(0,1), breaks = seq(0,1,0.1) ) + scale_y_continuous( limits=c(0,1), breaks = seq(0,1,0.1) )
    
    return(p)
}

## ----------------- ##
##
##  Function :: Function to plot raw vs normalised ( CPM and log-CPM ) RNA seq data
##
##  Input : 
##      1. Raw data object [ edgeR ]
##      2. Normalised and filtered data object [ edgeR ]
##      3. samplenames
##      4. file_path
##
plot_pre_post_norm = function( raw = NULL, normf = NULL, samplenames = NULL, file_path = "raw_vs_normalized.pdf" )
{
    library(RColorBrewer)
    require(edgeR)
    nsamples <- ncol(raw)

    col <- brewer.pal(nsamples, "Paired")
    pdf( file = file_path, width = 12, height = 8 );
    par(mfrow=c(1,2))

    # raw data 
    lcpm_raw = cpm( raw$counts, log = TRUE )
    plot(density(lcpm_raw[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, main="", xlab="")
    title(main="A. Raw data", xlab="Log-cpm")
    abline(v=0, lty=3)
    for (i in 2:nsamples){
     den <- density(lcpm_raw[,i])
     lines(den$x, den$y, col=col[i], lwd=2)
    }
    legend("topright", samplenames, text.col=col, bty="n")
    
    # raw data 
    lcpm_norm <- cpm(normf, log=TRUE)
    plot(density(lcpm_norm[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2, 
         main="", xlab="")
    title(main="B. Filtered data", xlab="Log-cpm")
    abline(v=0, lty=3)
    for (i in 2:nsamples){
       den <- density(lcpm_norm[,i])
       lines(den$x, den$y, col=col[i], lwd=2)
    }
    legend("topright", samplenames, text.col=col, bty="n")

    dev.off();

}




## ===================================================  ##
##    funtion to plot signature contribution per sample 

plot_sample_contribution = function( assigned_file_path = NULL, number_of_signatures = 2, id_file_path = NULL, clinical_df = NULL , label = NULL )
{
    require(reshape2)
    require(ggplot2)
    assigned = read.delim( assigned_file_path, sep = " ", header = F, stringsAsFactors = F );
    EmuReadyID = read.delim( id_file_path, sep = "\t", header = F, stringsAsFactors = F )

    sig_label = paste( "Signature ", 1:number_of_signatures, sep = "" ) 
    contribution = data.frame( SampleID = EmuReadyID[,1], assigned[ , 1:number_of_signatures] );
    rownames(contribution) = contribution$SampleID
    contribution_o = contribution[ clinical_df$SampleID, ]  
    colnames( contribution_o )[2:dim(contribution_o)[2]] = sig_label;

    contribution_o_m = melt(contribution_o)
    contribution_o_m$SampleID = factor( as.character( contribution_o_m$SampleID ), levels = unique( as.character( contribution_o_m$SampleID ) ), ordered = T )

    contrib_fig_file_path = paste( label, "_", number_of_signatures, "_signatures.pdf", sep = "" );
    pdf(file = contrib_fig_file_path, width=12, height=5 );
    contribution = ggplot(contribution_o_m,aes(x = SampleID, y = value,fill = variable,width=.85) ) + theme_bw() + geom_bar(colour="grey50",position = "fill",stat = "identity",  alpha = 0.95 ) + scale_y_continuous(labels = scales::percent ) + scale_fill_manual(values = c( "Signature 1" = "turquoise3", "Signature 2" = "khaki1", "Signature 3" = "yellowgreen") )
    contribution = contribution + theme(axis.text.x = element_text(size = 15, angle = 90, hjust = 0.5, vjust = 0.5) )
    contribution = contribution + theme(axis.text.y = element_text(size = 15, angle = 0, hjust = 0.5, vjust = 0.5) )
    print(contribution)
    dev.off()
}




my_lolliplot = function (SNP.gr, features = NULL, ranges = NULL, type = c("circle", 
    "pie", "pin", "pie.stack"), newpage = TRUE, ylab = TRUE, 
    yaxis = TRUE, xaxis = TRUE, legend = NULL, cex = 1, dashline.col = "gray80", 
    jitter = c("node", "label"), ...) 
{

    ## Loading all the relevent functions 
    require(trackViewer)    
    r_files = dir("/Users/rm8/Software/R_packages_Source/trackViewer/R/", full.names = T)
    for(i in 1:length(r_files))
    {
        source(r_files[i])
    }

    ## Actual lolliplot function ##

    stopifnot(inherits(SNP.gr, c("GRanges", "GRangesList", "list")))
    stopifnot(inherits(features, c("GRanges", "GRangesList", 
        "list")))
    jitter <- match.arg(jitter)
    if (type != "circle" && jitter == "label") {
        jitter <- "node"
        warning("if jitter set to label, type must be cirle.")
        message("jitter is set to node.")
    }
    SNP.gr.name <- deparse(substitute(SNP.gr))
    if (class(SNP.gr) == "GRanges") {
        SNP.gr <- GRangesList(SNP.gr)
        names(SNP.gr) <- SNP.gr.name
    }
    len <- length(SNP.gr)
    for (i in 1:len) {
        stopifnot(class(SNP.gr[[i]]) == "GRanges")
    }
    TYPES <- c("circle", "pie", "pin", "pie.stack")
    if (any(!type %in% TYPES)) {
        stop("Error in match argument: ", paste0("'type' should be one of '", 
            paste(TYPES, collapse = "', '"), "'."))
    }
    types <- rep(type, length = len)[1:len]
    rm(type)
    if (length(legend) > 0) {
        if (!is.list(legend)) {
            tmp <- legend
            legend <- vector(mode = "list", length = len)
            legend[[len]] <- tmp
            rm(tmp)
        }
        else {
            if (length(legend) == 1) {
                tmp <- legend[[1]]
                legend <- vector(mode = "list", length = len)
                legend[[len]] <- tmp
                rm(tmp)
            }
            else {
                if ("labels" %in% names(legend)) {
                  tmp <- legend
                  legend <- vector(mode = "list", length = len)
                  legend[[len]] <- tmp
                  rm(tmp)
                }
                else {
                  if (length(legend) < len) {
                    length(legend) <- len
                  }
                }
            }
        }
    }
    features.name <- deparse(substitute(features))
    if (length(ranges) > 0) {
        stopifnot(class(ranges) == "GRanges")
        ranges <- rep(ranges, length(SNP.gr))[1:length(SNP.gr)]
        stopifnot(length(ranges) == length(SNP.gr))
    }
    else {
        if (class(features) == "GRanges") {
            ranges <- range(features)[rep(1, len)]
        }
        else {
            if (length(features) != len) {
                stop("if both SNP.gr and features is GRangesList,", 
                  " the lengthes of them should be identical.")
            }
            ranges <- unlist(GRangesList(lapply(features, range)))
        }
    }
    if (class(ranges) == "GRanges") {
        for (i in len) {
            range <- ranges[i]
            stopifnot(all(width(SNP.gr[[i]]) == 1))
            ol <- findOverlaps(SNP.gr[[i]], range)
            SNP.gr[[i]] <- SNP.gr[[i]][queryHits(ol)]
        }
    }
    height <- 1/len
    height0 <- 0
    if (newpage) 
        grid.newpage()
    for (i in 1:len) {
        type <- match.arg(types[i], TYPES)
        if (type == "pin") {
            pinpath <- system.file("extdata", "map-pin-red.xml", 
                package = "trackViewer")
            pin <- readPicture(pinpath)
        }
        else {
            pin <- NULL
        }
        vp <- viewport(x = 0.5, y = height0 + height * 0.5, width = 1, 
            height = height)
        pushViewport(vp)
        LINEW <- as.numeric(convertX(unit(1, "line"), "npc"))
        LINEH <- as.numeric(convertY(unit(1, "line"), "npc"))
        GAP <- 0.2 * LINEH
        ratio.yx <- 1/as.numeric(convertX(unit(1, "snpc"), "npc"))
        if (inherits(features, c("GRangesList", "list"))) {
            feature <- features[[i]]
            stopifnot(class(feature) == "GRanges")
        }
        else {
            feature <- features
        }
        feature$height <- convertHeight2NPCnum(feature$height)
        if (length(feature$featureLayerID) != length(feature)) {
            feature$featureLayerID <- rep("1", length(feature))
        }
        feature <- feature[end(feature) >= start(ranges[i]) & 
            start(feature) <= end(ranges[i])]
        feature$featureLayerID <- as.character(feature$featureLayerID)
        start(feature)[start(feature) < start(ranges[i])] <- start(ranges[i])
        end(feature)[end(feature) > end(ranges[i])] <- end(ranges[i])
        feature.splited <- split(feature, feature$featureLayerID)
        bottomblank <- 4
        if (length(names(feature)) > 0) {
            feature.s <- feature[!duplicated(names(feature))]
            ncol <- getColNum(names(feature.s))
            bottomblank <- max(ceiling(length(names(feature.s))/ncol), 
                4)
            pushViewport(viewport(x = 0.5, y = bottomblank * 
                LINEH/2, width = 1, height = bottomblank * LINEH, 
                xscale = c(start(ranges[i]), end(ranges[i]))))
            color <- if (length(unlist(feature.s$color)) == length(feature.s)) 
                unlist(feature.s$color)
            else "black"
			
            fill <- if (length(unlist(feature.s$fill)) == length(feature.s)) 
                unlist(feature.s$fill)
            else "black"
            pch <- if (length(unlist(feature.s$pch)) == length(feature.s)) 
                unlist(feature.s$pch)
            else 22

            print(color)
            print(fill)
            print(ncol)
            
            grid.legend(label = names(feature.s), ncol = ncol, 
                byrow = TRUE, vgap = unit(0.2, "lines"), pch = pch, 
                gp = gpar(col = color, fill = fill))
            popViewport()
        }
        else {
            if (length(xaxis) > 1 || as.logical(xaxis[1])) {
                bottomblank <- 2
            }
            else {
                bottomblank <- 0
            }
        }
        SNPs <- SNP.gr[[i]]
        strand(SNPs) <- "*"
        SNPs <- sort(SNPs)
        scoreMax0 <- scoreMax <- if (length(SNPs$score) > 0) 
            ceiling(max(c(SNPs$score, 1), na.rm = TRUE))
        else 1
        if (type == "pie.stack") 
            scoreMax <- length(unique(SNPs$stack.factor))
        if (!type %in% c("pie", "pie.stack")) {
            if (length(yaxis) > 1 && is.numeric(yaxis)) {
                if (length(names(yaxis)) != length(yaxis)) {
                  names(yaxis) <- yaxis
                }
                scoreMax0 <- max(yaxis, scoreMax0)
            }
            if (scoreMax > 10) {
                SNPs$score <- 10 * SNPs$score/scoreMax
                scoreMax <- 10 * scoreMax0/scoreMax
            }
            else {
                scoreMax <- scoreMax0
            }
            scoreType <- if (length(SNPs$score) > 0) 
                all(floor(SNPs$score) == SNPs$score)
            else FALSE
        }
        else {
            scoreType <- FALSE
        }
        IsCaterpillar <- length(SNPs$SNPsideID) > 0
        if (IsCaterpillar) {
            if (any(is.na(SNPs$SNPsideID)) || !all(SNPs$SNPsideID %in% 
                c("top", "bottom"))) {
                warning("Not all SNPsideID is top or bottom")
                IsCaterpillar <- FALSE
            }
        }
        if (IsCaterpillar) {
            SNPs.top <- SNPs[SNPs$SNPsideID == "top"]
            SNPs.bottom <- SNPs[SNPs$SNPsideID == "bottom"]
        }
        else {
            SNPs.top <- SNPs
            SNPs.bottom <- GRanges()
        }
        if (length(SNPs.bottom) < 1) 
            IsCaterpillar <- FALSE
        if (!IsCaterpillar) {
            bottomblank <- bottomblank + 2
        }
        pushViewport(viewport(x = LINEW + 0.5, y = bottomblank * 
            LINEH/2 + 0.5, width = 1 - 7 * LINEW, height = 1 - 
            bottomblank * LINEH, xscale = c(start(ranges[i]), 
            end(ranges[i])), clip = "off"))
        plot.grid.xaxis <- function(col = "black") {

            print(col)
            if (length(xaxis) == 1 && as.logical(xaxis)) {
                grid.xaxis(gp = gpar(col = col))
            }
            if (length(xaxis) > 1 && is.numeric(xaxis)) {
                xaxisLabel <- names(xaxis)
                if (length(xaxisLabel) != length(xaxis)) 
                  xaxisLabel <- TRUE
                grid.xaxis(at = xaxis, label = xaxisLabel, gp = gpar(col = col))
            }
        }
        bottomHeight <- 0
        if (IsCaterpillar) {
            bottomHeight <- getHeight(SNPs = SNPs.bottom, ratio.yx = ratio.yx, 
                LINEW = LINEW, GAP = GAP, cex = cex, type = type, 
                scoreMax = scoreMax, level = "data&labels")
            vp <- viewport(y = bottomHeight, just = "bottom", 
                xscale = c(start(ranges[i]), end(ranges[i])))
            pushViewport(vp)
            plot.grid.xaxis("gray")
            popViewport()
        }
        else {
            plot.grid.xaxis()
        }
        baseline <- max(c(feature.splited[[1]]$height/2, 1e-04)) + 
            0.2 * LINEH
        baselineN <- max(c(feature.splited[[length(feature.splited)]]$height/2, 
            1e-04)) + 0.2 * LINEH
        feature.height <- plotFeatures(feature.splited, LINEH, 
            bottomHeight)
        if (length(SNPs.bottom) > 0) {
            custom_plotLollipops(SNPs.bottom, feature.height, bottomHeight, 
                baselineN, type, ranges[i], yaxis, scoreMax, 
                scoreMax0, scoreType, LINEW, cex, ratio.yx, GAP, 
                pin, dashline.col, side = "bottom", jitter = jitter)
        }
        feature.height <- feature.height + 2 * GAP
        if (length(SNPs.top) > 0) {
            custom_plotLollipops(SNPs.top, feature.height, bottomHeight, 
                baseline, type, ranges[i], yaxis, scoreMax, scoreMax0, 
                scoreType, LINEW, cex, ratio.yx, GAP, pin, dashline.col, 
                side = "top", jitter = jitter)
        }
        this.height <- getHeight(SNPs.top, ratio.yx, LINEW, GAP, 
            cex, type, scoreMax = scoreMax, level = "data&labels")
        this.height <- this.height + bottomHeight + feature.height
        this.height <- plotLegend(legend[[i]], this.height, LINEH)
        popViewport()
        this.height <- bottomblank * LINEH + this.height * (1 - 
            bottomblank * LINEH)
        vp <- viewport(x = 0.5, y = this.height * 0.5, width = 1, 
            height = this.height)
        pushViewport(vp)
        if (is.logical(ylab)) {
            if (ylab && length(names(SNP.gr)) > 0) {
                grid.text(names(SNP.gr)[i], x = LINEW, y = 0.5, 
                  rot = 90)
            }
        }
        if (is.character(ylab)) {
            if (length(ylab) == 1) 
                ylab <- rep(ylab, len)
            grid.text(ylab[i], x = LINEW, y = 0.5, rot = 90)
        }
        popViewport()
        popViewport()
        height0 <- height0 + this.height * height
    }
}


custom_plotLollipops = function(SNPs, feature.height, bottomHeight, baseline, 
                          type, ranges, yaxis, scoreMax, scoreMax0, scoreType,
                          LINEW, cex, ratio.yx, GAP, pin, dashline.col,
                          side=c("top", "bottom"), jitter=c("node", "label")){
    side <- match.arg(side)
    jitter <- match.arg(jitter)
    if(side=="top"){
        pushViewport(viewport(y=bottomHeight,
                              height=1,
                              just="bottom",
                              xscale=c(start(ranges), 
                                       end(ranges)),
                              clip="off"))
    }else{
        pushViewport(viewport(y=bottomHeight+feature.height,
                              height=1,
                              just="top",
                              xscale=c(start(ranges), 
                                       end(ranges)),
                              yscale=c(1, 0),
                              clip="off"))
    }
    if(type=="pie.stack" && length(SNPs$stack.factor)>0){
        stopifnot(is.vector(SNPs$stack.factor, mode="character"))
        if(length(SNPs$stack.factor.order)>0 || 
           length(SNPs$stack.factor.first)>0){
            warning("stack.factor.order and stack.factor.first are used by this function!",
                    "The values in these column will be removed.")
        }
        ## condense the SNPs
        stack.factors <- unique(as.character(SNPs$stack.factor))
        stack.factors <- sort(stack.factors)
        stack.factors.order <- 1:length(stack.factors)
        names(stack.factors.order) <- stack.factors
        SNPs <- SNPs[order(as.character(seqnames(SNPs)), start(SNPs), 
                           as.character(SNPs$stack.factor))]
        SNPs$stack.factor.order <- stack.factors.order[SNPs$stack.factor]
        SNPs$stack.factor.first <- !duplicated(SNPs)
        SNPs.condense <- SNPs
        SNPs.condense$oid <- 1:length(SNPs)
        SNPs.condense$factor <- paste(as.character(seqnames(SNPs)), start(SNPs), end(SNPs))
        SNPs.condense <- split(SNPs.condense, SNPs.condense$factor)
        SNPs.condense <- lapply(SNPs.condense, function(.ele){
            .oid <- .ele$oid
            .gr <- .ele[1]
            mcols(.gr) <- NULL
            .gr$oid <- NumericList(.oid)
            .gr
        })
        SNPs.condense <- unlist(GRangesList(SNPs.condense), use.names = FALSE)
        SNPs.condense <- sort(SNPs.condense)
        lab.pos.condense <- jitterLables(start(SNPs.condense), 
                                         xscale=c(start(ranges), end(ranges)), 
                                         lineW=LINEW*cex)
        lab.pos.condense <- reAdjustLabels(lab.pos.condense, 
                                           lineW=LINEW*cex)
        condense.ids <- SNPs.condense$oid
        lab.pos <- rep(lab.pos.condense, elementNROWS(condense.ids))
        lab.pos <- lab.pos[order(unlist(condense.ids))]
    }else{
        lab.pos <- jitterLables(start(SNPs), 
                                xscale=c(start(ranges), end(ranges)), 
                                lineW=LINEW*cex)
        lab.pos <- reAdjustLabels(lab.pos, 
                                  lineW=LINEW*cex)
    }
    
    if(length(SNPs)>0){
        yaxisat <- NULL
        yaxisLabel <- TRUE
        if(length(yaxis)>1 && is.numeric(yaxis)){
            yaxisat <- yaxis
            if(length(names(yaxis))==length(yaxis)) yaxisLabel <- names(yaxis)
            yaxis <- TRUE
        }
        if(yaxis && scoreMax>1 && !type %in% c("pie", "pie.stack")){
            if(side=="top"){
                grid.yaxis(at=yaxisat,
                           label=yaxisLabel,
                           vp=viewport(x=.5-LINEW,
                                       y=feature.height+5.25*GAP*cex+
                                           scoreMax*LINEW*ratio.yx/2*cex,
                                       width=1,
                                       height=scoreMax*LINEW*ratio.yx*cex,
                                       yscale=c(0, scoreMax0+.5)))
            }else{
                grid.yaxis(at=yaxisat,
                           label=yaxisLabel,
                           vp=viewport(x=.5-LINEW,
                                       y=1-(feature.height+5.25*GAP*cex+
                                           scoreMax*LINEW*ratio.yx/2*cex),
                                       width=1,
                                       height=scoreMax*LINEW*ratio.yx*cex,
                                       yscale=c(scoreMax0+.5, 0)))
            }
        }
        for(m in 1:length(SNPs)){
            this.dat <- SNPs[m]
            color <- if(is.list(this.dat$color)) this.dat$color[[1]] else this.dat$color
            border <- 
                if(is.list(this.dat$border)) this.dat$border[[1]] else this.dat$border
            fill <- if(is.list(this.dat$fill)) this.dat$fill[[1]] else this.dat$fill
            lwd <- if(is.list(this.dat$lwd)) this.dat$lwd[[1]] else this.dat$lwd
            id <- if(is.character(this.dat$label)) this.dat$label else NA
            id.col <- if(length(this.dat$label.col)>0) this.dat$label.col else "black"
            this.cex <- if(length(this.dat$cex)>0) this.dat$cex[[1]][1]*cex else cex
            this.dashline.col <- 
              if(length(this.dat$dashline.col)>0) this.dat$dashline.col[[1]][1] else dashline.col
            if(length(names(this.dat))<1) this.dashline.col <- NA
            this.dat.mcols <- mcols(this.dat)
            this.dat.mcols <- 
                this.dat.mcols[, 
                               !colnames(this.dat.mcols) %in% 
                                   c("color", "fill", "lwd", "id", 
                                     "cex", "dashline.col", 
                                     "id.col", "stack.factor", "SNPsideID"), 
                               drop=FALSE]
            if(type!="pie.stack"){
                this.dat.mcols <- 
                    this.dat.mcols[, !colnames(this.dat.mcols) %in% 
                                       c("stack.factor.order", 
                                         "stack.factor.first"), 
                                   drop=FALSE]
            }
            this.dat.mcols <- 
                this.dat.mcols[, !grepl("^label.parameter",
                                        colnames(this.dat.mcols)), 
                               drop=FALSE]

            grid.lollipop(x1=convertX(unit(start(this.dat), "native"), "npc", 
                                      valueOnly=TRUE),  
                          y1=baseline,
                          x2=convertX(unit(ifelse(jitter=="node", 
                                                  lab.pos[m], 
                                                  start(this.dat)), 
                                           "native"), "npc", valueOnly=TRUE), 
                          y2=feature.height,
                          y3=4*GAP*cex, y4=2.5*GAP*cex, 
                          radius=this.cex*LINEW/2,
                          col=color,
                          border=border,
                          percent=this.dat.mcols,
                          edges=100,
                          type=type,
                          ratio.yx=ratio.yx,
                          pin=pin,
                          scoreMax=(scoreMax-0.5) * LINEW * cex,
                          scoreType=scoreType,
                          id=id, id.col=id.col,
                          cex=this.cex, lwd=lwd, dashline.col=this.dashline.col,
                          side=side)

        }
        this.height <- getHeight(SNPs, 
                                 ratio.yx, LINEW, GAP, cex, type,
                                 scoreMax=scoreMax,
                                 level="data")
        labels.rot <- 90
        if(length(names(SNPs))>0){
            if(type=="pie.stack"){
                ## unique lab.pos and SNPs
                idx <- !duplicated(names(SNPs))
                lab.pos <- lab.pos[idx]
                SNPs <- SNPs[idx]
            }
            labels.x <- lab.pos
            labels.text <- names(SNPs)
            labels.just <- ifelse(side=="top", "left", "right")
            labels.hjust <- NULL
            labels.vjust <- NULL
            labels.check.overlap <- FALSE
            labels.default.units <- "native"
            labels.gp <- gpar(cex=cex)
            
            ## change the parameter by use definations.
            for(label.parameter in c("x", "y", "just", "hjust", "vjust",
                                     "rot", "check.overlap", "default.units",
                                     "gp")){
                label.para <- paste0("label.parameter.", label.parameter)
                if(label.para %in% colnames(mcols(SNPs))){
                    assign(paste0("labels.", label.parameter), 
                           mcols(SNPs)[, label.para])
                }
            }
            labels.gp <- c(labels.gp, cex=cex)
            labels.gp[duplicated(names(labels.gp))] <- NULL
            labels.gp <- do.call(gpar, labels.gp)
            if(jitter=="label"){
              ## add guide lines
              rased.height <- 4*GAP*cex
              guide.height <- 2.5*GAP*cex
              for(i in 1:length(SNPs)){
                this.dashline.col <- 
                  if(length(SNPs[i]$dashline.col)>0) 
                    SNPs[i]$dashline.col[[1]][1] else 
                      dashline.col
                if(length(names(SNPs[i]))<1) this.dashline.col <- NA
                grid.lines(x=c(start(SNPs[i]), labels.x[i]), 
                           y=c(this.height+feature.height-cex*LINEW, 
                               this.height+feature.height+rased.height),
                           default.units = labels.default.units,
                           gp=gpar(col=this.dashline.col, lty=3))
                grid.lines(x=c(labels.x[i], labels.x[i]),
                           y=c(this.height+rased.height+feature.height,
                               this.height+rased.height+
                                 guide.height+feature.height),
                           default.units = labels.default.units,
                           gp=gpar(col=this.dashline.col, lty=3))
              }
              ## add this height
              this.height <- this.height + rased.height + guide.height
            }
            grid.text(x=labels.x, y=this.height + feature.height, 
                      label = labels.text,  
                      just = labels.just, 
                      hjust = labels.hjust,
                      vjust = labels.vjust,
                      rot=labels.rot,
                      check.overlap = labels.check.overlap,
                      default.units = labels.default.units,
                      gp=labels.gp)
        }
    }
    popViewport()
}