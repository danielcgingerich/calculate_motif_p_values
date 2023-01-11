

convertGRangesToSeq <-
  function(ranges, genome, strand){
    
    chr <- as.character(seqnames(ranges))
    
    strt <- start(ranges)
    
    nd <- end(ranges)
    
    getSeq(genome, names = chr, start = strt, 
           end = nd, strand = strand)
  }
createSeqMatrix <-
  function(motif.ranges, mcol.seqs, add.snps, add.revmap){
    message('creating matrix')
    
    biostrings <- mcols(motif.ranges)[[mcol.seqs]]
    
    mtx <- as.character(biostrings) %>% 
      strsplit(split = '') %>% 
      unlist %>% matrix(nrow = length(biostrings), byrow = T)
    
    if (add.revmap){
      rownames(mtx) <- as.character(motif.ranges$revmap)
    }
    
    if (add.snps){
      message('adding snps')
      mtx[,median(1:ncol(mtx))] <- motif.ranges$alt
    }
    
    return(mtx)
  }
reverseSeqMatrix <-
  function(seq.matrix){
    
    mtx <- seq.matrix[, ncol(seq.matrix):1]
    
    comp.bases <- c('A' = 'T', 'T' = 'A', 'G' = 'C', 'C' = 'G')
    
    pos <- match(mtx, comp.bases)
    
    matrix(names(comp.bases)[pos], 
           nrow = nrow(mtx), 
           byrow = F, 
           dimnames = list(rownames(mtx), colnames(mtx)) )
  }
findAltMotifs <-
  function(motif.ranges, fwd.seq, rev.seq, pwm.list, 
           icm.list, scan.with, score.with){
    
    motif <- unique(motif.ranges$motif)
    
    message(motif)
    
    # preparation for scan
    
    pwm <- pwm.list[[motif]]@profileMatrix
    
    icm <- icm.list[[motif]]@profileMatrix
    
    colnames(pwm) <- colnames(icm) <- as.character(1:ncol(pwm))
    
    list.index <- unique(motif.ranges$motif.width) %>% as.character
    
    fwd.mtx <- fwd.seq[[list.index]][motif.ranges$revmap, ]
    
    rev.mtx <- rev.seq[[list.index]][motif.ranges$revmap, ]
    
    slider.start <- 1:ncol(pwm)
    
    slider.end <- slider.start + ncol(pwm) - 1
    
    # scan sequence matrix
    
    pwmscore.fwd <- pwmscore.rev <- 
      icmscore.fwd <- icmscore.rev <- 
      data.frame(row.names = motif.ranges$revmap)
    
    message('scanning seq matrices...')
    
    for (i in 1:ncol(pwm)){
      fwd.slide <- fwd.mtx[, slider.start[i]:slider.end[i]]
      rev.slide <- rev.mtx[, slider.start[i]:slider.end[i]]
      colnames(fwd.slide) <- colnames(rev.slide) <- 1:ncol(fwd.slide)
      fwd.slide <- melt(fwd.slide)
      rev.slide <- melt(rev.slide)
      fwd.index <- matrix(c(fwd.slide$value, fwd.slide$Var2), 
                          nrow = nrow(fwd.slide), byrow = F)
      rev.index <- matrix(c(rev.slide$value, rev.slide$Var2), 
                          nrow = nrow(rev.slide), byrow = F)
      pwmscore.fwd[,i] <- matrix(pwm[fwd.index], ncol = ncol(pwm), byrow = F) %>% rowSums
      pwmscore.rev[,i] <- matrix(pwm[rev.index], ncol = ncol(pwm), byrow = F) %>% rowSums
      icmscore.fwd[,i] <- matrix(icm[fwd.index], ncol = ncol(icm), byrow = F) %>% rowSums
      icmscore.rev[,i] <- matrix(icm[rev.index], ncol = ncol(icm), byrow = F) %>% rowSums
    }
    
    # extract best alt position
    
    colnames(pwmscore.fwd) <- colnames(icmscore.fwd) <- 
      paste0(1:ncol(pwm), '.fwd')
    
    colnames(pwmscore.rev) <- colnames(icmscore.rev) <- 
      paste0(ncol(pwm):1, '.rev')
    
    pwmscore <- cbind(pwmscore.fwd, pwmscore.rev)
    
    icmscore <- cbind(icmscore.fwd, icmscore.rev)
    
    pos <- switch(scan.with, 
                  'pwm' = apply(pwmscore, 1, which.max), 
                  'icm' = apply(icmscore, 1, which.max))
    
    pos.shift <- new.strnd <- colnames(pwmscore)[pos]
    
    new.strnd <- ifelse(grepl('fwd', new.strnd), '+', '-')
    
    pos.shift <- gsub('.fwd|.rev', '', pos.shift) %>% as.numeric
    
    mtf.start <- motif.ranges$seq.window.start + pos.shift - 1
    
    mtf.end <- mtf.start + ncol(pwm) - 1
    
    mtf.chr <- as.character(seqnames(motif.ranges))
    
    mtf.gr <- paste(mtf.chr, mtf.start, mtf.end, sep = '-')
    
    motif.ranges$old.position <- GRangesToString(motif.ranges)
    
    motif.ranges$new.position <- mtf.gr
    
    motif.ranges$old.strnd <- strand(motif.ranges)
    
    motif.ranges$new.strnd <- new.strnd
    
    # extract alt position scores
    
    score.index <- matrix(c(1:nrow(pwmscore), pos), ncol = 2, byrow = F)
    
    alt.score <- switch(score.with, 
                        'pwm' = pwmscore[score.index], 
                        'icm' = icmscore[score.index])
    
    motif.ranges$alt.score <- alt.score
    
    return(motif.ranges)
  }
scoreRefSeqs <-
  function(motif.ranges, ref.seq, score.with){
    
    motif <- unique(motif.ranges$motif)
    
    message(motif)
    
    pwm <- switch(score.with, 
                  'pwm' = pwm.list[[motif]]@profileMatrix, 
                  'icm' = icm.list[[motif]]@profileMatrix)
    
    colnames(pwm) <- 1:ncol(pwm)
    
    list.index <- as.character(unique(motif.ranges$motif.width))
    
    ref.mtx <- ref.seq[[list.index]][ motif.ranges$revmap, ]
    
    colnames(ref.mtx) <- as.character(1:ncol(ref.mtx))
    
    ref.index <- melt(ref.mtx)
    
    ref.index <- matrix(c(ref.index$value, ref.index$Var2), 
                        nrow = nrow(ref.index), byrow = F)
    
    ref.pssm <- matrix(pwm[ref.index], ncol = ncol(pwm), byrow = F, 
                       dimnames = list(rownames(ref.mtx), colnames(pwm)) )
    
    ref.scores <- rowSums(ref.pssm)
    
    motif.ranges$ref.score <- ref.scores
    
    return(motif.ranges)
  }
insertSNPs <-
  function(sequences, offset, snps){
    
    mtx <- matrix(F, ncol = width(sequences), nrow = length(sequences))
    
    mtx2 <- matrix(c(1:nrow(mtx), offset), ncol = 2)
    
    mtx[mtx2] <- T
    
    replaceLetterAt(x = sequences, at = mtx, letter = snps)
  }
scaleMotifScores <-
  function(motif.ranges, mcol.scores, score.type, pwm.list, icm.list){
    
    motif <- unique(motif.ranges$motif)
    
    message(motif)
    
    pwm <- switch(score.type, 
                  'pwm' = pwm.list[[motif]]@profileMatrix, 
                  'icm' = icm.list[[motif]]@profileMatrix)
    
    max.score <- sum(apply(pwm, 2, max))
    
    min.score <- sum(apply(pwm, 2, min))
    
    score <- mcols(motif.ranges)[, mcol.scores]
    
    scaled.score <- (score - min.score) / (max.score - min.score)
    
    return(scaled.score)
  }



splitByCols <- function(mtx){
  lapply(1:ncol(mtx), function(x){
    mtx[,x]
  })
}


# findThreshold -----------------------------------------------------------

findThreshold <- function(mtx.cols){
  max.score <- sapply(mtx.cols, max) %>% sum
  min.score <- sapply(mtx.cols, min) %>% sum
  0.8*(max.score-min.score)+min.score
}


# expandGridRowSums -------------------------------------------------------

expandGridRowSums <- function(min.possible.score, ATCG, previous.scores, threshold){
  A <- previous.scores[previous.scores >= min.possible.score - ATCG['A'] ] + ATCG['A']
  t <- previous.scores[previous.scores >= min.possible.score - ATCG['T'] ] + ATCG['T']
  C <- previous.scores[previous.scores >= min.possible.score - ATCG['C'] ] + ATCG['C']
  G <- previous.scores[previous.scores >= min.possible.score - ATCG['G'] ] + ATCG['G']
  unname(c(A, t, C, G))
}


# evenSplit ---------------------------------------------------------------

evenSplit <- function(obj, n){
  seq <- seq_along(obj)
  x <- length(obj)/n
  seq <- floor(seq/x)
  split(obj, f = seq)
}


# upperTail ---------------------------------------------------------------

# find right tail discrete probability distribution of PWM

upperTail <- function(pwm, thresh){
  seq.list <- splitByCols(pwm)
  max.vals <- sapply(seq.list, max)
  max.cum <- cumsum(max.vals[length(max.vals):1])
  max.cum <- max.cum[length(max.cum):1]
  ps.mins <- thresh - max.cum                                             # smallest value previous sequence can have and still reach threshold
  message('constructing upper tail distribution')
  upper.tail <- seq.list[[1]]
  for (i in 2:length(seq.list)){
    message(i)
    upper.tail <- expandGridRowSums(min.possible.score = ps.mins[[i]], 
                                    ATCG = seq.list[[i]], 
                                    previous.scores = upper.tail, 
                                    threshold = thresh)
  }
  upper.tail <- sort(upper.tail)
  return(upper.tail)
}

# closestMatch ------------------------------------------------------------

closestMatch <- function(val, vec){
  delta <- abs(vec - val)
  hit <- vec[which.min(delta)]
  d <- delta[which.min(delta)]
  message('difference: ', d)
  return(list('hit' = hit, 'difference' = d))
}
# motifPvalue -------------------------------------------------------------

motifPvalue <- 
  function(motif.ranges, pwm.list, mcol.scores){
    motif <- unique(motif.ranges$motif)
    message(motif)
    pwm <- pwm.list[[motif]]@profileMatrix
    scores <- mcols(motif.ranges)[, mcol.scores]
    thresh <- min(scores)
    
    upper.tail <- upperTail(pwm, thresh)                                  
    
    df <- data.frame('score' = scores)
    df$p.val <- 1
    
    message('account for decimal float point')
    df$score <- sapply(df$score, 
                       function(x){
                         closestMatch(x, vec = upper.tail)$hit
                       })                                               # decimal float point not accurate when evaluating scores.  must find closest match
    
    upper.tail <- evenSplit(upper.tail, n = 10)                           # most efficient to divide upper tail into chunks, then
    intervals <- sapply(upper.tail, min)                                  # use intervals to index each score to an upper tail chunk.
    intervals[1] <- -Inf                                                  # first interval must be open, or low scores may give error
    df$interval <- findInterval(x = df$score, vec = intervals)            # done this way, only chunk values need to be evaluated  by boolean operator
    top.shelf <- length(upper.tail)                                         # rathan than the whole upper tail. 
    
    message('finding p values')                                           # although confusing, this is fastest way to find p values
    for (i in 1:length(motif.ranges)){                                    # uppertail is split into chunks.  for each given score, the 
      index <- df$interval[i]                                             # total scores in its respective chunk >= score, and length of
      chunk <- upper.tail[[ index ]]                                      # the above chunks are summed.  This way, only the chunk values have to 
      if (index != top.shelf){                                            # be evaluated by boolean operator, not the whole upper tail
        chunks.above <- upper.tail[ (index+1):length(upper.tail)]
        n.above <- sum(sapply(chunks.above, length))
        n.chunk <- sum(chunk >= df$score[i] )
        df$p.val[i] <- n.above + n.chunk
      }                                                                   
      if (index == top.shelf){                                                
        n.chunk <- sum(chunk >= df$score[i] )                                       
        df$p.val[i] <- n.chunk                                                  
      }
    }
    
    n <- 4^ncol(pwm)
    df$p.val <- df$p.val/n
    mcols(motif.ranges)[, paste0(mcol.scores, '.p.value')] <- df$p.val
    
    return(motif.ranges)
  }


bitCrush <- function (v, N){ # v is the input vector (sorted), and keep every N sample
  seed <- c(TRUE, rep(FALSE, N-1) )
  cont <- rep(seed, ceiling( length(v)/N ) )[1:length(v)]
  return(v[which(cont)]) 
}

# combineTailScores ------------------------------------------------------
combineTailScores <- function(prefix, suffix, thresh){
  scores <- lapply(prefix, 
                   function(x){
                     y <-  x + suffix
                     y[y >= thresh - 0.001]
                   })
  unlist(scores)
}

# bigMotifPvalue ----------------------------------------------------------

bigMotifPvalue <- function(motif.ranges, pwm.list, downsample){
  
  # prepare pwm
  mstr <- motif.ranges
  motif <- unique(mstr$motif)
  message(motif)
  pwm <- pwm.list[[motif]]@profileMatrix
  ref.scores <- mstr$ref.score
  alt.scores <- mstr$alt.score
  thresh <- min(c(ref.scores, alt.scores))-0.01
  
  # split pwm in half                                                     # we divide the sequence into prefix and suffix
  n <- round(ncol(pwm)/2)
  pwm1 <- pwm[, 1:n]                                                      # prefix (first half)
  pwm2 <- pwm[, (n+1):ncol(pwm)]                                          # suffix (second half)
  
  # prefix pwm sums
  message('calculating prefix scores')
  thresh1 <- thresh - sum(apply(pwm2, 2, max))                            # calculate upper tail of prefix, assuming max possible score suffix
  upper.tail1 <- upperTail(pwm = pwm1, thresh = thresh1)                  
  
  # right pwm sums
  message('calculating suffix scores')
  thresh2 <- thresh - sum(apply(pwm1, 2, max))                            # calculate upper tail of suffix, assuming max possible score prefix
  upper.tail2 <- upperTail(pwm = pwm2, thresh = thresh2)
  
  # downsample each tail
  message('sorting upper tails for downsampling')                         # downsample each half before combining distributions
  upper.tail1 <- sort(upper.tail1)                                        # in order to downsample accurately, the tails must be ordered
  upper.tail2 <- sort(upper.tail2)
  message('downsampling upper tails by factor of ', downsample)
  upper.tail1 <- bitCrush(upper.tail1, downsample)
  upper.tail2 <- bitCrush(upper.tail2, downsample)
  
  
  # combine pwm sums for distribution upper tail
  upper.tail1 <- evenSplit(upper.tail1, n = 100)                          # step 1: split each tail into 100 equal parts
  upper.tail2 <- evenSplit(upper.tail2, n = 100)
  max1 <- sapply(upper.tail1, max)                                        # step 2: take the max of each part. When we look at max1 and max2,
  max2 <- sapply(upper.tail2, max)                                        # this tells us which combos are feasible for reaching threshold score, 
  pair.index <- expand.grid('max1' = names(max1), 'max2' = names(max2))   # because if max1 + max2 dont reach thresh, no other score will
  pair.index <- apply(pair.index, 2, as.character) %>% as.data.frame
  pairs <- expand.grid('max1' = max1, 'max2' = max2)
  keep <- which(rowSums(pairs) >= thresh)                                 # we only need to expand the grid of max pairs that reach threshold
  pair.index <- pair.index[keep, ]
  
  message('combining upper tails: ', nrow(pair.index), ' total')
  upper.tail <- list()
  for (j in 1:nrow(pair.index)){
    message(j)
    i1 <- pair.index$max1[j]
    i2 <- pair.index$max2[j]
    ut1 <- upper.tail1[[i1]]
    ut2 <- upper.tail2[[i2]]
    ut <- combineTailScores(prefix = ut1, suffix = ut2, thresh = thresh)
    ut <- ut[ut >= thresh]
    upper.tail[[j]] <- ut
  }
  message('done! sorting upper tails')
  upper.tail <- unlist(upper.tail)
  upper.tail <- sort(upper.tail)
  
  # obtain p values
  df <- data.frame('ref.score' = mstr$ref.score, 
                   'alt.score' = mstr$alt.score)
  df$ref.score.p.val <- 1
  df$alt.score.p.val <- 1
  message('splitting upper tail')
  upper.tail <- evenSplit(upper.tail, n = 10000)
  intervals <- sapply(upper.tail, min) %>% unname
  intervals[1] <- -Inf
  df$ref.interval <- findInterval(x = df$ref.score, vec = intervals)
  df$alt.interval <- findInterval(x = df$alt.score, vec = intervals)
  top.shelf <- length(upper.tail)
  
  
  message('finding ', length(motif.ranges), ' total p values')          
  for (i in 1:length(motif.ranges)){   
    # ref scores
    ref.index <- df$ref.interval[i]                                             
    ref.chunk <- upper.tail[[ ref.index ]]                                      
    if (ref.index != top.shelf){                                            
      ref.chunks.above <- upper.tail[ (ref.index+1):length(upper.tail)]
      n.ref.above <- sum(sapply(ref.chunks.above, length)) * downsample^2       # account for downsampling of both tails
      n.ref.chunk <- sum(ref.chunk >= df$ref.score[i] ) * downsample^2          # account for downsampling of both tails
      df$ref.p.val[i] <- n.ref.above + n.ref.chunk
    }                                                                   
    if (ref.index == top.shelf){                                                
      n.ref.chunk <- sum(ref.chunk >= df$ref.score[i] ) * downsample^2          # account for downsampling of both tails                                      
      df$ref.p.val[i] <- n.ref.chunk                                                  
    }    
    
    # alt scores
    alt.index <- df$alt.interval[i]                                             
    alt.chunk <- upper.tail[[ alt.index ]]                                      
    if (alt.index != top.shelf){                                            
      alt.chunks.above <- upper.tail[ (alt.index+1):length(upper.tail)]
      n.alt.above <- sum(sapply(alt.chunks.above, length)) * downsample^2       # account for downsampling of both tails
      n.alt.chunk <- sum(alt.chunk >= df$alt.score[i] ) * downsample^2          # account for downsampling of both tails
      df$alt.p.val[i] <- n.alt.above + n.alt.chunk
    }                                                                   
    if (alt.index == top.shelf){                                                
      n.alt.chunk <- sum(alt.chunk >= df$alt.score[i] ) * downsample^2          # account for downsampling of both tails                                      
      df$alt.p.val[i] <- n.alt.chunk                                                  
    }    
    
  }
  
  ###
  n <- 4^ncol(pwm)
  df$ref.score.p.val <- df$ref.p.val/n
  df$alt.score.p.val <- df$alt.p.val/n
  
  mstr$ref.score.p.value <- df$ref.score.p.val
  mstr$alt.score.p.value <- df$alt.score.p.val
  return(mstr)
}
