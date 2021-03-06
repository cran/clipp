pedigree_loglikelihood_g <- function(fam_penet, geno_freq, trans, target_id = NULL,
                                     sibships = NULL, sibedges = NULL) {
  # Calculate the log-likelihood for the family fam, based on the penetrance
  # matrix penet, genotype frequencies geno_freq and transmission
  # probabilities trans.  The function returns the log-likelihood if target_id is null,
  # or the updated penetrance matrix otherwise.
  # Notes:
  # *   The only possible genotypes are 1, ..., length(geno_freq)
  # *   geno_freq[j] > 0 is the population genotype frequency of genotype j,
  # *   penet = matrix of non-negative numbers, rows corr to people in fam and
  #     columns to genotypes, usually no row is 0 (otherwise the likelihood is 0)
  # *   trans[3*gm+gf-3,go] = probability that the offspring of a mother and
  #     father with genotypes gm and gf has genotype go
  # *   The sibling graph is assumed to be a connected tree
  # *   fam is assumed to be a dataframe with certain variable names and modes
  # *   If target_id is NULL (the default) then the family's log-likelihood is
  #     returned, otherwise the updated penetrance matrix is returned (updated by
  #     collapsing the pedigree onto the target person, in the absence of loops)




  # Some functions

  calc.sibgraph <- function(fam) {
    # Calculate all of the nuclear families (sibships) and a graph showing how
    # they overlap (with vertices 1:length(sibships) and edges given by sibedges[,1:2])
    # except remove edges where the overlapping person has multiple marriages.
    # We just keep track of the edges (pairs of vertices),
    # since these contain the names of all of the vertices.  The one exception to
    # this is isolated vertices, so include an edge from each isolated vertex
    # to itself just to remember the vertex.
    # Also, there is an edge between two sibships for each person in common
    # (can get 2 or more edges in common through incest).
    # Lastly, singletons (i.e. people in fam not connected to anyone else) are
    # ignored in the sibship graph, and treated separately later.
    # Each element of sibships lists the mother, father then children of the
    # corresponding nuclear family.

    # Calculate the sibships (nuclear families)
    parents <- unique(fam[, c("mother", "father")])
    parents <- parents[parents$mother != 0 & parents$father != 0, ]
    parents <- as.matrix(parents)
    rownames(parents) <- NULL
    colnames(parents) <- NULL
    if (nrow(parents)==0) {
      # Deal with the case where the data purely consists of singletons
      return(list(sibships = vector("list", 0),
                  sibedges = data.frame(start=c(), end=c(), person=c()) ))
    }
    # parents[1:10,]
    f <- function(i) {
      p <- parents[i, ]
      keep <- fam$mother == p[1] & fam$father == p[2]
      # Assumes mother and father are in that order
      return(c(p, fam$indiv[keep]))
    }
    # Returns a list always (not a matrix if all sibships are the same size)
    sibships <- lapply(1:nrow(parents), f)

    # Calculate the sibship graph
    # WARNING:  Gives edges with a random order (orientation) if i is a founder
    # WARNING:  If i has multiple marriages then this function selects a
    # (fairly arbitrary) subtree of the corr. complete graph
    f1 <- function(i)  {length(sibships[[i]])}
    sib.length <- sapply(1:length(sibships),f1)
    person <- unlist(sibships)
    sibnum <- rep(1:length(sibships),times=sib.length)
    dup.person <- unique(person[duplicated(person)]) # dup.person = all IDs that appear in 2+ sibships
    f2 <- function(dp)  {unique(sibnum[person==dp])}
    slist <- lapply(dup.person,f2) # slist[[i]] = all sibships containing dup.person[i]
    f3 <- function(i) {
      # i <- 0; i <- i+1
      s <- slist[[i]]
      #if (length(s) < 2) {return(c())}
      g2 <- function(x) {
        dup.person[i] %in% x[3:length(x)]
      }
      j <- which(sapply(sibships, g2))
      # j is of length 0 (if i is a founder) or 1 (otherwise)
      if (length(s) == 2) m <- matrix(c(j, setdiff(s, j)), 1, 2)
      if (length(s) > 2) {
        #m <- as.matrix(expand.grid(s, s))
        #g3 <- function(x) {
        #  x[1] < x[2]
        #}
        #m <- m[apply(m, 1, g3), ]
        if (length(j) == 1) {
          m <- matrix(cbind(j, setdiff(s, j)), length(s) - 1, 2)
        }
        if (length(j) == 0) {
          m <- matrix(cbind(s[1], s[-1]), length(s) - 1, 2)
        }
      }
      m <- cbind(m, dup.person[i])
      colnames(m) <- NULL
      m
      return(m)
    }
    sibedges <- matrix(numeric(0), 0, 3)
    if (length(slist) != 0) {
      st <- lapply(1:length(slist), f3)
      for (j in 1:length(st)) sibedges <- rbind(sibedges, st[[j]])
    }

    # Add a self-edge for each isolated sibship
    iso <- setdiff(1:length(sibships), c(sibedges[,1], sibedges[,2]))
    if (length(iso) > 0)  {
      x <- c(iso, iso, rep(-1,length(iso)))
      sibedges <- rbind(sibedges, matrix(x, length(iso), 3))
    }

    # Clean up and output
    sibedges <- as.data.frame(sibedges)
    names(sibedges) <- c("start", "end", "person")
    row.names(sibedges) <- 1:nrow(sibedges)
    return(list(sibships = sibships, sibedges = sibedges))
  }

  update.penet <- function(ll, attachID) {
    # Note that this function also depends on sibships, penet, geno_freq and trans
    ngeno <- length(geno_freq)
    kidIDs <- sibships[[ll]][-(1:2)]
    parentIDs <- sibships[[ll]][1:2]
    attach.is.parent <- (attachID %in% sibships[[ll]][1:2])
    if (attach.is.parent) {
      othparID <- setdiff(parentIDs, attachID)
      summand.p <- function(gatt, goth) {
        if (parentIDs[1] == attachID) {
          gm <- gatt
          gf <- goth
        } else {
          gf <- gatt
          gm <- goth
        }
        prob <- 1
        for (kID in kidIDs) {
          p <- penet[kID, ]
          prob <- prob * sum(p * trans[ngeno * gm + gf - ngeno, ] )
        }
        p <- penet[othparID, ]
        prob <- prob * p[goth] * geno_freq[goth]
        return(prob)
      }
      f1 <- function(gatt) {
        f2 <- function(goth) {
          summand.p(gatt, goth)
        }
        p <- penet[attachID, ]
        return(p[gatt] * sum(sapply(1:ngeno, f2)))
      }
      newpen <- sapply(1:ngeno, f1)
    }
    if (!attach.is.parent) {
      summand.c <- function(gatt, gm, gf) {
        prob <- 1
        for (kID in setdiff(kidIDs, attachID)) {
          p <- penet[kID, ]
          prob <- prob * sum(p * trans[ngeno * gm + gf - ngeno, ] )
        }
        pm <- penet[parentIDs[1], ]
        pf <- penet[parentIDs[2], ]
        prob <- prob * trans[ngeno * gm + gf - ngeno, gatt] * pm[gm] * pf[gf] *
          geno_freq[gm] * geno_freq[gf]
        return(prob)
      }
      f1 <- function(gatt) {
        f2 <- function(gpar) {
          summand.c(gatt, gpar[1], gpar[2])
        }
        EG <- expand.grid(1:ngeno, 1:ngeno)
        p <- penet[attachID, gatt]
        #if (!last.iter)  p <- p / geno_freq[gatt]
        p <- p / geno_freq[gatt]
        # return(p * sum(sapply(EG,f2)))
        return(p * sum(apply(EG, 1, f2)))
      }
      newpen <- sapply(1:ngeno, f1)
    }
    return(newpen)
  }

  choose.sibship <- function(sibships,sibedges) {
    # Choose a sibship to collapse, or return NA if there is no suitable one
    # If there are no others then choose target.sibship
    if (!is.null(target_id) & nrow(sibedges)==1 && (sibedges$start==target.sibship &
        sibedges$end==target.sibship))  return(target.sibship)
    # If there are any isolated leaves then choose one of these
    isoleaves <- unique(sibedges$end[sibedges$start==sibedges$end])
    if (!is.null(target_id))  isoleaves <- setdiff(isoleaves, target.sibship)
    if (length(isoleaves) > 0) {return(isoleaves[1])}
    # Otherwise
    s <- c(sibedges$start, sibedges$end)
    leaves <- setdiff(s, s[duplicated(s)])
    if (!is.null(target_id))  leaves <- setdiff(leaves, target.sibship)
    leavesdown <- intersect(leaves, sibedges$end)
    if (length(leavesdown) > 0) leaves <- leavesdown
    if (length(leaves) > 0) {ll <- leaves[1]} else {ll <- NA}
    return(ll)
  }



  # Perform the calculation

  # Unpack the data
  fam <- fam_penet$fam
  penet <- fam_penet$penet

  # To help prevent floating point underflow, divide each row of penet by a factor so
  # that the max of each row is 1, and calculate a corresponding adjustment to
  # the log-likelihood
  ll.max <- apply(penet, 1, max)
  if (any(ll.max == 0)) {
    if (is.null(target_id)) {return(-Inf)} else {return(rep(NA,length(geno_freq)))}
  }
  penet <- penet / ll.max
  loglik <- sum(log(ll.max))

  # Calculate the sibship graph
  if (is.null(sibships)) {
    cs <- calc.sibgraph(fam)
    sibships <- cs$sibships
    sibedges <- cs$sibedges
  }

  # Identify singletons (people without parents or offspring)
  # Change their penetrance so they don't contribute in each recursion
  singletons <- setdiff(fam$indiv, unlist(sibships))
  if (length(singletons) > 0) {
    f1 <- function(id)  {log(sum(geno_freq * penet[id, ]))}
    loglik <- loglik + sum(sapply(singletons, f1))
    #penet[singletons, ] <- rep(1, length(singletons) * length(geno_freq))
    penet[singletons, ] <- 1
  }

  # Choose a sibship that contains target_id, to be collapsed last
  if (!is.null(target_id)) {
    if (target_id %in% singletons) {
      # Deal with the case when target_id is a singleton
      p <- geno_freq * penet[target_id,]
      return(p/sum(p))
    } else {
      # Choose a sibship that contains target_id (this will be collapsed last)
      f2 <- function(x)  is.element(target_id, x)
      target.sibship <- which(sapply(sibships, f2))[1]
    }
  }

  # Collapse leaves one at a time until there are no more leaves left except
  # perhaps target.sibship (if this is non-null and not isolated), where
  # a leaf is a vertex that connects to 0 or 1 other vertices in the sibling graph.
  # Note that ll will always be chosen to be a leaf or NA.
  siborder <- c()
  while (!is.na(ll <- choose.sibship(sibships,sibedges))) {
    siborder <- c(siborder, ll)

    # Calculate the person attachID who we'll collapse sibship ll onto
    # and update sibedges and singletons
    if (any(sibedges$start == ll & sibedges$end == ll)) {
      # If sibship ll is an isolated sibship (indicated by a self-edge in sibedges)
      sibedges <- sibedges[!(sibedges$start == ll & sibedges$end == ll), ]
      attachID <- sibships[[ll]][3]
      if (!is.null(target_id) && target_id %in% sibships[[ll]]) {
        attachID <- target_id
      }
      singletons <- c(singletons, attachID)
    } else {
      # If sibship ll connects to one other sibship
      kk <- which(sibedges$start == ll | sibedges$end == ll)[1]
      mm <- setdiff(as.numeric(sibedges[kk,1:2]), ll)
      #attachID <- intersect(sibships[[ll]], sibships[[mm]])[1]
      attachID <- sibedges$person[kk]
      sibedges <- sibedges[-kk, ]
      if (!any(sibedges$start == mm | sibedges$end == mm)) {
        # If sibship mm is now isolated then add a self-edge
        sibedges <- rbind(sibedges, data.frame(start=mm, end=mm, person= -1))
      }
    }

    # Update the penetrance matrix penet
    newpen <- update.penet(ll, attachID)
    # Slow for large pedigrees, copy contents of update.penet above? ###########################################
    if (max(newpen) == 0) {
      if (is.null(target_id)) {return(-Inf)} else {return(rep(NA,length(geno_freq)))}
    }
    loglik <- loglik + log(max(newpen))
    newpen <- newpen / max(newpen)
    penet[attachID, ] <- newpen

  }

  # Deal with singletons
  if (is.null(target_id) & length(singletons) > 0) {
    f3 <- function(id)  {log(sum(geno_freq * penet[id, ]))}
    loglik <- loglik + sum(sapply(singletons, f3))
  }

  # If loops remain in the pedigree then use recursion
  if (nrow(sibedges)==0) {

    # The calculation has finished so output
    if (is.null(target_id)) {
      return(loglik)
    } else {
      p <- geno_freq * penet[target_id,]
      return(p/sum(p))
    }

  } else {

    # There are loops, so break the pedigree and use recursion

    # Deal with genotype probs first -- dependable but not the fastest approach
    # A faster version, though with bugs, is in clipp_20210324_genoprob_notwork
    if (!is.null(target_id)) {
      f4 <- function(i) {
        peneti <- penet
        peneti[fam$indiv==target_id,-i] <- 0
        fam_penet <- list(fam=fam, penet=peneti)
        pedigree_loglikelihood_g(fam_penet, geno_freq, trans, NULL, sibships, sibedges)
      }
      l <- sapply(1:length(geno_freq), f4)
      l <- l - max(l)
      p <- exp(l)
      return(p/sum(p))
    }

    # Choose the best person to split in two, to break a loop of the pedigree.
    # The sibling graph is loops, so cutting any edge will reduce the graph's
    # first Betti number.
    #    * Choose an  edge in the longest arc, since then the longest arc will
    #      be contracted the smallest number of times in the recursion.
    #    * Then choose an  edge in that arc whose corresponding person in the
    #      pedigree has the smallest number of possible genotypes, since we will
    #      need to sum over these in the recursion.

    # Make a copy of sibedges that We might alter slightly below
    # The rows of se will always correspond to the rows of sibedges
    # but some of the entries might be changed to disconnect edges in se
    se <- sibedges
    arclength <- rep(0,nrow(se))  # arclength[i] will be the length of the arc containing se[i,]

    # Split se at points where 3 or more edges join
    vert <- c(se$start,se$end)
    vert <- vert[duplicated(vert)]
    vert <- vert[duplicated(vert)]
    vert <- unique(vert)
    Ns <- sum(se$start %in% vert)
    Ne <- sum(se$end %in% vert)
    if (Ns > 0)  se$start[se$start %in% vert] <- seq(-1,-Ns,by=-1)
    if (Ne > 0)  se$end[se$end %in% vert] <- seq(-1-Ns,-Ne-Ns,by=-1)

    # Calculate lengths of the arcs, i.e. lengths of connected components of se
    while (any(arclength==0)) {
      arc <- which(arclength==0)[1]
      nold <- 0
      while (nold < length(arc)) {
        nold <- length(arc)
        vert <- unique(c(se$start[arc],se$end[arc]))
        arc <- which(se$start %in% vert | se$end %in% vert)
        arc <- unique(arc)
      }
      arclength[arc] <- length(arc)
    }
    # Swap back to sibedges now that we've calculated the arc lengths

    # Find the edge in the longest arc(s) whose person has the fewest possible genotypes
    f5 <- function(i) {
      breakID <- sibedges$person[i]
      return(sum(penet[breakID,] != 0))
    }
    npossgeno <- sapply(1:nrow(sibedges), f5)
    # We could save time by restricting the above calculation to the longest arc

    # Select an edge whose corresponding person is good for splitting the pedigree at
    good <- arclength==max(arclength)
    npossgeno[!good] <- max(npossgeno)
    good <- good & npossgeno==min(npossgeno)
    break.edge <- which(good)[1]
    breakID <- sibedges$person[break.edge]

    # Break a loop of fam by splitting breakID into 2 people
    # All edges corresponding to breakID will have exactly one vertex in common,
    # and this will correspond to the sibship where breakID is a child (if this exists).
    # Change the ID of breakID in any other sibship than the one corresponding to
    # this common vertex (since that would disconnect too much of the graph)
    keep <- sibedges$person==breakID
    sibs <- c(sibedges$start[keep],sibedges$end[keep])
    break.sib <- sibedges$end[break.edge]
    if (sum(sibs==break.sib) > 1)  break.sib <- sibedges$start[break.edge]
    # Change instances of breakID in break.sib, noting that
    # breakID is a parent in break.sib, by construction
    kids <- sibships[[break.sib]][-(1:2)]
    doppleganger <- max(fam$indiv) + 1
    fam$mother[fam$indiv %in% kids & fam$mother==breakID] <- doppleganger
    fam$father[fam$indiv %in% kids & fam$father==breakID] <- doppleganger
    # Add rows for breakID's doppleganger
    penet <- rbind(penet, 1/geno_freq)
    fam <- rbind(fam, fam[1,])
    fam$indiv[nrow(fam)] <- doppleganger
    fam$mother[nrow(fam)] <- 0
    fam$father[nrow(fam)] <- 0
    # Update sibships and sibedges
    s <- sibships[[break.sib]]
    s[s==breakID] <- doppleganger
    sibships[[break.sib]] <- s
    sibedges <- sibedges[-break.edge,]

    # Perform the recursion
    possgeno <- which(penet[fam$indiv==breakID,] != 0)
    if (length(possgeno) == 0) {return(-Inf)}
    f6 <- function(i) {
      peneti <- penet
      #peneti[fam$indiv %in% c(breakID,doppleganger),-i] <- 0
      peneti[fam$indiv==breakID,-possgeno[i]] <- 0
      peneti[fam$indiv==doppleganger,-possgeno[i]] <- 0
      fam_penet <- list(fam=fam, penet=peneti)
      pedigree_loglikelihood_g(fam_penet, geno_freq, trans, target_id, sibships, sibedges)
    }
    lli <- sapply(1:length(possgeno), f6)
    loglik <- loglik + log(sum(exp(lli)))
    return(loglik)

  }

}

