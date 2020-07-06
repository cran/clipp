pedigree_loglikelihood_g <- function(fam_penet, geno_freq, trans, cl, target_id = NULL) {
  # Calculate the log-likelihood for the family fam, based on the penetrance
  # matrix penet, genotype frequencies geno_freq and transmission
  # probabilities trans.
  # Notes:
  # *   The only possible genotypes are 1, ..., length(geno_freq)
  # *   geno_freq[j] > 0 is the population genotype frequency of genotype j,
  # *   penet = matrix of non-negative numbers, rows corr to people in fam and
  #     columns to genotypes,usually no row is 0 (otherwise the likelihood is 0)
  # *   trans[3*gm+gf-3,go] = probability that the offspring of a mother and
  #     father with genotypes gm and gf has genotype go
  # *   The sibling graph is assumed to be a connected tree
  # *   fam is assumed to be a dataframe with certain variable names and modes
  # *   If target_id is NULL (the default) then the family's log-likelihood is
  #     returned, otherwise a vector of genotype probabilities for the person with
  #     fam$ID==target_id is returned

  fam <- fam_penet$fam
  penet <- fam_penet$penet

  # Deal with singleton families first
  if (nrow(fam) == 1) {
    if (is.null(target_id)) {
      return(log(sum(geno_freq * penet[1, ])))
    } else {
      p <- geno_freq * penet[1, ]
      p <- p / sum(p)
      return(p)
    }
  }

  # Perform the penetrance calculations for everyone in the family.  To help
  # prevent floating point underflow, divide each row by a factor so that the
  # max of each row is 1, and calculate a corresponding adjustment to
  # the log-likelihood
  penmat <- penet
  ll.max <- apply(penmat, 1, max)
  if (any(ll.max == 0)) {
    return(-Inf)
  }
  penmat <- penmat / ll.max
  ll.adj <- sum(log(ll.max))

  calc.sibtree <- function(fam) {
    # Calculate all of the nuclear families (sibships) and a graph showing how
    # they overlap (sibtree) except remove edges where the overlapping person
    # has multiple marriages

    # Calculate the sibships (nuclear families)
    parents <- unique(fam[, c("mother", "father")])
    parents <- parents[parents$mother != 0 & parents$father != 0, ]
    parents <- as.matrix(parents)
    rownames(parents) <- NULL
    colnames(parents) <- NULL
    # parents[1:10,]
    f <- function(i) {
      p <- parents[i, ]
      keep <- fam$mother == p[1] & fam$father == p[2]
      # Assumes mother and father are in that order
      return(c(p, fam$indiv[keep]))
    }
    # It returns a list always (not a matrix if all sibships are the same size)
    sibships <- lapply(1:nrow(parents), f)

    # Calculate the sibship tree
    # WARNING:  Gives edges with a random order (orientation) if i is a founder
    # WARNING:  If i has multiple marriages then this function selects a
    # (fairly arbitrary) subtree of the corr. complete graph
    f1 <- function(i)  {length(sibships[[i]])}
    sib.length <- sapply(1:length(sibships),f1)
    person <- unlist(sibships)
    sibnum <- rep(1:length(sibships),times=sib.length)
    dup.person <- unique(person[duplicated(person)])
    f2 <- function(dp)  {unique(sibnum[person==dp])}
    slist <- lapply(dup.person,f2)
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
        m <- as.matrix(expand.grid(s, s))
        g3 <- function(x) {
          x[1] < x[2]
        }
        m <- m[apply(m, 1, g3), ]
        if (length(j) == 1) {
          m <- matrix(cbind(j, setdiff(s, j)), length(s) - 1, 2)
        }
        if (length(j) == 0) {
          m <- matrix(cbind(s[1], s[-1]), length(s) - 1, 2)
        }
      }
      colnames(m) <- NULL
      m
      return(m)
    }
    st <- lapply(1:length(slist), f3)
    sibtree <- c()
    for (j in 1:length(st)) sibtree <- rbind(sibtree, st[[j]])
    sibtree <- as.data.frame(sibtree)
    names(sibtree) <- c("start", "end")
    row.names(sibtree) <- 1:nrow(sibtree)
    return(list(sibships = sibships, sibtree = sibtree))
  }

  is.tree <- function(st) {
    # Return true if st lists the edges of a connected tree, false otherwise,
    # assuming vertices of st are 1:max(st)
    # st <- sibtree
    numv <- max(st)
    oldv <- c()
    v <- c(1)
    loop.detected <- F
    while (!loop.detected & (length(oldv) < length(v)) & (nrow(st) > 0)) {
      loop.detected <- any((st$end %in% v) & (st$start %in% v))
      boundary.edges <- which((st$end %in% v) | (st$start %in% v))
      oldv <- v
      if (length(boundary.edges) != 0) {
        boundary.edges <- boundary.edges[1]
        v <- unique(c(oldv, st$start[boundary.edges], st$end[boundary.edges]))
        st <- st[-boundary.edges, ]
      }
    }
    all.vertices <- length(v) == numv
    return(!loop.detected & all.vertices)
  }

  # Calculate the sibship graph then check that this is a connected tree --
  # if not then issue an error and stop the function
  # (or disconnect the pedigree?)
  cs <- calc.sibtree(fam)
  sibships <- cs$sibships
  sibtree <- cs$sibtree
  if (!is.tree(sibtree)) {
    stop("The sibship graph is not a connected tree",  call. = FALSE)
    #return(NA)
  }

  # Choose a sibship that contains target_id (this will be collapsed last)
  if (!is.null(target_id)) {
    f <- function(x)  is.element(target_id, x)
    target.sibship <- which(sapply(sibships, f))[1]
  }

  # Collapse sibships one at a time until there are no more
  sibtree.bak <- sibtree
  siborder <- c()
  last.iter <- F
  while (nrow(sibtree) > 0) {
    if (!last.iter) {
      # Choose a (fairly arbitrary) sibship to collapse
      s <- c(sibtree$start, sibtree$end)
      leaves <- setdiff(s, s[duplicated(s)])
      if (!is.null(target_id))  leaves <- setdiff(leaves, target.sibship)
      leavesdown <- intersect(leaves, sibtree$end)
      if (length(leavesdown) > 0) leaves <- leavesdown
      ll <- leaves[1] # We will collapse this sibship ll
      kk <- which(sibtree$start == ll | sibtree$end == ll)
      bad <- setdiff(unique(c(sibtree$start[kk], sibtree$end[kk])), ll)
      if (length(kk) != 1) {
        #warning("Expecting a leaf but sibship ", ll, " is joined to ", bad)
      }
      kk <- kk[1]

      # Find the ID of the person connecting the main tree to the
      # collapsed sibship
      mm <- setdiff(as.numeric(sibtree[kk, ]), ll)
      attachID <- intersect(sibships[[ll]], sibships[[mm]])
      #if (length(attachID) != 1)
      #  warning("A unique person doesn't connect sibships ", ll, " and ", mm)
      attachID <- attachID[1]
    }
    siborder <- c(siborder, ll)

    # Checks
    kk
    sibtree[kk, ]
    ll
    sibships[[ll]]
    attachID

    # Update the penetrance matrix penmat
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
          p <- penmat[kID, ]
          prob <- prob * sum(p * trans[ngeno * gm + gf - ngeno, ] )
        }
        p <- penmat[othparID, ]
        prob <- prob * p[goth] * geno_freq[goth]
        return(prob)
      }
      f1 <- function(gatt) {
        f2 <- function(goth) {
          summand.p(gatt, goth)
        }
        p <- penmat[attachID, ]
        return(p[gatt] * sum(sapply(1:ngeno, f2)))
      }
      newpen <- sapply(1:ngeno, f1)
    }
    if (!attach.is.parent) {
      summand.c <- function(gatt, gm, gf) {
        prob <- 1
        for (kID in setdiff(kidIDs, attachID)) {
          p <- penmat[kID, ]
          prob <- prob * sum(p * trans[ngeno * gm + gf - ngeno, ] )
        }
        pm <- penmat[parentIDs[1], ]
        pf <- penmat[parentIDs[2], ]
        prob <- prob * trans[ngeno * gm + gf - ngeno, gatt] * pm[gm] * pf[gf] *
          geno_freq[gm] * geno_freq[gf]
        return(prob)
      }
      f1 <- function(gatt) {
        f2 <- function(gpar) {
          summand.c(gatt, gpar[1], gpar[2])
        }
        EG <- expand.grid(1:ngeno, 1:ngeno)
        p <- penmat[attachID, gatt]
        #if (!last.iter)  p <- p / geno_freq[gatt]
        p <- p / geno_freq[gatt]
        # return(p * sum(sapply(EG,f2)))
        return(p * sum(apply(EG, 1, f2)))
      }
      newpen <- sapply(1:ngeno, f1)
    }
    if (max(newpen) == 0) {
      return(-Inf)
    }
    ll.adj <- ll.adj + log(max(newpen))
    newpen <- newpen / max(newpen)
    penmat[attachID, ] <- newpen

    # Remove from sibtree the edge connecting the collapsed sibship to the
    # rest of the tree
    if (nrow(sibtree) > 1 | last.iter) {
      sibtree <- sibtree[-kk, ]
    } else {
      # cat("Last iteration \n")
      last.iter <- T
      ll <- setdiff(as.numeric(sibtree[1, ]), ll)
      kk <- 1
      attachID <- sibships[[ll]][3]
      if (!is.null(target_id))  attachID <- target_id
    }
  }

  # Output
  if (is.null(target_id)) {
    # Output the log-likelihood
    ll <- log(sum(geno_freq * penmat[attachID, ])) + ll.adj
    return(ll)
  } else {
    # Output the carrier probabilities for the target person
    p <- geno_freq * penmat[attachID, ]
    p <- p / max(p)
    p <- p / sum(p)
    return(p)
  }


}

