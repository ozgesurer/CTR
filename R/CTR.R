CTR <- function(X, Y, Kmax = 20, cv = 0, warm.start = NULL){

  p <- ncol(X); n <- nrow(X); K <- Kmax; nfolds <- cv

  if(cv > 0){
    if (is.null(warm.start)){
      #Run CV without warm.start
      foldid <- sample(cut(seq(1, n), breaks = nfolds, labels = FALSE)) - 1; K <- 0; sumSqErr <- rep(0, Kmax);
      Groupsf <- rep(0, nfolds * p); Jf <- rep(0, nfolds * p);
      Zf <- rep(0, (Kmax + 1) * n * nfolds); ZZinvf <- rep(0, (Kmax + 1) * (Kmax + 1) * nfolds);

      cout <- .C("CTRwithCV", as.double(X), as.double(Y),
                 as.integer(p), as.integer(n),
                 as.integer(foldid), as.integer(nfolds),
                 as.integer(Kmax), K = as.integer(K),
                 sumSqErr = as.double(sumSqErr),
                 Groupsf = as.integer(Groupsf), Jf = as.integer(Jf),
                 Zf = as.double(Zf), ZZinvf = as.double(ZZinvf))

      K <- cout$K; sumSqErr <- cout$sumSqErr
      Groupsf <- cout$Groupsf; Jf <- cout$Jf
      Zf <- cout$Zf; ZZinvf <- cout$ZZinvf
    }else{
      #Run CV using warm.start
      foldid <- warm.start$foldid
      K <- ncol(warm.start$Zf)/nfolds - 1
      sumSqErr <- c(warm.start$sumSqErr, rep(0, Kmax - K))
      Jf <- matrix(0, nrow = p, ncol = nfolds)
      Groupsf <- matrix(0, nrow = p, ncol = nfolds)
      Zf <- matrix(0, nrow = n, ncol = nfolds*(Kmax + 1))
      ZZinvf <- matrix(0, nrow = Kmax + 1, ncol = (Kmax + 1) * nfolds);

      for(f in 1:nfolds){
        Zf[, (f - 1)*(Kmax + 1) + 1:(K + 1)] <- warm.start$Zf[, (f - 1)*(K + 1) + 1:(K + 1)]
        Groupsf[, f] <- warm.start$Groupsf[, f]; Jf[, f] <- warm.start$Jf[, f]
        ZZinvf[1:(K + 1), (f - 1)*(Kmax + 1) + 1:(K + 1)] <- warm.start$ZZinvf[, (f - 1)*(K + 1) + 1:(K + 1)]
      }

      cout <- .C("CTRwithCV", as.double(X), as.double(Y),
                 as.integer(p), as.integer(n),
                 as.integer(foldid), as.integer(nfolds),
                 as.integer(Kmax), K = as.integer(K),
                 sumSqErr = as.double(sumSqErr),
                 Groupsf = as.integer(Groupsf), Jf = as.integer(Jf),
                 Zf = as.double(Zf), ZZinvf = as.double(ZZinvf))

      K <- cout$K; sumSqErr <- cout$sumSqErr
      Groupsf <- cout$Groupsf; Jf <- cout$Jf
      Zf <- cout$Zf; ZZinvf <- cout$ZZinvf
    }
  }

  size.tree <- ((K + 1) * (K + 2)) / 2; rsq.tree <- rep(0, K); beta.tree <- rep(0, K*(p + 1))
  groups.split <- rep(0, K); group.tree <- rep(0, K*p); J.tree <- rep(0, K*p); nongr <- rep(0, 1)
  A <- matrix(0, nrow = p, ncol = K); Z <- rep(0, (K + 1)*n)

  CTR.model <- .C("CTRwithoutCV", as.double(X), as.double(Y), as.integer(p), as.integer(n), as.integer(K),
                  Z = as.double(Z), rsq.tree = as.double(rsq.tree), beta.tree = as.double(beta.tree),
                  groups.split = as.double(groups.split), group.tree = as.double(group.tree), J.tree = as.double(J.tree), nongr = as.integer(nongr))

  #Get the tree structure
  rsq.tree <- CTR.model$rsq.tree; beta.tree <- matrix(CTR.model$beta.tree, nrow = p + 1)
  group.tree <- matrix(CTR.model$group.tree, nrow = p); J.tree <- matrix(CTR.model$J.tree, nrow = p) + 1
  nongr <- CTR.model$nongr; groups.split <- matrix(CTR.model$groups.split, ncol = K)
  Z <- matrix(CTR.model$Z, nrow = n)
  groups <-  unique(group.tree[, K])
  for (j in 1:K)
    A[J.tree[, K][group.tree[, K] == groups[j]], j] = 1

  if(cv > 0){
    ZZinvf <- matrix(ZZinvf, ncol = (Kmax + 1)*nfolds)
    Zf <- matrix(Zf, ncol = (Kmax + 1)*nfolds)
    Jf <- matrix(Jf, ncol = nfolds)
    Groupsf <- matrix(Groupsf, ncol = nfolds)

    out <- list(A = A, Z = Z, rsq.tree = rsq.tree, beta.tree = beta.tree,
                group.tree = group.tree, J.tree = J.tree, nongr = nongr, groups.split = groups.split,
                sumSqErr = sumSqErr, foldid = foldid, Zf = Zf, Jf = Jf, Groupsf = Groupsf, ZZinvf = ZZinvf)
    class(out) <- "CTR"
  }else{
    out <- list(A = A, Z = Z, rsq.tree = rsq.tree, beta.tree = beta.tree,
                group.tree = group.tree, J.tree = J.tree, nongr = nongr, groups.split = groups.split)
    class(out) <- "CTR"
  }
  out
}

summary.CTR <- function(object, ...){

  A <- object$A; p <- nrow(A); K <- ncol(A); n <- nrow(object$Z)
  Rsq <- object$rsq.tree; group.tree <- object$group.tree; J.tree <- object$J.tree
  alpha <- matrix(coef(object)$alpha, nrow = K + 1)
  alpha <- c(alpha[2:(K+1)],0)
  groups <- unique(group.tree[, K])
  zerocol <- matrix(0, nrow = p)
  zerocol[J.tree[, K][group.tree[, K] == groups[K+1]],  1] = 1

  cat("Coefficient Tree Regression\n\n")
  cat("n =", n, "\n")
  cat("p =", p, "\n")
  cat("k =", K, "(number of derived predictors at termination)\n")
  cat(sum(A), "predictors included in", K, "groups at termination \n\n")

  A <- cbind(A, zerocol)

  if(is.null(object$sumSqErr))
    cat("Training r.squared (at termination) = ", round(Rsq[K], 4), "\n\n", sep = "")
  else
    cat("Training r.squared (at termination) = ", round(Rsq[K], 4), ", ", "CV r.squared = ", round(1 - object$sumSqErr[K], 4), "\n\n", sep = "")

  cat("Final group structure:\n\n")
  cat("(Group)", "(Size)", "(alpha)", "(Ids)\n", sep = "")

  for (k in 1: (K + 1)){
    if(is.null(object$sumSqErr))
      cat(k, "  ", sum(A[, k]), " ", sep = "")
    else
      cat(k, "  ", sum(A[, k]), " ", sep = "")

    iter <- 1
    cat(round(alpha[k], 2), " ", sep = "")
    for ( j in 1:nrow(A) ){
      if( A[j, k] == 1 ){
        if ( iter == 1 ){
          if( iter == sum(A[, k]) ){
            cat( "{", j, "}", sep = "" )
            iter = iter + 1
            }else{
              cat( "{", j, sep = "" )
              iter = iter + 1
              }
          }else if ( iter > 1 ){
            if( iter == sum(A[, k]) ){
              if( A[j - 1, k] == 0 ){
                cat(",", j, "}", sep = "")
                iter = iter + 1
              }else if( A[j - 1, k] == 1 ){
                cat("-", j, "}", sep = "")
                iter = iter + 1
              }
            }else if( A[j - 1, k] == 0 ){
              cat( ",", j, sep = "" )
              iter = iter + 1
            }else if(A[j - 1, k] == 1 & A[j + 1, k] == 0){
              cat("-", j, sep = "")
              iter = iter + 1
            }else if (A[j - 1, k] == 1 & A[j + 1, k] == 1){
              iter = iter + 1
            }
          }
        }
      }
      cat("\n")
    }

}
plot.CTR <- function(x, ...){
  relativeerror <- c(1, x$sumSqErr)
  y <- length(relativeerror)
  plot(0: (y - 1), relativeerror, xlab = "Model size", ylab = "Relative error", main = "Cross-validation result", xaxt = "n", xlim = c(0, y - 1))
  lines(0:(y - 1), relativeerror)
  axis(1, at = seq(0, (y - 1), by = 1), las = 2)
}
coef.CTR <- function(object, ...){
  A <- object$A; beta.tree <- object$beta.tree; K <- ncol(A)
  alpha <- matrix(0, nrow = K + 1); beta <- matrix(beta.tree[, K], nrow = nrow(A) + 1)
  alpha[1] <- beta[1]

  for(i in 1:K)
    alpha[i + 1] <- beta[which.max(A[, i]) + 1]

  rownames(alpha) <- c("intercept", paste0("z", 1:ncol(A)))
  rownames(beta) <- c("intercept", paste0("x", 1:nrow(A)))

  coefs <- list(alpha = alpha, beta = beta)
  coefs
}
predict.CTR <- function(object, newdata, ...){
  Xtest <- newdata

  coefs <- coef.CTR(object)
  beta <- coefs$beta

  predictions <- cbind(1, Xtest) %*% beta
  predictions
}
print.CTR <- function(x, option = c("text", "tree"), ...){
  groups.split <- x$groups.split; J.tree <- x$J.tree; group.tree <- x$group.tree; beta.tree <- x$beta.tree
  K <- ncol(group.tree); p <- nrow(group.tree); size.tree <- (K + 1)*(K + 2)/2
  coef.tree <- matrix(0, nrow = p, ncol = size.tree); coef.tree[, 1] = 1
  nongr <- x$nongr
  #Create a tree structure
  colno <- 2
  groups.split.tree <- rep(0, size.tree); groups.split.tree[1] <- 1
  for (i in 1:K){
    groups <- unique(group.tree[,i])
    if(length(groups) != i + 1)
      groups <- c(groups, nongr)
    for (j in 1:(i+1)){
      if(i < K)
        if(groups[j] == groups.split[i + 1])
          groups.split.tree[colno] <- 1
        coef.tree[J.tree[,i][group.tree[, i] == groups[j]], colno] <- 1
        colno = colno + 1
    }
  }
  #Create alpha tree
  colno <- 1
  alpha.tree <- c()
  for( i in  1:K){
    A <- matrix(coef.tree[, (colno + 1):(colno + i)], nrow = p)
    alpha <- matrix(0, nrow = i + 1)
    beta <- matrix(beta.tree[, i], nrow = p + 1)
    beta <- beta[-1]
    for(j in 1:i){
      alpha[j] <- beta[which.max(A[,j])]
      colno  <- colno + 1
    }
    alpha[i+1] <- 0
    colno <- colno + 1
    alpha.tree <- c(alpha.tree, alpha)
  }
  alpha.tree <- c(0, alpha.tree)

  if(option == "text"){
    #Create parent structure in the tree
    parent.tree <- c()
    parent.tree[1] <- 0
    c1 <- 1
    for (i in 1:K){
      level.id <- c1:(c1 + i - 1)
      for (j in 1:i){
        if(groups.split.tree[c1] == 1)
          Level = j
        c1 = c1 + 1
      }
      parent.id <- c(); c2 <- 1; c3 <- 1
      for (j in 1:(i+1)){
        if (j == Level){
          parent.id[c2] = level.id[c3]
          c2 = c2 + 1
          parent.id[c2] = level.id[c3]
          c2 = c2 + 1
          c3 = c3 + 1
        }else{
          if(length(parent.id) < i + 1){
            parent.id[c2] = level.id[c3]
            c2 = c2 + 1
          c3 = c3 + 1
          }
        }
      }
      parent.tree <- c(parent.tree, parent.id)
    }

    isPrint <- matrix(0, dim(coef.tree)[2], 1)
    K <- ncol(x$A)
    isFirstPrint <- TRUE

    parent1 <- 1
    parent2 <- 1

    level1 <- rep(1:(K + 1), 1:(K + 1))
    level2 <- matrix(1, nrow = K + 1)

    cat("(Group)(alpha){Ids}")
    printLevel <- function(coef.tree, k){
      iter <- 1
      for ( j in 1:nrow(coef.tree) ){
        if( coef.tree[j, k] == 1 ){
          if ( iter == 1 ){
            if( iter == sum(coef.tree[, k]) ){
              cat( "{", j, "}", sep = "" )
              iter = iter + 1
            }else{
              cat( "{", j, sep = "" )
              iter = iter + 1
            }
          }else if ( iter > 1 ){
            if( iter == sum(coef.tree[, k]) ){
              if( coef.tree[j - 1, k] == 0 ){
                cat(",", j, "}", sep = "")
                iter = iter + 1
              }else if( coef.tree[j - 1, k] == 1 ){
                cat("-", j, "}", sep = "")
                iter = iter + 1
              }
            }else if( coef.tree[j - 1, k] == 0 ){
              cat( ",", j, sep = "" )
              iter = iter + 1

            }else if(coef.tree[j - 1, k] == 1 & coef.tree[j + 1, k] == 0){
              cat("-", j, sep = "")
              iter = iter + 1
            }else if (coef.tree[j - 1, k] == 1 & coef.tree[j + 1, k] == 1){
              iter = iter + 1
            }
          }
        }
      }
    }

    parentList <- c()
    k <- 1
    cat("\n")
    cat("(", level1[k] - 1, ".", level2[level1[k]], ")", "(", round(alpha.tree[k], 4), ")", sep = "")
    printLevel(coef.tree, k)
    cat("(s)")
    isPrint[1] <- 1
    parent1 <- 1
    while (sum(isPrint) < size.tree)
    {
      isPrinted <- FALSE
      for (k in 1:size.tree){
        if (isPrint[k] == 0){
          if (parent.tree[k] == parent1){
            cat("\n")
            cat(character(level1[k] - 1), collapse = " ")
            cat("(", level1[k] - 1, ".", level2[level1[k]], ")", "(", round(alpha.tree[k], 4), ")", sep = "")
            printLevel(coef.tree, k)
            if(groups.split.tree[k] == 1)
              cat("(s)")
            isPrint[k] <- 1
            level2[level1[k]] <- level2[level1[k]] + 1
            if(k == size.tree - 1)
            {
              if(parent.tree[k + 1] == parent1){
                cat("\n")
                cat(character(level1[k + 1] - 1), collapse = " ")
                cat("(", level1[k + 1] - 1, ".", level2[level1[k + 1]], ")", "(", round(alpha.tree[k + 1], 4), ")", sep = "")
                printLevel(coef.tree, k + 1)
                isPrint[k + 1] = 1
                level2[level1[k + 1]] <- level2[level1[k + 1]] + 1
              }
            }else{
              parent1 <- k
              parentList <- union(parentList, k)
              parentList <- sort(parentList)
              isPrinted <- TRUE
            }
          }
        }
      }
      if(isPrinted == TRUE){
        parentList <- rev(parentList)
        parent1 <- parentList[1]
      }else{
        parentList <- parentList[-1]
        parent1 <- parentList[1]
      }
      if(length(parentList) == 0)
        parent1 <- 1
    }
  }else{
    #Print Tree Structure in Visual Form
    labels.tree <- c()
    printLevel <- function(coef.tree, k){
      iter <- 1
      lab0 <- c()
      for ( j in 1:nrow(coef.tree) ){
        if( coef.tree[j, k] == 1 ){
          if ( iter == 1 ){
            if( iter == sum(coef.tree[, k]) ){
              lab0 <- paste( "{", j, "}", sep = "" )
              iter = iter + 1
            }else{
              lab0 <- paste( "{", j, sep = "" )
              iter = iter + 1
            }
          }else if ( iter > 1 ){
            if( iter == sum(coef.tree[, k]) ){
              if( coef.tree[j - 1, k] == 0 ){
                lab0 <- paste(lab0, ",", j, "}", sep = "")
                iter = iter + 1
              }else if( coef.tree[j - 1, k] == 1 ){
                lab0 <- paste(lab0, "-", j, "}", sep = "")
                iter = iter + 1
              }
            }else if( coef.tree[j - 1, k] == 0 ){
              lab0 <- paste(lab0, ",", j, sep = "" )
              iter = iter + 1

            }else if(coef.tree[j - 1, k] == 1 & coef.tree[j + 1, k] == 0){
              lab0 <- paste(lab0, "-", j, sep = "")
              iter = iter + 1
            }else if (coef.tree[j - 1, k] == 1 & coef.tree[j + 1, k] == 1){
              iter = iter + 1
            }
          }
        }
      }
      if (length(lab0) > 0)
        lab0
      else{
        lab0 <- paste( " ", sep = "" )
        lab0
      }
    }

    label.tree <- c()
    for (k in 1:size.tree){
      label.tree[k] <- printLevel(coef.tree, k)
    }

    ndf <- create_node_df(
        n = size.tree,
        label = label.tree,
        shape = "rectangle")
    ndf

    fromvec <- rep(seq(1:(K*(K + 1)/2)), (groups.split.tree + 1)[1:(K*(K + 1)/2)])
    tovec <- 2:size.tree
    edf <- create_edge_df(
        from = fromvec,
        to   = tovec,
        rel = c("leading_to"),
        label = paste(round(alpha.tree[-1],2), sep=" "))
    edf

    the_graph <- create_graph(
        nodes_df = ndf,
        edges_df = edf, directed = TRUE)

    render_graph(graph = the_graph, layout = "tree", output = "graph")

  }
}
