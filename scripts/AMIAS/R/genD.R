#
# genD for 1d in trend filtering
#
genDtf1d=function(k, n = NULL, full = FALSE){

  # Check
  k <- as.integer(k)
  stopifnot(k>=1)

  if(full){

    # Check for full mode
    if(is.null(n))stop("n hasn't passed to genDtf1d in full mode")
    n <- as.integer(n)
    stopifnot(k<n-1)

    D1d <- matrix(c(-1,1),nrow=1,ncol=2)
    Dtf1 <- matrix(0,nrow=2,ncol=4)
    Dtf1[1,-4] = c(1,-2,1)
    Dtf1[2,-1] = c(1,-2,1)

    if(k == 1){
      result = lapply(1:(n-2),function(x){
        v = rep(0,n)
        v[x:(x+2)] = c(1,-2,1)
        return(v)
      })
    }else if(k == 2){
      vec = D1d %*% Dtf1
      result = lapply(1:(n-3),function(x){
        v = rep(0,n)
        v[x:(x+3)] = vec
        return(v)
      })
    }else{
      vec <- D1d %*% Dtf1
      for (i in 3:k){
        Dtf1 <- matrix(0,nrow=2,ncol=i+2)
        Dtf1[1,-(i+2)] <- vec
        Dtf1[2,-1] <- vec
        vec <- D1d %*% Dtf1
      }
      result = lapply(1:(n-k-1),function(x){
        v = rep(0,n)
        v[x:(x+k+1)] = vec
        return(v)
      })
    }
    return(do.call(rbind,result))

  }else{

    # only return the vec, nothing to do with n
    D1d = matrix(c(-1,1),nrow=1,ncol=2)
    Dtf1 = matrix(0,nrow=2,ncol=4)
    Dtf1[1,-4] = c(1,-2,1)
    Dtf1[2,-1] = c(1,-2,1)

    if(k == 1){
      return(c(1,-2,1))
    }else if(k == 2){
      return(D1d %*% Dtf1)
    }else{
      vec <- D1d %*% Dtf1
      for (i in 3:k){
        Dtf1 <- matrix(0,nrow=2,ncol=i+2)
        Dtf1[1,-(i+2)] <- vec
        Dtf1[2,-1] <- vec
        vec <- D1d %*% Dtf1
      }
      return(vec)
    }
  }

}
#
# genD for 2d in trend filtering
#
genDtf2d=function(k, dim1 = NULL, dim2 = NULL, full = FALSE){

	# Check
	k <- as.integer(k)
	stopifnot(k>=1)

	if(full){

		# Check for full mode
		if(is.null(dim2)|is.null(dim1))stop("dim2/dim1 hasn't passed to genDtf2d in full mode")
		dim2 <- as.integer(dim2)
		dim1 <- as.integer(dim1)
		stopifnot(k<dim2-1)
		stopifnot(k<dim1-1)

		D1d <- matrix(c(-1,1),nrow=1,ncol=2)
		Dtf1 <- matrix(0,nrow=2,ncol=4)
		Dtf1[1,-4] = c(1,-2,1)
		Dtf1[2,-1] = c(1,-2,1)

		if(k == 1){
			result.part1.item = lapply(1:(dim1-2),function(x){
				v = rep(0,dim1)
				v[x:(x+2)] = c(1,-2,1)
				return(v)
			})
			result.part1.item = do.call(rbind,result.part1.item)
			result.part1 = lapply(1:dim2,function(x){
				v = matrix(0,(dim1-2),dim1*dim2)
				v[,((x-1)*dim1+1):(x*dim1)] = result.part1.item
				return(v)
			})
			result.part1 = do.call(rbind,result.part1)
			result.part2 = lapply(1:((dim2-2)*dim1),function(x){
				v = rep(0,dim1*dim2)
				v[c(1,1+dim1,1+2*dim1)+(x-1)] = c(1,-2,1)
				return(v)
			})
			result.part2 = do.call(rbind,result.part2)
			result = rbind(result.part1,result.part2)
		}else if(k == 2){
			vec = D1d %*% Dtf1
			result.part1.item = lapply(1:(dim1-3),function(x){
				v = rep(0,dim1)
				v[x:(x+3)] = vec
				return(v)
			})
			result.part1.item = do.call(rbind,result.part1.item)
			result.part1 = lapply(1:dim2,function(x){
				v = matrix(0,(dim1-3),dim1*dim2)
				v[,((x-1)*dim1+1):(x*dim1)] = result.part1.item
				return(v)
			})
			result.part1 = do.call(rbind,result.part1)
			result.part2 = lapply(1:((dim2-3)*dim1),function(x){
				v = rep(0,dim1*dim2)
				v[c(1,1+dim1,1+2*dim1,1+3*dim1)+(x-1)] = vec
				return(v)
			})
			result.part2 = do.call(rbind,result.part2)
			result = rbind(result.part1,result.part2)
		}else{
			vec <- D1d %*% Dtf1
			for (i in 3:k){
				Dtf1 <- matrix(0,nrow=2,ncol=i+2)
				Dtf1[1,-(i+2)] <- vec
				Dtf1[2,-1] <- vec
				vec <- D1d %*% Dtf1
			}
			result.part1.item = lapply(1:(dim1-k-1),function(x){
				v = rep(0,dim1)
				v[x:(x+k+1)] = vec
				return(v)
			})
			result.part1.item = do.call(rbind,result.part1.item)
			result.part1 = lapply(1:dim2,function(x){
				v = matrix(0,(dim1-k-1),dim1*dim2)
				v[,((x-1)*dim1+1):(x*dim1)] = result.part1.item
				return(v)
			})
			result.part1 = do.call(rbind,result.part1)
			result.part2 = lapply(1:((dim2-k-1)*dim1),function(x){
				v = rep(0,dim1*dim2)
				v[seq(from=1,by=dim1,length.out=k+2)+(x-1)] = vec
				return(v)
			})
			result.part2 = do.call(rbind,result.part2)
			result = rbind(result.part1,result.part2)
		}
		return(result)

	}else{

		# only return the vec, nothing to do with dim1,dm2
		D1d = matrix(c(-1,1),nrow=1,ncol=2)
		Dtf1 = matrix(0,nrow=2,ncol=4)
		Dtf1[1,-4] = c(1,-2,1)
		Dtf1[2,-1] = c(1,-2,1)

		if(k == 1){
			return(c(1,-2,1))
		}else if(k == 2){
			return(D1d %*% Dtf1)
		}else{
			vec <- D1d %*% Dtf1
			for (i in 3:k){
				Dtf1 <- matrix(0,nrow=2,ncol=i+2)
				Dtf1[1,-(i+2)] <- vec
				Dtf1[2,-1] <- vec
				vec <- D1d %*% Dtf1
			}
			return(vec)
		}
	}

}
#
# genD for 2d
#
genD2d=function(dim1=NULL, dim2=NULL, full = FALSE){

	if(!full){
		# Return derectly
		return(c(-1,1))
	}else{

		# Check for full mode
		if(is.null(dim2)|is.null(dim1))stop("dim2/dim1 hasn't passed to genD2d in full mode")
		dim2 <- as.integer(dim2)
		dim1 <- as.integer(dim1)
		stopifnot(1<dim2)
		stopifnot(1<dim1)

		result.part1.item = lapply(1:(dim1-1),function(x){
			v = rep(0,dim1)
			v[x:(x+1)] = c(-1,1)
			return(v)
		})
		result.part1.item = do.call(rbind,result.part1.item)
		result.part1 = lapply(1:dim2,function(x){
			v = matrix(0,(dim1-1),dim1*dim2)
			v[,((x-1)*dim1+1):(x*dim1)] = result.part1.item
			return(v)
		})
		result.part1 = do.call(rbind,result.part1)
		result.part2 = lapply(1:((dim2-1)*dim1),function(x){
			v = rep(0,dim1*dim2)
			v[c(1,1+dim1)+(x-1)] = c(-1,1)
			return(v)
			})
		result.part2 = do.call(rbind,result.part2)
		result = rbind(result.part1,result.part2)
		return(result)
	}

}
