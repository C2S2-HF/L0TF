invtf2d=function(dim1,dim2,k){

	# Check
	stopifnot(dim1>k+1)
	stopifnot(dim2>k+1)
	
	D1 <- getDtfSparse(dim1, k)
	diags <- lapply(1:dim2,function(x){D1})
	D.part1 <- Matrix::bdiag(diags)

	if(is.null(dim(D1))){
		diags <- rep(D1, each = dim1*dim2)
	}else{
		diags <- rep(D1[1,1:(k+2)], each = dim1*dim2)
	}
	diags <- matrix(diags, nrow=dim1*dim2)
	D.part2 <- Matrix::bandSparse(dim1*(dim2-k-1), dim1*dim2, k = seq(0, dim1*(k+1), dim1), diag = diags, symm=FALSE)
	D <- rbind(D.part1, D.part2)
	DTD <- D%*%Matrix::t(D);rm(D);gc()

	if(dim1*(dim2-k-1)+dim2*(dim1-k-1)<dim1*dim2){
		if(rcond(as.matrix(DTD))<1e-12){
			Matrix::diag(DTD)<-Matrix::diag(DTD)+1e-6
		}
	}else{
		Matrix::diag(DTD)<-Matrix::diag(DTD)+1e-6
	}
	return(Matrix::solve(sparse = TRUE, DTD))

}

inv2d=function(dim1,dim2){
	
	D <- getD2dSparse(dim1=dim1,dim2=dim2)
	DTD <- D%*%	Matrix::t(D);rm(D);gc()

	if(dim1*(dim2-1)+dim2*(dim1-1)<dim1*dim2){
		if(rcond(as.matrix(DTD))<1e-12){
			Matrix::diag(DTD)<-Matrix::diag(DTD)+1e-6
		}
	}else{
		Matrix::diag(DTD)<-Matrix::diag(DTD)+1e-6
	}
	return(Matrix::solve(sparse = TRUE, DTD))

}

#
# function that is the same with the one in genlasso
#
getDtfSparse=function (n, k)
{
  D = Matrix::bandSparse(n, m = n, c(0, 1), diagonals = list(rep(-1, n), rep(1, n - 1)))
  D0 = D
  for (i in seq(1, k)) D = D0 %*% D
  return(D[seq(1, n - k - 1), ])
}
getD2dSparse=function(dim1, dim2){
  D1 = bandSparse(dim1 * dim2, m = dim1 * dim2, k = c(0, 1), diagonals = list(rep(-1, dim1 * dim2), rep(1, dim1 * dim2 - 1)))
  D1 = D1[(seq(1, dim1 * dim2)%%dim1) != 0, ]
  D2 = bandSparse(dim1 * dim2 - dim1, m = dim1 * dim2, k = c(0, dim1), diagonals = list(rep(-1, dim1 * dim2), rep(1, dim1 * dim2 - 1)))
  return(rbind(D1, D2))
}
