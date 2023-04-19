getM = function(type="identity",k=NULL,n){
	if(type ==  "fused.1d"){
		M = matrix(0,n-1,n)
		diag(M[,-n]) = -1
		diag(M[,-1]) = 1
	}else if(type == "fused.tfk"){
		M = genDtf1d(k,n,T)
	}else if(type == "identity"){
		M = diag(n)
	}
	M
}

getM2d = function(type="identity",k=NULL,dim1,dim2){
	if(type ==  "fused.1d"){
		M = genD2d(dim1,dim2,T)
	}else if(type == "fused.tfk"){
		M = genDtf2d(k,dim1,dim2,T)
	}else if(type == "identity"){
		M = diag(dim1*dim2)
	}
	M
}
