!*****************************************************************************80
!
!! fusedl0(for 1d) return the fitting coefficients of given k with a small loss
!! function in fusedl0 via primal dual active set algorithm.
!
!  Licensing:
!
!    This code is distributed under LGPL-3.
!
!  Author:
!
!    Original R version by Wen canhong, replica fortran90 version by Hong Zelin.
!
!  Modified:
!
!    29 August 2017
!
!  Parameters:
!
!    Input, integer (kind = 4) N, The length of y.
!    Input, real (kind = 8) y(N), The observe variable.
!    Input, real (kind = 8) beta(N), The fitting coefficients of given 
!         T with a small loss function in fusedl0.
!    Input, real (kind = 8) z(N-1), the split variable of the argumented 
!         lagrange form.
!    Input, real (kind = 8) u(N-1), the lagrange operator of the argumented 
!         lagrange form for linear item.
!    Input, integer (kind = 4) T, the number of discover point by fusedl0.
!    Input, real (kind = 8) rho, the lagrange operator of the argumented 
!         lagrange form for split item.
!    Input, interger (kind = 4) iter, The iterations algorithm need to converge.
!    Input, integer (kind = 4) itermax, Maximum allowed number of iterations 
!         over the data for given parameters.
!
!
subroutine fusedl0(n, y, beta, z, u, A, I, Anull, T, rho, iter, itermax, adjust, adjust_max, delta)

  implicit none

  real(kind = 8), dimension (:), allocatable :: tempvec,tempDTD
  integer, dimension (:), allocatable :: ordervec,Alast
  
  integer, intent(in) :: adjust
  integer, intent(in) :: adjust_max
  integer, intent(in) :: delta
  logical endmark
  integer n, iter, itermax, j, k, T, Anull
  real(kind = 8) y(n), beta(n), z(n-1), u(n-1), rho
  integer A(T), I(n-1-T)

  allocate(tempvec(1:n))
  allocate(tempDTD(1:n-1))
  allocate(ordervec(1:n-1))
  allocate(Alast(1:T))

  tempvec = 0
  tempDTD = 0
  Alast = 0
  ! open(2, file='D:/data1.txt', status='old')
  if(Anull == 1) then
      tempvec(1:(n-1)) = y(2:n) - y(1:(n-1))
      Do j = 1, n-1
        call itergen(j,n-1,tempDTD(1:n-1))
        u(j) = dot_product(tempDTD(1:(n-1)),tempvec(1:(n-1)))
      end do

      u(n-1) = dot_product(tempDTD(1:(n-2)),tempvec(1:(n-2)))
      Do j = 1,n-1
        ordervec(j) = j
      end do
      call Qsortorder(abs(z + u/rho),n-1,ordervec)
      A = ordervec(1:T)
      I = ordervec((T+1):n-1)
  
      if(adjust == 1 .and. T /= 1) then
          call adjustA(A, I, abs(z + u/rho), adjust_max, T, delta, n-1)
      else
          if( T>1 )then
              call qSortvec(A,T)
          end if
      end if
  else
      if( T>1 )then
          call qSortvec(A,T)
      end if
  end if  

  iter = 1
  Do while(iter < itermax)

    if(T .eq. 1)then
        if(A(1) .eq. 1 )then
            tempvec(1:(n-2)) = y(3:n) - y(2:(n-1))
            Do j = 1, n-2
                call itergen(j,n-2,tempDTD(1:n-2))
                u(j+1) = dot_product(tempDTD(1:(n-2)),tempvec(1:(n-2)))
            end do
        else if(A(1) .eq. n-1)then
            tempvec(1:(n-2)) = y(2:(n-1)) - y(1:(n-2))
            Do j = 1, n-2
                call itergen(j,n-2,tempDTD(1:n-2))
                u(j) = dot_product(tempDTD(1:(n-2)),tempvec(1:(n-2)))
            end do
        else
            tempvec(1:(A(1)-1)) = y(2:A(1)) - y(1:(A(1)-1))
            Do j = 1, A(1)-1
                call itergen(j,A(1)-1,tempDTD(1:A(1)-1))
                u(j) = dot_product(tempDTD(1:A(1)-1),tempvec(1:A(1)-1))
            end do
            j = 1
            tempvec(1:(n-1-A(1))) = y(A(1)+2:n) - y(A(1)+1:(n-1))
            Do j = 1, n-1-A(1)
                call itergen(j,n-1-A(1),tempDTD(1:n-1-A(1)))
                u(j+A(1)) = dot_product(tempDTD(1:n-1-A(1)),tempvec(1:n-1-A(1)))
            end do
        end if
    else
        if( A(1) > 1 )then
            tempvec(1:(A(1)-1)) = y(2:A(1)) - y(1:(A(1)-1))
            Do j = 1, A(1)-1
                call itergen(j,A(1)-1,tempDTD(1:A(1)-1))
                u(j) = dot_product(tempDTD(1:A(1)-1),tempvec(1:A(1)-1))
            end do
        end if
        Do k = 1,T-1
            if( A(k+1)-A(k)>1 )then
                j = 1
                tempvec(1:(A(k+1)-A(k)-1)) = y(A(k)+2:A(k+1)) - y(A(k)+1:(A(k+1)-1))
                Do j = 1, A(k+1)-A(k)-1
                    call itergen(j,A(k+1)-A(k)-1,tempDTD(1:A(k+1)-A(k)-1))
                    u(j+A(k)) = dot_product(tempDTD(1:A(k+1)-A(k)-1),tempvec(1:A(k+1)-A(k)-1))
                end do
            end if
        end do
        if( A(T) < n-1 )then
            j = 1
            tempvec(1:(n-1-A(T))) = y(A(T)+2:n) - y(A(T)+1:(n-1))
            Do j = 1, n-1-A(T)
                call itergen(j,n-1-A(T),tempDTD(1:n-1-A(T)))
                u(j+A(T)) = dot_product(tempDTD(1:n-1-A(T)),tempvec(1:n-1-A(T)))
            end do
        end if
    end if
    u(A(1:T)) = 0
    tempvec(1) = -u(1)
    tempvec(2:n-1) = u(1:n-2) - u(2:n-1)
    tempvec(n) = u(n-1)
    tempvec = y - tempvec
    z(A(1:T)) = tempvec(A(1:T)+1)-tempvec(A(1:T))
    z(I(1:n-T-1)) = 0

    j = 1
    Do j = 1,n-1
        ordervec(j) = j
    end do
    call Qsortorder(abs(z + u/rho),n-1,ordervec)
    A = ordervec(1:T)
    I = ordervec((T+1):n-1)

	if(adjust == 1 .and. T /= 1) then
        call adjustA(A, I, abs(z + u/rho), adjust_max, T, delta, n-1)
    else
        if( T>1 )then
            call qSortvec(A,T)
        end if
    end if
    endmark = .true.
    j = 1
    Do j = 1,T
        if(Alast(j) .ne. A(j))then
            endmark = .false.
            EXIT
        end if
    end do

    if(endmark)then
        iter = iter + 1
        EXIT
    else
        iter = iter + 1
    end if

    Alast = A

  end do
  ! close(2)
  beta = tempvec
  iter = iter - 1
  deallocate(tempvec,tempDTD,ordervec,Alast)
  return
end
subroutine afusedl0(n, y, beta, z, u, tao, Kmax, L, eps, rho, miter, itermax, adjust, adjust_max, delta)

  implicit none

  real(kind = 8), dimension (:), allocatable :: tempvec,tempDTD
  integer, dimension (:), allocatable :: A,I,ordervec,Alast

  integer, intent(in) :: adjust
  integer, intent(in) :: adjust_max
  integer, intent(in) :: delta
  logical endmark
  integer n, iter, itermax, tao, Kmax, L, j, k, T, kk
  integer miter(Kmax)
  real(kind = 8) y(n), beta(n,Kmax), z(n-1, Kmax), u(n-1, Kmax), rho, eps

  allocate(tempvec(1:n))
  allocate(tempDTD(1:n-1))
  allocate(A(1:L))
  allocate(I(1:n-1))
  allocate(ordervec(1:n-1))
  allocate(Alast(1:L))

  kk = 1
  tempvec = 0
  tempDTD = 0
  Alast = 0

  Do while(kk < Kmax)

    T = tao * kk
    if(kk == 1) then
        tempvec(1:(n-1)) = y(2:n) - y(1:(n-1))
        Do j = 1, n-1
          call itergen(j,n-1,tempDTD(1:n-1))
          u(j, kk) = dot_product(tempDTD(1:(n-1)),tempvec(1:(n-1)))
        end do
        Do j = 1,n-1
          ordervec(j) = j
        end do
        call Qsortorder(abs(z(:,kk) + u(:,kk)/rho),n-1,ordervec)
    end if

    A(1:T) = ordervec(1:T)
    I(1:n-1-T) = ordervec((T+1):n-1)
    if(kk == 1) then
        if(adjust == 1 .and. T /= 1) then
            call adjustA(A, I, abs(z(:,kk) + u(:,kk)/rho), adjust_max, T, delta, n-1)
        else
            if( T>1 )then
                call qSortvec(A,T)
            end if
        end if
    else
        if( T>1 )then
            call qSortvec(A,T)
        end if    
    end if
    iter = 1
    Do while(iter < itermax)

      if(T .eq. 1)then
        if(A(1) .eq. 1 )then
            tempvec(1:(n-2)) = y(3:n) - y(2:(n-1))
            Do j = 1, n-2
                call itergen(j,n-2,tempDTD(1:n-2))
                u(j+1,kk) = dot_product(tempDTD(1:(n-2)),tempvec(1:(n-2)))
            end do
        else if(A(1) .eq. n-1)then
            tempvec(1:(n-2)) = y(2:(n-1)) - y(1:(n-2))
            Do j = 1, n-2
                call itergen(j,n-2,tempDTD(1:n-2))
                u(j,kk) = dot_product(tempDTD(1:(n-2)),tempvec(1:(n-2)))
            end do
        else
            tempvec(1:(A(1)-1)) = y(2:A(1)) - y(1:(A(1)-1))
            Do j = 1, A(1)-1
                call itergen(j,A(1)-1,tempDTD(1:A(1)-1))
                u(j,kk) = dot_product(tempDTD(1:A(1)-1),tempvec(1:A(1)-1))
            end do
            j = 1
            tempvec(1:(n-1-A(1))) = y(A(1)+2:n) - y(A(1)+1:(n-1))
            Do j = 1, n-1-A(1)
                call itergen(j,n-1-A(1),tempDTD(1:n-1-A(1)))
                u(j+A(1),kk) = dot_product(tempDTD(1:n-1-A(1)),tempvec(1:n-1-A(1)))
            end do
        end if
      else
        if( A(1) > 1 )then
            tempvec(1:(A(1)-1)) = y(2:A(1)) - y(1:(A(1)-1))
            Do j = 1, A(1)-1
                call itergen(j,A(1)-1,tempDTD(1:A(1)-1))
                u(j,kk) = dot_product(tempDTD(1:A(1)-1),tempvec(1:A(1)-1))
            end do
        end if
        Do k = 1,T-1
            if( A(k+1)-A(k)>1 )then
                j = 1
                tempvec(1:(A(k+1)-A(k)-1)) = y(A(k)+2:A(k+1)) - y(A(k)+1:(A(k+1)-1))
                Do j = 1, A(k+1)-A(k)-1
                    call itergen(j,A(k+1)-A(k)-1,tempDTD(1:A(k+1)-A(k)-1))
                    u(j+A(k),kk) = dot_product(tempDTD(1:A(k+1)-A(k)-1),tempvec(1:A(k+1)-A(k)-1))
                end do
            end if
        end do
        if( A(T) < n-1 )then
            j = 1
            tempvec(1:(n-1-A(T))) = y(A(T)+2:n) - y(A(T)+1:(n-1))
            Do j = 1, n-1-A(T)
                call itergen(j,n-1-A(T),tempDTD(1:n-1-A(T)))
                u(j+A(T),kk) = dot_product(tempDTD(1:n-1-A(T)),tempvec(1:n-1-A(T)))
            end do
        end if
      end if
      u(A(1:T),kk) = 0
      tempvec(1) = -u(1,kk)
      tempvec(2:n-1) = u(1:n-2,kk) - u(2:n-1,kk)
      tempvec(n) = u(n-1,kk)
      tempvec = y - tempvec
      z(A(1:T),kk) = tempvec(A(1:T)+1)-tempvec(A(1:T))
      z(I(1:n-T-1),kk) = 0

      j = 1
      Do j = 1,n-1
        ordervec(j) = j
      end do
      call Qsortorder(abs(z(:,kk) + u(:,kk)/rho),n-1,ordervec)
      A(1:T) = ordervec(1:T)
      I(1:n-1-T)= ordervec((T+1):n-1)
      if(adjust == 1 .and. T /= 1) then
          call adjustA(A, I, abs(z(:,kk) + u(:,kk)/rho), adjust_max, T, delta, n-1)
      else
          if( T>1 )then
              call qSortvec(A,T)
          end if
      end if
	  
      endmark = .true.
      j = 1
      Do j = 1,T
        if(Alast(j) .ne. A(j))then
            endmark = .false.
            EXIT
        end if
      end do

      if(endmark)then
        iter = iter + 1
        EXIT
      else
        iter = iter + 1
      end if

      Alast(1:T) = A(1:T)
    end do
    beta(:,kk) = tempvec
    miter(kk) = iter - 1
    kk = kk + 1
    if( dot_product(y-tempvec,y-tempvec)<n*eps .or. T>=L )EXIT

  end do

  Kmax = kk-1
  deallocate(tempvec,tempDTD,A,I,ordervec,Alast)
  return
end
subroutine btfusedl0(n, y, beta, z, u, A, I, Anull, T, rho, iter, itermax, inv, vec, veck, adjust, adjust_max, delta)

  implicit none

  real(kind = 8), dimension (:), allocatable :: tempvec, tempITI
  real(kind = 8), dimension (:,:), allocatable :: tempDTD,tempATA
  integer, dimension (:), allocatable :: ordervec,Alast, Asort

  integer, intent(in) :: adjust
  integer, intent(in) :: adjust_max
  integer, intent(in) :: delta

  logical setequal
  integer n,iter,itermax,j,T,veck,Anull
  real(kind = 8) y(n),beta(n),z(n-veck+1),u(n-veck+1),rho,inv(2*veck-1,n-veck+1),vec(veck)
  integer A(T),I(n-veck+1-T)

  external DGBSV   ! solve banded linear system (AX=B) function in LAPACK

  allocate(tempvec(1:n-veck+1))
  allocate(tempDTD(1:2*veck-1,1:n-veck+1))
  allocate(tempATA(1:3*veck-2,1:n-veck+1))
  allocate(tempITI(1:n-veck+1))
  allocate(ordervec(1:n-veck+1))
  allocate(Alast(1:T))
  allocate(Asort(1:T))

  tempDTD = 0
  tempITI = 0
  Alast = 0
  Asort = 0
  tempvec = 0
  call dymul(vec,n,veck,y,tempvec)
  !open(2, file='D:/data1.txt', status='old')

  if(Anull == 1) then
      call DTDsub(inv,tempDTD,I(1:n-veck+1),n-veck+1,n-veck+1,veck)
      tempATA = 0
      tempATA(veck:(3*veck-2),:) = tempDTD

      tempITI(1:n-veck+1) = tempvec
      call DGBSV(n-veck+1,veck-1,veck-1,1,tempATA,3*veck-2,ordervec(1:n-veck+1-T),tempITI,n-veck+1,j)
      u= tempITI
      Do j = 1,n-veck+1
          ordervec(j) = j
      end do
      call Qsortorder(abs(z + u/rho),n-veck+1,ordervec)
      A = ordervec(1:T)
      I = ordervec((T+1):n-veck+1)

      if(adjust == 1 .and. T /= 1) then
          call adjustA(A, I, abs(z + u/rho), adjust_max, T, delta, n-veck+1)
      else
          if( T>1 )then
              call qSortvec(A,T)
          end if
      end if 
  else
      if( T>1 )then
          call qSortvec(A,T)
      end if
  end if

  iter = 1

  Do while(iter < itermax)
    if(n-veck+1-T .eq. 1)then
        u(I(1)) = tempvec(I(1))/inv(veck,1)
    else
        call qSortvec(I(1:n-veck+1-T),n-veck+1-T)
        call DTDsub(inv,tempDTD(:,1:n-veck+1-T),I(1:n-veck+1-T),n-veck+1-T,n-veck+1,veck)
        tempATA = 0
        tempATA(veck:(3*veck-2),1:n-veck+1-T) = tempDTD(:,1:n-veck+1-T)

        tempITI(1:n-veck+1-T) = tempvec(I(1:(n-veck+1-T)))
        call DGBSV(n-veck+1-T,veck-1,veck-1,1,tempATA(:,1:n-veck+1-T),3*veck-2,ordervec(1:n-veck+1-T),tempITI,n-veck+1-T,j)
        u(I(1:(n-veck+1-T))) = tempITI(1:n-veck+1-T)

        !write(2, *) u(I(1:(n-veck+1-T)))
    end if
	
    u(A(1:T)) = 0
    call tdymul(vec,n,veck,u,beta)
    beta = y - beta
    call sdymul(vec,n,veck,beta,A,T,z)

    j = 1
    Do j = 1,n-veck+1
        ordervec(j) = j
    end do
    call Qsortorder(abs(z + u/rho),n-veck+1,ordervec)
    A = ordervec(1:T)
    I = ordervec((T+1):n-veck+1)

    if(adjust == 1 .and. T /= 1) then
        call adjustA(A, I, abs(z + u/rho), adjust_max, T, delta, n-veck+1)
    else
        if( T>1 )then
            call qSortvec(A,T)
        end if  
    end if
    

    setequal = .TRUE.
    j = 1
    Do j = 1,T
        if(A(j) .NE. Alast(j))then
            setequal = .FALSE.
            EXIT
        end if
    end do
    if(setequal)then
        iter = iter + 1
        EXIT
    else
        Alast = A
        iter = iter + 1
    end if
    
  end do
  iter = iter - 1

  deallocate(tempvec,tempDTD,tempATA,tempITI,Alast,Asort,ordervec)
  !close(2)
  return
end
subroutine batfusedl0(n, y, beta, z, u, tao, Kmax, L, eps, rho, miter, itermax, inv, vec, veck, adjust, adjust_max, delta)

  implicit none

  real(kind = 8), dimension (:), allocatable :: tempvec,tempbeta,tempITI
  real(kind = 8), dimension (:,:), allocatable :: tempDTD,tempATA
  integer, dimension (:), allocatable :: A,I,ordervec,Alast,Asort

  integer, intent(in) :: adjust
  integer, intent(in) :: adjust_max
  integer, intent(in) :: delta

  logical setequal
  integer n,iter,itermax,tao,Kmax,L,j,T,k,veck
  integer miter(Kmax)
  real(kind = 8) y(n),beta(n,Kmax),z(n-veck+1,Kmax),u(n-veck+1,Kmax),rho,eps,inv(2*veck-1,n-veck+1),vec(veck)

  external DGBSV   ! solve banded linear system (AX=B) function in LAPACK

  allocate(tempvec(1:n-veck+1))
  allocate(tempbeta(1:n))
  allocate(tempDTD(1:2*veck-1,1:n-veck+1))
  allocate(tempATA(1:3*veck-2,1:n-veck+1))
  allocate(tempITI(1:n-veck+1))
  allocate(A(1:L))
  allocate(I(1:n-veck+1))
  allocate(ordervec(1:n-veck+1))
  allocate(Alast(1:L))
  allocate(Asort(1:L))

  tempDTD = 0
  Alast = 0
  tempvec = 0
  tempbeta = 0
  call dymul(vec,n,veck,y,tempvec)
  z = 0
  k = 1
  !open(2, file='D:/data1.txt', status='old')
  Do while( k<Kmax )

    T = tao*k
    j = 1
    if(k == 1) then
        Do j = 1,n-veck+1
          I(j) = j
        end do
        call DTDsub(inv,tempDTD,I(1:n-veck+1),n-veck+1,n-veck+1,veck)
        tempATA = 0
        tempATA(veck:(3*veck-2),:) = tempDTD

        tempITI(1:n-veck+1) = tempvec
        call DGBSV(n-veck+1,veck-1,veck-1,1,tempATA,3*veck-2,ordervec(1:n-veck+1),tempITI,n-veck+1,j)
        u(:,k)= tempITI
        Do j = 1,n-veck+1
          ordervec(j) = j
        end do
        call Qsortorder(abs(z(:,k) + u(:,k)/rho),n-veck+1,ordervec)
    end if
    A(1:T) = ordervec(1:T)
    !write(2, *) A(1:T)
    I(1:n-veck-T+1) = ordervec((T+1):n-veck+1)
    if(k == 1) then
        if(adjust == 1 .and. T /= 1) then
            call adjustA(A, I, abs(z(:,k) + u(:,k)/rho), adjust_max, T, delta, n-veck+1)
        else
            if( T>1 )then
                call qSortvec(A,T)
            end if
        end if 
    else
        if( T>1 )then
            call qSortvec(A,T)
        end if
    end if
	
    iter = 1
    Do while(iter < itermax)

      if(n-veck+1-T .eq. 1)then
        u(I(1),k) = tempvec(I(1))/inv(veck,1)
      else
        tempITI = 0
        j = 1
        call qSortvec(I(1:n-veck+1-T),n-veck+1-T)
        call DTDsub(inv,tempDTD(1:2*veck-1,1:n-veck+1-T),I(1:n-veck+1-T),n-veck+1-T,n-veck+1,veck)
        tempATA = 0
        tempATA(veck:(3*veck-2),1:n-veck+1-T)=tempDTD(1:2*veck-1,1:n-veck+1-T)
        tempITI(1:n-veck+1-T) = tempvec(I(1:(n-veck+1-T)))
        call DGBSV(n-veck+1-T,veck-1,veck-1,1,tempATA(:,1:n-veck+1-T),3*veck-2,ordervec(1:n-veck+1-T), &
        tempITI(1:n-veck+1-T),n-veck+1-T,j)
        u(I(1:(n-veck+1-T)),k) = tempITI(1:n-veck+1-T)
      end if

      u(A(1:T),k) = 0
      call tdymul(vec,n,veck,u(:,k),tempbeta)
      tempbeta = y - tempbeta
      call sdymul(vec,n,veck,tempbeta,A(1:T),T,z(:,k))

      j = 1
      Do j = 1,n-veck+1
          ordervec(j) = j
      end do
      call Qsortorder(abs(z(:,k) + u(:,k)/rho),n-veck+1,ordervec)
      A(1:T) = ordervec(1:T)
      Asort(1:T) = A(1:T)
      I(1:n-veck+1-T) = ordervec((T+1):n-veck+1)
      if(adjust == 1 .and. T /= 1) then
          call adjustA(Asort, I, abs(z(:,k) + u(:,k)/rho), adjust_max, T, delta, n-veck+1)
          A = Asort
      else
          if( T>1 )then
              call qSortvec(Asort,T)
          end if
      end if 
      setequal = .TRUE.
      j = 1
      Do j = 1,T
          if(Asort(j) .NE. Alast(j))then
            setequal = .FALSE.
            EXIT
          end if
      end do

      if(setequal)then
        iter = iter + 1
        EXIT
      else
        Alast(1:T) = Asort(1:T)
        iter = iter + 1
      end if

    end do
      beta(:,k) = tempbeta
      miter(k) = iter - 1
      k = k + 1
      if( dot_product(y-tempbeta,y-tempbeta)<n*eps .or. T>=L )EXIT
    End do

  Kmax = k-1
  !close(2)
  deallocate(tempvec,tempbeta,tempDTD,tempATA,tempITI,A,I,ordervec,Alast,Asort)
  return
end
subroutine gfusedl0(n, y, beta, z, u, A, I, Anull, T, rho, iter, itermax, inv, D, m, adjust, adjust_max, delta)

  implicit none

  real(kind = 8), dimension (:), allocatable :: tempvec, tempITI
  real(kind = 8), dimension (:,:), allocatable :: tempATA
  integer, dimension (:), allocatable :: ordervec,Alast,Asort

  integer, intent(in) :: adjust
  integer, intent(in) :: adjust_max
  integer, intent(in) :: delta
  
  logical setequal
  integer n,iter,itermax,j,T,m,Anull
  real(kind = 8) y(n),beta(n),z(m),u(m),rho,inv(m,m),D(m,n)
  integer A(T), I(m-T)

  external DGESV   ! solve general linear system (AX=B) function in LAPACK

  allocate(tempvec(1:m))
  allocate(tempITI(1:m))
  allocate(tempATA(1:m,1:m))
  allocate(ordervec(1:m))
  allocate(Alast(1:T))
  allocate(Asort(1:T))

  tempATA = 0
  Alast = 0
  tempvec = 0
  tempvec = matmul(D,y)
  !open(2, file='D:/data2.txt', status='old')

  if(Anull == 1) then
      tempITI(1:m) = tempvec(1:m)
      tempATA = inv
      !call DGBSV(m,3,3,1,tempATA,m,ordervec,tempITI,m,j)
      call DGESV(m,1,tempATA,m,ordervec,tempITI,m,j)
      u = tempITI(1:m)
      Do j = 1,m
        ordervec(j) = j
      end do
      call Qsortorder(abs(z + u/rho),m,ordervec)
      A = ordervec(1:T)
      I = ordervec((T+1):m)
      if(adjust == 1 .and. T /= 1) then
          call adjustA(A, I, abs(z + u/rho), adjust_max, T, delta, m)
      else
          if( T>1 )then
              call qSortvec(A,T)
          end if
      end if 
  else
      if( T>1 )then
          call qSortvec(A,T)
      end if  
  end if

  iter = 1
  Do while(iter < itermax)
    if(m-T .eq. 1)then
        u(I(1)) = tempvec(I(1))/inv(I(1),I(1))
    else
        j = 1
        tempITI(1:m-T) = tempvec(I(1:(m-T)))
        tempATA(1:m-T,1:m-T) = inv(I(1:m-T),I(1:m-T))
        !call DGBSV(m-T,3,3,1,tempATA(1:m-T,1:m-T),m-T,ordervec(1:m-T),tempITI,m-T,j)
        call DGESV(m-T,1,tempATA(1:m-T,1:m-T),m-T,ordervec(1:m-T),tempITI(1:m-T),m-T,j)
        u(I(1:(m-T))) = tempITI(1:m-T)
        !write(2, *) u
    end if
    u(A(1:T)) = 0
    call gtdymul(D,n,m,u,beta)
    beta = y - beta
    z = 0
    if(T .eq. 1)then
        z(A(1)) = dot_product(D(A(1),:),beta)
    else
        z(A(1:T)) = matmul(D(A(1:T),:),beta)
    end if

    j = 1
    Do j = 1,m
        ordervec(j) = j
    end do
    call Qsortorder(abs(z + u/rho),m,ordervec)
    A = ordervec(1:T)
    Asort = A
    I = ordervec((T+1):m)
    if(adjust == 1 .and. T /= 1) then
        call adjustA(Asort, I, abs(z + u/rho), adjust_max, T, delta, m)
        A = Asort
    else
        if( T>1 )then
            call qSortvec(Asort,T)
        end if
    end if 
    if(T>1)call qSortvec(Asort(1:T),T)
    setequal = .TRUE.
    j = 1
    Do j = 1,T
        if(Asort(j) .NE. Alast(j))then
            setequal = .FALSE.
            EXIT
        end if
    end do

    if(setequal)then
        iter = iter + 1
        EXIT
    else
        Alast = Asort
        iter = iter + 1
    end if

  end do
  !close(2)
  iter = iter - 1
  deallocate(tempvec,tempITI,tempATA,ordervec,Alast,Asort)
  return
end
subroutine agfusedl0(n, y, beta, z, u, tao, Kmax, L, eps, rho, miter, itermax, inv, D, m, adjust, adjust_max, delta)

  implicit none

  real(kind = 8), dimension (:), allocatable :: tempvec,tempbeta,tempITI
  real(kind = 8), dimension (:,:), allocatable :: tempATA
  integer, dimension (:), allocatable :: A,I,ordervec,Alast,Asort

  integer, intent(in) :: adjust
  integer, intent(in) :: adjust_max
  integer, intent(in) :: delta
  
  logical setequal
  integer n,iter,itermax,j,T,m,tao,Kmax,L,k
  integer miter(Kmax)
  real(kind = 8) y(n),beta(n,Kmax),z(m,Kmax),u(m,Kmax),rho,inv(m,m),D(m,n),eps

  external DGESV   ! solve general linear system (AX=B) function in LAPACK

  allocate(tempvec(1:m))
  allocate(tempbeta(1:n))
  allocate(tempATA(1:m,1:m))
  allocate(tempITI(1:m))
  allocate(A(1:L))
  allocate(I(1:m))
  allocate(ordervec(1:m))
  allocate(Alast(1:L))
  allocate(Asort(1:L))

  tempbeta = 0
  tempATA = 0
  tempITI = 0
  Alast = 0
  tempvec = 0
  tempvec = matmul(D,y)
  z = 0
  !open(2, file='D:/data1.txt', status='old')
  k = 1
  Do while(k < Kmax)

    T=tao*k
    if(k == 1) then
		tempITI(1:m) = tempvec(1:m)
        tempATA = inv
        call DGESV(m,1,tempATA,m,ordervec,tempITI,m,j)
		u(:,k)= tempITI
        Do j = 1,m
          ordervec(j) = j
        end do
        call Qsortorder(abs(z(:,k) + u(:,k)/rho),m,ordervec)
    end if
    A(1:T) = ordervec(1:T)
    I(1:m-T) = ordervec((T+1):m)
	!write(2,*) A(1:T)
    if(k == 1) then
        if(adjust == 1 .and. T /= 1) then
            call adjustA(A, I, abs(z(:,k) + u(:,k)/rho), adjust_max, T, delta, m)
        else
            if( T>1 )then
                call qSortvec(A,T)
            end if
        end if 
    else
        if( T>1 )then
            call qSortvec(A,T)
        end if
    end if
    iter = 1
    Do while(iter < itermax)

      if(m-T .eq. 1)then
        u(I(1),k) = tempvec(I(1))/inv(I(1),I(1))
      else
	    j = 1
        tempATA(1:m-T,1:m-T) = inv(I(1:m-T),I(1:m-T))
		tempITI(1:m-T) = tempvec(I(1:(m-T)))
        call DGESV(m-T,1,tempATA(1:m-T,1:m-T),m-T,ordervec(1:m-T),tempITI(1:m-T),m-T,j)
        u(I(1:(m-T)),k) = tempITI(1:m-T)
      end if
      u(A(1:T),k) = 0
      call gtdymul(D,n,m,u(:,k),tempbeta)
      tempbeta = y - tempbeta
      z(:,k) = 0
      if(T .eq. 1)then
        z(A(1),k) = dot_product(D(A(1),:),tempbeta)
      else
        z(A(1:T),k) = matmul(D(A(1:T),:),tempbeta)
      end if

      j = 1
      Do j = 1,m
        ordervec(j) = j
      end do
      call Qsortorder(abs(z(:,k) + u(:,k)/rho),m,ordervec)
      A(1:T) = ordervec(1:T)
      Asort(1:T) = A(1:T)
      I(1:m-T) = ordervec((T+1):m)
      if(adjust == 1 .and. T /= 1) then
          call adjustA(Asort, I, abs(z(:,k) + u(:,k)/rho), adjust_max, T, delta, m)
          A = Asort
      else
          if( T>1 )then
              call qSortvec(Asort,T)
          end if
      end if 

      setequal = .TRUE.
      j = 1
      Do j = 1,T
        if(Asort(j) .NE. Alast(j))then
            setequal = .FALSE.
            EXIT
        end if
      end do

      if(setequal)then
        iter = iter + 1
        EXIT
      else
        Alast = Asort
        iter = iter + 1
      end if

    end do

    beta(:,k) = tempbeta
    miter(k) = iter - 1
    k = k + 1
    if( dot_product(y-tempbeta,y-tempbeta)<n*eps .or. T>=L )EXIT
  end do

  Kmax = k - 1
  !close(2)
  deallocate(tempvec,tempbeta,tempATA,tempITI,A,I,ordervec,Alast,Asort)
  return
end
subroutine comggfusedl0(n,y,beta,alpha,u,T,rho,iter,itermax,inv,gamma,v,D,m,invw,W,mw,Tw,rhow)

  implicit none

  real(kind = 8), dimension (:), allocatable :: tempvec,tdu,twv
  real(kind = 8), dimension (:,:), allocatable :: tempATA,tempITI
  integer, dimension (:), allocatable :: A,Aw,I,Iw,ordervec,Alast,Alastw,Asort,Asortw

  logical setequal,setequalw
  integer n,iter,itermax,j,T,Tw,m,mw,lm,lmt
  real(kind = 8) y(n),beta(n),alpha(m),gamma(mw),u(m),v(mw),rho,rhow,inv(m,m),invw(mw,mw),D(m,n),W(mw,n)

  external DGESV   ! solve general linear system (AX=B) function in LAPACK

  lmt = max(m-T,mw-Tw)
  lm = max(m,mw)
  allocate(tempvec(1:lmt))
  allocate(tdu(1:n))
  allocate(twv(1:n))
  allocate(tempATA(1:lmt,1:lmt))
  allocate(tempITI(1:lmt,1:lmt))
  allocate(A(1:T))
  allocate(Aw(1:Tw))
  allocate(I(1:m-T))
  allocate(Iw(1:mw-Tw))
  allocate(ordervec(1:lm))
  allocate(Alast(1:T))
  allocate(Asort(1:T))
  allocate(Alastw(1:Tw))
  allocate(Asortw(1:Tw))

  tempvec = 0
  tdu = 0 
  twv = 0
  tempATA = 0
  tempITI = 0
  A = 0
  Aw = 0
  I = 0
  Iw = 0
  ordervec = 0
  Alast = 0
  Alastw = 0
  Asort = 0
  Asortw = 0
  alpha = matmul(D,y)
  gamma = matmul(W,y)
  u = 0
  v = 0

  Do j = 1,m
    ordervec(j) = j
  end do
  call Qsortorder(abs(alpha + u/rho),m,ordervec(1:m))
  A = ordervec(1:T)
  I = ordervec((T+1):m)
  j = 1
  Do j = 1,mw
    ordervec(j) = j
  end do
  call Qsortorder(abs(gamma + v/rhow),mw,ordervec(1:mw))
  Aw = ordervec(1:Tw)
  Iw = ordervec((Tw+1):mw)

  iter = 1

  DO WHILE(iter<itermax)
    ! the update of coff related to D
    tempvec(1:m-T) = matmul(D(I(1:(m-T)),:),(y-twv))
    if(m-T .eq. 1)then
        u(I(1)) = tempvec(1)/inv(I(1),I(1))
    else
        tempITI = 0
        j = 1
        DO j = 1,m-T
            tempITI(j,j) = 1
        END DO
        tempATA(1:m-T,1:m-T) = inv(I(1:m-T),I(1:m-T))
        call DGESV(m-T,m-T,tempATA(1:m-T,1:m-T),m-T,ordervec(1:m-T),tempITI(1:m-T,1:m-T),m-T,j)
        u(I(1:(m-T))) = matmul(tempITI(1:m-T,1:m-T),tempvec(1:m-T))
    end if
    u(A(1:T)) = 0
    call gtdymul(D,n,m,u,tdu)
    beta = y-tdu-twv
    alpha = 0
    if(T .eq. 1)then
        alpha(A(1)) = dot_product(D(A(1),:),beta)
    else
        alpha(A(1:T)) = matmul(D(A(1:T),:),beta)
    end if
    ! the update of coff related to W
    tempvec(1:mw-Tw) = matmul(W(Iw(1:(mw-Tw)),:),(y-tdu))
    if(mw-Tw .eq. 1)then
        v(Iw(1:(mw-Tw))) = tempvec(1)/invw(Iw(1),Iw(1))
    else
        tempITI = 0
        j = 1
        DO j = 1,mw-Tw
            tempITI(j,j) = 1
        END DO
        tempATA(1:mw-Tw,1:mw-Tw) = invw(Iw(1:mw-Tw),Iw(1:mw-Tw))
        call DGESV(mw-Tw,mw-Tw,tempATA(1:mw-Tw,1:mw-Tw),mw-Tw,ordervec(1:mw-Tw),tempITI(1:mw-Tw,1:mw-Tw),mw-Tw,j)
        v(Iw(1:(mw-Tw))) = matmul(tempITI(1:mw-Tw,1:mw-Tw),tempvec(1:mw-Tw))
    end if
    v(Aw(1:Tw)) = 0
    call gtdymul(W,n,mw,v,twv)
    beta = y-tdu-twv
    gamma = 0
    if(Tw .eq. 1)then
        gamma(Aw(1)) = dot_product(W(Aw(1),:),beta)
    else
        gamma(Aw(1:Tw)) = matmul(W(Aw(1:Tw),:),beta)
    end if

    j = 1
    Do j = 1,m
      ordervec(j) = j
    end do
    call Qsortorder(abs(alpha + u/rho),m,ordervec(1:m))
    A(1:T) = ordervec(1:T)
    Asort(1:T) = A(1:T)
    I(1:m-T) = ordervec((T+1):m)
    j = 1
    Do j = 1,mw
      ordervec(j) = j
    end do
    call Qsortorder(abs(gamma + v/rhow),mw,ordervec(1:mw))
    Aw(1:Tw) = ordervec(1:Tw)
    Asortw(1:Tw) = Aw(1:Tw)
    Iw(1:mw-Tw) = ordervec((Tw+1):mw)

    if(T>1)call qSortvec(Asort(1:T),T)
    if(Tw>1)call qSortvec(Asortw(1:Tw),Tw)
    setequal = .TRUE.
    j = 1
    Do j = 1,T
        if(Asort(j) .NE. Alast(j))then
            setequal = .FALSE.
            EXIT
        end if
    end do
    setequalw = .TRUE.
    j = 1
    Do j = 1,Tw
        if(Asortw(j) .NE. Alastw(j))then
            setequalw = .FALSE.
            EXIT
        end if
    end do
    if(setequal .and. setequalw)then
        iter = iter + 1
        EXIT
    else
        Alast(1:T) = Asort(1:T)
        Alastw(1:Tw) = Asortw(1:Tw)
        iter = iter + 1
    end if
  end do

  iter = iter - 1
  deallocate(tempvec,tdu,twv,tempATA,tempITI,A,Aw,I,Iw,ordervec,Alast,Alastw,Asort,Asortw)
  return
end
subroutine acomggfusedl0(n,y,beta,alpha,u,tao,Kmax,L,eps,rho,miter,itermax,gamma,v,inv,D,m,invw,W,mw,Tw,rhow)

  implicit none

  real(kind = 8), dimension (:), allocatable :: tempvec,tdu,twv,tempbeta
  real(kind = 8), dimension (:,:), allocatable :: tempATA,tempITI
  integer, dimension (:), allocatable :: A,Aw,I,Iw,ordervec,Alast,Alastw,Asort,Asortw

  logical setequal,setequalw
  integer n,iter,itermax,j,T,Tw,m,mw,lm,lmt,tao,Kmax,L,miter(Kmax),k
  real(kind = 8) y(n),beta(n,Kmax),alpha(m),gamma(mw),u(m),v(mw),rho,rhow,inv(m,m),invw(mw,mw),D(m,n),W(mw,n),eps

  external DGESV   ! solve general linear system (AX=B) function in LAPACK

  lmt = max(m-1,mw-Tw)
  lm = max(m,mw)
  allocate(tempbeta(1:n))
  allocate(tempvec(1:lmt))
  allocate(tdu(1:n))
  allocate(twv(1:n))
  allocate(tempATA(1:lmt,1:lmt))
  allocate(tempITI(1:lmt,1:lmt))
  allocate(A(1:L))
  allocate(Aw(1:Tw))
  allocate(I(1:m-1))
  allocate(Iw(1:mw-Tw))
  allocate(ordervec(1:lm))
  allocate(Alast(1:L))
  allocate(Asort(1:L))
  allocate(Alastw(1:Tw))
  allocate(Asortw(1:Tw))

  tempbeta = 0
  tempvec = 0
  tdu = 0 
  twv = 0
  tempATA = 0
  tempITI = 0
  A = 0
  Aw = 0
  I = 0
  Iw = 0
  ordervec = 0
  Alast = 0
  Alastw = 0
  Asort = 0
  Asortw = 0
  alpha = matmul(D,y)
  gamma = matmul(W,y)
  u = 0
  v = 0

  k = 1
  Do while(k < Kmax)

    T=tao*k
    Do j = 1,m
      ordervec(j) = j
    end do
    call Qsortorder(abs(alpha + u/rho),m,ordervec(1:m))
    A(1:T) = ordervec(1:T)
    I(1:m-T) = ordervec((T+1):m)
    j = 1
    Do j = 1,mw
      ordervec(j) = j
    end do
    call Qsortorder(abs(gamma + v/rhow),mw,ordervec(1:mw))
    Aw(1:Tw) = ordervec(1:Tw)
    Iw(1:mw-Tw) = ordervec((Tw+1):mw)
  
    iter = 1

    DO WHILE(iter<itermax)
      ! the update of coff related to D
      tempvec(1:m-T) = matmul(D(I(1:(m-T)),:),(y-twv))
      if(m-T .eq. 1)then
        u(I(1)) = tempvec(1)/inv(I(1),I(1))
      else
        tempITI = 0
        j = 1
        DO j = 1,m-T
            tempITI(j,j) = 1
        END DO
        tempATA(1:m-T,1:m-T) = inv(I(1:m-T),I(1:m-T))
        call DGESV(m-T,m-T,tempATA(1:m-T,1:m-T),m-T,ordervec(1:m-T),tempITI(1:m-T,1:m-T),m-T,j)
        u(I(1:(m-T))) = matmul(tempITI(1:m-T,1:m-T),tempvec(1:m-T))
      end if
      u(A(1:T)) = 0
      call gtdymul(D,n,m,u,tdu)
      tempbeta = y-tdu-twv
      alpha = 0
      if(T .eq. 1)then
        alpha(A(1)) = dot_product(D(A(1),:),tempbeta)
      else
        alpha(A(1:T)) = matmul(D(A(1:T),:),tempbeta)
      end if
      ! the update of coff related to W
      tempvec(1:mw-Tw) = matmul(W(Iw(1:(mw-Tw)),:),(y-tdu))
      if(mw-Tw .eq. 1)then
        v(Iw(1)) = tempvec(1)/invw(Iw(1),Iw(1))
      else
        tempITI = 0
        j = 1
        DO j = 1,mw-Tw
            tempITI(j,j) = 1
        END DO
        tempATA(1:mw-Tw,1:mw-Tw) = invw(Iw(1:mw-Tw),Iw(1:mw-Tw))
        call DGESV(mw-Tw,mw-Tw,tempATA(1:mw-Tw,1:mw-Tw),mw-Tw,ordervec(1:mw-Tw),tempITI(1:mw-Tw,1:mw-Tw),mw-Tw,j)
        v(Iw(1:(mw-Tw))) = matmul(tempITI(1:mw-Tw,1:mw-Tw),tempvec(1:mw-Tw))
      end if
      v(Aw(1:Tw)) = 0
      call gtdymul(W,n,mw,v,twv)
      tempbeta = y-tdu-twv
      gamma = 0
      if(Tw .eq. 1)then
        gamma(Aw(1)) = dot_product(W(Aw(1),:),tempbeta)
      else
        gamma(Aw(1:Tw)) = matmul(W(Aw(1:Tw),:),tempbeta)
      end if

      j = 1
      Do j = 1,m
        ordervec(j) = j
      end do
      call Qsortorder(abs(alpha + u/rho),m,ordervec(1:m))
      A(1:T) = ordervec(1:T)
      Asort(1:T) = A(1:T)
      I(1:m-T) = ordervec((T+1):m)
      j = 1
      Do j = 1,mw
        ordervec(j) = j
      end do
      call Qsortorder(abs(gamma + v/rhow),mw,ordervec(1:mw))
      Aw(1:Tw) = ordervec(1:Tw)
      Asortw(1:Tw) = Aw(1:Tw)
      Iw(1:mw-Tw) = ordervec((Tw+1):mw)

      if(T>1)call qSortvec(Asort(1:T),T)
      if(Tw>1)call qSortvec(Asortw(1:Tw),Tw)
      setequal = .TRUE.
      j = 1
      Do j = 1,T
        if(Asort(j) .NE. Alast(j))then
            setequal = .FALSE.
            EXIT
        end if
      end do
      setequalw = .TRUE.
      j = 1
      Do j = 1,Tw
        if(Asortw(j) .NE. Alastw(j))then
            setequalw = .FALSE.
            EXIT
        end if
      end do
        if(setequal .and. setequalw)then
        iter = iter + 1
        EXIT
      else
        Alast(1:T) = Asort(1:T)
        Alastw(1:Tw) = Asortw(1:Tw)
        iter = iter + 1
      end if
    end do

    beta(:,k) = tempbeta
    miter(k) = iter - 1
    k = k + 1
    if( dot_product(y-tempbeta,y-tempbeta)<n*eps .or. T>=L )EXIT
  End do

  Kmax = k-1
  deallocate(tempbeta,tempvec,tdu,twv,tempATA,tempITI,A,Aw,I,Iw,ordervec,Alast,Alastw,Asort,Asortw)
  return
end
subroutine comdifusedl0(n,y,beta,alpha,u,T,rho,iter,itermax,gamma,v,Tw,rhow)

  implicit none

  real(kind = 8), dimension (:), allocatable :: tempvec,tdu,twv,tempDTD
  integer, dimension (:), allocatable :: A,Alast,Aw,Alastw,I,Iw,ordervec

  logical setequal,setequalw
  integer n,iter,itermax,j,T,Tw,k
  real(kind = 8) y(n),beta(n),alpha(n-1),gamma(n),u(n-1),v(n),rho,rhow

  allocate(tempvec(1:n))
  allocate(tdu(1:n))
  allocate(twv(1:n))
  allocate(A(1:T))
  allocate(Alast(1:T))
  allocate(Aw(1:Tw))
  allocate(Alastw(1:Tw))
  allocate(I(1:n-1-T))
  allocate(Iw(1:n-Tw))
  allocate(ordervec(1:n))
  allocate(tempDTD(1:n-1))

  tempvec = 0
  tdu = 0 
  twv = 0
  A = 0
  Alast = 0
  Aw = 0
  Alastw = 0
  I = 0
  Iw = 0
  ordervec = 0
  alpha = y(2:n) - y(1:n-1)
  gamma = y
  u = 0
  v = 0

  Do j = 1,n-1
    ordervec(j) = j
  end do
  call Qsortorder(abs(alpha + u/rho),n-1,ordervec(1:n-1))
  A = ordervec(1:T)
  I = ordervec((T+1):n-1)
  j = 1
  Do j = 1,n
    ordervec(j) = j
  end do
  call Qsortorder(abs(gamma + v/rhow),n,ordervec(1:n))
  Aw = ordervec(1:Tw)
  Iw = ordervec((Tw+1):n)
  if( T>1 )call qSortvec(A,T)
  if( Tw>1 )call qSortvec(Aw,Tw)

  iter = 1
  DO WHILE(iter<itermax)
    ! the update of coff related to D
    if(T .eq. 1)then
        if(A(1) .eq. 1 )then
            tempvec(1:(n-2)) = y(3:n) - twv(3:n) - y(2:(n-1)) + twv(2:(n-1))
            Do j = 1, n-2
                call itergen(j,n-2,tempDTD(1:n-2))
                u(j+1) = dot_product(tempDTD(1:(n-2)),tempvec(1:(n-2)))
            end do
        else if(A(1) .eq. n-1)then
            tempvec(1:(n-2)) = y(2:(n-1)) - twv(2:(n-1)) - y(1:(n-2)) + twv(1:(n-2))
            Do j = 1, n-2
                call itergen(j,n-2,tempDTD(1:n-2))
                u(j) = dot_product(tempDTD(1:(n-2)),tempvec(1:(n-2)))
            end do
        else
            tempvec(1:(A(1)-1)) = y(2:A(1)) - twv(2:A(1)) - y(1:(A(1)-1)) + twv(1:(A(1)-1))
            Do j = 1, A(1)-1
                call itergen(j,A(1)-1,tempDTD(1:A(1)-1))
                u(j) = dot_product(tempDTD(1:A(1)-1),tempvec(1:A(1)-1))
            end do
            j = 1
            tempvec(1:(n-1-A(1))) = y(A(1)+2:n) - twv(A(1)+2:n) - y(A(1)+1:(n-1)) + twv(A(1)+1:(n-1))
            Do j = 1, n-1-A(1)
                call itergen(j,n-1-A(1),tempDTD(1:n-1-A(1)))
                u(j+A(1)) = dot_product(tempDTD(1:n-1-A(1)),tempvec(1:n-1-A(1)))
            end do
        end if
    else
        if( A(1) > 1 )then
            tempvec(1:(A(1)-1)) = y(2:A(1)) - twv(2:A(1)) - y(1:(A(1)-1)) + twv(1:(A(1)-1))
            Do j = 1, A(1)-1
                call itergen(j,A(1)-1,tempDTD(1:A(1)-1))
                u(j) = dot_product(tempDTD(1:A(1)-1),tempvec(1:A(1)-1))
            end do
        end if
        Do k = 1,T-1
            if( A(k+1)-A(k)>1 )then
                j = 1
                tempvec(1:(A(k+1)-A(k)-1)) = y(A(k)+2:A(k+1)) - twv(A(k)+2:A(k+1)) - y(A(k)+1:(A(k+1)-1)) &
                + twv(A(k)+1:(A(k+1)-1))
                Do j = 1, A(k+1)-A(k)-1
                    call itergen(j,A(k+1)-A(k)-1,tempDTD(1:A(k+1)-A(k)-1))
                    u(j+A(k)) = dot_product(tempDTD(1:A(k+1)-A(k)-1),tempvec(1:A(k+1)-A(k)-1))
                end do
            end if
        end do
        if( A(T) < n-1 )then
            j = 1
            tempvec(1:(n-1-A(T))) = y(A(T)+2:n) - twv(A(T)+2:n) - y(A(T)+1:(n-1)) + twv(A(T)+1:(n-1))
            Do j = 1, n-1-A(T)
                call itergen(j,n-1-A(T),tempDTD(1:n-1-A(T)))
                u(j+A(T)) = dot_product(tempDTD(1:n-1-A(T)),tempvec(1:n-1-A(T)))
            end do
        end if
    end if
    u(A(1:T)) = 0
    tdu(1) = -u(1)
    tdu(2:n-1) = u(1:n-2) - u(2:n-1)
    tdu(n) = u(n-1)
    beta = y-tdu-twv
    alpha = 0
    alpha(A(1:T)) = beta(A(1:T)+1)-beta(A(1:T))
    ! the update of coff related to W
    tempvec(1:n-Tw) = y(Iw(1:(n-Tw))) - tdu(Iw(1:(n-Tw)))
    v = 0
    v(Iw(1:(n-Tw))) = tempvec(1:n-Tw)
    twv = v
    beta = y-tdu-twv
    gamma = 0
    if(Tw .eq. 1)then
        gamma(Aw(1)) = beta(Aw(1))
    else
        gamma(Aw(1:Tw)) = beta(Aw(1:Tw))
    end if

    j = 1
    Do j = 1,n-1
      ordervec(j) = j
    end do
    call Qsortorder(abs(alpha + u/rho),n-1,ordervec(1:n-1))
    A(1:T) = ordervec(1:T)
    I(1:n-1-T) = ordervec((T+1):n-1)
    j = 1
    Do j = 1,n
      ordervec(j) = j
    end do
    call Qsortorder(abs(gamma + v/rhow),n,ordervec(1:n))
    Aw(1:Tw) = ordervec(1:Tw)
    Iw(1:n-Tw) = ordervec((Tw+1):n)

    if(T>1)call qSortvec(A(1:T),T)
    if(Tw>1)call qSortvec(Aw(1:Tw),Tw)
    setequal = .TRUE.
    j = 1
    Do j = 1,T
        if(A(j) .NE. Alast(j))then
            setequal = .FALSE.
            EXIT
        end if
    end do
    setequalw = .TRUE.
    j = 1
    Do j = 1,Tw
        if(Aw(j) .NE. Alastw(j))then
            setequalw = .FALSE.
            EXIT
        end if
    end do
    if(setequal .and. setequalw)then
        iter = iter + 1
        EXIT
    else
        Alast(1:T) = A(1:T)
        Alastw(1:Tw) = Aw(1:Tw)
        iter = iter + 1
    end if
  end do

  iter = iter - 1
  deallocate(tempvec,tdu,twv,A,Aw,I,Iw,ordervec,Alast,Alastw,tempDTD)
  return
end
subroutine acomdifusedl0(n,y,beta,alpha,u,tao,Kmax,L,eps,rho,miter,itermax,gamma,v,Tw,rhow)

  implicit none

  real(kind = 8), dimension (:), allocatable :: tempvec,tdu,twv,tempbeta,tempDTD
  integer, dimension (:), allocatable :: A,Aw,I,Iw,ordervec,Alast,Alastw

  logical setequal,setequalw
  integer n,iter,itermax,j,T,Tw,tao,Kmax,L,miter(Kmax),k,kk
  real(kind = 8) y(n),beta(n,Kmax),alpha(n-1),gamma(n),u(n-1),v(n),rho,rhow,eps

  external DGESV   ! solve general linear system (AX=B) function in LAPACK

  allocate(tempDTD(1:n-1))
  allocate(tempbeta(1:n))
  allocate(tempvec(1:n))
  allocate(tdu(1:n))
  allocate(twv(1:n))
  allocate(A(1:L))
  allocate(Aw(1:Tw))
  allocate(I(1:n-1))
  allocate(Iw(1:n-Tw))
  allocate(ordervec(1:n))
  allocate(Alast(1:L))
  allocate(Alastw(1:Tw))

  tempbeta = 0
  tempvec = 0
  tdu = 0 
  twv = 0
  A = 0
  Aw = 0
  I = 0
  Iw = 0
  ordervec = 0
  Alast = 0
  Alastw = 0
  alpha = y(2:n) - y(1:n-1)
  gamma = y
  u = 0
  v = 0

  kk = 1
  Do while(kk < Kmax)

    T=tao*kk
    Do j = 1,n-1
      ordervec(j) = j
    end do
    call Qsortorder(abs(alpha + u/rho),n-1,ordervec(1:n-1))
    A(1:T) = ordervec(1:T)
    I(1:n-1-T) = ordervec((T+1):n-1)
    j = 1
    Do j = 1,n
      ordervec(j) = j
    end do
    call Qsortorder(abs(gamma + v/rhow),n,ordervec(1:n))
    Aw(1:Tw) = ordervec(1:Tw)
    Iw(1:n-Tw) = ordervec((Tw+1):n)
    if(T>1)call qSortvec(A(1:T),T)
    if(Tw>1)call qSortvec(Aw(1:Tw),Tw)
  
    iter = 1
    DO WHILE(iter<itermax)
      ! the update of coff related to D
      if(T .eq. 1)then
        if(A(1) .eq. 1 )then
            tempvec(1:(n-2)) = y(3:n) - twv(3:n) - y(2:(n-1)) + twv(2:(n-1))
            Do j = 1, n-2
                call itergen(j,n-2,tempDTD(1:n-2))
                u(j+1) = dot_product(tempDTD(1:(n-2)),tempvec(1:(n-2)))
            end do
        else if(A(1) .eq. n-1)then
            tempvec(1:(n-2)) = y(2:(n-1)) - twv(2:(n-1)) - y(1:(n-2)) + twv(1:(n-2))
            Do j = 1, n-2
                call itergen(j,n-2,tempDTD(1:n-2))
                u(j) = dot_product(tempDTD(1:(n-2)),tempvec(1:(n-2)))
            end do
        else
            tempvec(1:(A(1)-1)) = y(2:A(1)) - twv(2:A(1)) - y(1:(A(1)-1)) + twv(1:(A(1)-1))
            Do j = 1, A(1)-1
                call itergen(j,A(1)-1,tempDTD(1:A(1)-1))
                u(j) = dot_product(tempDTD(1:A(1)-1),tempvec(1:A(1)-1))
            end do
            j = 1
            tempvec(1:(n-1-A(1))) = y(A(1)+2:n) - twv(A(1)+2:n) - y(A(1)+1:(n-1)) + twv(A(1)+1:(n-1))
            Do j = 1, n-1-A(1)
                call itergen(j,n-1-A(1),tempDTD(1:n-1-A(1)))
                u(j+A(1)) = dot_product(tempDTD(1:n-1-A(1)),tempvec(1:n-1-A(1)))
            end do
        end if
      else
        if( A(1) > 1 )then
            tempvec(1:(A(1)-1)) = y(2:A(1)) - twv(2:A(1)) - y(1:(A(1)-1)) + twv(1:(A(1)-1))
            Do j = 1, A(1)-1
                call itergen(j,A(1)-1,tempDTD(1:A(1)-1))
                u(j) = dot_product(tempDTD(1:A(1)-1),tempvec(1:A(1)-1))
            end do
        end if
        Do k = 1,T-1
            if( A(k+1)-A(k)>1 )then
                j = 1
                tempvec(1:(A(k+1)-A(k)-1)) = y(A(k)+2:A(k+1)) - twv(A(k)+2:A(k+1)) - y(A(k)+1:(A(k+1)-1)) &
                + twv(A(k)+1:(A(k+1)-1))
                Do j = 1, A(k+1)-A(k)-1
                    call itergen(j,A(k+1)-A(k)-1,tempDTD(1:A(k+1)-A(k)-1))
                    u(j+A(k)) = dot_product(tempDTD(1:A(k+1)-A(k)-1),tempvec(1:A(k+1)-A(k)-1))
                end do
            end if
        end do
        if( A(T) < n-1 )then
            j = 1
            tempvec(1:(n-1-A(T))) = y(A(T)+2:n) - twv(A(T)+2:n) - y(A(T)+1:(n-1)) + twv(A(T)+1:(n-1))
            Do j = 1, n-1-A(T)
                call itergen(j,n-1-A(T),tempDTD(1:n-1-A(T)))
                u(j+A(T)) = dot_product(tempDTD(1:n-1-A(T)),tempvec(1:n-1-A(T)))
            end do
        end if
      end if
      u(A(1:T)) = 0
      tdu(1) = -u(1)
      tdu(2:n-1) = u(1:n-2) - u(2:n-1)
      tdu(n) = u(n-1)
      tempbeta = y-tdu-twv
      alpha = 0
      alpha(A(1:T)) = tempbeta(A(1:T)+1)-tempbeta(A(1:T))
      ! the update of coff related to W
      tempvec(1:n-Tw) = y(Iw(1:(n-Tw))) - tdu(Iw(1:(n-Tw)))
      v = 0
      v(Iw(1:(n-Tw))) = tempvec(1:n-Tw)
      twv = v
      tempbeta = y-tdu-twv
      gamma = 0
      if(Tw .eq. 1)then
        gamma(Aw(1)) = tempbeta(Aw(1))
      else
        gamma(Aw(1:Tw)) = tempbeta(Aw(1:Tw))
      end if

      j = 1
      Do j = 1,n-1
        ordervec(j) = j
      end do
      call Qsortorder(abs(alpha + u/rho),n-1,ordervec(1:n-1))
      A(1:T) = ordervec(1:T)
      I(1:n-1-T) = ordervec((T+1):n-1)
      j = 1
      Do j = 1,n
        ordervec(j) = j
      end do
      call Qsortorder(abs(gamma + v/rhow),n,ordervec(1:n))
      Aw(1:Tw) = ordervec(1:Tw)
      Iw(1:n-Tw) = ordervec((Tw+1):n)

      if(T>1)call qSortvec(A(1:T),T)
      if(Tw>1)call qSortvec(Aw(1:Tw),Tw)
      setequal = .TRUE.
      j = 1
      Do j = 1,T
        if(A(j) .NE. Alast(j))then
            setequal = .FALSE.
            EXIT
        end if
      end do
      setequalw = .TRUE.
      j = 1
      Do j = 1,Tw
        if(Aw(j) .NE. Alastw(j))then
            setequalw = .FALSE.
            EXIT
        end if
      end do
       if(setequal .and. setequalw)then
        iter = iter + 1
        EXIT
      else
        Alast(1:T) = A(1:T)
        Alastw(1:Tw) = Aw(1:Tw)
        iter = iter + 1
      end if
    end do

    beta(:,kk) = tempbeta
    miter(kk) = iter - 1
    kk = kk + 1
    if( dot_product(y-tempbeta,y-tempbeta)<n*eps .or. T>=L )EXIT
  End do

  Kmax = kk-1
  deallocate(tempbeta,tempvec,tdu,twv,A,Aw,I,Iw,ordervec,Alast,Alastw,tempDTD)
  return
end
recursive subroutine qSort(a,vec,na)

    integer, intent(in) :: nA
    real(kind = 8), dimension(nA), intent(in out) :: A
    integer, dimension(nA), intent(in out) :: vec
 
    integer :: left, right
    real(kind = 8) :: random
    real(kind = 8) :: pivot
    real(kind = 8) :: temp
    integer :: itemp
    integer :: marker
 
    if (nA > 1) then
 
        call random_number(random)
        pivot = A(int(random*real(nA-1))+1) ! random pivor
        left = 0
        right = nA + 1
 
        do while (left < right)
            right = right - 1
            do while (A(right) < pivot)
                right = right - 1
            end do
            left = left + 1
            do while (A(left) > pivot)
                left = left + 1
            end do
            if (left < right) then
                temp = A(left)
                A(left) = A(right)
                A(right) = temp
                itemp = vec(left)
                vec(left) = vec(right)
                vec(right) = itemp
            end if
        end do
 
        if (left == right) then
            marker = left + 1
        else
            marker = left
        end if
 
        call qSort(A(:marker-1),vec(:marker-1),marker-1)
        call qSort(A(marker:),vec(marker:),nA-marker+1)
 
    end if

end subroutine qSort
subroutine qsortorder(vec,nvec,ordervec)

    implicit none

    real(kind = 8), dimension (:), allocatable :: tempvec

    integer, intent(in) :: nvec
    real(kind = 8), dimension(nvec), intent(in) :: vec
    integer, dimension(nvec), intent(in out) :: ordervec

    allocate(tempvec(1:nvec))

    tempvec = vec
    call qSort(tempvec, ordervec, nvec)

    deallocate(tempvec)
    return
end
subroutine adjustA(A, I, vu, adjust_max, T, delta, nvec)

  implicit none

	integer, intent(in) :: adjust_max, T, delta, nvec
	integer, dimension(T), intent(in out) :: A
	real(kind = 8), dimension(nvec), intent(in) :: vu
	integer, dimension(nvec - T), intent(in out) :: I

  logical endmark
	integer :: aj, sumA, sumI, k, j, tmp, tmp2, n_rm, mark, max_index, index
	real(kind = 8) :: max_num, random

  integer, dimension (:), allocatable :: Alast, dA, Akeep, existance, candidate

  allocate(Alast(1:T))
  allocate(dA(1:T-1))
  allocate(Akeep(1:T))
  allocate(existance(1:nvec))
  allocate(candidate(1:nvec))
	
	Alast = A
	I = 0
	dA = 0
	Akeep = 0
	existance = 0
	candidate = 0
	sumA = 1
	call qSortvec(A(1:T), T)
	Akeep(1) = A(1)

	Do k = 1, T-1
      dA(k) = A(k+1) - A(k)
		  if(dA(k) > delta) then
		      sumA = sumA + 1
		      Akeep(sumA) = A(k+1)
	end if
  end do
	aj = 0
 ! open(1, file='D:/data.txt', status='old') 

	Do while(T-sumA /= 0 .and. aj < adjust_max)
	    Do k = 1, sumA
			    max_num = -1.0
			    max_index = Akeep(k)
			    tmp = Akeep(k) - delta
          Do j = 1, 2*delta+1
              tmp2 = tmp
		          if(tmp < 1) then
		              tmp2 = 1
			        else if(tmp > nvec) then
				          tmp2 = nvec
			        end if
			        if(vu(tmp2) > max_num) then
				          max_num = vu(tmp2)
					        max_index = tmp2
			        end if
			        tmp = tmp + 1
		      end do
			    Akeep(k) = max_index
         ! write(1, *) Akeep(k)
		    end do     
 
		  Do k = 1, sumA
		      tmp = Akeep(k) - delta
		      Do j = 1, 2*delta+1
              tmp2 = tmp
			        if(tmp < 1) then
		              tmp2 = 1
		          else if(tmp > nvec) then
				          tmp2 = nvec
				      end if
			        existance(tmp2) = 1		
			        tmp = tmp + 1
			    end do
		  end do
		  sumI = 0
		  candidate = 0
		  Do k = 1, nvec 
		      if(existance(k) == 0) then
			        sumI = sumI + 1
			        candidate(sumI) = k		  
			    end if
		  end do

		  n_rm = T - sumA
		  mark = 0
		  if(n_rm>0 .and. sumI>0) then
		      Do k = 1, n_rm
			        sumA = sumA + 1
			        max_num = -1
              Do j = 1, sumI
                  if(vu(candidate(j))>max_num) then
                	    max_num = vu(candidate(j))
					            max_index = candidate(j)
                      index = j
				        	end if
				      end do
              ! write(1, *) index, sumI,candidate(index),candidate(sumI), sumA
              if(index /= sumI) then
                  candidate(index) = candidate(sumI)
                 ! write(1, *) 1111, candidate(index),candidate(sumI), max_index, index
              end if
			        Akeep(sumA) = max_index
			        sumI = sumI - 1
			        if(sumI == 0)then
                  EXIT
              end if
              tmp = max_index - delta
		          Do j = 1, 2*delta+1
                  tmp2 = tmp
			            if(tmp < 1) then
		                  tmp2 = 1
		              else if(tmp > nvec) then
				              tmp2 = nvec
				          end if
			            existance(tmp2) = 1		
			            tmp = tmp + 1
		  	      end do
		          sumI = 0
		          candidate = 0
		          Do j = 1, nvec 
		              if(existance(j) == 0) then
			                sumI = sumI + 1
			                candidate(sumI) = j	  
			            end if
		          end do		 
		          if(sumI == 0)then
                  EXIT
              end if
			    end do
		    	if(sumA == T) then
			        Do k = 1, T
		              A(k) = Akeep(k)
		          end do
		      end if
	  	end if
		  call qSortvec(A(1:T), T)

	  	endmark = .true.
		  Do k = 1,T
          if(Alast(k) /= A(k))then
              endmark = .false.
              EXIT
          end if
      end do
      if(endmark)then
          EXIT
      else
          Alast = A
      end if
      sumA = 1
	    Akeep = 0
	    Akeep(1) = A(1)
	    Do k = 1,T-1
          dA(k) = A(k+1) - A(k)
		      if(dA(k) > delta) then
		          sumA = sumA + 1
		          Akeep(sumA) = A(k+1)
		      end if
      end do        
      aj = aj + 1
	end do
	Do k = 1, nvec
      existance(k) = 0
  end do
  Do k = 1, T
      existance(A(k)) = 1	
  end do
	sumI = 0
	Do k = 1, nvec
	    if(existance(k) == 0) then
		      sumI = sumI + 1
          I(sumI)	= k
      end if
  end do
 ! close(1)
  deallocate(Alast, dA, Akeep, existance, candidate)
  return
end
recursive subroutine qSortvec(vec,nvec)

    integer, intent(in) :: nvec
    integer, dimension(nvec), intent(in out) :: vec
 
    integer :: left, right
    real(kind = 8) :: random
    integer :: pivot
    integer :: itemp
    integer :: marker

    if (nvec > 1) then
 
        call random_number(random)
        pivot = vec(int(random*real(nvec-1))+1) ! random pivor
        left = 0
        right = nvec + 1

        do while (left < right)
            right = right - 1
            do while (vec(right) > pivot)
                right = right - 1
            end do
            left = left + 1
            do while (vec(left) < pivot)
                left = left + 1
            end do
            if (left < right) then
                itemp = vec(left)
                vec(left) = vec(right)
                vec(right) = itemp
            end if
        end do

        if (left == right) then
            marker = left + 1
        else
            marker = left
        end if
 
        call qSortvec(vec(:marker-1),marker-1)
        call qSortvec(vec(marker:),nvec-marker+1)

    end if

end subroutine qSortvec
subroutine itergen(rows,ni,tempDTD) ! caculate (D%*%t(D))^-1 in fused.1d by row

    implicit none

    integer, intent(in) :: rows
    integer, intent(in) :: ni
    real(kind = 8), dimension(ni), intent(in out) :: tempDTD
    integer :: i

    if(rows .eq. 1)then
        Do i = 1, ni
            tempDTD(i) = ni+1-i
        end do
        tempDTD = tempDTD/(ni+1)
    else if(rows .eq. ni)then
        Do i = 1, ni
            tempDTD(i) = i
        end do
        tempDTD = tempDTD/(ni+1)
    else
        Do i = 1, rows-1
            tempDTD(i) = i*(ni+1-rows)
        end do
        Do i = rows, ni
            tempDTD(i) = rows*(ni+1-i)
        end do
        tempDTD = tempDTD/(ni+1)
    end if

    return
end
subroutine dymul(vec,n,k,y,resultvec) ! caculate D%*%y in banded case

    integer, intent(in) :: k
    integer, intent(in) :: n
    real(kind = 8), dimension(k), intent(in) :: vec
    real(kind = 8), dimension(n), intent(in) :: y
    real(kind = 8), dimension(n-k+1), intent(in out) :: resultvec
    integer :: i

    Do i = 1,n-k+1
        resultvec(i) = dot_product(vec,y(i:k+i-1))
    end do

    return
end
subroutine sdymul(vec,n,k,y,sel,t,resultvec) ! caculate D%*%y for a select index in banded case

    integer, intent(in) :: k
    integer, intent(in) :: n
    integer, intent(in) :: t
    integer, dimension(t), intent(in) :: sel
    real(kind = 8), dimension(k), intent(in) :: vec
    real(kind = 8), dimension(n), intent(in) :: y
    real(kind = 8), dimension(n-k+1), intent(in out) :: resultvec
    integer :: i

    resultvec = 0
    Do i = 1,t
        resultvec(sel(i)) = dot_product(vec,y(sel(i):k+sel(i)-1))
    end do

    return
end
subroutine tdymul(vec,n,k,y,resultvec) ! caculate t(D)%*%y in banded case

    integer, intent(in) :: k
    integer, intent(in) :: n
    real(kind = 8), dimension(k), intent(in) :: vec
    real(kind = 8), dimension(n-k+1), intent(in) :: y
    real(kind = 8), dimension(n), intent(in out) :: resultvec
    integer :: i

    if(n-k+1>=k)then
        Do i = 1,k-1
            resultvec(i) = dot_product(vec(i:1:-1),y(1:i))
        end do
        i = 1
        Do i = k,n-k+1
            resultvec(i) = dot_product(vec(k:1:-1),y(i-k+1:i))
        end do
        i = 1
        Do i = n-k+2,n
            resultvec(i) = dot_product(vec(k:(i+k-n):-1),y(i-k+1:n-k+1))
        end do
    else
        Do i = 1,n-k+1
            resultvec(i) = dot_product(vec(i:1:-1),y(1:i))
        end do
        i = 1
        Do i = n-k+2,k
            resultvec(i) = dot_product(vec(i:i-n+k:-1),y(1:n-k+1))
        end do
        i = 1
        Do i = k+1,n
            resultvec(i) = dot_product(vec(k:(i-n+k):-1),y(i-k+1:n-k+1))
        end do
    end if

    return
end
subroutine dtdmul(vec,n,k,mat) ! caculate D%*%t(D) and return in a banded form in banded case

    real(kind = 8), dimension (:), allocatable :: doublevec

    integer, intent(in) :: k
    integer, intent(in) :: n
    real(kind = 8), dimension(k), intent(in) :: vec
    real(kind = 8), dimension(2*k-1,n-k+1), intent(in out) :: mat
    integer :: i

    allocate(doublevec(1:(2*k)))

    doublevec = 0
    doublevec((k+1):(2*k)) = vec

    if(n-k+1>=k)then
        Do i = 1,k
            mat(i,(k-i+1):(n-k+1)) = dot_product(vec,doublevec(i+1:i+k))
            mat(2*k - i,1:(n-2*k+i+1)) = mat(i,k)
        end do
    else
        Do i = (2*k-n),k
            mat(i,(k-i+1):(n-k+1)) = dot_product(vec,doublevec(i+1:i+k))
            mat(2*k - i,1:(n-2*k+i+1)) = mat(i,n-k+1)
        end do
    end if

    deallocate(doublevec)
    return
end
subroutine gtdymul(D,n,m,y,resultvec) ! caculate t(D)%*%y in general case

    integer, intent(in) :: m
    integer, intent(in) :: n
    real(kind = 8), dimension(n), intent(in out) :: resultvec
    real(kind = 8), dimension(m,n), intent(in) :: D
    real(kind = 8), dimension(m), intent(in) :: y
    integer :: i

    Do i = 1,n
       resultvec(i) = dot_product(D(:,i),y)
    end do

    return
end
subroutine DTDsub(DTD,sDTD,I,Ilength,Dn,veclen) ! caculate the submatrix of D%*%t(D) for fused.tfk

    implicit none

    real(kind = 8), dimension (:), allocatable :: tempvec

    integer, intent(in) :: veclen
    integer, intent(in) :: Ilength
    integer, intent(in) :: Dn
    integer, dimension(Ilength), intent(in) :: I
    real(kind = 8), dimension(2*veclen-1,Dn),intent(in) :: DTD
    real(kind = 8), dimension(2*veclen-1,Ilength), intent(in out) :: sDTD
    integer j,k,sm

    allocate(tempvec(1:(veclen-1)))

    sDTD = DTD(:,I(1:Ilength))
    Do j = 1,(Ilength-1)
        if(I(j+1)>I(j)+veclen-1)then
            sDTD((veclen+1):(2*veclen-1),j) = 0
        else
            sm = 0
            tempvec = sDTD((veclen+1):(2*veclen-1),j)
            sDTD((veclen+1):(2*veclen-1),j) = 0
            Do k = 1, min(veclen-1,Ilength-j)
                if(I(j+k)-I(j)<=veclen-1)then
                    sm = sm + 1
                    sDTD(veclen+sm,j) = tempvec(I(j+k)-I(j))
                end if
            end Do
        end if
    end Do
    sDTD((veclen+1):(2*veclen-1),Ilength) = 0
    sDTD(1:(veclen-1),:) = 0
    j = 1
    DO j = 1,(veclen-1)
        sDTD(j,(veclen-j+1):Ilength) = sDTD((2*veclen-j),1:(Ilength+j-veclen))
    end do

    deallocate(tempvec)
    return
end

!*****************************************************************************80
!
!! tfusedl02d return the fitting coefficients of given k with a small loss
!! function in tfusedl02d via primal dual active set algorithm.
!
!  Licensing:
!
!    This code is distributed under LGPL-3.
!
!  Author:
!
!    Original R version by Wen canhong, replica fortran90 version by Hong Zelin.
!
!  Modified:
!
!    29 August 2017
!
!  Parameters:
!
!    Input, integer (kind = 4) dim1, The number of rows for 2d X.
!    Input, integer (kind = 4) dim1, The length of columns for 2d X.
!    Input, real (kind = 8) y(dim1*dim2), The observe variable(geted by X).
!    Input, real (kind = 8) beta(dim1*dim2), The fitting coefficients of given 
!         T with a small loss function in fusedl0.
!    Input, real (kind = 8) z(dim1*(dim2-k+1)+dim2*(dim1-k+1)), the split 
!         variable of the argumented lagrange form.
!    Input, real (kind = 8) u(dim1*(dim2-k+1)+dim2*(dim1-k+1)), the lagrange 
!         operator of the argumented lagrange form for linear item.
!    Input, integer (kind = 4) T, the number of discover point by fusedl0.
!    Input, real (kind = 8) rho, the lagrange operator of the argumented 
!         lagrange form for split item.
!    Input, interger (kind = 4) iter, The iterations algorithm need to converge.
!    Input, integer (kind = 4) itermax, Maximum allowed number of iterations 
!         over the data for given parameters.
!    Input, integer (kind = 4) inv, Inverse of D%*%t(D),where D if the penalty 
!         Matrix.
!    Input, integer (kind = 4) vec, nozero elements in penalty Matrix D.
!    Input, integer (kind = 4) veck, Length of vec.
!
!
subroutine tfusedl02d(dim1, dim2, y, beta, z, u, T, rho, iter, itermax, inv, vec, veck)

  implicit none

  real(kind = 8), dimension (:), allocatable :: tempvec
  real(kind = 8), dimension (:,:), allocatable :: tempDTD,tempATA,tempITI
  integer, dimension (:), allocatable :: A,I,ordervec,Alast,Asort

  logical setequal
  integer dim1,dim2,iter,itermax,j,T,veck, rowlength
  real(kind = 8) y(dim1*dim2),beta(dim1*dim2),z(dim1*(dim2-veck+1)+dim2*(dim1-veck+1)), &
  u(dim1*(dim2-veck+1)+dim2*(dim1-veck+1))
  real(kind = 8) rho,inv(dim1*(dim2-veck+1)+dim2*(dim1-veck+1),dim1*(dim2-veck+1)+dim2*(dim1-veck+1)),vec(veck)

  external DGESV   ! solve general linear system (AX=B) function in LAPACK

  rowlength = dim1*(dim2-veck+1)+dim2*(dim1-veck+1)
  allocate(tempvec(1:rowlength))
  allocate(tempDTD(1:T,1:T))
  allocate(tempATA(1:T,1:T))
  allocate(tempITI(1:rowlength-T,1:rowlength-T))
  allocate(A(1:T))
  allocate(I(1:rowlength-T))
  allocate(ordervec(1:rowlength))
  allocate(Alast(1:T))
  allocate(Asort(1:T))

  tempDTD = 0
  tempATA = 0
  tempITI = 0
  Alast = 0
  tempvec = 0
  call dymul2d(vec,dim1,dim2,veck,y,tempvec)
  z = tempvec

  Do j = 1,rowlength
    ordervec(j) = j
  end do
  call Qsortorder(abs(z + u/rho),rowlength,ordervec)
  A = ordervec(1:T)
  I = ordervec((T+1):rowlength)

  iter = 1
  Do while(iter < itermax)

    if(T .eq. 1)then
        j = 1
        Do j = 1,rowlength-T
            tempITI(j,:) = inv(I(j),A(1)) * inv(A(1),I(1:(rowlength-T)))
        end Do
        tempITI = inv(I(1:(rowlength-T)),I(1:(rowlength-T)))&
        - tempITI/inv(A(1),A(1))
        u(I(1:(rowlength-T))) = matmul(tempITI,tempvec(I(1:(rowlength-T))))
    else
        tempDTD = 0
        j = 1
        DO j = 1,T
            tempDTD(j,j) = 1
        END DO
        tempATA = inv(A(1:T),A(1:T))
        call DGESV(T,T,tempATA,T,ordervec(1:T),tempDTD,T,j)
        tempITI = inv(I(1:(rowlength-T)),I(1:(rowlength-T))) - &
        matmul(matmul(inv(I(1:(rowlength-T)),A(1:T)),tempDTD),inv(A(1:T),I(1:(rowlength-T))))
        u(I(1:(rowlength-T))) = matmul(tempITI,tempvec(I(1:(rowlength-T))))
    end if
    u(A(1:T)) = 0
    call tdymul2d(vec,dim1,dim2,veck,u,beta)
    beta = y - beta
    call sdymul2d(vec,dim1,dim2,veck,beta,A,T,z)

    j = 1
    Do j = 1,rowlength
        ordervec(j) = j
    end do
    call Qsortorder(abs(z + u/rho),rowlength,ordervec)
    A = ordervec(1:T)
    Asort = A
    I = ordervec((T+1):rowlength)

    if(T>1)call qSortvec(Asort(1:T),T)
    setequal = .TRUE.
    j = 1
    Do j = 1,T
        if(Asort(j) .NE. Alast(j))then
            setequal = .FALSE.
            EXIT
        end if
    end do

    if(setequal)then
        iter = iter + 1
        EXIT
    else
        Alast = Asort
        iter = iter + 1
    end if

  end do

  iter = iter - 1
  deallocate(tempvec,tempDTD,tempATA,tempITI,A,I,ordervec,Alast,Asort)
  return
end
subroutine atfusedl02d(dim1, dim2, y, beta, z, u, tao, Kmax, L, eps, rho, miter, itermax, inv, vec, veck)

  implicit none

  real(kind = 8), dimension (:), allocatable :: tempvec,tempbeta
  real(kind = 8), dimension (:,:), allocatable :: tempDTD,tempATA,tempITI
  integer, dimension (:), allocatable :: A,I,ordervec,Alast,Asort

  logical setequal
  integer dim1,dim2,iter,itermax,tao,Kmax,L,j,T,k,veck,rowlength
  integer miter(Kmax)
  real(kind = 8) y(dim1*dim2),beta(dim1*dim2,Kmax),z(dim1*(dim2-veck+1)+dim2*(dim1-veck+1)),u(dim1*(dim2-veck+1)+dim2*(dim1-veck+1))
  real(kind = 8) rho,eps,inv(dim1*(dim2-veck+1)+dim2*(dim1-veck+1),dim1*(dim2-veck+1)+dim2*(dim1-veck+1)),vec(veck)

  external DGESV   ! solve general linear system (AX=B) function in LAPACK

  rowlength = dim1*(dim2-veck+1)+dim2*(dim1-veck+1)
  allocate(tempvec(1:rowlength))
  allocate(tempbeta(1:dim1*dim2))
  allocate(tempDTD(1:L,1:L))
  allocate(tempATA(1:L,1:L))
  allocate(tempITI(1:rowlength-tao,1:rowlength-tao))
  allocate(A(1:L))
  allocate(I(1:rowlength-tao))
  allocate(ordervec(1:rowlength))
  allocate(Alast(1:L))
  allocate(Asort(1:L))

  tempDTD = 0
  tempATA = 0
  tempITI = 0
  Alast = 0
  tempvec = 0
  tempbeta = 0
  call dymul2d(vec,dim1,dim2,veck,y,tempvec)
  z = tempvec

  k = 1
  Do while( k<Kmax )

    T = tao*k
    j = 1
    Do j = 1,rowlength
      ordervec(j) = j
    end do
    call Qsortorder(abs(z + u/rho),rowlength,ordervec)
    A(1:T) = ordervec(1:T)
    I(1:rowlength-T) = ordervec((T+1):rowlength)

    iter = 1
    Do while(iter < itermax)

      if(T .eq. 1)then
          j = 1
          Do j = 1,rowlength-T
            tempITI(j,1:rowlength-T) = inv(I(j),A(1)) * inv(A(1),I(1:(rowlength-T)))
          end Do
          tempITI(1:rowlength-T,1:rowlength-T) = inv(I(1:(rowlength-T)),I(1:(rowlength-T)))  &
          - tempITI(1:rowlength-T,1:rowlength-T)/inv(A(1),A(1))
          u(I(1:(rowlength-T))) = matmul(tempITI(1:rowlength-T,1:rowlength-T),tempvec(I(1:(rowlength-T))))
      else
          tempDTD = 0
          j = 1
          DO j = 1,T
            tempDTD(j,j) = 1
          END DO
          tempATA(1:T,1:T) = inv(A(1:T),A(1:T))
          call DGESV(T,T,tempATA(1:T,1:T),T,ordervec(1:T),tempDTD(1:T,1:T),T,j)
          tempITI(1:rowlength-T,1:rowlength-T) = inv(I(1:(rowlength-T)),I(1:(rowlength-T))) - &
          matmul(matmul(inv(I(1:(rowlength-T)),A(1:T)),tempDTD(1:T,1:T)),inv(A(1:T),I(1:(rowlength-T))))
          u(I(1:(rowlength-T))) = matmul(tempITI(1:rowlength-T,1:rowlength-T),tempvec(I(1:(rowlength-T))))
      end if
      u(A(1:T)) = 0
      call tdymul2d(vec,dim1,dim2,veck,u,tempbeta)
      tempbeta = y - tempbeta
      call sdymul2d(vec,dim1,dim2,veck,tempbeta,A(1:T),T,z)

      j = 1
      Do j = 1,rowlength
          ordervec(j) = j
      end do
      call Qsortorder(abs(z + u/rho),rowlength,ordervec)
      A(1:T) = ordervec(1:T)
      Asort(1:T) = A(1:T)
      I(1:rowlength-T) = ordervec((T+1):rowlength)

      if(T>1)call qSortvec(Asort(1:T),T)
      setequal = .TRUE.
      j = 1
      Do j = 1,T
          if(Asort(j) .NE. Alast(j))then
            setequal = .FALSE.
            EXIT
          end if
      end do

      if(setequal)then
        iter = iter + 1
        EXIT
      else
        Alast(1:T) = Asort(1:T)
        iter = iter + 1
      end if

    end do
      beta(:,k) = tempbeta
      miter(k) = iter - 1
      k = k + 1
      if( dot_product(y-tempbeta,y-tempbeta)<dim1*dim2*eps .or. T>=L )EXIT
    End do

  Kmax = k-1
  deallocate(tempvec,tempbeta,tempDTD,tempATA,tempITI,A,I,ordervec,Alast,Asort)
  return
end
subroutine dymul2d(vec,dim1,dim2,k,y,resultvec) ! caculate D%*%y in banded case for 2d

    real(kind = 8), dimension (:), allocatable :: tempvec
    integer, intent(in) :: k
    integer, intent(in) :: dim1
    integer, intent(in) :: dim2
    real(kind = 8), dimension(k), intent(in) :: vec
    real(kind = 8), dimension(dim1*dim2), intent(in) :: y
    real(kind = 8), dimension(dim1*(dim2-k+1)+dim2*(dim1-k+1)), intent(in out) :: resultvec
    integer :: i
    allocate(tempvec(dim1-k+1))
    tempvec = 0

    Do i = 1,dim2
        call dymul(vec,dim1,k,y((1+(i-1)*dim1):(i*dim1)),tempvec)
        resultvec((1+(i-1)*(dim1-k+1)):(i*(dim1-k+1))) = tempvec
    End Do
    i=1
    Do i = 1,dim1*(dim2-k+1)
        resultvec(dim2*(dim1-k+1)+i) = dot_product(vec,y(i:(i+dim1*(k-1)):dim1))
    End Do

    deallocate(tempvec)
    return
end
subroutine sdymul2d(vec,dim1,dim2,k,y,sel,t,resultvec) ! caculate D%*%y with a select index in banded case for 2d

    integer, intent(in) :: k
    integer, intent(in) :: dim1
    integer, intent(in) :: dim2
    integer, intent(in) :: t
    integer, dimension(t), intent(in) :: sel
    real(kind = 8), dimension(k), intent(in) :: vec
    real(kind = 8), dimension(dim1*dim2), intent(in) :: y
    real(kind = 8), dimension(dim1*(dim2-k+1)+dim2*(dim1-k+1)), intent(in out) :: resultvec
    integer :: i
    integer :: j
    integer :: mm

    resultvec = 0
    Do i = 1,t
        if(sel(i)<=dim2*(dim1-k+1))then
            j = int(sel(i)/(dim1-k+1))
            mm = sel(i) - j*(dim1-k+1)
            if(mm .eq. 0)then
                j = j-1
                mm = mm + (dim1-k+1)
            end if
            resultvec(sel(i)) = dot_product(vec,y((mm+j*dim1):(mm+j*dim1+k-1)))
        else
            resultvec(sel(i)) = dot_product(vec,y((sel(i)-dim2*(dim1-k+1)): &
            ((sel(i)-dim2*(dim1-k+1))+dim1*(k-1)):dim1))
        end if
    end do

    return
end
subroutine tdymul2d(vec,dim1,dim2,k,y,resultvec) ! caculate t(D)%*%y in banded case for 2d

    real(kind = 8), dimension (:), allocatable :: tempvec
    integer, intent(in) :: k
    integer, intent(in) :: dim1
    integer, intent(in) :: dim2
    real(kind = 8), dimension(k), intent(in) :: vec
    real(kind = 8), dimension(dim1*(dim2-k+1)+dim2*(dim1-k+1)), intent(in) :: y
    real(kind = 8), dimension(dim1*dim2), intent(in out) :: resultvec
    integer :: i
    integer :: j
    allocate(tempvec(1:dim1))
    tempvec=0

    if(dim2-k+1>k)then
        Do i = 1,dim2-k+1
            call tdymul(vec,dim1,k,y((1+(i-1)*(dim1-k+1)):(i*(dim1-k+1))),tempvec)
            j = 1
            Do j = 1,dim1
                if(i<=k)then
                    tempvec(j) = tempvec(j) + dot_product(vec(i:1:-1),y((j+dim2*(dim1-k+1)):&
                    (j+dim2*(dim1-k+1)+dim1*(i-1)):dim1))
                else
                    tempvec(j) = tempvec(j) + dot_product(vec(k:1:-1),y((j+dim1*(i-k)+dim2*(dim1-k+1)):&
                    (j+dim1*(i-k)+dim2*(dim1-k+1)+dim1*(k-1)):dim1))
                end if
            end Do
            resultvec((1+dim1*(i-1)):(dim1*i)) = tempvec
        end do
        i = 1
        Do i = dim2-k+2,dim2
            call tdymul(vec,dim1,k,y((1+(i-1)*(dim1-k+1)):(i*(dim1-k+1))),tempvec)
            j = 1
            Do j = 1,dim1
                tempvec(j) = tempvec(j) + dot_product(vec(k:(i+k-dim2):-1),y((j+dim1*(i-k)+dim2*(dim1-k+1)):&
                (j+dim1*(i-k)+dim2*(dim1-k+1)+dim1*(dim2-i)):dim1))
            end Do
            resultvec((1+dim1*(i-1)):(dim1*i)) = tempvec
        end do
    else
        Do i = 1,dim2-k+1
            call tdymul(vec,dim1,k,y((1+(i-1)*(dim1-k+1)):(i*(dim1-k+1))),tempvec)
            j = 1
            Do j = 1,dim1
                tempvec(j) = tempvec(j) + dot_product(vec(i:1:-1),y((j+dim2*(dim1-k+1)):&
                (j+dim2*(dim1-k+1)+dim1*(i-1)):dim1))
            end Do
            resultvec((1+dim1*(i-1)):(dim1*i)) = tempvec
        end Do
        i = 1
        Do i = dim2-k+2,k
            call tdymul(vec,dim1,k,y((1+(i-1)*(dim1-k+1)):(i*(dim1-k+1))),tempvec)
            j = 1
            Do j = 1,dim1
                tempvec(j) = tempvec(j) + dot_product(vec(i:(i-dim2+k):-1),y((j+dim2*(dim1-k+1)):&
                (j+dim2*(dim1-k+1)+dim1*(dim2-k)):dim1))
            end Do
            resultvec((1+dim1*(i-1)):(dim1*i)) = tempvec
        end do
        if(dim2 > k)then
            i = 1
            Do i = k+1,dim2
                call tdymul(vec,dim1,k,y((1+(i-1)*(dim1-k+1)):(i*(dim1-k+1))),tempvec)
                j = 1
                Do j = 1,dim1
                    tempvec(j) = tempvec(j) + dot_product(vec(k:(i-dim2+k):-1),y((j+(i-k)*dim1+dim2*(dim1-k+1)):&
                    (j+(i-k)*dim1+dim2*(dim1-k+1)+dim1*(dim2-i)):dim1))
                end Do
                resultvec((1+dim1*(i-1)):(dim1*i)) = tempvec
            end do
        end if
    end if

    deallocate(tempvec)
    return
end