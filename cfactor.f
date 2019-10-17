        implicit double precision (a-h,o-z)
        double precision, allocatable :: w(:)
        double complex sum,ima
        double complex, allocatable :: a(:,:),r(:)
        double complex, allocatable :: x(:,:),y(:,:),z(:,:),z2(:,:)
        integer, allocatable :: icols(:)
c
        data ima / (0.0d0,1.0d0) /
c
c       Form a matrix A(n,m) of lowish rank.
c
        n = 210
        m = 200
c
        dd1 = n
        dd1 = 1.0d0/dd1
c
        dd2 = m
        dd2 = 1.0d0/dd2       
c
        allocate(a(n,m),icols(m),r(n*m))
        allocate(w(2*(3*n+1)*m+1))
c       
        do 1000 i=1,n
        do 1100 j=1,m
        xx = -0.001d0-dd1*(i-1)
        yy =  0.001d0+dd2*(j-1)
        dd = sqrt(abs(xx-yy))
        a(i,j) = log(dd)*cos(13*dd)
 1100 continue
 1000 continue
c
        eps = 1.0d-15
        call factor_right(eps,n,m,a,krank,icols,r,w)      
        call prinf("after factor_right, krank = *",krank,1)
c
c       Form a matrix X of random entries.
c
        k  = 1010
        allocate(x(m,k),y(krank,k),z(n,k),z2(n,k))
        do 2000 i=1,m
        do 2100 j=1,k
        ione=1
        call corrand3(ione,dd1)
        call corrand3(ione,dd2)
        x(i,j)=dd1+ima*dd2
 2100 continue
 2000 continue
c
c       Apply the matrix C from the left to X to form Y.
c
        call apply_left(m,k,krank,icols,r,x,y)
c
c       Apply the subset of the columns of A to Y to form z.
c
        do 2200 i=1,n
        do 2300 j=1,k
        sum=0
        do 2400 l=1,krank
        sum = sum + a(i,icols(l))*y(l,j)
 2400 continue
        z(i,j)=sum
 2300 continue
 2200 continue
c
c       Now apply A to X directly to form z2.
c
        do 2500 i=1,n
        do 2600 j=1,k
        sum=0
        do 2700 l=1,m
        sum = sum + a(i,l)*x(l,j)
 2700 continue
        z2(i,j)=sum
 2600 continue
 2500 continue
c
c       Check the error.
c     
        errl2=0
        do 2800 i=1,n
        do 2900 j=1,k
        dd = abs(z(i,j)-z2(i,j))
        errl2 = errl2+dd**2
 2900 continue
 2800 continue
        errl2=sqrt(errl2)
        call prin2("errl2=*",errl2,1)
c
c       Repeat the experiment, but for apply_rightt.
c
        k  = 103
        deallocate(x,y,z,z2)
        allocate(x(k,m),y(k,krank),z(k,n),z2(k,n))
c
        do 3000 i=1,k
        do 3100 j=1,m
        ione=1
        call corrand3(ione,dd1)
        call corrand3(ione,dd2)
        x(i,j)=dd1+ima*dd2
 3100 continue
 3000 continue
c
c       Form X*C^t in Y.
c
        call apply_rightt(m,k,krank,icols,r,x,y)
c
c       Multiply Y on right by the tranpose of the subset of A.
c
        do 3200 i=1,k
        do 3300 j=1,n
        sum=0
        do 3400 l=1,krank
        sum = sum + y(i,l)*a(j,icols(l))
 3400 continue
        z(i,j)=sum
 3300 continue
 3200 continue
c
c       Now compute X A^t directly.
c
        do 3500 i=1,k
        do 3600 j=1,n
        sum=0
        do 3700 l=1,m
        sum = sum + x(i,l)*a(j,l)
 3700 continue
        z2(i,j)=sum
 3600 continue
 3500 continue
c
c       Check the error.
c     
        errl2=0
        do 3800 i=1,k
        do 3900 j=1,n
        dd = abs(z(i,j)-z2(i,j))
        errl2 = errl2+dd**2
 3900 continue
 3800 continue
        errl2=sqrt(errl2)
        call prin2("errl2=*",errl2,1)

        end



        subroutine apply(n,m,a,x,y)
        implicit double precision (a-h,o-z)
        double complex a(n,m),x(1),y(1),sum
        do 1000 i=1,n
        sum=0
        do 1100 j=1,m
        sum=sum+a(i,j)*x(j)
 1100 continue
        y(i)=sum
 1000 continue
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and begining of the 
c       "skeletonization" code proper.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains code for computing interpolative 
c       decompositions of arbitrary complex-valued matrices.  The 
c       following routines are user-callable:
c
c   factor_right - factor an input matrix A on the right; that is, as
c
c           A = B * C                                                    (1)
c
c       where B is a subset of the columns of A.  Note that the 
c       matrix C has internal structure which makes it more efficient 
c       to store and apply than a general matrix.
c       
c   apply_left - given the data returned by factor_right apply the 
c       matrix C in (1) from the left; that is, compute the matrix
c       product
c
c            X = C * Y
c
c       for a user-specified matrix Y.
c
c   apply_rightt - given the data returned by factor_right, apply the
c       transpose of the matrix C in (1) from the right; that is,
c       compute
c
c            X = Y * C^t
c
c       IMPORTANT NOTE: this subroutine applies the transpose and NOT 
c       the conjugate of C
c       
c   form_right - given the data returned by factor_right, explicitly
c       form the matrix C in (1)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



        subroutine factor_right(eps,n,m,a,krank,icols,r,w)
        implicit double precision (a-h,o-z)
        double complex a(n,m),r(1)
        dimension icols(1),w(1)
c
c       Factor the input matrix a as
c
c         A(n,m) = B(n,krank) * C(krank,m)
c
c       where krank is a proxy for the numerical rank of A to precision
c       eps and B is a subset of the columns of A.
c
c       More specifically, B is a matrix consisting of the subset of the 
c       columns of A defined by
c
c          B = A(:,icols(1:krank))                                      (1)
c
c       where icols is a permutation vector returned by the routine
c       and R is a matrix of the form
c
c         [I(krank,krank) R(krank,m-krank)] \Pi(m,m)^-1
c
c       with \Pi the permutation matrix defined by icols and I the
c       identity matrix.
c
c       The requested precision eps for the factorization determines
c       the value of krank, which is a proxy for the numerical rank of
c       the matrix a to precision eps.
c
c                           Input Parameters:
c
c   eps - precision for the computation
c   (n,m) - dimensions of the input matrix A 
c   a - the (n,m) input matrix, which is NOT destroyed by this routine
c
c   w - a work array whose length (in double works) is at least
c
c               2*(n+1)*m+1+2*krank*m.  
c
c       To be safe in cases when krank is not known, the length should 
c       be at least
c
c               2*(n+1)*m + 1 + 2*min(n,m)*m.
c
c       More simply, 
c
c               2*(3*n+1)*m + 1 
c
c       will always suffice.
c           
c                         Output Parameters:
c
c   krank - approximately the numerical rank of the input matrix a to 
c       precision eps
c   icols - a permutation vector which specifies the permutation \Pi
c   r - the (krank,m-krank) matrix r
c
        krank = 0
        if (n .eq. 0 .OR. m .eq. 0) return
c
c       Make a copy of the matrix in the work array.
c
        iq = 1
        lq = n*m*2
c
        irnorms = iq+lq
        lrnorms = m+1
c
        call move(2*n*m,a,w(iq))
c
c       Gram Schmidt the input matrix.  Let A_1 denote the matrix
c       consisting of the pivot columns of A in the order in which
c       they were chosen and A_2 the remaining portion of the matrix.
c
        call cgspiv(w(iq),n,m,eps,w(irnorms),icols,krank)
        call quicksorti(m-krank,icols(krank+1))
c
        call prinf("in factor_right, krank = *",krank,1)
        call prin2("in factor_right, rnorms = *",w(irnorms),krank)
c        call prinf("in factor_right, icols = *",icols,m)
c
c       Allocate space for the matrices Q^* A_1 and Q^* A_2 from the work 
c       array.
c
        ir1 = irnorms+lrnorms
        lr1 = 2*krank*krank
c
        ir2 = ir1+lr1
        lr2 = 2*krank*(m-krank)
c
c       Call the auxillary routine to form R = R_1 ^(-1) Q^* R_2
c
        call factor_right0(n,m,krank,icols,a,w(iq),w(ir1),w(ir2),r)
c
        end


        subroutine factor_right0(n,m,krank,icols,a,q,r1,r2,r)
        implicit double precision (a-h,o-z)
        double complex q(n,krank),r1(krank,krank),r2(krank,m-krank)
        double complex r(krank,m-krank),a(n,m),sum
        dimension icols(1)
c
c       Form the product R_1 = Q^* A_1
c
        do 1000 i=1,krank
        do 1100 j=1,krank
        sum = 0
        do 1200 l=1,n
        sum = sum + conjg(q(l,i))*a(l,icols(j))
 1200 continue
        r1(i,j)=sum
 1100 continue
 1000 continue
c
c       Form the product R_2 = Q^* A_2
c
        do 2000 i=1,krank
        do 2100 j=1,m-krank
        sum = 0
        do 2200 l=1,n
        sum = sum + conjg(q(l,i))*a(l,icols(krank+j))
 2200 continue
        r2(i,j)=sum
 2100 continue
 2000 continue
c
c       Apply the inverse of R_1 to R_2 to form R.
c
        do 3000 j=1,m-krank
        r(krank,j) = r2(krank,j) / r1(krank,krank)
        do 3100 i=krank-1,1,-1
        r(i,j) = r2(i,j)
        do 3200 l=krank,i+1,-1
        r(i,j) = r(i,j)-r1(i,l)*r(l,j)
 3200 continue      
        r(i,j) = r(i,j) / r1(i,i)
 3100 continue
 3000 continue
c      
        end


        subroutine apply_left(m,n,krank,icols,r,x,y)
        implicit double precision (a-h,o-z)
        double complex r(krank,m-krank),x(m,n),y(krank,n),sum
        dimension icols(1)
c
c       Apply the (krank,m) matrix C appearing in (1) above to a user-
c       specified matrix (m,n) matrix X.
c
c                         Input Parameters:
c
c   m - the number of columns in the orignal input matrix
c   n - the number of columns of the matrix X to which C is to be
c       applied
c   krank - the numerical rank of the matrix
c   icols - the permutation vector returned by factor_right
c   r - the (krank,m-krank) matrix r returned by factor_right
c
c   x - the (m,n) matrix to which
c
c
c                        Output Parameters:
c
c   y - the result of applying the matrix C to X
c
        do 1000 l=1,n
        do 1100 i=1,krank
        sum=0
        do 1200 j=1,m-krank
        sum=sum+r(i,j)*x(icols(j+krank),l)
 1200 continue
        y(i,l)=sum+x(icols(i),l)
 1100 continue
 1000 continue
        end



        subroutine apply_rightt(m,n,krank,icols,r,x,y)
        implicit double precision (a-h,o-z)
        double complex r(krank,m-krank),x(n,m),y(n,krank),sum
        dimension icols(1)
c
c       Apply a user-supplied matrix X(n,m) to the *transpose* (not
c       conjugate) of the (krank,m) matrix C appearing in (1).
c
c                       Input Parameters:
c
c   m - the number of columns in the orignal input matrix
c   n - the number of rows of the matrix X
c   krank - the numerical rank of the matrix
c   icols - the permutation vector returned by factor_right
c   r - the (krank,m-krank) matrix r returned by factor_right
c
c   x - the matrix which is to be applied to the transpose of C
c
c
c                      Output Parameters:
c
c   y - the result of applying the matrix to C
c
        do 1000 l=1,n
        do 1100 i=1,krank
        sum=0
        do 1200 j=1,m-krank
        sum=sum+r(i,j)*x(l,icols(j+krank))
 1200 continue
        y(l,i)=sum+x(l,icols(i))
 1100 continue
 1000 continue

        end



        subroutine form_right(m,krank,icols,r,c)
        implicit double precision (a-h,o-z)
        double complex r(krank,m-krank),c(krank,m)
        dimension icols(1)
c
c       Explicitly form the matrix c in the factorization (1).
c      
c                       Input Parameters:
c
c   m - the number of columns in the original input matrix
c   krank - rank of the input matrix
c   icols - the permutation vector returned by factor_right
c   r2 - the matrix r2 returned by factor_right
c
c                      Output Parameters:
c  
c   r - the matrix r appearing in the factorization (1)
c
        do 1000 i=1,krank
        do 1100 j=1,m
        c(i,icols(j))=0
 1100 continue
 1000 continue
c
        do 1200 i=1,krank
        c(i,icols(i)) = 1
 1200 continue
c
        do 1300 i=1,krank
        do 1400 j=1,m-krank
        c(i,icols(j+krank)) = r(i,j)
 1400 continue
 1300 continue
c
        end



c
c
        subroutine cgspiv(b,n,m,eps,rnorms,ipivots,ncols)
        implicit double precision (a-h,o-z)
        save
        dimension rnorms(1),ipivots(1)
        double complex  b(n,m),cd
c 
c        . . . initialize the array of pivots
c 
        do 1100 i=1,m
        ipivots(i)=i
 1100 continue
c 
c       . . . prepare the array of values of norms
c             of columns
c 
        done=1
        dtot=0
        do 1400 i=1,m
c 
        d=0
        do 1200 j=1,n
        d=d+b(j,i)*conjg(b(j,i))
 1200 continue
        rnorms(i)=sqrt(d)
        dtot=dtot+d
 1400 continue
c
        dtot=sqrt(dtot)
        if(dtot .lt. eps) then
        ncols=0
        return
        endif
        
c 
        thresh=dtot*eps
c 
c       . . . conduct gram-schmidt iterations
c 
cccc        do 4000 i=1,n
        do 4000 i=1,m
c 
c       find the pivot
c 
        ipivot=i
        rn=rnorms(i)
c 
        do 2200 j=i+1,m
        if(rnorms(j) .le. rn) goto 2200
        rn=rnorms(j)
        ipivot=j
 2200 continue
 2400 continue
c
c       put the column number ipivot in the i-th place
c 
        do 2600 j=1,n
        cd=b(j,i)
        b(j,i)=b(j,ipivot)
        b(j,ipivot)=cd
 2600 continue
c 
        iijj=ipivots(i)
        ipivots(i)=ipivots(ipivot)
        ipivots(ipivot)=iijj
c 
        d=rnorms(ipivot)
        rnorms(ipivot)=rnorms(i)
        rnorms(i)=d
c 
c       orthogonalize the i-th column to all preceding ones
c 
        if(i .eq. 1) goto 2790
        do 2780 j=1,i-1
c 
        call gspiv_cleascap(b(1,i),b(1,j),n,cd)
c 
        do 2770 l=1,n
        b(l,i)=b(l,i)-b(l,j)*cd
 2770 continue
 2780 continue
 2790 continue
c 
c       normalize the i-th column
c 
        call gspiv_cleascap(b(1,i),b(1,i),n,cd)
c 
        d=cd
        if(d .lt. thresh**2 ) return
c 
        ncols=i
c 
        d=done/sqrt(d)
        do 2800 j=1,n
        b(j,i)=b(j,i)*d
 2800 continue
c 
        if(i .eq. m) goto 3400
c 
c        orthogonalize everything else to it
c 
        do 3200 j=i+1,m
c 
        if(rnorms(j) .lt. thresh) goto 3200
c 
        call gspiv_cleascap(b(1,i),b(1,j),n,cd)
c 
        cd=conjg(cd)
c 
        rrn=0
        do 3000 l=1,n
        b(l,j)=b(l,j)-b(l,i)*cd
        rrn=rrn+b(l,j)*conjg(b(l,j))
 3000 continue
        rnorms(j)=sqrt(rrn)
 3200 continue
 3400 continue
c 
 4000 continue
c 
        return
        end
c
c
c
        subroutine gspiv_cleascap(x,y,n,prod)
        implicit double complex  (a-h,o-z)
        save
        dimension x(1),y(1)
c 
        prod=0
        do 1200 i=1,n
        prod=prod+x(i)*dconjg(y(i))
 1200 continue
        return
        end
c 
