        implicit real *8 (a-h,o-z)
ccc        dimension filt(100),t(1100)
        complex *16 ima,a(1100),c(1000),b(1000),
     1       work(1000)
        data ima/(0.0d0,1.0d0)/
c 
        call prini(6,13)
c 
c       SET ALL PARAMETERS
C 
        PRINT *, 'ENTER n '
        READ *,n
        CALL PRINF('n=*',n,1 )
C 
        call invtes(a,b,n,c,work)
c 
  
  
        stop
        end
c 
c 
c 
c 
c 
        subroutine matpr(a,n,mes)
        save
        complex *16 a(n,n)
        character *1 mes(1)
c 
c        print the message
c 
        call prin(mes,a,0)
c 
c       print the matrix
c 
        do 2000 i=1,n
 1200 format(20x)
 1400 format(2x,' row number',i5)
        write(6,1200)
        write(6,1400) i
c 
        write(13,1200)
        write(13,1400) i
 1600 format(3(1x,e10.4,1x,e10.4,1x))
        write(6,1600)(a(j,i),j=1,n)
        write(13,1600)(a(j,i),j=1,n)
 2000 continue
        return
        end
c 
c 
c 
c 
c 
         subroutine invtes(a,b,n,c,work)
         implicit real *8 (a-h,o-z)
        save
         complex *16 a(n,n),b(n,n),ima,cd,
     1        c(n,n),work(1),x(100),y(100),z(100)
        data ima/(0.0d0,1.0d0)/
c 
c       construct the test matrix
c 
         eps=1.0d-3
         eps=10
        do 1400 i=1,n
        do 1200 j=1,n
        a(i,j)=i**2-j**2*ima
        a(i,j)=a(i,j)*eps
 1200 continue
        a(i,i)=i+a(i,i)
 1400 continue
c 
        do 1800 i=1,n
        do 1600 j=1,n
        c(j,i)=a(j,i)
 1600 continue
 1800 continue
c 
         call prin2('before cortho1, a=*',a,n*n*2)
c 
c        invert the matrix a
c 
        call corthom(a,n,work,cond)
         call prin2('after corthom, cond=*',cond,1)
c 
cccc        call prin2('after cortho1, a=*',a,n*n*2)
c 
c        multiply a and its putative inverse
c 
        call matmul(c,a,b,n)
        call prin2('after matmul, b=*',b,n*n*2)
        call matpr(b,n, 'and be printed politely*')
c 
c        verify the whole nonsense by applying it to the vector
c 
c      . . . create the test vector
c 
        do 2200 i=1,n
        x(i)=i+i**2*ima
 2200 continue
         call prin2('x as created *',x,n*2)
c 
c       apply a to x
c 
ccc        call matmul(a,n,x,y)
        call corthom_matvec(c,n,x,y)
          call prin2('a x =*',y,n*2)
c 
c        apply inverse of a to x
c 
        call corthom_matvec(a,n,y,z)
         call prin2('z=*',z,n*2)
  
  
        return
        end
c 
c 
c 
c 
c 
        subroutine matmul(a,b,c,n)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n),b(n,n),c(n,n),cd
c 
c        multiply matrices a,b, getting c
c 
        do 2000 i=1,n
        do 1800 j=1,n
        cd=0
        do 1200 k=1,n
        cd=cd+a(i,k)*b(k,j)
 1200 continue
        c(i,j)=cd
 1800 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c          this is the end of the debugging code and the beginning
c          of the actual matrix inversion routine
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
c 
c 
        subroutine corthom(a,n,work,cond)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(1),work(1)
c 
c        this subroutine inverts the user-supplied matrix
c        a by means of the gram-schmidt process. in fact, this
c        is a memory management routine for the routine cortho1,
c        which performs the actual inversion
c 
c               input parameters:
c 
c  a - the n*n - matrix to be inverted (destroyed by the
c         subroutine)
c  n - the dimensionality of the matrix a
c 
c               output parameters:
c 
c  a - the inverse of a on input
c  cond - the estimate of the condition number of the matrix a
c         that has been inverted
c 
c               work array
c 
c  work - must be at least 2 * n * (n+1)+2 real *8 locations long
c 
c 
c        . . . allocate memory for the matrix inverter
c 
        ib=1
        lb=n*n
c 
        iw=ib+lb+1
c 
c       invert the matrix
c 
        call cortho1(a,work(ib),n,work(iw),cond)
        return
        end
c 
c 
c 
c 
c 
        subroutine cortho1(a,b,n,work,cond)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n),b(n,n),cd,zero,
     1      work(1)
c 
c       set the matrix b to unity
c 
        done=1
        zero=0
        do 1400 i=1,n
        do 1200 j=1,n
        b(j,i)=zero
 1200 continue
        b(i,i)=1
 1400 continue
c 
c        conduct the gram-schmidt process
c 
        dnormax=-1
        dnormin=1.0d23
c 
        do 4000 i=1,n
c 
c       orthogonalize the i-th column to all preceeding ones
c 
        if(i .eq. 1) goto 2190
        do 2180 j=1,i-1
c 
        call corthom_scap(a,n,i,j,cd)
c 
        do 2170 k=1,n
        a(i,k)=a(i,k)-a(j,k)*cd
        b(i,k)=b(i,k)-b(j,k)*cd
 2170 continue
 2180 continue
 2190 continue
c 
c       normalize the i-th row
c 
        call corthom_scap(a,n,i,i,cd)
c 
        d=cd
        d=done/dsqrt(d)
        if(dnormax .lt. d) dnormax=d
        if(dnormin .gt. d) dnormin=d
        do 2400 j=1,n
        a(i,j)=a(i,j)*d
        b(i,j)=b(i,j)*d
 2400 continue
c 
c       orthogonalize all subsequent rows to the i-th one
c 
        if(i .eq. n) goto 4000
c 
        do 3000 j=i+1,n
c 
        cd=0
        do 2600 k=1,n
        cd=cd+a(i,k)*dconjg(a(j,k))
 2600 continue
        cd=dconjg(cd)
c 
        do 2800 k=1,n
        a(j,k)=a(j,k)-cd*a(i,k)
        b(j,k)=b(j,k)-cd*b(i,k)
 2800 continue
 3000 continue
c 
 4000 continue
        cond=dnormax/dnormin
c 
c       now, multiply the abjoint of the resulting
c       orthogonal matrix by the triangular one,
c       obtaining the inverse of the original a
c 
        do 5000 i=1,n
        do 4800 j=1,n
        cd=0
        do 4200 k=1,n
        cd=cd+dconjg(a(k,i))*b(k,j)
 4200 continue
        work(j)=cd
 4800 continue
        do 4900 j=1,n
        a(j,i)=work(j)
 4900 continue
 5000 continue
c 
c        now, transpose a
c 
        do 5400 i=2,n
        do 5200 j=1,i
        cd=a(i,j)
        a(i,j)=a(j,i)
        a(j,i)=cd
 5200 continue
 5400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine corthom_matvec(a,n,x,y)
        implicit real *8 (a-h,o-z)
        save
        complex *16 a(n,n),x(1),y(1),cd
c 
c        apply the matrix a to the vector x obtaining y
c 
        do 1400 i=1,n
        cd=0
        do 1200 j=1,n
        cd=cd+a(i,j)*x(j)
 1200 continue
        y(i)=cd
 1400 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine corthom_scap(a,n,i,j,prod)
        implicit complex *16 (a-h,o-z)
        save
        dimension a(n,n)
c 
        prod=0
        do 1200 k=1,n
        prod=prod+a(i,k)*conjg(a(j,k))
 1200 continue
        return
        end
  
  
