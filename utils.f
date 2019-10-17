cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains a number of basic utility routines which are used 
c       throughout the solver.  The following subroutines are user-callable:
c
c   quicksort - an implementation of the quicksort algorithm which
c       sorts a list of reals and keeps track of their indices 
c
c   insort - an implementation of the insertion sort algorithm which
c       sorts a list of reals and keeps track of their indices
c
c   reverse - reverse the order of a list of double precision reals 
c       together with their indices
c
c   quicksort0 - a version of quicksort which does not keep track of
c       indices
c
c   insort0 - a version of insertion sort which does not keep track of
c       indices
c
c   reverse0 - a version of reverse which does not keep track of the
c       indices
c
c   isort - sort a list of integers using the quicksort algorithm
c
c   insorti - sort a list of integers using the insertion sort algorithm
c
c   insorti2 - sort a list of integer using the insertion sort algorithm;
c        keep track of their indices
c 
c   quicksorti - sort a list of integers using the quicksort algorithm
c
c   quicksorti2 - sort a list of integers using the quicksort algorithm;
c       keep track of their indices
c
c   iduplicates - sort a list of integers and remove any duplicates
c       which occur from the list
c
c   iremove - remove from a sorted list of integers all integers 
c       appearing in a second sorted list.  The resulting list will
c       be sorted.
c
c   imerge - merge two sorted lists of integers
c
c   move - move a block of double precision reals
c
c   imove - move a block of integers
c  
c   izero - zero a block of integers
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      
        subroutine iremove(n,ia,m,ib)
        implicit double precision (a-h,o-z)
        dimension ia(n),ib(m)
c
c       Remove from the list ia of length n all integers appearing in 
c       the list ib of length m.  Both the list ia and the list ib
c       must be sorted before this call is made.  The results will
c       also be sorted.
c
        isrc = 1
        itar = 1
        ii   = 1
 1000 continue
c
        if (ii .gt. m)   goto 2000
        if (isrc .gt. n) goto 3000
c
        if (ia(isrc) .gt. ib(ii)) then
        ii=ii+1
        goto 1000
        endif
c
        if (ia(isrc) .lt. ib(ii)) then         
        ia(itar) = ia(isrc)
        itar=itar+1
        isrc=isrc+1
        goto 1000
        endif
c
        isrc=isrc+1
        goto 1000
c
 2000 continue
        if (isrc .gt. n) goto 3000
        ia(itar) = ia(isrc)
        itar=itar+1
        isrc=isrc+1
        goto 2000
c
 3000 continue
        n = itar-1
        return        
        
        end


        subroutine quicksort(n,vals,idxs)
        implicit double precision (a-h,o-z)
        dimension istack(2,20 000)
        dimension vals(1),idxs(1)
c
c       Sorts the list of double precision numbers in vals and keep track of
c       indices.
c
        if (n .lt. 100) then
        call insort(n,vals,idxs)
        return
        endif
c
        maxstack = 10 000
        k        = 60
c
        m = 1
        istack(1,1) = 1
        istack(2,1) = n
c
 1000 continue
        if (m .eq. 0) goto 1100
        i1 = istack(1,m)
        i2 = istack(2,m)
        m=m-1
c
        l = i2-i1+1
        if (l .le. k) then
        call insort(l,vals(i1),idxs(i1))
        goto 1000
        endif
c
c       Otherwise perform quicksort.
c
        call quicksort1(vals,idxs,i1,i2,i3)
c
c       This should never happen, but just in case ...
c
        if (m+2 .ge. maxstack) then
        print *,"quicksort out of memory"
        stop
        endif
c
c       Make sure the smaller half is processed first to ensure storage
c       requirements are O(log(n))
c             
        n1 = i3-i1+1
        n2 = i2-i3
c
        if (n2 .lt. n1) then
c
        m = m+1
        istack(1,m) = i1
        istack(2,m) = i3
c
        m = m+1
        istack(1,m) = i3+1
        istack(2,m) = i2
c
        else
c
        m = m+1
        istack(1,m) = i3+1
        istack(2,m) = i2
c
        m = m+1
        istack(1,m) = i1
        istack(2,m) = i3
c
        endif
c
        goto 1000
 1100 continue
        end


        subroutine quicksort1(vals,idxs,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension vals(1),idxs(1)
c
c       Randomly choose a pivot index.
c
        call corrand3(1,r)
        ipiv = i1+floor((i2-i1)*r)
c
        ipiv = i1+(i2-i1)/2
c
        val  = vals(ipiv)
        ival = idxs(ipiv)
c
c       Swap the pivot element and the last element.
c
        vals(ipiv) = vals(i2)
        vals(i2)   = val
c
        idxs(ipiv) = idxs(i2)
        idxs(i2)   = ival
c
        i3 = i1
c
        do 1000 i=i1,i2-1
c
        if( vals(i) .lt. val) then
        d  = vals(i)
        id = idxs(i)
c
        vals(i)  = vals(i3)
        vals(i3) = d       
c
        idxs(i)  = idxs(i3)
        idxs(i3) = id
        i3=i3+1
        endif
 1000 continue
c
        dd = vals(i3)
        vals(i3) = vals(i2)
        vals(i2) = dd
c
        idd = idxs(i3)
        idxs(i3) = idxs(i2)
        idxs(i2) = idd
c
        end


        subroutine insort(k,a,ib)
        implicit double precision (a-h,o-z)
        dimension a(1),ib(1)
c
        if (k .le. 1) return
c
        do 1000 i=2,k
        val=a(i)
        ival=ib(i)
        j=i-1
        do 1100 while (j .ge. 1 .AND. a(j) .gt. val) 
        a(j+1)=a(j)
        ib(j+1)=ib(j)
        j=j-1
 1100 continue
        a(j+1)=val
        ib(j+1)=ival
 1000 continue
        end


        subroutine reverse(n,vals,idxs)
        implicit double precision (a-h,o-z)
        dimension vals(1),idxs(1)
c
c       Reverse the order of a list of reals.
c
        do 1000 j=1,n/2
        d = vals(j)
        vals(j)     = vals(n-j+1)
        vals(n-j+1) = d
        i = idxs(j)
        idxs(j)     = idxs(n-j+1)
        idxs(n-j+1) = i
 1000 continue
c
        end


        subroutine insort0(k,a)
        implicit double precision (a-h,o-z)
        dimension a(1)
c
        if (k .le. 1) return
c
        do 1000 i=2,k
        val=a(i)
        j=i-1
        do 1100 while (j .ge. 1 .AND. a(j) .gt. val) 
        a(j+1)=a(j)
        j=j-1
 1100 continue
        a(j+1)=val
 1000 continue
        end


        subroutine quicksort0(n,vals)
        implicit double precision (a-h,o-z)
        dimension istack(2,20 000)
        dimension vals(1),idxs(1)
c
c       Sort a list of double precision numbers.
c
        if (n .lt. 100) then
        call insort0(n,vals)
        return
        endif
c
        maxstack = 10 000
        k        = 60
c
        m = 1
        istack(1,1) = 1
        istack(2,1) = n
c
 1000 continue
        if (m .eq. 0) goto 1100
        i1 = istack(1,m)
        i2 = istack(2,m)
        m=m-1
c
        l = i2-i1+1
        if (l .le. k) then
        call insort0(l,vals(i1))
        goto 1000
        endif
c
c       Otherwise perform quicksort.
c
        call quicksort01(vals,i1,i2,i3)
c
c       This should never happen, but just in case ...
c
        if (m+2 .ge. maxstack) then
        print *,"quicksort out of memory"
        stop
        endif
c
c       Make sure the smaller half is processed first to reduce storage
c       to O(logn).
c             
        n1 = i3-i1+1
        n2 = i2-i3
c
        if (n2 .lt. n1) then
c
        m = m+1
        istack(1,m) = i1
        istack(2,m) = i3
c
        m = m+1
        istack(1,m) = i3+1
        istack(2,m) = i2
c
        else
c
        m = m+1
        istack(1,m) = i3+1
        istack(2,m) = i2
c
        m = m+1
        istack(1,m) = i1
        istack(2,m) = i3
c
        endif
c
        goto 1000
 1100 continue
        end


        subroutine quicksort01(vals,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension vals(1)
c
c       Randomly choose a pivot index.
c
        call corrand3(1,r)
        ipiv = i1+floor((i2-i1)*r)
c
c        ipiv = i1+(i2-i1)/2
c
        val  = vals(ipiv)
c
c       Swap the pivot element and the last element.
c
        vals(ipiv) = vals(i2)
        vals(i2)   = val
c
        i3 = i1
c
        do 1000 i=i1,i2-1
c
        if( vals(i) .lt. val) then
        d  = vals(i)
c
        vals(i)  = vals(i3)
        vals(i3) = d       
c
        i3=i3+1
        endif
c
 1000 continue
c
        dd = vals(i3)
        vals(i3) = vals(i2)
        vals(i2) = dd
c
        end


        subroutine reverse0(n,vals)
        implicit double precision (a-h,o-z)
        dimension vals(1)
c
c       Reverse the order of a list of reals.
c
        do 1000 j=1,n/2
        d = vals(j)
        vals(j)     = vals(n-j+1)
        vals(n-j+1) = d
 1000 continue
c
        end


        subroutine insorti(k,ia)
        implicit double precision (a-h,o-z)
        dimension ia(1)
c
        if (k .le. 1) return
c
        do 1000 i=2,k
        ival=ia(i)
        j=i-1
        do 1100 while (j .ge. 1 .AND. ia(j) .gt. ival) 
        ia(j+1)=ia(j)
        j=j-1
 1100 continue
        ia(j+1)=ival
 1000 continue
        end



        subroutine insorti2(k,ia,ib)
        implicit double precision (a-h,o-z)
        dimension ia(1),ib(1)
c
        if (k .le. 1) return
c
        do 1000 i=2,k
        ival=ia(i)
        ival2=ib(i)
        j=i-1
        do 1100 while (j .ge. 1 .AND. ia(j) .gt. ival) 
        ia(j+1)=ia(j)
        ib(j+1)=ib(j)
        j=j-1
 1100 continue
        ia(j+1)=ival
        ib(j+1)=ival2
 1000 continue
        end


        subroutine quicksorti(n,ivals)
        implicit double precision (a-h,o-z)
        dimension istack(2,20 000),ivals(1)
c
c       Sort a list of integers.
c
        if (n .lt. 60) then
        call insorti(n,ivals)
        return
        endif
c
        maxstack = 10 000
        k        = 60
c
        nstack      = 1
        istack(1,1) = 1
        istack(2,1) = n
c
 1000 continue
        if (nstack .eq. 0) goto 1100
        i1 = istack(1,nstack)
        i2 = istack(2,nstack)
        nstack=nstack-1
c
c
        l = i2-i1+1
        if (l .le. k) then
        call insorti(l,ivals(i1))
        goto 1000
        endif
c
c       Otherwise perform quicksort.
c
        call quicksorti0(ivals,i1,i2,i3)
c
c       This should never happen, but just in case ...
c
        if (nstack+2 .ge. maxstack) then
        print *,"quicksort out of memory"
        stop
        endif
c
c       Make sure the smaller half is processed first to reduce storage
c       to O(logn).
c             
        n1 = i3-i1+1
        n2 = i2-(i3+1)+1
c
        if (n2 .lt. n1) then
c
        nstack = nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
c
        nstack = nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
c
        else
c
        nstack=nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
c
        nstack=nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
c
        endif
c
        goto 1000
 1100 continue
        end


        subroutine quicksorti0(ivals,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension ivals(1)
c
c       Randomly choose a pivot index.
c
        call corrand3(1,r)
        ipiv = i1+floor((i2-i1)*r)
c
c        ipiv = i1+(i2-i1)/2
c
        ival  = ivals(ipiv)
c
c       Swap the pivot element and the last element.
c
        ivals(ipiv) = ivals(i2)
        ivals(i2)   = ival
c
        i3 = i1
c
        do 1000 i=i1,i2-1
        if( ivals(i) .lt. ival) then
        id  = ivals(i)
c
        ivals(i)  = ivals(i3)
        ivals(i3) = id
c
        i3=i3+1
        endif
 1000 continue
c
        id        = ivals(i3)
        ivals(i3) = ivals(i2)
        ivals(i2) = id
c
        end


        subroutine quicksorti2(n,ivals,ivals2)
        implicit double precision (a-h,o-z)
        dimension istack(2,20 000),ivals(1),ivals2(1)
c
c       Sort a list of integers.
c
        if (n .lt. 60) then
        call insorti2(n,ivals,ivals2)
        return
        endif
c
        maxstack = 10 000
        k        = 60
c
        nstack      = 1
        istack(1,1) = 1
        istack(2,1) = n
c
 1000 continue
        if (nstack .eq. 0) goto 1100
        i1 = istack(1,nstack)
        i2 = istack(2,nstack)
        nstack=nstack-1
c
        l = i2-i1+1
        if (l .le. k) then
        call insorti2(l,ivals(i1),ivals2(i1))
        goto 1000
        endif
c
c       Otherwise perform quicksort.
c
        call quicksorti20(ivals,ivals2,i1,i2,i3)
c
c       This should never happen, but just in case ...
c
        if (nstack+2 .ge. maxstack) then
        print *,"quicksort out of memory"
        stop
        endif
c
c       Make sure the smaller half is processed first to reduce storage
c       to O(logn).
c             
        n1 = i3-i1+1
        n2 = i2-(i3+1)+1
c
        if (n2 .lt. n1) then
c
        nstack = nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
c
        nstack = nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
c
        else
c
        nstack=nstack+1
        istack(1,nstack) = i3+1
        istack(2,nstack) = i2
c
        nstack=nstack+1
        istack(1,nstack) = i1
        istack(2,nstack) = i3
c
        endif
c
        goto 1000
 1100 continue
        end


        subroutine quicksorti20(ivals,ivals2,i1,i2,i3)
        implicit double precision (a-h,o-z)
        dimension ivals(1),ivals2(1)
c
c       Randomly choose a pivot index.
c
        call corrand3(1,r)
        ipiv = i1+floor((i2-i1)*r)
c
c        ipiv = i1+(i2-i1)/2
c
        ival  = ivals(ipiv)
        ival2 = ivals2(ipiv)
c
c       Swap the pivot element and the last element.
c
        ivals(ipiv) = ivals(i2)
        ivals(i2)   = ival
c
        ivals2(ipiv) = ivals2(i2)
        ivals2(i2)   = ival2
c
        i3 = i1
c
        do 1000 i=i1,i2-1
        if( ivals(i) .lt. ival) then
        id  = ivals(i)
        id2 = ivals2(i)
c
        ivals(i)  = ivals(i3)
        ivals(i3) = id
c
        ivals2(i)  = ivals2(i3)
        ivals2(i3) = id2

        i3=i3+1
        endif
 1000 continue
c
        id        = ivals(i3)
        ivals(i3) = ivals(i2)
        ivals(i2) = id
c
        id2        = ivals2(i3)
        ivals2(i3) = ivals2(i2)
        ivals2(i2) = id2
c
        end


        subroutine move(n,a,b)
        implicit double precision (a-h,o-z)
        dimension a(1),b(1)
        do 1000 j=1,n
        b(j)=a(j)
 1000 continue
        end




        subroutine iduplicates(n,idxs)
        implicit double precision (a-h,o-z)
        dimension idxs(1)
c
        if ( n.le. 1) return
c
c       First sort the list.
c
        if ( n .gt. 60 ) then
        call quicksorti(n,idxs)
        else
        call insorti(n,idxs)        
        endif
c
c       Now remove duplicates.
c
        j=2
        do 1000 i=2,n
        if (idxs(j-1) .ne. idxs(i)) then
        idxs(j) = idxs(i)
        j=j+1
        endif        
 1000 continue
c
        n=j-1
        end


        subroutine imove(n,ia,ib)
        implicit double precision (a-h,o-z)
        dimension ia(1),ib(1)
        do 1000 j=1,n
        ib(j)=ia(j)
 1000 continue
        end


        SUBROUTINE corrand3(N,Y)
        IMPLICIT double precision (A-H,O-Z)
        save
        DIMENSION Y(1)
        DATA IFCALL /0/
c
c        This subroutine returns to the user a collection of
c        reasonably random numbers uniformly distributed on
c        the interval [0,1]
c
C       . . . CONDUCT PRELIMINARY RANDOMIZATION
C 
        IF (IFCALL.EQ.1 ) GOTO  1300
        IFCALL =1
c
c        construct parameters for the first process
c
        LAMBDA=3661
        MU=30809
        IP=145800
        M=41
c
c        construct parameters for the second process
c
        LAMBDAq=8121
        MUq=28411
        IPq=134456
        Mq=43
c
c        construct parameters for the third process
c
        LAMBDAqq=3613
        MUqq=45289
        IPqq=214326
        Mqq=53
c
        DO 1200 I=1,100
        M1=M*LAMBDA  +MU
        J=M1/IP
        M=M1-J*IP
c
        M1=Mq*LAMBDAq  +MUq
        J=M1/IPq
        Mq=M1-J*IPq
c
        M1=Mqq*LAMBDAqq  +MUqq
        J=M1/IPqq
        Mqq=M1-J*IPqq
c
 1200 CONTINUE
c
 1300 CONTINUE
c
        D=1
        D=D/IP
        Dq=1
        Dq=Dq/IPq
        Dqq=1
        Dqq=Dqq/IPqq
C 
C       GENERATE PSEUDO-RANDOM NUMBERS
C 
        DO 1400 I=1,N
c
        M1=M*LAMBDA +MU
        J=M1/IP
        M=M1-J*IP
c
        M1=Mq*LAMBDAq +MUq
        J=M1/IPq
        Mq=M1-J*IPq
c
        M1=Mqq*LAMBDAqq +MUqq
        J=M1/IPqq
        Mqq=M1-J*IPqq
c
        Y(I)=M*D+Mq*Dq+mqq*dqq
c
        call corrand_comp(y(i),rint)
c
        y(i)=rint
 1400 CONTINUE
c
        RETURN
c
c
c
c
        entry corrand3_init(dummy)
c
        ifcall=0
        return
        END
c


        SUBROUTINE corrand_comp(x,rint)
        IMPLICIT double precision (A-H,O-Z)
        save
c
        if (x .lt. 1) then
            rint=x**3/6
            return
        endif
c      
        if ( (x .ge. 1) .and. (x .le. 2) ) then
            rint= 0.75d0*x-(x-1.5d0)**3/3 - 0.625d0
            return
        endif
c
        rint= (x-3)**3/6 +1
c
        return
        end


        subroutine izero(k,ia)
        implicit double precision (A-h,o-z)
        dimension ia(1)
        do 1000 j=1,k
        ia(j)=0
 1000 continue
        end





        subroutine zero(k,a)
        implicit double precision (a-h,o-z)
        dimension a(1)
        do 1000 i=1,k
        a(i)=0
 1000 continue
        end


        SUBROUTINE corrand_norm(N,Y1,y2)
        IMPLICIT double precision (A-H,O-Z)
        save
        DIMENSION Y1(1),y2(1)
c
c        This subroutine returns to the user two random
c        vectors y1, y2, each of which is distiibuted
c        normally with the distribution density
c
c        1/sqrt(2 * \pi) * e^(y^2/2)                         (1)
c   
c
c              Input parameters:
c
c  n - the number of elements to be returned in each of the 
c        arrays y1, y2
c
c              Output parameters:
c  y1, y2 - two pseudo-random arrays distributed normally 
c        with the distribution density (1)
c
c
c        . . . construct vectors of variables distributed
c              uniformly on the interval [0,1]
c
        call corrand3(N,Y1)
        call corrand3(N,Y2)
c
c       combine the variables y1, y2 converting them
c       into variables distributed normally (Box-Muller 
c       algorithm)
c
        done=1
        pi=atan(done)*4
        do 1400 i=1,n
c
        z1=sqrt(-2*log(y1(i)))*cos(2*pi*y2(i))
        z2=sqrt(-2*log(y1(i)))*sin(2*pi*y2(i))
c
        y1(i)=z1
        y2(i)=z2
 1400 continue
c
        return
        end


        subroutine corrand_integer_knuth(n,ixs)
        implicit double precision (a-h,o-z)
        save
        dimension ixs(n),tt(12)
c
c        This subroutine returns to the user a random permutation 
c        of n integer numbers 1,2,...,n.
c
c              Input parameters:
c
c  n - random numbers in the array ixs will be the distributed
c        uniformly on the interval [1,n]
c
c              Output parameters:
c
c  ixs - a pseudo-random permutation of length n
c
        call corrand3(11,tt)
        call corrand3(11,tt)
        call corrand3(11,tt)
        do 1200 i=1,n
c
        ixs(i)=i
 1200 continue
c
        done=1
        do 1400 i=1,n-1
c
        call corrand3(1,tt)
c
        k=n-i+1

cccc        call prinf('k=*',k,1)

        h=done/k
        j=tt(1)/h+1

cccc        call prinf('and j=*',j,1)
c
        jj=ixs(k)
        ixs(k)=ixs(j)
        ixs(j)=jj
 1400 continue
c
        return
        end
