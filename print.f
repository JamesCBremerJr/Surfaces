        subroutine prini()
c
c       This is a now defunct initialization routine for the print 
c       routines found below.  
c
        end


        subroutine prin2(msg,a,n)
        implicit double precision (a-h,o-z)
        dimension a(1)
        character*1 msg(1),AST
        data AST / '*' /
        data iw1  / 13 /, iw2 / 6 /
c
c       Output an array of double words as a formatted column; show the 
c       most significant 5 digits of each value.
c
 0100 format(1X,80A)
 0200 format(8(2X,E12.5))
c
        len = 0
        do 1000 j=1,10 000
        if (msg(j) .eq. AST) goto 1100
        len=len+1
 1000 continue
 1100 continue
        write (iw1,0100) (msg(j),j=1,len)
        write (iw2,0100) (msg(j),j=1,len)
        write (iw1,0200) (a(j),j=1,n)
        write (iw2,0200) (a(j),j=1,n)
c
        end



        subroutine prin3(msg,a,n)
        implicit double precision (a-h,o-z)
        dimension a(1)
        character*1 msg(1),AST
        data AST / '*' /
        data iw1  / 13 /, iw2 / 6 /
c
c       Output an array of double words as a formatted column; show the 
c       most significant 15 digits of each value.
c
 0100 format(1X,80A)
 0200 format(4(2X,E22.15))
c
        len = 0
        do 1000 j=1,10 000
        if (msg(j) .eq. AST) goto 1100
        len=len+1
 1000 continue
 1100 continue
        write (iw1,0100) (msg(j),j=1,len)
        write (iw2,0100) (msg(j),j=1,len)
        write (iw1,0200) (a(j),j=1,n)
        write (iw2,0200) (a(j),j=1,n)
c
        end


        subroutine prinf(msg,ia,n)
        implicit double precision (a-h,o-z)
        dimension ia(1)
        character*1 msg(1),AST
        data AST / '*' /
        data iw1  / 13 /, iw2 / 6 /
c
c       Output an array of integers in a formatted column.
c
 0100 format (1X,80A1)
 0200 format (10(1X,I10))
c
        len = 0
        do 1000 j=1,10 000
c           print *,j,msg(j)
        if (msg(j) .eq. AST) goto 1100
        len=len+1
 1000 continue
 1100 continue
c
        write (iw1,0100) (msg(j),j=1,len)
        write (iw2,0100) (msg(j),j=1,len)
        write (iw1,0200) (ia(j),j=1,n)
        write (iw2,0200) (ia(j),j=1,n)
c        
        end


        subroutine prina(msg)
        implicit double precision (a-h,o-z)
        character*1 msg(1),AST
        data AST / '*' /
        data iw1  / 13 /, iw2 / 6 /
c
c       Output an asterisk terminated string.
c
 0100 format (1X,80A1)
c
        len = 0
        do 1000 j=1,10 000
        if (msg(j) .eq. AST) goto 1100
        len=len+1
 1000 continue
 1100 continue
c
        write (iw1,0100) (msg(j),j=1,len)
        write (iw2,0100) (msg(j),j=1,len)
c
        end




        subroutine print_sing(msg,n0,m0,amatr)
        implicit double precision (a-h,o-z)
        double precision amatr(n0,m0)
        double precision, allocatable :: bmatr(:,:),w(:),sing(:)
        integer*4, allocatable :: iw(:)
        integer*4 n,m,lw,info,l
c
        character*1 msg(1)
c
c       Print the singular values of a real-valued matrix; also, 
c       compute and print the condition number.  The standard LAPACK
c       singular value routine dgesdd is used to compute the singular 
c       values.
c
c                          Input Parameters:
c
c   msg - an asterisk terminated string to print along with the 
c       singular values of the matrix
c
c   (n,m) - dimensions of the input matrix
c   amatr - the (n,m) real-valued input matrix
c
        call mach_zero(eps0)
        if (eps0 .lt. 1.0d-20) then
        return
        endif
c
        n = n0
        m = m0
        l = max(n,m)        
        allocate(bmatr(n,m),sing(l),iw(8*(l)))
c
        do 1000 i=1,n
        do 1100 j=1,m
        bmatr(i,j)=amatr(i,j)
 1100 continue
 1000 continue
c      
c       Allocate workspace for the SVD.
c
        lw = -1
        call dgesdd('N',n,m,bmatr,n,sing,u,n,vt,n,ww,lw,iw,info)
        lw = ww
        allocate(w(lw))
c
c       Perform the SVD.
c
        call dgesdd('N',n,m,bmatr,n,sing,u,n,vt,n,w,lw,iw,info)
c
c
        if (info .ne. 0) then
        info0 = info
        call prinf("in print_sing, ddegsdd failed with info=*",info0,1)
        return
        endif
c
c       Report to the user as promised.
c
        call prin2(msg,sing,min(n,m))
        cond = sing(1)/sing(min(n,m))
        call prin3("cond = *",cond,1)
c
        end


        subroutine print_singz(msg,n,m,amatr)
        implicit double precision (a-h,o-z)
        double complex amatr(n,m)
        double complex, allocatable :: bmatr(:,:),w(:)
        double precision, allocatable :: sing(:),rw(:)
        integer, allocatable :: iw(:)
c
        character*1 msg(1)
c
c       Print the singular values of a complex-valued matrix; also, 
c       compute and print the condition number.  The standard LAPACK
c       singular value routine dgesdd is used to compute the singular 
c       values.
c
c                          Input Parameters:
c
c   msg - an asterisk terminated string to print along with the 
c       singular values of the matrix
c
c   (n,m) - dimensions of the input matrix
c   amatr - the (n,m) real-valued input matrix
c
        allocate(bmatr(n,m),sing(n+m),iw(10*(n+m)),rw(10*(n+m)))
c
        do 0100 i=1,n
        do 0200 j=1,m
        bmatr(i,j)=amatr(i,j)
 0200 continue
 0100 continue
c
c
c
c
c       Use LAPACK SVD zgesvd.
c      
        lw = -1
        call zgesvd('N','N',n,m,bmatr,n,sing,u,n,vt,n,ww,lw,rw,info)
        lw = ww
        lw = lw + 5 000 000
        if (allocated(w)) deallocate(w)
        allocate(w(lw))
        call zgesvd('N','N',n,m,bmatr,n,sing,u,n,vt,n,w,lw,rw,info)
c
        if (info .ne. 0) then
        call prinf("in print_singz, zgesvd failed with info = *",info,1)
        return
        endif
c
        goto 2000
c
c       Perform the SVD via lapack divide and conquer.
c$$$c
c$$$ 1000 continue
c$$$        lw = -1
c$$$        call zgesdd('N',n,m,bmatr,n,sing,u,n,vt,n,ww,lw,rw,iw,info)
c$$$        lw = ww
c$$$        lw = lw + 50 000 000
c$$$        if (allocated(w)) deallocate(w)
c$$$        allocate(w(lw))
c$$$
c$$$        call zgesdd('N',n,m,bmatr,n,sing,u,n,vt,n,w,lw,rw,iw,info)
c$$$c
c$$$        if (info .ne. 0) then
c$$$        call prinf("in print_singz, zegsdd failed with info=*",info,1)
c$$$        return
c$$$        return
c$$$        endif
c$$$        goto 2000
c$$$c
c$$$ 1100 continue
c$$$c
c$$$c       Compute the SVD using d2mcudv.
c$$$c
c$$$        ifuv=0
c$$$        eps=1.0d-12
c$$$        call d2mcudv(bmatr,n,u,v,eps,ifuv,numrows)
c$$$c
c$$$        do 1500 j=1,n
c$$$        sing(j)=bmatr(j,1)
c$$$ 15000 continue
c
c       Report to the user as promised.
c
 2000 continue
        call prin2(msg,sing,min(n,m))
        cond = sing(1)/sing(min(n,m))
        call prin3("cond = *",cond,1)
c
        end


        subroutine mach_zero(zero_mach)
        implicit double precision (a-h,o-z)
        data eps0 / -1.0d0 /
        save
c
c       Only perform the computation once and save the results.
c
        if (eps0 .ne. -1.0d0) then
        zero_mach=eps0
        return
        endif
c
c       Approximate machine zero in some reasonable way.
c
        zero_mach=1
        d = 1
c
 1000 continue
        if (d+zero_mach/2 .eq. d) goto 1100
        zero_mach=zero_mach/2
        goto 1000
 1100 continue
        return
        
        end


        subroutine elapsed(t)
        implicit double precision (a-h,o-z)
        integer*8 i,irate
        real t1
c
c       Return the elapsed time from the beginning of program execution
c       in seconds as a double word.  
c
c       NOTE: the resolution of the timer which is queried is operating
c       system and compiler dependent.  This routines seems to produce
c       reasonable results on a wide range of modern systems, however.
c
        call system_clock(i,irate)
c
        dd = i
        dd = dd/irate
        t = dd
c
c        call cpu_time(t1)
c        t = t1
c
        return
        end
