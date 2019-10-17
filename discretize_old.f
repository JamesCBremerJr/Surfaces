        implicit double precision (a-h,o-z)
        double precision, allocatable :: disc(:),eval(:,:),zs(:,:)
        double complex, allocatable   :: x(:),y(:),sums(:),vals(:,:)
        double complex, allocatable   :: vals2(:),sums2(:,:)
        double precision, allocatable :: tris(:,:),errs(:)
        integer, allocatable          :: itars(:)
c
        dimension errmaxs(1000),errareas(1000),errl2s(1000)
        dimension ntris00(1000)
c
        double complex val,sum,val2
        external surface,kernel
c
        integer omp_get_thread_num
c
        pi = acos(-1.0d0)
c
        norder  = 8
        nself   = norder
        dnear   = 8.00d0
        alpha   = 0.0125d0
c
        isplit  = 1
c
        call prinf("nself  = *",nself,1)
        call prin2("alpha = *",alpha,1)
c
        call surface_init(alpha)
c
        ldisc   = 200 000 000
        maxtris = 10 000
c
        allocate(disc(ldisc),tris(6,maxtris))
c
        niters=7
        do 0100 iter=1,niters
c
        call disc_init(ier,norder,nself,dnear,disc,ldisc)
        if (ier .ne. 0) then
        call prinf("after disc_init, ier = *",ier,1)
        stop
        endif
c
c       Add a decomposition of the ellipsoid to the decomposition.
c
        x1 = -1.0d0
        x2 =  1.0d0
c
        y1 = -1.0d0
        y2 =  1.0d0
c
        do 1000 iparam=1,6
c
        if (isplit .eq. 1) then
        nlevels = iter-1
        call triangulate2(ier,iparam,x1,y1,x2,y2,nlevels,
     -    maxtris,ntris,tris)
        else
c
        epsarea = 16.0d0 / 4.0d0**iter
        call triangulate1(ier,iparam,x1,y1,x2,y2,epsarea,
     -    maxtris,ntris,tris)
c
        endif
c
        if (ier .ne. 0) then
        call prinf("after triangulate, ier = *",ier,1)
        stop
        endif
c
        call disc_add(ier,disc,iparam,ntris,tris)
        if (ier .ne. 0) then
        call prinf("after disc_add, ier = *",ier,1)
        stop
        endif
 1000 continue
c
c       Fetch and display information about the decomposition stored
c       in the disc structure.
c
        call disc_info(disc,ntris,nquad,n,lkeep)
c
        call prinf("total # of triangles ntris = *",ntris,1)
        call prinf("total # of nodes n = *",n,1)
c
        dd = lkeep*8
        dd = dd/(1024d0**2)
c
        ntris00(iter) = ntris
c
        call prin2("disc structure size (in MB) = *",dd,1)
c
c       Fetch evaluation data and plot the discretization nodes.
c
        if (allocated(eval)) deallocate(eval)
        allocate(eval(13,n))
        call disc_data(disc,eval)
c
        dd = 0
        do 0200 ii=1,n
        dd = dd + eval(1,ii)**2
 0200 continue
c     
        errarea=abs(dd-4*pi)
        call prin2("area error = *",errarea,1)
c
c       Plot the discretization nodes.
c
c$$$        if (allocated(zs)) deallocate(zs)
c$$$        allocate(zs(3,n))
c$$$        call disc_points(disc,zs)
c
c        iplot=1
c        call plot3d_points("discretization nodes *",iplot,n,zs)
c
c       Form the vector to which the matrix will be applied and choose a
c       set of random targets.
c
        ntars = 10
c
        if (allocated(x)) deallocate(x,y,itars,errs,sums)
        allocate(x(n),y(n),itars(n),errs(n),sums(n))
c
        do 2000 i=1,n
        x(i)     = eval(1,i)
        y(i)     = 0
        itars(i) = i
        sums(i)  = 0
 2000 continue
c
        call corrand_integer_knuth(n,itars)
c
c       Test the evaluation code.
c
        call elapsed(t1)
c
        if (allocated(sums2)) deallocate(sums2)
        allocate(sums2(ntars,ntris))
c
c
!$omp   parallel default (shared)
!$omp!  private(jtri,i,j,vals,sum)
        allocate(vals(nquad,ntars))
!$omp   do 
        do 3100 jtri=1,ntris
c
        call disc_eval(disc,jtri,ntars,itars,vals,
     -     kernel,par1,par2,par3,par4)
c
        do 3200 i=1,ntars
        sum = 0
        do 3300 j=1,nquad
        sum = sum + vals(j,i)*x(j+(jtri-1)*nquad)
 3300 continue
        sums2(i,jtri) = sum
 3200 continue

 3100 continue
!$omp   end do
!$omp   end parallel
c
        call elapsed(t2)
c
        do 3400 i=1,ntars
        sum = 0
        do 3500 jtri=1,ntris
        sum = sum + sums2(i,jtri)
 3500 continue
        sums(i) = sum
 3400 continue
c
        errmax = 0
        errl2  = 0
        do 3600 i=1,ntars
        errabs  = abs(sums(i)-eval(1,itars(i))/2)
        errl2   = errl2+errabs**2
        errmax  = max(errabs,errmax)
        errs(i) = errabs
 3600 continue
c
        errl2       = sqrt(errl2)
        errmaxs(iter)  = errmax
        errl2s(iter)   = errl2
        errareas(iter) = errarea
c
        call prin2("errs = *",errs,ntars)
        call prin2("errl2 = *",errl2,1)
        call prin2("errmax = *",errmax,1)
c
        call disc_stats(disc,avgself)
        call prin2("avgself = *",avgself,1)
c
        call prin2("time elapsed = *",t2-t1,1)
 0100 continue
c
        call prin2("errmaxs  = *",errmaxs,niters)
        call prin2("errareas = *",errareas,niters)
        call prin2("errl2s   = *",errl2s,niters)
c
        call prina("*")
c
 9999 format(I6.6," ",D9.2,"  (",I12.12,")")
c
        call prina("---------------------*")
        call prina("*")        
c
        call prinf("norder = *",norder,1)
        call prinf("nself  = *",nself,1)
        call prin2("dnear  = *",dnear,1)
        call prin2("alpha  = *",alpha,1)
        call prinf("isplit = *",isplit,1)
        call prina("*")
c
        do 1111 i=1,niters
        dratio = 0
        if (i .gt. 1) dratio = errl2s(i-1)/errl2s(i)
c
        ii = ceiling(dratio)
        write(*,9999)  ntris00(i),errl2s(i),ii
        write(13,9999) ntris00(i),errl2s(i),ii
 1111 continue
c
        call prina("*")
        call prina("---------------------*")

c$$$        do 9999 i=1,niters-1
c$$$        print *,errmaxs(i)/errmaxs(i+1)
c$$$ 9999 continue
c
        end


        subroutine triangle_list(ier,x1,y1,x2,y2,nlevels,
     -    maxtris,ntris,tris)
        implicit double precision (a-h,o-z)
        integer, allocatable :: istack(:,:)
c
c        dimension istack(2,maxtris+100),tris(6,1)
        dimension children(6,4),tris(6,1)
c
        allocate(istack(2,maxtris+100))
c
c       Subdivide the parameterization domain.
c
        ier = 0
        pi  = acos(-1.0d0)
c
        if (maxtris .lt. 10) then
        ier = 4
        return
        endif
c
c       Construct the top-level triangles.
c 
        nstack        = 0
        ntris         = 0
c
        nstack = nstack+1
        ntris  = ntris+1
        istack(1,nstack) = ntris
        istack(2,nstack) = 0
        tris(1,ntris)    = x1
        tris(2,ntris)    = y1
        tris(3,ntris)    = x2
        tris(4,ntris)    = y1
        tris(5,ntris)    = x2
        tris(6,ntris)    = y2
c
        nstack = nstack+1
        ntris  = ntris+1
        istack(1,nstack) = ntris
        istack(2,nstack) = 0
        tris(1,ntris)    = x1
        tris(2,ntris)    = y1
        tris(3,ntris)    = x1
        tris(4,ntris)    = y2
        tris(5,ntris)    = x2
        tris(6,ntris)    = y2
c
 1000 continue
        if (nstack .eq. 0) goto 2000
c
        itri   = istack(1,nstack)
        level  = istack(2,nstack)
        nstack = nstack-1
c
        if (level .ge. nlevels) goto 1000
c
c       Otherwise, split the triangle.
c
        if (ntris+5 .gt. maxtris) then
        ier = 4
        return
        endif
c
        call split_tri(tris(1,itri),children)
c
        do 1100 i=1,6
        tris(i,itri) = children(i,1)
 1100 continue
        nstack=nstack+1
        istack(1,nstack)=itri
        istack(2,nstack)=level+1
c        
        do 1200 l=2,4
        ntris=ntris+1
        nstack=nstack+1
        istack(1,nstack)=ntris
        istack(2,nstack)=level+1
        do 1300 i=1,6
        tris(i,ntris) = children(i,l)
 1300 continue
 1200 continue
c
        goto 1000
 2000 continue
c
        end


        subroutine torus_triangles(M,N,ntris,tris)
        implicit double precision (a-h,o-z)
        dimension tris(6,1)
c
        pi = acos(-1.0d0)
c
        ntris = 0
c
        dd1 = 2*pi/M
        dd2 = 2*pi/N
c
        do 1000 i=1,M
        do 1100 j=1,N
c
        x1 = (i-1)*dd1
        x2 = i*dd1
c
        y1 = (j-1)*dd2
        y2 = j*dd2
c
        ntris = ntris + 1
        tris(1,ntris) = x1
        tris(2,ntris) = y1
        tris(3,ntris) = x2
        tris(4,ntris) = y1
        tris(5,ntris) = x1
        tris(6,ntris) = y2
c
        ntris = ntris + 1
        tris(1,ntris) = x2
        tris(2,ntris) = y1
        tris(3,ntris) = x2
        tris(4,ntris) = y2
        tris(5,ntris) = x1
        tris(6,ntris) = y2
c
 1100 continue
 1000 continue
c
        end


        subroutine split_tri(tri,children)
        implicit double precision (a-h,o-z)
        dimension tri(2,3),children(2,3,4)
c
        x1 = tri(1,1)
        y1 = tri(2,1)
        x2 = tri(1,2)
        y2 = tri(2,2)
        x3 = tri(1,3)
        y3 = tri(2,3)
c
        x4 = (x1+x2)/2
        y4 = (y1+y2)/2
c
        x5 = (x1+x3)/2
        y5 = (y1+y3)/2
c
        x6 = (x2+x3)/2
        y6 = (y2+y3)/2
c
        children(1,1,1) = x1
        children(2,1,1) = y1
        children(1,2,1) = x4
        children(2,2,1) = y4
        children(1,3,1) = x5
        children(2,3,1) = y5
c
        children(1,1,2) = x6
        children(2,1,2) = y6
        children(1,2,2) = x4
        children(2,2,2) = y4
        children(1,3,2) = x5
        children(2,3,2) = y5
c
        children(1,1,3) = x3
        children(2,1,3) = y3
        children(1,2,3) = x5
        children(2,2,3) = y5
        children(1,3,3) = x6
        children(2,3,3) = y6
c
        children(1,1,4) = x4
        children(2,1,4) = y4
        children(1,2,4) = x2
        children(2,2,4) = y2
        children(1,3,4) = x6
        children(2,3,4) = y6
c        
        end



        subroutine kernel(x,y,val,par1,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension x(13),y(13)
        double complex val
c
        call kdouble0(x,y,val)
c
c        val = x(1)*y(1)
c
        end
        

        subroutine surface(iparam,s,t,r,dr)
        implicit double precision (a-h,o-z)
        dimension r(3),dr(3,2)
        save
c
        pi = acos(-1.0d0)
c
        if (iparam .eq. 1) then
        dd   = sqrt(1+s**2+t**2)
c
        r(1)    = alpha*s/dd
        r(2)    = t/dd
        r(3)    = 1/dd
c
        dr(1,1) = -alpha*s**2/dd**3+alpha/dd
        dr(2,1) = -s*t/dd**3
        dr(3,1) = -s/dd**3
c
        dr(1,2) = -alpha*s*t/dd**3
        dr(2,2) = -t**2/dd**3+1/dd
        dr(3,2) = -t/dd**3
c
        return
        endif
c
        if (iparam .eq. 2) then
        dd      = sqrt(1+s**2+t**2)
c
        r(1)    = alpha*t/dd
        r(2)    = s/dd
        r(3)    = -1/dd
c
        dr(1,1) = -alpha*s*t/dd**3
        dr(2,1) = -s**2/dd**3+1/dd
        dr(3,1) = s/dd**3
c
        dr(1,2) = -alpha*t**2/dd**3+alpha/dd
        dr(2,2) = -s*t/dd**3
        dr(3,2) = t/dd**3
        return
        endif
c
        if (iparam .eq. 3) then
        dd      = sqrt(1+s**2+t**2)
c
        r(1)    = alpha/dd
        r(2)    = s/dd
        r(3)    = t/dd
c
        dr(1,1) = -alpha*s/dd**3
        dr(2,1) = -s**2/dd**3+1/dd
        dr(3,1) = -s*t/dd**3
c
        dr(1,2) = -alpha*t/dd**3
        dr(2,2) = -s*t/dd**3
        dr(3,2) = -t**2/dd**3 + 1/dd
        return
        endif
c
        if (iparam .eq. 4) then
        dd      = sqrt(1+s**2+t**2)
c
        r(1)    = -alpha/dd
        r(2)    = t/dd
        r(3)    = s/dd
c
        dr(1,1) = alpha*s/dd**3
        dr(2,1) = -s*t/dd**3
        dr(3,1) = -s**2/dd**3+1/dd
c
        dr(1,2) = alpha*t/dd**3
        dr(2,2) =-t**2/dd**3+1/dd
        dr(3,2) = -s*t/dd**3
        return
        endif
c
        if (iparam .eq. 5) then
c
        dd      = sqrt(1+s**2+t**2)
c
        r(1)    = alpha*s/dd
        r(2)    = -1/dd
        r(3)    = t/dd
c
        dr(1,1) = -alpha*s**2/dd**3+alpha/dd
        dr(2,1) = s/dd**3
        dr(3,1) = -s*t/dd**3
c
        dr(1,2) = -alpha*s*t/dd**3
        dr(2,2) = t/dd**3
        dr(3,2) = -t**2/dd**3+1/dd
c
        return
        endif
c
        if (iparam .eq. 6) then
c
        dd      = sqrt(1+s**2+t**2)
c
        r(1)    = alpha*t/dd
        r(2)    = 1/dd
        r(3)    = s/dd
c
        dr(1,1) = -alpha*s*t/dd**3
        dr(2,1) = -s/dd**3
        dr(3,1) = -s**2/dd**3 + 1/dd
c
        dr(1,2) = -alpha*t**2/dd**3+alpha/dd
        dr(2,2) = -t/dd**3
        dr(3,2) = -s*t/dd**3
        return
        endif
c
        print *,"!!!!!!!!"
        stop
        return
c
        entry surface_init(alpha0)
        alpha=alpha0
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the 
c       discretization code proper.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       !!EXPLAIN WHAT EVALUATION DATA IS!!!
c
c       !! EXPLAIN WHAT AN INTERACTION BETWEEN A SOURCE TRIANGLE
c          AND A POINT IS !!!!!!!!!!!!!!!!!!!!!!!!!!
c
c       REBUILD QUADRATURES FOR WEAKLY SINGULAR CASE SO AS TO INCREASE
c       EFFICIENCY
c
c       This file contains code for discretizing the L^2 action of
c       weakly singular integral operators given on surfaces.
c       The discretization procedure is described in some detail in the
c       paper ``A Nystr\"om method for the discretization of weakly
c       singular integral operators on surfaces."  In brief, the
c       user inputs a decomposition of a surface which gives rise
c       to a scheme for discretizing integral operators with weakly
c       singular kernels.
c
c       A decomposition of a surface \Sigma is a collection of smooth
c       parameterizations
c
c         \rho_j: T_j \to \mathbb{R}^3 
c
c       given on triangles such that \rho_j(T_j) form a disjoint union
c       of \Sigma together with a positive integer N.  Associated
c       with a decomposition of the surface \Sigma is the subspace 
c       S of L^2(\Sigma) which is the span of all functions of the 
c       form
c
c          \rho_j(p) \sqrt{|d\rho_j| },
c
c       where p is a polynomial of order N and |d\rho_j| denotes the
c       determinate of the Jacobian derviative d\rho_j of the 
c       parameterization \rho_j, and the mapping 
c
c          \Phi: S \to C^M
c
c       defined by
c
c          \Phi(f) = (f(x_1)\sqrt{w_1}, \ldots, f(x_N)\sqrt{w_N})
c
c       where x1,..,x_N,w_1,..,w_N is a quadrature compo
c       the triangle T_j.
c
c       The sets \rho_j(T_j) are referred to as surface elements; by
c       abuse of terminology we will often conflate the surface elements
c       with their corresponding triangles T_j in the parameterization
c       domain.
c
c       The user specifies the decomposition of a surface by providing
c       two inputs:
c
c         (1) a GLOBAL external subroutine called surface which supplies 
c             the parameterizations \rho_j
c
c         (2) the triangles T_j which comprise the parameterization
c             domains
c
c       Descriptions of the triangles T_j are stored in a discretization
c       structure along with other information about the decomposition.
c       Such a structure is constructed by the user by calling the
c       disc_init routine (below).  Information about the T_j is
c       added to the disc structure through calls to the routine
c       disc_add (below).
c
c
c       A decomposition of a surface \Sigma gives rise to a scheme for
c       discretizing integral operators of the form
c
c          Tf(x) = \int_       K(x,y) f(y) ds(y)
c                       \Sigma
c
c       where K(x,y) is a weakly singular kernel.  In particular,
c       the discretization A of the operator T associated with the
c       decomposition is the matrix which makes the diagram 
c
c                                
c                                   T
c                L^2(\Sigma) --------------->  L^2(\Sigma)
c                     |                             |
c                \Phi |                             | \Phi
c                     |                             |
c                    \_/            A              \_/
c                    C^n     --------------->      C^n
c       
c
c       commute.  
c
c       Each element of the matrix A discretizing an integral operator
c       with kernel K corresponds to a pair consisting of one source 
c       point and one target point.
c
c       This code provides routines for evaluating entries of the matrix 
c       A.  IT is most convenient
c

c       
c
c
c       The following subroutines are user-callable:
c
c   disc_init - initialize a discretization structure which will hold
c       data related to a decomposition
c
c   disc_add - add a collection of surface elements to the decomposition
c       stored in a discretization structure
c
c   disc_info - return some basic information about a discretization
c       structure
c
c   disc_data - return the evaluation data for the discretization nodes
c        associated with the decomposition stored in a discretization
c        structure
c
c   disc_points - return the coordinates of the discretization nodes
c        associated with the decomposition stored in a discretization
c        structure
c
c   disc_tri - return information about one of the triangles in the
c        discretization structure (for instance, coordinates of a 
c        bounding triangle)
c
c   disc_eval - evaluate a portion of the matrix discretizing a user-
c        supplied integral operator.  The operator is specified
c        via  arising from a decomposition described by a disc
c        structure.  More specifically, this routine evaluates 
c        a block of the discretization matrix corresonding to a 
c        specified source triangle and a collection of user-specified
c        target points
c
c   disc_eval_far - evaluate an entry of the matrix discretizing
c        a user-specfied integral operator under the assumption
c        the the target and source points are well-separated.
c
c   disc_eval2 - 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c$$$c        subroutine disc_eval(disc,jtri,ntars,itars,vals,
c$$$c     -     kernel,par1,par2,par3,par4)
c$$$
c$$$        subroutine disc_eval2(ier,disc,ntars,itars,nsrc,isrcs,vals,
c$$$     -   kernel,par1,par2,par3,par4)
c$$$        implicit double precision (a-h,o-z)
c$$$        dimension disc(1)
c$$$        dimension itars(1),isrcs(1)
c$$$        double complex vals(ntars,nsrcs)
c$$$c
c$$$        integer, allocatable :: jtris(:)
c$$$c
c$$$        ier = 0
c$$$c
c$$$c
c$$$c
c$$$
c$$$c
c$$$c       Make a list of the triangles
c$$$c        
c$$$        end




        subroutine disc_stats(disc,avgself)
        implicit double precision (a-h,o-z)
        dimension disc(1)
        dd1 = disc(50)
        dd2 = disc(51)
        avgself = dd2/dd1
c
c        dd1 = disc(52)
c        dd2 = disc(53)
c        avgnear = dd2/dd1
c
        end


        subroutine disc_eval_self(disc,i,jtri,vals,
     -    kernel,par1,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension disc(1)
c
        double complex vals(1)
c
        external kernel
        integer omp_get_thread_num
c
c       Evaluate a portion of the matrix discretizing an integral operator
c       of the form
c
c          Tf(x) = \int K(x,y) f(y) ds(y)
c
c       where K(x,y) is a weakly singular integral kernel.
c
c     
c                           Input Parameters:
c
c   disc - a discretization structure
c   i -
c   jtri -
c
c                           Output Parameters:
c
c   vals - the 
c
        ier = 0
c
        ithread = omp_get_thread_num()
c
c       Fetch data from the discretization structure.
c
        norder   = disc(1)
        nself    = disc(2)
        maxtris  = disc(3)
        ntris    = disc(4)
        len      = disc(5)
        itris    = disc(6)
        maxquad  = disc(7)
        lw       = disc(8)
        iw       = disc(9)+ithread*lw
        eps0     = disc(10)
c
        nquad    = disc(20)
        ixs      = disc(21)
        iys      = disc(22)
        iwhts    = disc(23)
        iu       = disc(24)
        iv       = disc(25)
c
        irad     = disc(40)
c
c       Determine the index of the triangle containing the target point
c       and the index of the target point in that triangle.
c
        itri = jtri-1
        ipt  = i-(itri*nquad)
        itri = itri+1
c
        iptr = itris+(itri-1)*len+11+(ipt-1)*13
        jptr = itris+(jtri-1)*len
c
c       Handle a self interaction.
c
        ixsself   = iw
        lxsself   = maxquad
c
        iysself   = ixsself+lxsself
        lysself   = maxquad
c
        iwhtsself = iysself+lysself
        lwhtsself = maxquad
c        
        ivals0    = iwhtsself+lwhtsself
        lvals0    = nquad
c
        ivals1    = ivals0+lvals0
        lvals1    = 2*maxquad
c
        ivals2    = ivals1+lvals1
        lvals2    = 2*nquad
c
        call disc_eval_self1(norder,nself,nquad,disc(ixs),
     -    disc(iys),disc(iwhts),disc(iu),disc(iptr),len,disc(itris),
     -    jtri,ipt,maxquad,disc(ixsself),disc(iysself),disc(iwhtsself),
     -    disc(ivals0),disc(ivals1),disc(ivals2),disc(irad),vals,
     -    kernel,par1,par2,par3,par4,nquad00)
c
        disc(50) = disc(50) + 1
        disc(51) = disc(51) + nquad00
        return
        end



        subroutine disc_eval_self1(norder,nself0,nquad,xs,ys,
     -    whts,u,x,len,tris,jtri,ipt,maxquad,xsself,ysself,whtsself,
     -    vals0,vals1,vals2,rad,vals,kernel,par1,par2,par3,par4,nquad00)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1),u(nquad,nquad),x(13),y(13)
        dimension tris(len,1),xsself(1),ysself(1),whtsself(1)
        dimension amap(2,3),ainv(2,3),r(3),dr(3,2),dn(3)
c
        double complex vals1(maxquad),vals2(nquad)
c
        double complex val,sum,vals(1)
        dimension vals0(nquad)
        external kernel
c
c       Form the map which takes the simplex to the integration
c       domain and its inverse.
c
        ier = 0
c
        iparam=tris(7,jtri)
        call disc_trimap2(tris(1,jtri),amap,ainv,det)
c
c       First construct the quadrature.
c
        s0 = amap(1,1)*xs(ipt)+amap(1,2)*ys(ipt)+amap(1,3)
        t0 = amap(2,1)*xs(ipt)+amap(2,2)*ys(ipt)+amap(2,3)
c
        if (nself0 .eq. 0) then
        n            = 300
        norder0      = 16
c
        call selfquad0(ier,n,norder0,tris(1,jtri),s0,t0,
     -    nself,xsself,ysself,whtsself)
c
        else
        call selfquad(ier,rad,tris(1,jtri),s0,t0,x(5),
     -    nself,xsself,ysself,whtsself)
        endif
c
        if (ier .ne. 0) then
        call prinf("in disc_eval, selfquad failed with ier = *",
     -    ier,1)
        stop
        endif
c
        nquad00 = nself
c
c       Evaluate the kernel at the nodes of the quadrature.
c
        do 1000 i=1,nself
        s   = xsself(i)
        t   = ysself(i)
        wht = whtsself(i)
c
        call surface(iparam,s,t,r,dr)
c
c       Construct the normal vector.
c
        dn(1) = dr(2,1)*dr(3,2) - dr(3,1)*dr(2,2)
        dn(2) = dr(3,1)*dr(1,2) - dr(1,1)*dr(3,2)
        dn(3) = dr(1,1)*dr(2,2) - dr(2,1)*dr(1,2)
        da    = sqrt(dn(1)**2+dn(2)**2+dn(3)**2)
c
        dn(1) = dn(1)/da
        dn(2) = dn(2)/da
        dn(3) = dn(3)/da
c
        wht = wht * da
c
        y(1)  = sqrt(wht)
        y(2)  = r(1)
        y(3)  = r(2)
        y(4)  = r(3)
        y(5)  = dr(1,1)
        y(6)  = dr(2,1)
        y(7)  = dr(3,1)
        y(8)  = dr(1,2)
        y(9)  = dr(2,2)
        y(10) = dr(3,2)
        y(11) = dn(1)
        y(12) = dn(2)
        y(13) = dn(3)
c
        call kernel(x,y,vals1(i),par1,par2,par3,par4)
 1000 continue
c
c
c       Form the product of the values vector with the matrix of 
c       polynomials values.
c
        do 1100 i=1,nquad
        vals2(i)=0
 1100 continue
c
        do 1200 i=1,nself
        s   = ainv(1,1)*xsself(i)+ainv(1,2)*ysself(i)+ainv(1,3)
        t   = ainv(2,1)*xsself(i)+ainv(2,2)*ysself(i)+ainv(2,3)
        wht = whtsself(i)/det
c
        call koorn(norder,s,t,vals0)
c
        do 1300 j=1,nquad
        vals2(j) = vals2(j) + vals0(j)*vals1(i)*sqrt(wht)
 1300 continue
 1200 continue
c
        do 1400 j=1,nquad
        sum = 0
        do 1500 i=1,nquad
        sum = sum + vals2(i)*u(i,j)
 1500 continue
        vals(j) = sum
 1400 continue
c
        end



        subroutine disc_eval(disc,jtri,ntars,itars,vals,
     -     kernel,par1,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension disc(1),itars(1)
        double complex vals(1)
        external kernel
c
c       Evaluate a block of the matrix discretizing an integral 
c       operator.  The operator is specified via an external subroutine
c       for evaluating its kernel and the 
c
c                              Input Parameters:
c
c                              Output Parameters:
c
c
c
c       Fetch data from the disc structure.
c
        norder     = disc(1)
        nself      = disc(2)
        maxtris    = disc(3)
        ntris      = disc(4)
        len        = disc(5)
        itris      = disc(6)
        maxquad    = disc(7)
        lw         = disc(8)
        iw         = disc(9)
c
        nquad      = disc(20)
        ixs        = disc(21)
        iys        = disc(22)
        iwhts      = disc(23)
        iu         = disc(24)
        iv         = disc(25)
c
        nquadnear  = disc(30)
        ixsnear    = disc(31)
        iysnear    = disc(32)
        iwhtsnear  = disc(33)
        iamatrin   = disc(34)
c
        irad       = disc(40)
c
c       Call an auxilliary routine to shape arrays.
c
        call disc_eval0(disc,norder,nquad,disc(ixs),disc(iys),
     -   disc(iwhts),disc(iu),disc(iamatrin),ntars,itars,len,
     -   disc(itris),jtri,vals,
     -   kernel,par1,par2,par3,par4)
c
        end



        subroutine disc_eval0(disc,norder,nquad,xs,ys,whts,u,
     -    amatrin,ntars,itars,len,tris,jtri,vals,
     -    kernel,par1,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension disc(1)
        dimension xs(1),ys(1),whts(1),u(nquad,nquad),itars(1)
        dimension amatrin(4*nquad,nquad),tris(len,1)
        double complex vals(nquad,ntars)
c
        dimension children(6,4),xeval(13),verts(6),yeval(13)
        double complex val

        double precision, allocatable :: tree(:,:)
        integer, allocatable          :: levels(:,:),istack(:)
        integer, allocatable          :: ifprocessed(:)
        double complex, allocatable   :: vals0(:,:),vals1(:)
c
        external kernel
c
        dnear = disc(99)
c
c       Allocate memory for the tree structure.
c
        iparam = tris(7,jtri)
c
        minlevels = 0
        maxlevels = 30
        maxtris   = 10 000
c
        allocate(tree(16+13*nquad,maxtris),levels(2,maxlevels))
        allocate(vals1(4*nquad))
c
c       tree(1-6,j)  - coordinates of the triangle vertices
c       tree(7,j)    - bx
c       tree(8,j)    - by
c       tree(9,j)    - bz
c       tree(10,j)   - br
c       tree(11,j)   - ichild1
c       tree(12,j)   - ichild2
c       tree(13,j)   - ichild3
c       tree(14,j)   - ichild4
c       tree(15,j)   - level of the node
c       tree(16,j)   - first entry of the evaluation data of length 13*nquad
c
c       Construct the first level.
c
        nnodes      = 1
        levels(1,1) = 1
        levels(2,1) = 1
c
        call disc_eval1(nquad,xs,ys,whts,tris(1,jtri),tree(1,1),
     -    iparam)
c
        tree(15,1) = 1
c
c       Recursively decompose the triangle until each target point is in 
c       the far field.
c
        do 1000 level = 1,maxlevels-1
c
        levels(1,level+1) = nnodes+1
c
        inode1 = levels(1,level)
        nn     = levels(2,level)
c
        do 1100 inode=inode1,inode1+nn-1
c
        bx = tree(7,inode)
        by = tree(8,inode)
        bz = tree(9,inode)
        br = tree(10,inode)*dnear
c
c       Check to see if the target nodes are in the far field.
c
        ifsplit = 0
        do 1200 ii=1,ntars
        i    = itars(ii)
c
        itri = (i-1)/nquad
        ipt  = i-(itri*nquad)
        itri = itri+1
c
        if (itri .eq. jtri) goto 1200
c
        iptr = 12 + (ipt-1)*13
        x    = tris(iptr+1,itri)
        y    = tris(iptr+2,itri)
        z    = tris(iptr+3,itri)
c
        dd   = (x-bx)**2+(y-by)**2+(z-bz)**2
c
        if (dd .lt. br .OR. level .lt. minlevels) then
        ifsplit = 1
        goto 1300
        endif
 1200 continue
 1300 continue
c
        if (ifsplit .eq. 0) goto 1100
c
        if (nnodes+4 .ge. maxtris) then
        print *,
     -  "in disc_eval, maximum number of triangles exceeded"
c
        stop
        endif
c
c       If not, split it and by so doing create 4 nodes on the next
c       level.
c
        call disc_split_tri(tree(1,inode),children)
c
        do 1400 ii=1,4
        nnodes            = nnodes+1
        tree(10+ii,inode) = nnodes
        tree(15,nnodes)   = level+1
c
        call disc_eval1(nquad,xs,ys,whts,children(1,ii),
     -    tree(1,nnodes),iparam)
c
 1400 continue
 1100 continue
c
       newnodes = nnodes - levels(1,level+1)+1
       levels(2,level+1) = newnodes
       if (newnodes .eq. 0) goto 1999
 1000 continue
c
        print *,"in disc_eval0, maximum number of levels exceeded"
        stop
        return
 1999 continue
c
        nlevels = level
c
        call prinf("nlevels = *",nlevels,1)
        call prinf("nnodes = *",nnodes,1)
c
c       Allocate space for the evaluation of the kernel and for a
c       stack.
c
        allocate(vals0(nquad,nnodes),istack(nnodes))
        allocate(ifprocessed(nnodes))
c
c       Each near interaction is evaluated by traversing the tree.
c
        do 2000 ii=1,ntars
        i = itars(ii)
c
        ifprocessed = 0
c
c       Copy over the evaluation data and find the coordinates of the
c       target node.
c
        itri = (i-1)/nquad
        ipt  = i-(itri*nquad)
        itri = itri+1
c
c       Handle self interactions.
c
        if (itri .eq. jtri) then
        call disc_eval_self(disc,i,jtri,vals(1,ii),
     -    kernel,par1,par2,par3,par4)
        goto 2000
        endif
c
c        disc(52) = disc(52) + ntars
c
        iptr = 12 + (ipt-1)*13
        x    = tris(iptr+1,itri)
        y    = tris(iptr+2,itri)
        z    = tris(iptr+3,itri)
c
        call disc_move(13,tris(iptr,itri),xeval)
c
c       Initialize the stack used to traverse the tree.
c
        nstack    = 1
        istack(1) = 1
 2100 continue
c
c       Get the index of the top node on the stack (but do not pop
c       it off the stack).
c
        if (nstack .le. 0) goto 2900
        inode  = istack(nstack)
c
c       Determine whether or not the target node is in the far
c       field.
c
        bx    = tree(7,inode)
        by    = tree(8,inode)
        bz    = tree(9,inode)
        br    = tree(10,inode)
        level = tree(15,inode)
c
        dd = (bx-x)**2+(by-y)**2+(bz-z)**2
c     
c       If the node is distant, evaluate the kernel at the
c       appropriate nodes and mark the node as processed.
c
        if (dd .gt. br*dnear .AND. level .ge. minlevels) then
        iptr = 16
c
        do 2400 i=1,nquad
        call kernel(xeval,tree(iptr,inode),val,par1,par2,par3,par4)
        vals0(i,inode) = val
        iptr=iptr+13
 2400 continue
c
        ifprocessed(inode) = 1
        nstack             = nstack-1
        goto 2100
        endif
c
c       Check to see if the node's children haven been processed;
c       push any children onto the stack.
c
        ifproc = 1
c
        ichild1 = tree(11,inode)
        ichild2 = tree(12,inode)
        ichild3 = tree(13,inode)
        ichild4 = tree(14,inode)
c
        if (ifprocessed(ichild1) .eq. 0) then
        nstack         = nstack + 1
        istack(nstack) = ichild1
        ifproc         = 0
        endif
c
        if (ifprocessed(ichild2) .eq. 0) then
        nstack         = nstack + 1
        istack(nstack) = ichild2
        ifproc         = 0
        endif
c
        if (ifprocessed(ichild3) .eq. 0) then
        nstack         = nstack + 1
        istack(nstack) = ichild3
        ifproc         = 0
        endif
c
        if (ifprocessed(ichild4) .eq. 0) then
        nstack         = nstack + 1
        istack(nstack) = ichild4
        ifproc         = 0
        endif
c
        if (ifproc .eq. 0) goto 2100
c
c       If all children have been processed, pop the node off the
c       stack and apply the interpolation matrix.
c
        call disc_move(2*nquad,vals0(1,ichild1),vals1(1))
        call disc_move(2*nquad,vals0(1,ichild2),vals1(1+nquad))
        call disc_move(2*nquad,vals0(1,ichild3),vals1(1+2*nquad))
        call disc_move(2*nquad,vals0(1,ichild4),vals1(1+3*nquad))
c
c       Apply the interpolation matrix to compute the values for this
c       node.
c
        call disc_apply(4*nquad,nquad,amatrin,vals1,vals0(1,inode))
c
c       Mark the node as processed, pop it off the stack and proceed.
c
        ifprocessed(inode) = 1
        nstack             = nstack-1
        goto 2100
 2900 continue
c
c       Copy the values out of the top node.
c
        do 2500 i=1,nquad
        vals(i,ii) = vals0(i,1)
 2500 continue
 2000 continue
c
c
        end



        subroutine disc_eval1(nquad,xs,ys,whts,verts,tree,iparam)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1),verts(6)
        dimension tree(1),amap(2,3),dn(3),r(3),dr(3,2)
        dimension zs(3,nquad)
c
c       Construct evaluation data for a triangle.
c
        tree(1) = verts(1)
        tree(2) = verts(2)
        tree(3) = verts(3)
        tree(4) = verts(4)
        tree(5) = verts(5)
        tree(6) = verts(6)
c
        call disc_trimap(verts,amap,det)
c
c       Construct the evaluation data for the surface element.
c
        do 1000 i=1,nquad
        s   = amap(1,1)*xs(i)+amap(1,2)*ys(i)+amap(1,3)
        t   = amap(2,1)*xs(i)+amap(2,2)*ys(i)+amap(2,3)
        wht = whts(i)*det
c
        call surface(iparam,s,t,r,dr)
c
        zs(1,i) = r(1)
        zs(2,i) = r(2)
        zs(3,i) = r(3)
c
c       Construct the normal vector.
c
        dn(1) = dr(2,1)*dr(3,2) - dr(3,1)*dr(2,2)
        dn(2) = dr(3,1)*dr(1,2) - dr(1,1)*dr(3,2)
        dn(3) = dr(1,1)*dr(2,2) - dr(2,1)*dr(1,2)
        da    = sqrt(dn(1)**2+dn(2)**2+dn(3)**2)
        dn(1) = dn(1)/da
        dn(2) = dn(2)/da
        dn(3) = dn(3)/da
c
        wht = wht * da
c
        iptr = 16 + (i-1)*13
c
        tree(iptr)    = sqrt(wht)
        tree(iptr+1)  = r(1)
        tree(iptr+2)  = r(2)
        tree(iptr+3)  = r(3)
        tree(iptr+4)  = dr(1,1)
        tree(iptr+5)  = dr(2,1)
        tree(iptr+6)  = dr(3,1)
        tree(iptr+7)  = dr(1,2)
        tree(iptr+8)  = dr(2,2)
        tree(iptr+9)  = dr(3,2)
        tree(iptr+10) = dn(1)
        tree(iptr+11) = dn(2)
        tree(iptr+12) = dn(3)
 1000 continue
c
c
c       Construct a bounding sphere for the triangle.  Note: br is the
c       square of the radius of the bounding ball.
c
        call bounding_ball(nquad,zs,bx,by,bz,br)
c
        tree(7)  = bx
        tree(8)  = by
        tree(9)  = bz
        tree(10) = br
c
        end



        subroutine disc_eval_far(disc,i,j,val,kernel,par1,par2,par3,
     -    par4)
        implicit double precision (a-h,o-z)
        dimension disc(1)
        double complex val
        external kernel
c
c       Evaluate an entry of a discretization matrix corresponding to
c       a single source point and a single target in the far regime.
c
c       NOTE: 
c
c
c                          Input Parameters:
c
c   disc - the discretization structure
c   (i,j) - entry of the matrix to evaluate
c
c   kernel - 
c
c                         Output Parameters:
c
c   val - 
c
c
c       Fetch data from the disc structure.
c
        norder   = disc(1)
        nself    = disc(2)
        maxtris  = disc(3)
        ntris    = disc(4)
        len      = disc(5)
        itris    = disc(6)
c
        nquad    = disc(20)
        ixs      = disc(21)
        iys      = disc(22)
        iwhts    = disc(23)
        iu       = disc(24)
        iv       = disc(25)
c
        nquadnear  = disc(30)
        ixsnear    = disc(31)
        iysnear    = disc(32)
        iwhtsnear  = disc(33)
c
        irad     = disc(40)
c
c       Determine the index of the triangle and point within the
c       triangle.
c        
        itri = (i-1)/nquad
        ipt  = i-(itri*nquad)
        itri = itri+1
c
        jtri = (j-1)/nquad
        jpt  = j-(jtri*nquad)
        jtri = jtri+1
c
c       Find the index of the evaluation data for the target node and
c       source node.
c
        iptr = itris+11+(itri-1)*len+(ipt-1)*13
        jptr = itris+11+(jtri-1)*len+(jpt-1)*13
c
c       Call the kernel subroutine to evaluate.
c
        call kernel(disc(iptr),disc(jptr),val,par1,par2,par3,par4)
c
        end



        subroutine disc_init(ier,norder,nself,dnear,disc,ldisc)
        implicit double precision (a-h,o-z)
        dimension disc(1),simp(2,3),rnodes(2,1000)
c
        integer omp_get_max_threads
c
c       This subroutine intializes a discretization structure.  That 
c       is, it constructs the structure header and stores the
c       necessary quadrature data in the structure.  
c
c                            Input Parameters:
c
c   norder - order of polynomials used to represent solutions;
c       see discquadfor a table specifying the allowable orders
c
c   nself - order for the self-interaction quadratures; see selfquad.f
c       for a table specifying the allowable orders
c           
c       NOTE: if nself is set to 0, then self-interactions will be
c       evaluated using a massively oversampled quadrature formula
c
c   dnear - 
c
c
c   ldisc - the length of the user-supplied array disc
c
c                            Output Parameters:
c
c   ier - an error return code;
c       ier = 0     indicates successful execution
c       ier = 2     means that the disc array is of woefully insufficient 
c                   length
c       ier = 4     means that the disc array is of insufficient length
c       ier = 32    means that the parameter norder was invalid
c       ier = 128   means that the parameter nself was invalid
c       ier = 1024  means that an error was encountered while 
c                   reading the table of self-interaction quadrature
c                   rules from the disk
c
c   disc - upon return this user-supplied array will contian a properly
c       initialized discretization structure formatted as described below
c
c   
c                           Discretization structure:
c
c   The header of the discretization structure is formatted as follows:
c
c   disc(1)  - the order of the polynomials used to represent solutions   (norder)
c   disc(2)  - the order of the self-interaction quadratures              (nself)
c   disc(3)  - the maximum number of triangles which can be accomodate    (maxtris)
c              by the tris array
c   disc(4)  - the current number of triangles in the tris array          (ntris)
c   disc(5)  - the length of each entry in the tri array                  (len)
c   disc(6)  - a pointer to the tris array                                (itris)
c   disc(7)  - the maximum size for the quadratures used to evaluate      (maxquad
c              self interactions
c   disc(8)  - the length of the work arrays needed by the evaluation     (lw)
c              subroutines (this is the length of the work array per
c              thread)
c   disc(9)  - a pointer to a work array of length lw                     (iw)
c   disc(10) - machine epsilon                                            (eps0)
c
c
c   disc(20) - the number of nodes in the discretization quadrature       (nquad)
c   disc(21) - a pointer to the x-coords of the discretization quadrature (ixs)
c              nodes
c   disc(22) - a pointer to the y-coords of the discretization quadrature (iys)
c              nodes
c   disc(23) - a pointer to the weights of the discretization quadrature  (iwhts)
c              nodes
c   disc(24) - a pointer  to the matrix which takes scaled values at the  (iu)
c              discretization nodes to coefficient expansions
c   disc(25) - a pointer to the inverse of the matrix u                   (iv)
c
c   disc(30) - the number of nodes in the quadrature used for near        (nquadnear)
c              interactions 
c   disc(31) - a pointer to an array containing the x-coordinates of the  (ixsnear)
c              near quadrature nodes
c   disc(32) - a pointer to an array containing the y-coordinates of the  (iysnear)
c              near quadrature nodes
c   disc(33) - a pointer to an array containing the near quadrature       (iwhtsnear)
c              weights
c   disc(34) - a pointer to the (4*nquad,nquad) interpolation matrix      (iamatrin)
c              used in the computation of near interactions
c   disc(40) - a pointer to the rad structure which contains the tables   (irad)
c              used to construct self interaction quadratures
c
c
c
c   The tris array is a (len,ntris) matrix which is formatted as a structure
c   with ntris entries.  Each entries described one surface element:
c
c   tris(1,j)  - the x-coordinates of the first vertex of the triangle
c   tris(2,j)  - the y-coordinates of the first vertex of the triangle
c   tris(3,j)  - the x-coordinates of the second vertex of the triangle
c   tris(4,j)  - the y-coordinates of the second vertex of the triangle
c   tris(5,j)  - the x-coordinates of the third vertex of the triangle
c   tris(6,j)  - the y-coordinates of the third vertex of the triangle
c   tris(7,j)  - the index of the parameterization 
c   tris(8,j)  - x-coordinate of the center of the associated bounding    (bx)
c                ball 
c   tris(9,j)  - y-coordinate of the center of the associated bounding    (by)
c                ball
c   tris(10,j) - z-coordinate of the center of the associated bounding    (bz)
c                ball
c   tris(11,j) - the square of the radius of the bounding ball for the    (br)
c                surface element
c   tris(12,j) - first entry of the evaluation data for the discretization
c                nodes
c
c   Evaluation data for the discretization nodes on the surface element are
c   stored
c
c    <<<<<<<<<DESCRIBE EVALUATION DATA>>>>>>>>>>>
c
        ier = 0
c
        call mach_zero(eps0)
        maxthreads = omp_get_max_threads()
c
        call prin2("in disc_init, eps0 = *",eps0,1)
        call prinf("in disc_init, maxthreads = *",maxthreads,1)
c
c       Perform some basic sanity checks.
c
        if (ldisc .lt. 100 000) then
        ier = 2
        return
        endif
c
c       Allocate memory for the quadratures.
c
        nquad     = (norder+1)*(norder+2)/2
        maxquad   = 100 000
c
        ixs       = 100
        lxs       = nquad
c
        iys       = ixs+lxs
        lys       = nquad
c
        iwhts     = iys+lys
        lwhts     = nquad
c
        iu        = iwhts+lwhts
        lu        = nquad*nquad
c
        iv        = iu+lu
        lv        = nquad*nquad
c
        iamatrin  = iv+lv
        lamatrin  = 4*nquad*nquad
c
        iw        = iamatrin+lamatrin
        lw        = 3*nquad+nquad**2+5*maxquad
c        
        irad      = iw+lw*maxthreads
        lrad      = ldisc-irad
c
        if (lrad .le. 0) then 
        ier = 4
        return
        endif
c
c       Fetch the discretization quadrature.
c
        call discquad(jer,norder,disc(ixs),disc(iys),disc(iwhts))
        if (jer .ne. 0) then
        ier = 32
        return
        endif
c
c       Construct the matrices u and v.
c
        if (lrad .lt. nquad*(nquad+1)+10 000) then
        ier = 4
        return
        endif
c
        call disc_init0(norder,nquad,disc(ixs),disc(iys),disc(iwhts),
     -    disc(iu),disc(iv),disc(irad))
c
c       Construct amatrin.
c
        call disc_amatrin(norder,nquad,disc(ixs),disc(iys),disc(iwhts),
     -    disc(iu),disc(iamatrin))
c
c       Initialize the self-interaction quadratures.
c
        if (nself .ne. 0) then
        call selfquad_init(ier,nself,disc(irad),lrad,lkeep)
        if (ier .ne. 0) return
        else
        lkeep = 0
        endif
c
        lrad = lkeep
c
        itris   = irad+lrad
        ltris   = ldisc-itris
c
        len     = 11+13*nquad
        maxtris = ltris/len
        ntris   = 0
c
        call prinf("in disc_init, norder = *",norder,1)
        call prinf("in disc_init, nquad = *",nquad,1)
        call prinf("in disc_init, len = *",len,1)
        call prinf("in disc_init, maxtris = *",maxtris,1)
c
        disc(1)  = norder       
        disc(2)  = nself
        disc(3)  = maxtris
        disc(4)  = ntris
        disc(5)  = len
        disc(6)  = itris
        disc(7)  = maxquad
        disc(8)  = lw
        disc(9)  = iw
        disc(10) = eps0
c
        disc(20) = nquad
        disc(21) = ixs
        disc(22) = iys
        disc(23) = iwhts
        disc(24) = iu
        disc(25) = iv
c
        disc(34) = iamatrin
c
        disc(40) = irad
c
        disc(50) = 0
        disc(51) = 0
c        disc(52) = 0
c        disc(53) = 0
c     
        disc(99)=dnear
        end



        subroutine disc_amatrin(norder,nquad,xs,ys,whts,u,amatrin)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1),u(nquad,nquad)
        dimension amatrin(4*nquad,nquad)
        dimension simp(6)
c
        double precision, allocatable :: vals(:,:)
c
        dimension vals0(nquad)
        dimension tris(6,4),amap(2,3)
c
        dimension xs0(4*nquad),ys0(4*nquad),whts0(4*nquad)
        dimension vals1(4*nquad),vals2(nquad)
c
c       Coordinates of the subtriangles.
c
        simp(1) = 0.0d0
        simp(2) = 0.0d0
        simp(3) = 1.0d0
        simp(4) = 0.d00
        simp(5) = 0.0d0
        simp(6) = 1.0d0
c
        call disc_split_tri(simp,tris)
c
        allocate(vals(4*nquad,nquad))
c
c       Construct the matrix of values of the koornwinder polynomials
c       at the discretization nodes on the subtriangles.
c
        do 1000 itri=1,4
        call disc_trimap(tris(1,itri),amap,det)
c
        do 1100 j=1,nquad
        x   = amap(1,1)*xs(j)+amap(1,2)*ys(j)+amap(1,3)
        y   = amap(2,1)*xs(j)+amap(2,2)*ys(j)+amap(2,3)
        wht = whts(j)*det
c
        xs0(nquad*(itri-1)+j)   = x
        ys0(nquad*(itri-1)+j)   = y
        whts0(nquad*(itri-1)+j) = wht
c
        call koorn(norder,x,y,vals0)
        idx = (itri-1)*nquad+j
        do 1200 i=1,nquad
        vals(idx,i) = vals0(i)*sqrt(wht)
 1200 continue
 1100 continue
 1000 continue
c
c       Multiply on the right by u in order to construct amatrin.
c
        do 2000 i=1,4*nquad
        do 2100 j=1,nquad
        sum = 0
        do 2200 l=1,nquad
        sum = sum + vals(i,l)*u(l,j)
 2200 continue
        amatrin(i,j) = sum
 2100 continue
 2000 continue
c
c$$$c
c$$$c       Test it.
c$$$c
c$$$        do 3000 i=1,nquad
c$$$        x   = xs(i)
c$$$        y   = ys(i)
c$$$        wht = whts(i)
c$$$c
c$$$        vals2(i) = (x**2*y+1)*sqrt(wht)
c$$$ 3000 continue
c$$$c
c$$$        do 3100 i=1,4*nquad
c$$$        sum = 0
c$$$        do 3200 j=1,nquad
c$$$        sum = sum + amatrin(i,j)*vals2(j)
c$$$ 3200 continue
c$$$        vals1(i) = sum
c$$$c
c$$$        x   = xs0(i)
c$$$        y   = ys0(i)
c$$$        wht = whts0(i)
c$$$        print *,sum - (x**2*y+1)*sqrt(wht)
c$$$ 3100 continue
c$$$c
c$$$        stop
c
        end


        subroutine disc_init0(norder,nquad,xs,ys,whts,u,v,w)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1),u(nquad,nquad)
        dimension v(nquad,nquad),w(1)
        dimension vals0(nquad)
c
c       Explicity reconstruct u and v.
c
        do 0200 i=1,nquad
c
        x   = xs(i)
        y   = ys(i)
        wht = whts(i)
c
        call koorn(norder,x,y,vals0)
c
        do 0300 j=1,nquad
        v(i,j) = vals0(j)*sqrt(wht)
        u(i,j) = vals0(j)*sqrt(wht)
 0300 continue
 0200 continue
c
        call orthom(u,nquad,w,cond)
c
c        call print_sing("in disc_init, u sings = *",nquad,nquad,u)
c
        end



        subroutine disc_add(ier,disc,iparam,ntris0,tris0)
        implicit double precision (a-h,o-z)
        dimension disc(1),trisnew(6,1)
c
c       Add a collection of surface elements to the decomposition
c       contained in the disc structure.
c
c                           Input Parameters:
c
c   disc - the discretization structure 
c   iparam - the index of the parameterization
c   ntrisnew - the number of triangles to add the the decomposition
c   trisnew - a list of the new triangles
c
c                          Output Parameters:
c 
c   ier - an error return code
c
        ier = 0
c
c       Fetch data from the disc structure.
c
        norder     = disc(1)
        nself      = disc(2)
        maxtris    = disc(3)
        ntris      = disc(4)
        len        = disc(5)
        itris      = disc(6)
        maxquad    = disc(7)
        lw         = disc(8)
        iw         = disc(9)
c
        nquad      = disc(20)
        ixs        = disc(21)
        iys        = disc(22)
        iwhts      = disc(23)
        iu         = disc(24)
        iv         = disc(25)
c
        irad       = disc(40)
c
c       Check that the maximum number of triangles is not exceeded.
c
        if (ntris+ntris0 .gt. maxtris) then
        ier = 4
        return
        endif
c
c       Call an auxillary routine to shape arrays.
c
        call disc_add0(nquad,disc(ixs),disc(iys),disc(iwhts),ntris,
     -    len,disc(itris),iparam,ntris0,tris0,disc(iw))
c
c       Update the header.
c
        disc(4) = ntris
c
        end

       
        subroutine disc_add0(nquad,xs,ys,whts,ntris,len,tris,
     -    iparam,ntris0,tris0,zs)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1),tris(len,1),dn(3)
        dimension tris0(6,1),r(3),dr(3,2),amap(2,3),zs(3,1)
c
        do 1000 j=1,ntris0
c
        itri = ntris+j
c
c       Store the vertices of the parameterization domain.
c
        tris(1,itri) = tris0(1,j)
        tris(2,itri) = tris0(2,j)
        tris(3,itri) = tris0(3,j)
        tris(4,itri) = tris0(4,j)
        tris(5,itri) = tris0(5,j)
        tris(6,itri) = tris0(6,j)
        tris(7,itri) = iparam
c
c       Construct the mapping taking the simplex to the triangle
c       which constitutes the parameterization domain.
c
        call disc_trimap(tris(1,itri),amap,det)
c
c       Construct the evaluation data for the surface element.
c
        iptr = 12
c
        do 1100 i=1,nquad
        s   = amap(1,1)*xs(i)+amap(1,2)*ys(i)+amap(1,3)
        t   = amap(2,1)*xs(i)+amap(2,2)*ys(i)+amap(2,3)
        wht = whts(i)*det
c
        call surface(iparam,s,t,r,dr)
c
        zs(1,i) = r(1)
        zs(2,i) = r(2)
        zs(3,i) = r(3)
c
c       Construct the normal vector.
c
        dn(1) = dr(2,1)*dr(3,2) - dr(3,1)*dr(2,2)
        dn(2) = dr(3,1)*dr(1,2) - dr(1,1)*dr(3,2)
        dn(3) = dr(1,1)*dr(2,2) - dr(2,1)*dr(1,2)
        da    = sqrt(dn(1)**2+dn(2)**2+dn(3)**2)
        dn(1) = dn(1)/da
        dn(2) = dn(2)/da
        dn(3) = dn(3)/da
c
        wht = wht * da
c
        tris(iptr,    itri)  = sqrt(wht)
        tris(iptr+1,  itri)  = r(1)
        tris(iptr+2,  itri)  = r(2)
        tris(iptr+3,  itri)  = r(3)
        tris(iptr+4,  itri)  = dr(1,1)
        tris(iptr+5,  itri)  = dr(2,1)
        tris(iptr+6,  itri)  = dr(3,1)
        tris(iptr+7,  itri)  = dr(1,2)
        tris(iptr+8,  itri)  = dr(2,2)
        tris(iptr+9,  itri)  = dr(3,2)
        tris(iptr+10, itri)  = dn(1)
        tris(iptr+11, itri)  = dn(2)
        tris(iptr+12, itri)  = dn(3)
c
        iptr=iptr+13
 1100 continue
c
c       Construct a bounding sphere for the triangle.  NOTE: br is 
c       the SQUARE of the radius of the ball.
c
        call bounding_ball(nquad,zs,bx,by,bz,br)
c
        tris(8, itri) = bx
        tris(9, itri) = by
        tris(10,itri) = bz
        tris(11,itri) = br
c
 1000 continue
        ntris = ntris+ntris0
c
c
        end


        subroutine disc_info(disc,ntris,nquad,n,lkeep)
        implicit double precision (a-h,o-z)
        dimension disc(1)
c
c       Return some rudimentry information from the disc structure.
c
c                          Input Parameters:
c
c   disc - the discretization structure to interogate
c
c                          Output Parameters;
c 
c   ntris - the number of surface elements or "triangles" 
c   nquad - the number of discretization nodes per surfac element
c   n - the total number of discretization odes
c   lkeep - the length of the data stored in the disc array
c
c
c       Fetch data from the disc structure.
c
        norder     = disc(1)
        nself      = disc(2)
        maxtris    = disc(3)
        ntris      = disc(4)
        len        = disc(5)
        itris      = disc(6)
        maxquad    = disc(7)
        lw         = disc(8)
        iw         = disc(9)
c
        nquad      = disc(20)
        ixs        = disc(21)
        iys        = disc(22)
        iwhts      = disc(23)
        iu         = disc(24)
        iv         = disc(25)
c
        nquadnear  = disc(30)
        ixsnear    = disc(31)
        iysnear    = disc(32)
        iwhtsnear  = disc(33)
c
        irad       = disc(40)
c
c       Return some of it.
c
        n     = nquad*ntris
        lkeep = itris+ntris*len
c
        end


        subroutine disc_data(disc,eval)
        implicit double precision (a-h,o-z)
        dimension disc(1),eval(13,1)
c
c       Return evaluation data for the discretization nodes of the 
c       decomposition.
c
c
c                       Input Parameters:
c
c   disc - the discretization structure
c
c                       Output Parameters:
c
c   eval - a (13,n) 
c        

c
c       Fetch data from the disc structure.
c
        norder     = disc(1)
        nself      = disc(2)
        maxtris    = disc(3)
        ntris      = disc(4)
        len        = disc(5)
        itris      = disc(6)
        maxquad    = disc(7)
        lw         = disc(8)
        iw         = disc(9)
c
        nquad      = disc(20)
        ixs        = disc(21)
        iys        = disc(22)
        iwhts      = disc(23)
        iu         = disc(24)
        iv         = disc(25)
c
        nquadnear  = disc(30)
        ixsnear    = disc(31)
        iysnear    = disc(32)
        iwhtsnear  = disc(33)
c
        irad       = disc(40)
c
c
c
        n     = nquad*ntris
        lkeep = itris+ntris*len
c
c       Call an auxillary subroutine to shape arrays.
c
        call disc_data0(nquad,len,ntris,disc(itris),eval)
        
        end



        subroutine disc_data0(nquad,len,ntris,tris,eval)
        implicit double precision (a-h,o-z)
        dimension eval(13,1),tris(len,1)
c
c       Return evaluation data for the discretization nodes of the 
c       decomposition.
c
        ipt = 1
c
        do 1000 itri=1,ntris
c
        iptr = 12
c
        do 1100 i=1,nquad
c
        eval(1, ipt) = tris(iptr,   itri)
        eval(2, ipt) = tris(iptr+1, itri)
        eval(3, ipt) = tris(iptr+2, itri)
        eval(4, ipt) = tris(iptr+3, itri)
        eval(5, ipt) = tris(iptr+4, itri)
        eval(6, ipt) = tris(iptr+5, itri)
        eval(7, ipt) = tris(iptr+6, itri)
        eval(8, ipt) = tris(iptr+7, itri)
        eval(9, ipt) = tris(iptr+8, itri)
        eval(10,ipt) = tris(iptr+9, itri)
        eval(11,ipt) = tris(iptr+10,itri)
        eval(12,ipt) = tris(iptr+11,itri)
        eval(13,ipt) = tris(iptr+12,itri)
c
        ipt  = ipt+1
        iptr = iptr+13
 1100 continue
c
 1000 continue
c
        end


        subroutine disc_points(disc,zs)
        implicit double precision (a-h,o-z)
        dimension disc(1),zs(3,1)
c
c       Return evaluation data for the discretization nodes of the 
c       decomposition.
c
c
c                       Input Parameters:
c
c   disc - the discretization structure
c
c                       Output Parameters:
c
c   zs - a (3,n) array specifying the coordinates of the discretization
c       nodes on the surface

c
c       Fetch data from the disc structure.
c
        norder   = disc(1)
        nself    = disc(2)
        maxtris  = disc(3)
        ntris    = disc(4)
        len      = disc(5)
        itris    = disc(6)
        maxquad  = disc(7)
        lw       = disc(8)
        iw       = disc(9)
        eps0     = disc(10)
c
        nquad    = disc(20)
        ixs      = disc(21)
        iys      = disc(22)
        iwhts    = disc(23)
        iu       = disc(24)
        iv       = disc(25)
c
        nquadnear  = disc(30)
        ixsnear    = disc(31)
        iysnear    = disc(32)
        iwhtsnear  = disc(33)
c
        irad     = disc(40)
c
        n     = nquad*ntris
        lkeep = itris+ntris*len
c
c       Call an auxillary subroutine to shape arrays.
c
        call disc_points0(nquad,len,ntris,disc(itris),zs)
c        
        end



        subroutine disc_points0(nquad,len,ntris,tris,zs)
        implicit double precision (a-h,o-z)
        dimension zs(3,1),tris(len,1)
c
c       Return evaluation data for the discretization nodes of the 
c       decomposition.
c
        ipt = 1
c
        do 1000 itri=1,ntris
c
        iptr = 12
c
        do 1100 i=1,nquad
c
        zs(1,ipt)= tris(iptr+1, itri)
        zs(2,ipt)= tris(iptr+2, itri)
        zs(3,ipt)= tris(iptr+3, itri)
c
        ipt  = ipt+1
        iptr = iptr+13
 1100 continue
c
 1000 continue
c
        end



        subroutine disc_trimap(verts,amap,det)
        implicit double precision (a-h,o-z)
        dimension verts(2,3),amap(2,3)
c
c       Return the affine mapping (amap) taking the simplex to a user-
c       specified triangle.
c
        amap(1,1) = verts(1,2) - verts(1,1)
        amap(2,1) = verts(2,2) - verts(2,1)
        amap(1,2) = verts(1,3) - verts(1,1)
        amap(2,2) = verts(2,3) - verts(2,1)
        amap(1,3) = verts(1,1)
        amap(2,3) = verts(2,1)
c
        det = abs(amap(1,1)*amap(2,2)-amap(1,2)*amap(2,1))
c
        end


        subroutine disc_trimap2(verts,amap,ainv,det)
        implicit double precision (a-h,o-z)
        dimension verts(2,3),amap(2,3),ainv(2,3)
c
c       Return the affine mapping (amap) taking the simplex to a user-
c       specified triangle and its inverse.
c
        amap(1,1) = verts(1,2) - verts(1,1)
        amap(2,1) = verts(2,2) - verts(2,1)
        amap(1,2) = verts(1,3) - verts(1,1)
        amap(2,2) = verts(2,3) - verts(2,1)
        amap(1,3) = verts(1,1)
        amap(2,3) = verts(2,1)
c
        det = amap(1,1)*amap(2,2)-amap(1,2)*amap(2,1)
c
        ainv(1,1) = amap(2,2)/det
        ainv(2,2) = amap(1,1)/det
c
        ainv(1,2) =-amap(1,2)/det
        ainv(2,1) =-amap(2,1)/det
c
        ainv(1,3) = - ( ainv(1,1)*verts(1,1) + ainv(1,2)*verts(2,1))
        ainv(2,3) = - ( ainv(2,1)*verts(1,1) + ainv(2,2)*verts(2,1))
c
        det = abs(det)
c
        end


        subroutine disc_split_tri(tri,children)
        implicit double precision (a-h,o-z)
        dimension tri(2,3),children(2,3,4)
c
        x1 = tri(1,1)
        y1 = tri(2,1)
        x2 = tri(1,2)
        y2 = tri(2,2)
        x3 = tri(1,3)
        y3 = tri(2,3)
c
        x4 = (x1+x2)/2
        y4 = (y1+y2)/2
c
        x5 = (x1+x3)/2
        y5 = (y1+y3)/2
c
        x6 = (x2+x3)/2
        y6 = (y2+y3)/2
c
        children(1,1,1) = x1
        children(2,1,1) = y1
        children(1,2,1) = x4
        children(2,2,1) = y4
        children(1,3,1) = x5
        children(2,3,1) = y5
c
        children(1,1,2) = x6
        children(2,1,2) = y6
        children(1,2,2) = x4
        children(2,2,2) = y4
        children(1,3,2) = x5
        children(2,3,2) = y5
c
        children(1,1,3) = x3
        children(2,1,3) = y3
        children(1,2,3) = x5
        children(2,2,3) = y5
        children(1,3,3) = x6
        children(2,3,3) = y6
c
        children(1,1,4) = x4
        children(2,1,4) = y4
        children(1,2,4) = x2
        children(2,2,4) = y2
        children(1,3,4) = x6
        children(2,3,4) = y6
c        
        end



        subroutine disc_move(k,a,b)
        implicit double precision (a-h,o-z)
        dimension a(1),b(1)
        do 1000 i=1,k
        b(i)=a(i)
 1000 continue
        end



        subroutine disc_apply(n,m,a,x,y)
        implicit double precision (a-h,o-z)
        double precision a(n,m)
        double complex x(n),y(m),sum
c
c       Compute y(m) = x(n) a(n,m)
c
        do 1000 i=1,m
        sum = 0
        do 1100 j=1,n
        sum = sum + x(j)*a(j,i)
 1100 continue
        y(i) = sum
 1000 continue
c
        end




        subroutine koorn(n,x,y,pols)
        implicit double precision (a-h,o-z)
        dimension pols(1)
        dimension plege(n+1),pjac(n+1)
        data sq32 / 5.65685424949238019520675489683879231d0/
c
c       Evaluate the orthonormal polynomials of degree 0 through n
c       at the point (x,y) using the Koornwinder representation.  There 
c       are (n+1)(n+2)/2 such polynomials.
c
c       Note that this routine normalizes the polynomials so that their
c       L^2 norms over the simplex are 1.
c
c       The n+1 polynomials of degree n are given by 
c
c             (2k+1,0)         (0,0)
c         \phi (2x-1)  *   \phi (2y/(1-x)-1) * (1-x)^k 
c             n-k              k
c
c       k=0,...,n, where 
c
c             (a,b)
c         \phi (x)
c             k
c
c       denotes the Jacobi polynomial of degree k with parameters 
c       \alpha=a and \beta=b.       
c
c       The values of the polynomials are sorted first by degree and
c       then by the value of k.        
c
c       
c                          Input Parameters:
c
c   n - maximum order of the polynomials to evaluate
c   (x,y) - point at which to evaluate the polynomials
c
c                         Output Parameters:
c
c   vals - the (n+1)*(n+2)/2 values of the polynomials
c             
c     
c       First, evaluate the scaled Legendre polynomials.
c
        a=0
        b=0
c
        z = 2*y/(1-x)-1
        call legen(n,z,plege)
c
        alpha=sq32
        do 1000 j=0,n
        plege(j+1)=plege(j+1)*alpha
        alpha=alpha*(1-x)
 1000 continue     
c
        z = 2*x-1
c
        do 1100 k=0,n
c
c       Evaluate the (other) Jacobi polynomials
c
        a = 2*k+1
        call jacobi(n-k,a,z,pjac)
c
        do 1200 l=0,n-k
        m = l+k
        idx = ((m+1)*(m))/2+1+k
c
        pols(idx)=plege(k+1)*pjac(l+1)
c
 1200 continue
 1100 continue
c
        end
c
c
c
        subroutine koorn2(n,x,y,pols,dersx,dersy)
        implicit double precision (a-h,o-z)
        dimension pols(1),dersx(1),dersy(1)
        dimension plege(n+1),pjac(n+1),dlege(n+1),djac(n+1)       
c
        data sq32 / 5.65685424949238019520675489683879231d0/
c
c       Evaluate the orthonormal polynomials of degree 0 through n
c       and their derivatives at the point (x,y) using the Koornwinder 
c       representation.  There are (n+1)(n+2)/2 such polynomials.
c
c       Note that this routine normalizes the polynomials so that their
c       L^2 norms over the simplex are 1.
c
c       The values of the polynomials are sorted first by degree and
c       then by the value of k.        
c
c                          Input Parameters:
c
c   n - maximum order of the polynomials to evaluate
c   (x,y) - point at which to evaluate the polynomials
c
c                         Output Parameters:
c
c   vals - the (n+1)*(n+2)/2 values of the polynomials at (x,y)
c   dersx - the derivatives of said polynomials w.r.t. x at the (x,y)
c   dersy - the derivatives of said polynomials w.r.t. y at the (x,y) 
c
c       Evaluate the scaled Legendre polynomials.
c
        z1    = 2*y/(1-x)-1
        call legedersn(n,z1,plege,dlege)
c
        z2 = 2*x-1
c
        nn=0
c
        do 1100 k=0,n
c
c       Evaluate the (other) Jacobi polynomials
c
        a = 2*k+1
        call jacobiders(n-k,a,z2,pjac,djac)
c
        d = sq32*(1-x)**k
c
        do 1200 l=0,n-k
c
        m = l+k
        idx = (m+1)*(m)/2+1+k
c
        pols(idx)=pjac(l+1)*plege(k+1)*d
c
        dd1=djac(l+1)*2*plege(k+1)*d
        dd2=pjac(l+1)*dlege(k+1)*2*y/(1-x)**2*d
c
        if (k .eq. 0) then
        dd3=0
        else
        dd3=-sq32*k*pjac(l+1)*plege(k+1)*(1-x)**(k-1)
        endif
c
        dersx(idx)=(dd1+dd2+dd3)
        dersy(idx)=pjac(l+1)*d*dlege(k+1)*2/(1-x)
c        
 1200 continue
 1100 continue
c
c
        end
c
c
c
        subroutine legen(n,x,vals)
        implicit double precision (a-h,o-z)
        dimension vals(1)
c     
c       Evaluate the normalized Legendre polynomials of degree 0 
c       through n at the point x.  
c
c       The polynomials are normalized so as to make their L^2[-1,1]
c       norms 1.
c
        vals(1) = 1
        if (n .eq. 0) goto 2000
c
        vals(2) = x
        if (n .eq. 1) goto 2000
c
        do 1000 j=2,n
        vals(j+1)=((2*j-1)*x*vals(j)-(j-1)*vals(j-1))/j
 1000 continue
c
c
c       Normalize the polynomials.
c
 2000 continue
c
        do 2100 j=0,n
        vals(j+1)=vals(j+1)*sqrt((2*j+1.0d0)/2.0d0)
 2100 continue
        end
c
c
c
        subroutine legedersn(n,x,pols,ders)
        implicit double precision (a-h,o-z)
        dimension pols(1),ders(1)
c     
c       Evaluate the normalized Legendre polynomials of degree 0 
c       through n and their derivatives at the point x.
c
c       The polynomials are normalized so as to make their L^2[-1,1]
c       norms 1.
c
        pols(1) = 1
        ders(1) = 0
        if (n .eq. 0) goto 2000
c
        pols(2) = x
        ders(2) = 1
        if (n .eq. 1) goto 2000        
c
        do 1000 j=2,n
        pols(j+1)=((2*j-1)*x*pols(j)-(j-1)*pols(j-1))/j
 1000 continue
 2000 continue
c
c       Compute the derivatives.
c
        d=x**2-1
        do 3000 j=3,n+1
        ders(j)=(j-1)*(x*pols(j)-pols(j-1))/d
 3000 continue
c
c       Normalize the polynomials and derivatives.
c
        do 2100 j=0,n
        d=sqrt( (2.0d0*j+1.0d0)/2.0d0 ) 
        pols(j+1)=pols(j+1)*d
        ders(j+1)=ders(j+1)*d
 2100 continue
c
        end
c
c
c
        subroutine jacobi(n,a,x,pols)
        implicit double precision (a-h,o-z)
        dimension pols(1)
c     
c       Evaluate the normalized Jacobi polynomials (with parameters
c       \alpha = a and \beta = 0) of degrees 0 through n at the point
c       x.
c
c       The polynomials are normalized so as to make their L^2 norms
c       with respect to the measure (1-x)^a dx equal to 1.
c
        pols(1) = 1
        if (n .eq. 0) goto 2000
c
        pols(2) = (a+(2+a)*x)/2
        if (n .eq. 1) goto 2000
c
        do 1000 l=2,n
c
        a1=2*(l)*(l+a)*(2*l-2+a)
        a2=(2*l-1+a)*a**2
        a3=(2*l-2+a)*(2*l+a-1)*(2*l+a)
        a4=2*(l+a-1)*(l-1)*(2*l+a)
c
        pols(l+1)=((a2+x*a3)*pols(l)-a4*pols(l-1))/a1
 1000 continue
c
c       Normalize the polynomials.
c
 2000 continue
c
        alpha=0
        d = 2**(alpha+2)
        do 2100 j=0,n
        pols(j+1)=pols(j+1)*sqrt(2*j+a+1.0d0)/(d)
 2100 continue
        end
c
c
c
        subroutine jacobiders(n,a,x,pols,ders)
        implicit double precision (a-h,o-z)
        dimension pols(1),ders(1)
c     
c       Evaluate the normalized Jacobi polynomials (with parameters
c       \alpha = a and \beta = 0) of degrees 0 through n at the point
c       x and their derivatives.
c
c       The polynomials are normalized so as to make their L^2 norms
c       with respect to the measure (1-x)^a dx equal to 1.
c
        pols(1) = 1
        ders(1) = 0
        if (n .eq. 0) goto 2000
c
        pols(2) = (a+(2+a)*x)/2
        ders(2) = 1+a/2
        if (n .eq. 1) goto 2000
c
        do 1000 l=2,n
c
        a1=2*(l)*(l+a)*(2*l-2+a)
        a2=(2*l-1+a)*a**2
        a3=(2*l-2+a)*(2*l+a-1)*(2*l+a)
        a4=2*(l+a-1)*(l-1)*(2*l+a)
c
        pols(l+1)=((a2+x*a3)*pols(l)-a4*pols(l-1))/a1
        ders(l+1)=((a2+x*a3)*ders(l)+a3*pols(l)-a4*ders(l-1))/a1
 1000 continue
c
c       Normalize the polynomials.
c
 2000 continue
c
c        call prin2("pols = *",pols,n)
c
        alpha=0
        d = 2.0d0**(alpha+2)
c
        do 2100 j=0,n
        pols(j+1)=pols(j+1)*sqrt(2.0d0*j+a+1.0d0)/(d)
        ders(j+1)=ders(j+1)*sqrt(2.0d0*j+a+1.0d0)/(d)
 2100 continue
c
c

        end


        subroutine discquad(ier,norder,xsdisc,ysdisc,whtsdisc)
        implicit double precision (a-h,o-z)
        dimension xsdisc(1),ysdisc(1),whtsdisc(1)
c
c       Return one of the discretization quadratures for integrating
c       polynomials over the standard simplex.  The returned quadrature
c       will be of length
c
c               (norder+1)*(norder+2)
c               --------------------- ,
c                         2
c
c       which is the dimension of the space of polynomials of order 
c       norder, but it will accurately integrate polynomials of
c       somewhat higher order.
c
c
c       The following tables gives the valid values of norder,
c       the size of the resulting quadrature and the order of
c       polynomials it integrates.
c
c
c                 disc order    nodes     poly order
c                 ----------    -----     ----------
c                    4            15         7
c                    8            45        14
c                   12            91        21
c                   16           153        28 (?)
c
c       All quadrature formulas returned by this routine are accurate
c       to around 30 digits, have positive weights and interior nodes.
c
c
c                            Input Parameters:
c
c   norder - the order of the quadrature to return
c
c                            Output Parameters:
c
c   ier - an error return code;
c      ier = 0     indicates successful execution
c      ier = 1024  means the parameter norder was invalid
c
c   xsdisc - this user-supplied array will contain the x-coordinates of
c      the requested quadrature
c   ysdisc - this user-supplied array will contain the y-coordinates of
c      the requested quadrature
c   whtsdisc - this user-supplied array will contain the weights of
c      the quadrature rule
c
c       4th order
c
        parameter ( nquad4  =  15 , m4 =   7)
        dimension xs4( 15)
        dimension ys4( 15)
        dimension whts4( 15)

cccccccccccccccccccccccc GENERALIZED GAUSSIAN
        data xs4 / 
     -    0.649324151964922343073344376958858851D-01,
     -    0.198394775029309042060805502983842700D+00,
     -    0.517025481808406434532985798726530320D+00,
     -    0.517034323752486790056943273087229677D+00,
     -    0.284575348788002070011479595459183315D+00,
     -    0.284564382149961645361253651217656216D+00,
     -    0.198358239433028602067512937265594915D+00,
     -    0.438545449708890472439860188814423925D-01,
     -    0.313566869370115942992699092079369253D+00,
     -    0.642586774478727964581729963445351702D+00,
     -    0.438694718033534620388981591892551633D-01,
     -    0.313556583778533638153154691394046038D+00,
     -    0.642564648250388813401576637675356063D+00,
     -    0.649296185707683681457988043272763069D-01,
     -    0.870139662314929268514296796457279655D+00 /
        data ys4 / 
     -    0.870138296511631838009526464100738549D+00,
     -    0.284575999635906221842302444861830831D+00,
     -    0.284586492306177274277580877669790364D+00,
     -    0.198411019704622710308592095978531310D+00,
     -    0.198373684968258137616778558821978919D+00,
     -    0.517054470688654194849719812369102657D+00,
     -    0.517045161142257044625505226551689201D+00,
     -    0.642590048910570913317540800617133017D+00,
     -    0.642567849579552102576519210999789459D+00,
     -    0.313551514311710042040074310273103313D+00,
     -    0.313561847386493219155553912022242630D+00,
     -    0.438573535388267305547591591385213403D-01,
     -    0.438724674370882008136126959943077118D-01,
     -    0.649314190420260343309686495981383691D-01,
     -    0.649286155350458270104800877431443352D-01 /
        data whts4 /
     -    0.265391657928760840005171139929677959D-01,
     -    0.354151357454423214493437682872669937D-01,
     -    0.354299884607202743291229982574722383D-01,
     -    0.354199814358862819650116903193868595D-01,
     -    0.354380021884711887864384270725737352D-01,
     -    0.354231285309258776874741345392141622D-01,
     -    0.354330134427217266597175429681126316D-01,
     -    0.346309653500952028276659554051694402D-01,
     -    0.346387569539455452284168923120143896D-01,
     -    0.346359579409124104555270871893193739D-01,
     -    0.346417498655874444309181471739052989D-01,
     -    0.346328496152912078179703879921096685D-01,
     -    0.346437678836335400244616422862426219D-01,
     -    0.265389052844368669995383560524879865D-01,
     -    0.265386315090540273378758561517574180D-01 /


cccccccccccccccccccccccccc VR
c$$$          data xs4/
c$$$     -      0.242178877429365427415842317545529884D+00,
c$$$     -      0.515642245141269145168315364908936909D+00,
c$$$     -      0.242178877429365427415842317545530029D+00,
c$$$     -      0.474101304688892224560294156365469921D+00,
c$$$     -      0.517973906222155508794116872690563784D-01,
c$$$     -      0.474101304688892224560294156365470258D+00,
c$$$     -      0.233916149991311356180523963867127575D+00,
c$$$     -      0.717926712780444024794562208158499202D+00,
c$$$     -      0.481571372282446190249138279743716583D-01,
c$$$     -      0.717926712780444024794562208158501705D+00,
c$$$     -      0.233916149991311356180523963867129212D+00,
c$$$     -      0.481571372282446190249138279743708638D-01,
c$$$     -      0.514351415804550152570918459875889853D-01,
c$$$     -      0.897129716839089969485816308024824437D+00,
c$$$     -      0.514351415804550152570918459875865478D-01/
c$$$          data ys4/
c$$$     -      0.515642245141269145168315364908936813D+00,
c$$$     -      0.242178877429365427415842317545529812D+00,
c$$$     -      0.242178877429365427415842317545530005D+00,
c$$$     -      0.474101304688892224560294156365470114D+00,
c$$$     -      0.474101304688892224560294156365470354D+00,
c$$$     -      0.517973906222155508794116872690563303D-01,
c$$$     -      0.717926712780444024794562208158499202D+00,
c$$$     -      0.233916149991311356180523963867127431D+00,
c$$$     -      0.717926712780444024794562208158501705D+00,
c$$$     -      0.481571372282446190249138279743715620D-01,
c$$$     -      0.481571372282446190249138279743708277D-01,
c$$$     -      0.233916149991311356180523963867129188D+00,
c$$$     -      0.897129716839089969485816308024824437D+00,
c$$$     -      0.514351415804550152570918459875888408D-01,
c$$$     -      0.514351415804550152570918459875868005D-01/
c$$$          data whts4/
c$$$     -      0.635086547941006332905328779420276919D-01,
c$$$     -      0.635086547941006332905328779420276919D-01,
c$$$     -      0.635086547941006332905328779420276919D-01,
c$$$     -      0.323254868211278648566249984545284920D-01,
c$$$     -      0.323254868211278648566249984545277156D-01,
c$$$     -      0.323254868211278648566249984545277156D-01,
c$$$     -      0.273576739148287444345776008262900919D-01,
c$$$     -      0.273576739148287444345776008262900919D-01,
c$$$     -      0.273576739148287444345776008262897067D-01,
c$$$     -      0.273576739148287444345776008262897067D-01,
c$$$     -      0.273576739148287444345776008262893215D-01,
c$$$     -      0.273576739148287444345776008262889393D-01,
c$$$     -      0.161171772217806796503535886175348952D-01,
c$$$     -      0.161171772217806796503535886175348952D-01,
c$$$     -      0.161171772217806796503535886175348952D-01/
c
c       8th order
c
        parameter ( nquad8  =  45 , m8 =  14)
        dimension xs8(45)
        dimension ys8(45)
        dimension whts8(45)
          data xs8/
     -      0.558243537563118223900322345008054564D+00,
     -      0.386028294454679857386767989870977846D+00,
     -      0.386028294454679857386767989870874423D+00,
     -      0.220878231218440888049838827496039548D+00,
     -      0.220878231218440888049838827495911353D+00,
     -      0.227943411090640285226464020258152281D+00,
     -      0.976856167428592844390855286593242875D-01,
     -      0.369152488188135087815437016812445146D+00,
     -      0.976856167428592844390855286593572690D-01,
     -      0.533161895069005627745477454528240750D+00,
     -      0.369152488188135087815437016812707169D+00,
     -      0.533161895069005627745477454527881468D+00,
     -      0.199351725530086044110629058082517698D-01,
     -      0.490032413723495697794468547096121405D+00,
     -      0.490032413723495697794468547095570251D+00,
     -      0.924206243631347164496734987249412486D-01,
     -      0.221384619787592680628127744079283071D+00,
     -      0.182698666764278740464419880939257108D-01,
     -      0.348360091360732930340708152478132391D+00,
     -      0.686194755849272602922198757195870051D+00,
     -      0.924206243631347164496734987249975219D-01,
     -      0.182698666764278740464419880939250608D-01,
     -      0.221384619787592680628127744079842602D+00,
     -      0.686194755849272602922198757195041420D+00,
     -      0.633370041962839195612849859427932798D+00,
     -      0.348360091360732930340708152478531300D+00,
     -      0.633370041962839195612849859427491182D+00,
     -      0.184371498704658312332689028556717367D-01,
     -      0.992485769679237492889912597299867063D-01,
     -      0.199943975938913669799180174371952787D+00,
     -      0.184371498704658312332689028556898284D-01,
     -      0.781618874190620498967550922771945849D+00,
     -      0.801502846064152501422017480539678957D+00,
     -      0.781618874190620498967550922772302820D+00,
     -      0.992485769679237492889912597303208551D-01,
     -      0.199943975938913669799180174372423533D+00,
     -      0.899630497824555747567092213056902787D+00,
     -      0.899630497824555747567092213056160245D+00,
     -      0.194270438632592782377755233646741943D-01,
     -      0.809424583121849741951322635786297672D-01,
     -      0.809424583121849741951322635788864816D-01,
     -      0.194270438632592782377755233647542528D-01,
     -      0.972665693372719704666628843226822511D+00,
     -      0.136671533136401476666855783865507556D-01,
     -      0.136671533136401476666855783866338204D-01/
          data ys8/
     -      0.220878231218440888049838827496032037D+00,
     -      0.227943411090640285226464020258195543D+00,
     -      0.386028294454679857386767989871018916D+00,
     -      0.220878231218440888049838827495960464D+00,
     -      0.558243537563118223900322345008095394D+00,
     -      0.386028294454679857386767989870934464D+00,
     -      0.533161895069005627745477454528062409D+00,
     -      0.533161895069005627745477454528234876D+00,
     -      0.369152488188135087815437016812582754D+00,
     -      0.369152488188135087815437016812439272D+00,
     -      0.976856167428592844390855286593661524D-01,
     -      0.976856167428592844390855286593625894D-01,
     -      0.490032413723495697794468547095454936D+00,
     -      0.490032413723495697794468547095638670D+00,
     -      0.199351725530086044110629058082631087D-01,
     -      0.686194755849272602922198757195478798D+00,
     -      0.686194755849272602922198757195723680D+00,
     -      0.633370041962839195612849859427627730D+00,
     -      0.633370041962839195612849859427941079D+00,
     -      0.221384619787592680628127744079209934D+00,
     -      0.221384619787592680628127744079602486D+00,
     -      0.348360091360732930340708152478314970D+00,
     -      0.924206243631347164496734987249373727D-01,
     -      0.924206243631347164496734987250118821D-01,
     -      0.348360091360732930340708152478147125D+00,
     -      0.182698666764278740464419880939220034D-01,
     -      0.182698666764278740464419880939275404D-01,
     -      0.781618874190620498967550922772341628D+00,
     -      0.801502846064152501422017480539853350D+00,
     -      0.781618874190620498967550922772360598D+00,
     -      0.199943975938913669799180174372436870D+00,
     -      0.199943975938913669799180174372380777D+00,
     -      0.992485769679237492889912597299276765D-01,
     -      0.184371498704658312332689028556904664D-01,
     -      0.992485769679237492889912597301857752D-01,
     -      0.184371498704658312332689028556764553D-01,
     -      0.194270438632592782377755233646571979D-01,
     -      0.809424583121849741951322635790730681D-01,
     -      0.899630497824555747567092213056792624D+00,
     -      0.899630497824555747567092213056663587D+00,
     -      0.194270438632592782377755233647268715D-01,
     -      0.809424583121849741951322635789923957D-01,
     -      0.136671533136401476666855783867173485D-01,
     -      0.972665693372719704666628843226955978D+00,
     -      0.136671533136401476666855783866795431D-01/
          data whts8/
     -      0.226843308437213620015618859369815019D-01,
     -      0.259623713402991498139240506829434283D-01,
     -      0.259623713402991498139240506829434283D-01,
     -      0.226843308437213620015618859369861212D-01,
     -      0.226843308437213620015618859369764975D-01,
     -      0.259623713402991498139240506829434283D-01,
     -      0.171294207506684469666459370206662877D-01,
     -      0.171294207506684469666459370206928445D-01,
     -      0.171294207506684469666459370206662877D-01,
     -      0.171294207506684469666459370206990044D-01,
     -      0.171294207506684469666459370206435798D-01,
     -      0.171294207506684469666459370206378050D-01,
     -      0.688859943626271717032432190540610911D-02,
     -      0.688859943626271717032432190542535411D-02,
     -      0.688859943626271717032432190539225294D-02,
     -      0.139119010181500454704467661290159098D-01,
     -      0.139119010181500454704467661290186031D-01,
     -      0.704118755793328342802064630902593490D-02,
     -      0.704118755793328342802064630901746834D-02,
     -      0.139119010181500454704467661290170639D-01,
     -      0.139119010181500454704467661290205260D-01,
     -      0.704118755793328342802064630898359680D-02,
     -      0.139119010181500454704467661290174475D-01,
     -      0.139119010181500454704467661290201423D-01,
     -      0.704118755793328342802064630898513604D-02,
     -      0.704118755793328342802064630899013969D-02,
     -      0.704118755793328342802064630902632084D-02,
     -      0.651166218169209675165003504125206928D-02,
     -      0.111317640131603406629743748298020702D-01,
     -      0.651166218169209675165003504124552564D-02,
     -      0.651166218169209675165003504124706563D-02,
     -      0.651166218169209675165003504123089985D-02,
     -      0.111317640131603406629743748298190063D-01,
     -      0.651166218169209675165003504126631063D-02,
     -      0.111317640131603406629743748298351706D-01,
     -      0.651166218169209675165003504125052929D-02,
     -      0.472640077673590400124718231947651779D-02,
     -      0.472640077673590400124718231950731010D-02,
     -      0.472640077673590400124718231947805778D-02,
     -      0.472640077673590400124718231948498586D-02,
     -      0.472640077673590400124718231949999758D-02,
     -      0.472640077673590400124718231950577086D-02,
     -      0.135845646286354378186089967265265773D-02,
     -      0.135845646286354378186089967264284263D-02,
     -      0.135845646286354378186089967266301148D-02/
c
c       12th order quadrature
c
        parameter ( nquad12  =  91 , m12 =  21)
        dimension xs12(91),ys12(91),whts12(91)

          data xs12/
     -      0.333333333333333329625344425380629097D+00,
     -      0.453485958400438141752080165796997293D+00,
     -      0.423376866886948430894533053247046668D+00,
     -      0.318681771755455104040715999793795735D+00,
     -      0.336268018099417166348302742135652120D+00,
     -      0.123137174712613383899096485983947555D+00,
     -      0.576578269784060891109753055815250344D+00,
     -      0.216312025283207195932848496066213507D+00,
     -      0.220645477532419357383797975902767332D+00,
     -      0.202668324715234788672947915819174906D+00,
     -      0.220753405500704308576207902086405623D+00,
     -      0.447419956617375624450937115283897907D+00,
     -      0.460672750712125525314219774913796527D+00,
     -      0.346547144301823793464971600385073969D+00,
     -      0.312593585179935895222292861065545548D+00,
     -      0.117079315746674577670660806012604893D+00,
     -      0.530575302132329047366005629445571682D+00,
     -      0.570327099073389493841391630087212796D+00,
     -      0.122877553565847132166374093159788673D+00,
     -      0.526117713932697428855478723830700035D+00,
     -      0.427040465172698394267540113116710540D+00,
     -      0.422796314366814428459646842871517612D+00,
     -      0.521535392157938777007801456270350991D+00,
     -      0.514241426693628175886562346380505374D-01,
     -      0.510859717004881286150755701524625829D-01,
     -      0.491015539538148134053433303364163890D+00,
     -      0.498994318340720903625091105487098087D+00,
     -      0.637273922968470467777762275913892607D+00,
     -      0.630141169090691958531176172547828597D+00,
     -      0.638309471793628178116678136188041190D+00,
     -      0.683687206963930285447546957361475358D+00,
     -      0.596776462552731370201022101135655521D+00,
     -      0.596981389696067658197670997842785481D+00,
     -      0.317964549010314339105917471351492631D+00,
     -      0.238144744452491644408509084585633631D+00,
     -      0.110621517328064472596710377024873598D+00,
     -      0.393388209889084133145426167011647053D+00,
     -      0.476921546525338066300916354321954650D-01,
     -      0.313998373553838020109912619488036713D+00,
     -      0.999014212113090048090847904078203275D-02,
     -      0.205691275708005208359225156020206236D+00,
     -      0.124581332579037883293366207071866214D+00,
     -      0.393426231829825299649964896017299947D+00,
     -      0.959237847410703784073978092950682885D-02,
     -      0.518942818989937114194496878288752253D-01,
     -      0.983532755818445044205730968984147333D-02,
     -      0.743248343984407048572923601603841951D+00,
     -      0.709495811881816818934590463205674517D+00,
     -      0.750051716886630805112402853962666184D+00,
     -      0.743404748389881844482487077196064308D+00,
     -      0.717278189947551883953036646456091284D+00,
     -      0.280438471054996379709403345058988586D+00,
     -      0.205174969358066699348514688024303383D+00,
     -      0.137200762947225428000878278752850096D+00,
     -      0.469458253122042592041177119337163517D-01,
     -      0.891267454191634696506687992859014688D-02,
     -      0.273809135510531779901095234661318495D+00,
     -      0.203002457801164942646045649272051534D+00,
     -      0.515766866575262614592790223154877054D-01,
     -      0.119394488662892722063653328091131181D+00,
     -      0.100657170631867980410356446137301012D-01,
     -      0.816488069180110039441213287686157423D+00,
     -      0.824974361710174346121561212542945134D+00,
     -      0.834348293717586124544912390001573559D+00,
     -      0.843988834191637596025263836669565460D+00,
     -      0.173796909009967971619619977313971966D+00,
     -      0.165374898318327731504849089004733087D+00,
     -      0.106607031789178734625176559800896321D+00,
     -      0.965073997149793015351298410831596544D-02,
     -      0.562472075517419931691259535016920453D-01,
     -      0.109404498730671887524555400562639721D+00,
     -      0.971502180992200693158633041817407519D-02,
     -      0.494041340191836790574051769106381964D-01,
     -      0.903686716557178626676179625820052778D+00,
     -      0.890786345929322891319301273309844703D+00,
     -      0.924078454781401327285736857498249921D+00,
     -      0.987645360305592567968881871738099886D-01,
     -      0.374467902950908360389512891664162952D-01,
     -      0.111559284893352789387092572911913226D-01,
     -      0.104491180401179062415919346065715562D-01,
     -      0.851573549534860853105908123047697951D-01,
     -      0.384747549235078402866192976584749204D-01,
     -      0.964573477007098585000154765657181513D+00,
     -      0.937182472201716171104980741525171560D+00,
     -      0.177168316118745478598494262977789148D-02,
     -      0.567607399982020155281925837494285659D-01,
     -      0.605678780008183561677193571048268137D-02,
     -      0.336548398317139361284814847076681681D-01,
     -      0.745448371186314990636081730402619155D-02,
     -      0.979479656502450198980185101050269768D+00,
     -      0.130658597856866427728633282171096293D-01/
          data ys12/
     -      0.333333333333333338926511306387301352D+00,
     -      0.123137174712613389796284484355989413D+00,
     -      0.453485958400438183053346446157253315D+00,
     -      0.220645477532419365026898892204410917D+00,
     -      0.447419956617375636109755148982114570D+00,
     -      0.423376866886948481426582060297305796D+00,
     -      0.202668324715234796932306410944216206D+00,
     -      0.336268018099417170540644503119893171D+00,
     -      0.460672750712125530519436549330124016D+00,
     -      0.220753405500704304289863906167151755D+00,
     -      0.576578269784060894510387637023801050D+00,
     -      0.216312025283207203470348734906773651D+00,
     -      0.318681771755455115257891003381090174D+00,
     -      0.122877553565847136911391517534639384D+00,
     -      0.570327099073389531600070609140107842D+00,
     -      0.312593585179935921362494007016748491D+00,
     -      0.346547144301823820365148669531048315D+00,
     -      0.117079315746674578675598292595261930D+00,
     -      0.530575302132329097040292556130853363D+00,
     -      0.422796314366814442108088792559266698D+00,
     -      0.521535392157938786424693291029149574D+00,
     -      0.510859717004881323335107774423356730D-01,
     -      0.514241426693628196525008889822289954D-01,
     -      0.427040465172698325266849760348589465D+00,
     -      0.526117713932697428565579492768084165D+00,
     -      0.999014212113090200494671513174945898D-02,
     -      0.491015539538148194956845426908695302D+00,
     -      0.238144744452491647663034043638474035D+00,
     -      0.317964549010314330711643698097164726D+00,
     -      0.476921546525338090389833259096528104D-01,
     -      0.110621517328064484914354330943304550D+00,
     -      0.393388209889084179659191082726404404D+00,
     -      0.959237847410703744034437823293852118D-02,
     -      0.518942818989937105161482903292540267D-01,
     -      0.124581332579037879966443593048999709D+00,
     -      0.205691275708005222429957844724013038D+00,
     -      0.983532755818445029236344773597293803D-02,
     -      0.313998373553837959848518394209317887D+00,
     -      0.638309471793628172472004688258676635D+00,
     -      0.498994318340720859368102914084577991D+00,
     -      0.683687206963930301194114686978280527D+00,
     -      0.637273922968470500015978504615167254D+00,
     -      0.596981389696067663028891354025626804D+00,
     -      0.393426231829825262341989607197846522D+00,
     -      0.630141169090692006787945536755572145D+00,
     -      0.596776462552731389824669675025927225D+00,
     -      0.205174969358066687350640292871338057D+00,
     -      0.280438471054996383085088808600925463D+00,
     -      0.469458253122042662463262463472804273D-01,
     -      0.137200762947225423913916586988929075D+00,
     -      0.891267454191634810367102068506602066D-02,
     -      0.100657170631867979550404453170246996D-01,
     -      0.515766866575262619533851560687099077D-01,
     -      0.119394488662892729429674377926756581D+00,
     -      0.203002457801164903677927469506823680D+00,
     -      0.273809135510531773985945227933536589D+00,
     -      0.717278189947551872122078515295834880D+00,
     -      0.750051716886630786777141404650095196D+00,
     -      0.743248343984407091835331155747005729D+00,
     -      0.743404748389881869842925775034762250D+00,
     -      0.709495811881816854575469908479398689D+00,
     -      0.173796909009967952144992819659049230D+00,
     -      0.965073997149793160467515103965484761D-02,
     -      0.562472075517419904070637701739729941D-01,
     -      0.106607031789178723721131119661634524D+00,
     -      0.971502180992200755994669137315468616D-02,
     -      0.824974361710174335723880653687367049D+00,
     -      0.494041340191836809269848713426103775D-01,
     -      0.165374898318327742279251644177757360D+00,
     -      0.109404498730671878598156425315514568D+00,
     -      0.834348293717586133878606298655573638D+00,
     -      0.816488069180110064945397987067566016D+00,
     -      0.843988834191637631363087771554255257D+00,
     -      0.111559284893352783907849072506271645D-01,
     -      0.987645360305592112763792649106752310D-01,
     -      0.374467902950908285279999881273433413D-01,
     -      0.104491180401178991512496411547300236D-01,
     -      0.384747549235078436383193132432503648D-01,
     -      0.851573549534861118798358070245020654D-01,
     -      0.890786345929322824987162415232208480D+00,
     -      0.903686716557178641916311394686810332D+00,
     -      0.924078454781401350832900364354651790D+00,
     -      0.177168316118744458114996798241174717D-02,
     -      0.567607399982020008759362396976951041D-01,
     -      0.336548398317139700745307342683687804D-01,
     -      0.605678780008183092033215218077137419D-02,
     -      0.937182472201716073248430965903544141D+00,
     -      0.964573477007098629952387592828151968D+00,
     -      0.979479656502450191921120269695053927D+00,
     -      0.130658597856866423956882845120568920D-01,
     -      0.745448371186315475756935984930375185D-02/
          data whts12/
     -      0.144337970552756069678532581553471415D-01,
     -      0.932731826914862573744300647801222265D-02,
     -      0.932731826914862572317846176470315303D-02,
     -      0.133414502592278536276382634845548973D-01,
     -      0.132742614018645820320490082945317632D-01,
     -      0.932731826914862716963742726399427236D-02,
     -      0.112924526214995025903210306435216326D-01,
     -      0.132742614018645824323888943878849106D-01,
     -      0.133414502592278533245217177342119775D-01,
     -      0.112924526214995027458958797535833730D-01,
     -      0.112924526214995018795457453764226075D-01,
     -      0.132742614018645824892200647889677453D-01,
     -      0.133414502592278536326511672416151888D-01,
     -      0.908703195954351706452414546876551217D-02,
     -      0.935157125515536388802285926730309762D-02,
     -      0.935157125515536548626683003335347303D-02,
     -      0.908703195954351882568334291363922717D-02,
     -      0.935157125515536504616712843549578078D-02,
     -      0.908703195954351743224408698199646630D-02,
     -      0.547457668520058607113533337585169630D-02,
     -      0.613972585967859798199335110619528480D-02,
     -      0.547457668520058534948129702145195758D-02,
     -      0.613972585967859904686344070744502382D-02,
     -      0.613972585967859870655957091976995533D-02,
     -      0.547457668520059151999619163442746915D-02,
     -      0.242044605701187826886006244428228969D-02,
     -      0.242044605701187700752145609770442757D-02,
     -      0.926522340175687688728137443667639231D-02,
     -      0.649799900408112598144092103362485034D-02,
     -      0.625035529062928911790244470404942024D-02,
     -      0.764443844999929615109626688771881622D-02,
     -      0.268633037625562631532042921933237854D-02,
     -      0.287120763916888119870166252638653908D-02,
     -      0.649799900408112525431771053494992097D-02,
     -      0.926522340175687533342856096753138190D-02,
     -      0.764443844999929596362635834632657365D-02,
     -      0.268633037625562502273050152418603527D-02,
     -      0.625035529062928718653951068427955469D-02,
     -      0.625035529062928827259870084807301341D-02,
     -      0.242044605701187790831719280987930658D-02,
     -      0.764443844999929484347412420495158903D-02,
     -      0.926522340175687534502296104478897483D-02,
     -      0.287120763916888021974041979815342651D-02,
     -      0.287120763916887947536051526522180407D-02,
     -      0.649799900408112620304452518593696035D-02,
     -      0.268633037625562757828847857832358650D-02,
     -      0.626168281551884286997165836612699734D-02,
     -      0.293390633207459151832955734910682988D-02,
     -      0.534889320349777347256104311846305886D-02,
     -      0.756148376260634250407847837555738681D-02,
     -      0.269310521628365013404043790226052484D-02,
     -      0.293390633207459076037758922156821539D-02,
     -      0.626168281551884271502036584503426714D-02,
     -      0.756148376260634220771502502127544489D-02,
     -      0.534889320349777150651477987515455999D-02,
     -      0.269310521628364923552611323523998294D-02,
     -      0.269310521628365003561122948210356452D-02,
     -      0.534889320349777394268147953143437781D-02,
     -      0.626168281551884209977956039312305823D-02,
     -      0.756148376260634079880516834832094621D-02,
     -      0.293390633207459144664108233404699655D-02,
     -      0.237899791896875977650412068963234244D-02,
     -      0.235561888924248391206500830125995255D-02,
     -      0.522958752630325134979443393750609435D-02,
     -      0.464375501038344592818763522385505103D-02,
     -      0.237899791896875903190464295118641807D-02,
     -      0.235561888924248441819813956082152834D-02,
     -      0.464375501038344592712461698378333248D-02,
     -      0.235561888924248340213283223932878332D-02,
     -      0.522958752630325094695988573188260163D-02,
     -      0.522958752630325102269169024696529860D-02,
     -      0.237899791896875884828277915592515253D-02,
     -      0.464375501038344502821066760119404255D-02,
     -      0.185123842534245045642152424523734433D-02,
     -      0.133602633352073711936052192239925326D-02,
     -      0.253983745289682117104431550561110727D-02,
     -      0.133602633352073706721550666005480491D-02,
     -      0.253983745289682141266395853190635430D-02,
     -      0.185123842534245093614988352760757316D-02,
     -      0.133602633352073439722171701863745652D-02,
     -      0.185123842534245025454024657012139802D-02,
     -      0.253983745289682046646751921952733996D-02,
     -      0.416673666706645018964671763530887448D-03,
     -      0.785647715964213819702516503481224355D-03,
     -      0.416673666706645627599941307437333452D-03,
     -      0.785647715964214725113788618324958244D-03,
     -      0.785647715964216159196092971657402644D-03,
     -      0.416673666706644495756873072324502098D-03,
     -      0.594558182043188711788356789297320778D-03,
     -      0.594558182043188613542907998875130600D-03,
     -      0.594558182043188401643956376365038440D-03/
c
        parameter ( nquad16  = 153 , m16 =  28)
        dimension xs16(153)
        dimension ys16(153)
        dimension whts16(153)
          data xs16/
     -      0.302431316268704796042908925414429991D+00,
     -      0.302529607273387195840090942884628791D+00,
     -      0.395063931591685766580790370315143449D+00,
     -      0.214454624832134983293520189795954943D+00,
     -      0.392800214661270243073631605301961402D+00,
     -      0.392707801300601402856617542651978303D+00,
     -      0.215337681019017224560746971385909230D+00,
     -      0.215459457248522571277767806649661870D+00,
     -      0.313866827106886098836378411046253892D+00,
     -      0.313680385257148910094559617289742490D+00,
     -      0.136014642972966562425872284425551616D+00,
     -      0.136122014369109587307019286345147308D+00,
     -      0.470874637716072828708269395971127229D+00,
     -      0.470742049552093874784092312816953938D+00,
     -      0.411842922722062449627975514066590870D+00,
     -      0.412101270891199374868676075450660076D+00,
     -      0.451861994471425207511387308829479835D+00,
     -      0.452006092444166009354045989524006783D+00,
     -      0.704236186042401044205665719780452612D-01,
     -      0.142798091556413265575712728592364645D+00,
     -      0.229435176839617811539241107931836866D+00,
     -      0.790979204216395370415979228785669766D-01,
     -      0.319046810417511801280426860253490161D+00,
     -      0.142755289398359432143929059080280715D+00,
     -      0.229580614381211444862966608085057818D+00,
     -      0.790747093586052062773249077171236731D-01,
     -      0.319051113419125136499004548798746625D+00,
     -      0.302092334804568866756321171252597931D-01,
     -      0.301999179529065382726625305332961063D-01,
     -      0.464692522537515681093527973426992026D+00,
     -      0.379096314684441817016010293812640848D+00,
     -      0.464877676276147431037631419513818260D+00,
     -      0.379309363514795621564214527103716526D+00,
     -      0.540955735628189996549776000823291549D+00,
     -      0.538077058005980497174920185452171794D+00,
     -      0.538248894381397390104474899872679727D+00,
     -      0.541666122754701756803449020522569051D+00,
     -      0.541727056588530671363995033933215990D+00,
     -      0.436312825084017629484148418855896117D+00,
     -      0.808807722015413097593870694685197958D-01,
     -      0.148277622038385319135942908979684882D+00,
     -      0.343683357911596655578125628959215728D-01,
     -      0.222981777289783111623034382301376524D+00,
     -      0.436120365567610864581161952467669709D+00,
     -      0.808798331723407056866056378960210539D-01,
     -      0.148351474939894356259964230692998276D+00,
     -      0.343714210911214422478052892029922188D-01,
     -      0.223001185448292710350772858819606479D+00,
     -      0.570629496803482582331161353659122729D-02,
     -      0.282588438564370454780751891038166941D+00,
     -      0.593238615216649726220714810038395926D-02,
     -      0.533727165799185248336524981234301003D+00,
     -      0.593490617588697134572074466200151541D-02,
     -      0.342063761518813358415613921181516747D+00,
     -      0.282392878512121888372618444554984747D+00,
     -      0.533506699039126054321742932575199881D+00,
     -      0.342129986680270623924724086362298011D+00,
     -      0.325440594206227473703889550807903073D-01,
     -      0.636670451161602809395319825718889836D+00,
     -      0.810186712115616221779897200860145583D-01,
     -      0.628718994082772383001095165476771117D+00,
     -      0.628641043064132488364403400693794356D+00,
     -      0.678950312264433475002648429802974174D-02,
     -      0.636547229891117437361497703144190920D+00,
     -      0.623505358703628084303469882949964460D+00,
     -      0.497312604908278420474244671680644687D+00,
     -      0.410110188893687161685741019900360188D+00,
     -      0.139771563743534697386438422832818416D+00,
     -      0.195497094483522524279545095037057925D+00,
     -      0.256034463426759594648698558433339034D+00,
     -      0.583924576490433769420604538909497739D+00,
     -      0.331730655223533413939525902653474592D+00,
     -      0.623580589410633476872821605277651295D+00,
     -      0.497019311879614643128890083439809192D+00,
     -      0.325369488384764471446406593559536406D-01,
     -      0.679788439682447416074613022434303129D-02,
     -      0.809539807646663610346484130294200715D-01,
     -      0.139745700249707659703867808328520299D+00,
     -      0.410114160725183169568424747779369546D+00,
     -      0.195382170809990366377774747364760428D+00,
     -      0.583975772597892459870844398468411805D+00,
     -      0.256572605061291280154468836722930107D+00,
     -      0.332039480616851747279058170079814858D+00,
     -      0.710976355567288538122899630658194868D+00,
     -      0.723566361791087094873375757557110548D+00,
     -      0.661147225451088997308206370121085565D+00,
     -      0.720475120669841255285854263497732895D+00,
     -      0.723496593909361049700640460004546768D+00,
     -      0.711296815925722935268565536529695027D+00,
     -      0.661427979010520012744181105212079590D+00,
     -      0.346330349392944591990797267909754967D-01,
     -      0.608498245336996135565922807261259583D-02,
     -      0.735299673011199805957818158576469860D-01,
     -      0.189835770123654892944939229266498793D+00,
     -      0.120670382592678134809782370486853114D+00,
     -      0.256643550191060746094826903026243876D+00,
     -      0.607473138257149765615981052144018049D-02,
     -      0.346543468339740953746565367285376925D-01,
     -      0.735671485913032494842963124946850544D-01,
     -      0.189971510345967543605787386292094616D+00,
     -      0.120670383898476665141025834179849337D+00,
     -      0.256423400907761779828832629022521124D+00,
     -      0.737507802811185601329991748857187262D+00,
     -      0.775463562856272680570907925638621352D+00,
     -      0.805729308090636442822674720257585150D+00,
     -      0.775501819349450326203952869084154400D+00,
     -      0.805728051460471777406311997001320774D+00,
     -      0.737302338909043658230321351086818047D+00,
     -      0.659189817872052775481348615809581159D-02,
     -      0.658576706402559977850605805077527984D-02,
     -      0.295944438229564915292844316507010762D-01,
     -      0.181820601962446774734848948159423705D+00,
     -      0.631481449194973626672948770181718869D-01,
     -      0.131673645884505720909166058451635457D+00,
     -      0.296669794377902856156850691202255738D-01,
     -      0.181372697758542400101559789792474659D+00,
     -      0.631622223884014726494127337101266170D-01,
     -      0.130738223783696236467367080418320221D+00,
     -      0.514130897896413517091633778227846815D-02,
     -      0.812005602644382361947724604742832728D+00,
     -      0.112522934301934473892299084087267144D+00,
     -      0.839547898752369226496979459102329322D+00,
     -      0.277228855877048647708636644258354339D-01,
     -      0.754562257706843572911137672849185701D-01,
     -      0.838970100731769884041696583085674329D+00,
     -      0.811678904133436929663380399687301820D+00,
     -      0.873653971675868590589563887358446063D+00,
     -      0.517208922751558675302322326435178891D-02,
     -      0.112473372415410432552997875289198339D+00,
     -      0.277450302699471893196474991747660898D-01,
     -      0.746382767337317840935809240806939676D-01,
     -      0.599417935348845665213930193583089614D-02,
     -      0.882351468893980981594998196921061486D+00,
     -      0.283127918702113744410707280311096603D-01,
     -      0.586750472274116602191316119243568267D-01,
     -      0.897515201642183193553820542520593072D+00,
     -      0.882357950464212301923939093934990858D+00,
     -      0.897000811817480910219103274334164871D+00,
     -      0.594972831709593512676842787964971349D-02,
     -      0.596779279106414468812379941919327535D-01,
     -      0.279934967976983222407162133959849913D-01,
     -      0.467386635060392000991445300819857411D-02,
     -      0.255674849837572119225141607664514332D-01,
     -      0.934403792731877932046717796275418114D+00,
     -      0.448570049906792238720874795935560360D-02,
     -      0.270887529332885512293658220422164004D-01,
     -      0.935206741317988752650562001571016261D+00,
     -      0.943728547283893839247936886394643257D+00,
     -      0.566483367989507478744732640768804684D-02,
     -      0.968403158649102667023666426325496207D+00,
     -      0.969635566346902796594897153857305494D+00,
     -      0.604718232198114385285917248211208580D-02,
     -      0.988243812623773650098104457318075309D+00/
          data ys16/
     -      0.302506639757888189013497704164165343D+00,
     -      0.395023014355676166275666030979515299D+00,
     -      0.302405996986061132992224601068261933D+00,
     -      0.392724814541112431328685073458221097D+00,
     -      0.214449448003336246221337185458778626D+00,
     -      0.392818490000195511817778236458288368D+00,
     -      0.313703112464622089080780519846717467D+00,
     -      0.470783977926051083866072868460406771D+00,
     -      0.215431228216478041841876166229139431D+00,
     -      0.470968517913461039082400273277917608D+00,
     -      0.412172586648968035347741620293481726D+00,
     -      0.452017177439447533969453397140966127D+00,
     -      0.215306845578487967066960775878442769D+00,
     -      0.313816971885511836602196402368492742D+00,
     -      0.136092961137842142823008309158415790D+00,
     -      0.451810764236290959740827108404993301D+00,
     -      0.136029782544618383691365546596706839D+00,
     -      0.411922268527803597951174716793526567D+00,
     -      0.464724918652927127298945992206154761D+00,
     -      0.538095044838347680064568392111584771D+00,
     -      0.541060947019147681897201703966523168D+00,
     -      0.541681588682656427162388773154869800D+00,
     -      0.538213869071652947618556070816756479D+00,
     -      0.319058906889390376175616124096419301D+00,
     -      0.229432000472273582840362931936235530D+00,
     -      0.379087070593918008312594519813443010D+00,
     -      0.142811380691040796989472935729613393D+00,
     -      0.533480426790031455877099114475920354D+00,
     -      0.436077568198375159377711037662870615D+00,
     -      0.464876474824343383504633098177494644D+00,
     -      0.541822970659543677254135639595132502D+00,
     -      0.703731021601924599502464361990750640D-01,
     -      0.791084530757112954533094035504834557D-01,
     -      0.229530378001768189198713303929281805D+00,
     -      0.319108677726218944823632069664368171D+00,
     -      0.142748636032904484572535932161216948D+00,
     -      0.379235965456933009428768732373850994D+00,
     -      0.790827386355190834510825990755846345D-01,
     -      0.302026365176680066939431704636438858D-01,
     -      0.636587233574720824834180693943197575D+00,
     -      0.628676401895567326091464592434181796D+00,
     -      0.623576364706813956687066047331697863D+00,
     -      0.628697636820986106786270948047996953D+00,
     -      0.533684375074890656596805194269090409D+00,
     -      0.282408830134979313498308287540765642D+00,
     -      0.222993821442986748412522059086861716D+00,
     -      0.342135405218172268996662583393054439D+00,
     -      0.148269842150567330120426806623653546D+00,
     -      0.496962472826357519401084600127588770D+00,
     -      0.808782926371964177666504435447663702D-01,
     -      0.410114933105515248089098351305762761D+00,
     -      0.301868201156704933515873983154655834D-01,
     -      0.583990124589469130542476420933565214D+00,
     -      0.343744755415233284514050500643565623D-01,
     -      0.636740676880711625593146746151735334D+00,
     -      0.436272413049082663311148230272533390D+00,
     -      0.623495955707887405761703772569121973D+00,
     -      0.256666017481659167720198765656216357D+00,
     -      0.808686295412089851681278108870748312D-01,
     -      0.195426732061812465955074140773945481D+00,
     -      0.148331042830191477608926875476962181D+00,
     -      0.223052104529295399765296516505212326D+00,
     -      0.332084040491790397180703616330650315D+00,
     -      0.282563764973116401435206767943775641D+00,
     -      0.343799986968225864826255129365268020D-01,
     -      0.570598973835654113038022323270193797D-02,
     -      0.593298221447241169112544765601974802D-02,
     -      0.139751006650496131688994826586614911D+00,
     -      0.809483380355830476372411786106487408D-01,
     -      0.325303704551144566003975463133009808D-01,
     -      0.592749394412468488535711874646897473D-02,
     -      0.680068349740502027488739017066542050D-02,
     -      0.342054288934002137894823640449908733D+00,
     -      0.497273739021581327085798544328304140D+00,
     -      0.711342921714536800830358074558413137D+00,
     -      0.661456929394697327323333448016323113D+00,
     -      0.723555098947398454555413718568346648D+00,
     -      0.720438507970551145780696711402273123D+00,
     -      0.583955840996834646181292302407530482D+00,
     -      0.723613422581703055080849327492366272D+00,
     -      0.410085605416883119080658599521281476D+00,
     -      0.710895225957279083616844676813283244D+00,
     -      0.661166752122906905507080872220305180D+00,
     -      0.325291288420212115282898305788117248D-01,
     -      0.810100986544337710535529181659827351D-01,
     -      0.679600214364638479504962098110219113D-02,
     -      0.139810910353030919359673096419107166D+00,
     -      0.195532571417567819014366038013559889D+00,
     -      0.256158485738086570874783166626244994D+00,
     -      0.331778657044049100638639341197860310D+00,
     -      0.189957824052571577178524253418237842D+00,
     -      0.256377595325605179991965367234527409D+00,
     -      0.120726876191067716986855944278299710D+00,
     -      0.346686197087038294372357793967345150D-01,
     -      0.736000543527183188912282150808868403D-01,
     -      0.607032781611058988631753759206280090D-02,
     -      0.737267012017499438581445512214053181D+00,
     -      0.775451182219839736795478582150875765D+00,
     -      0.805707744570418641094355017521021707D+00,
     -      0.775387763482114600382249918078445426D+00,
     -      0.805731140149949025168941227639397365D+00,
     -      0.737497915704729386592511711008432240D+00,
     -      0.607616583064277159708280403445932297D-02,
     -      0.346564068196919889137421052860318973D-01,
     -      0.736100288467653936786082773570347493D-01,
     -      0.189840692548069573293472154149878905D+00,
     -      0.120754104754303069436128081961219185D+00,
     -      0.256617886154587488882576378713328595D+00,
     -      0.181307699544105378876863594762710867D+00,
     -      0.811633186444285223636652136511531428D+00,
     -      0.130567287998519356380558549019232250D+00,
     -      0.658908867670266920787433362328961100D-02,
     -      0.631902214720203394630491725497080701D-01,
     -      0.296916723233359053913078616928311408D-01,
     -      0.838813687821253997373240933935807822D+00,
     -      0.812032546773107255671761965785738726D+00,
     -      0.873641360547783384796064491219506864D+00,
     -      0.839609784817584699785771886277906482D+00,
     -      0.112472906303878340183125107314468163D+00,
     -      0.659786000905907757035073902036725596D-02,
     -      0.517751591171327666870949540320616152D-02,
     -      0.296509061670442601117380724752973331D-01,
     -      0.745193086448869926508886706900978463D-01,
     -      0.277570770874521871535922393538214212D-01,
     -      0.131411188814850913563986670067505631D+00,
     -      0.181737061408330033530409093175911464D+00,
     -      0.631927474312769659199911506244472057D-01,
     -      0.882336302533193072878741217553961682D+00,
     -      0.882369702557596172606941555803321742D+00,
     -      0.896943547314205755386059337013847410D+00,
     -      0.897609276668989874800432263004019415D+00,
     -      0.598467358876075174803549580344357859D-01,
     -      0.515522280757245403765741792705181386D-02,
     -      0.279477132316105768634962508981903040D-01,
     -      0.594334741484913691075806183384319236D-02,
     -      0.277482794633089236659259802468935096D-01,
     -      0.112481709281704536646865016953522450D+00,
     -      0.752611688353871003104010799564950707D-01,
     -      0.935287722641008374735228501540757032D+00,
     -      0.934332564496973307177003641356347279D+00,
     -      0.943753403753687139807420114838495031D+00,
     -      0.273119569951263209338965375875747241D-01,
     -      0.446410052264522508460080116863518489D-02,
     -      0.598661006119535102952189120791516211D-02,
     -      0.969769515827914555657248027915457853D+00,
     -      0.968270867185850458462276565111977140D+00,
     -      0.588364717728768611579237262995369845D-01,
     -      0.282339601726270622906286674485433530D-01,
     -      0.609906567195781611116655535736528786D-02,
     -      0.462826336151628344382652276196379392D-02,
     -      0.258585197952659479565873342862277801D-01,
     -      0.988239976128234044425749131411634614D+00,
     -      0.574041957657484045498441022803137394D-02/
          data whts16/
     -      0.833543564879251463932355962815366172D-02,
     -      0.833571822870689721045968127156860249D-02,
     -      0.833616010380554738390408691519172149D-02,
     -      0.671087035218141377914110488407572244D-02,
     -      0.669959946586777032676551211754413124D-02,
     -      0.670871709642214336375940608474668652D-02,
     -      0.634948376070186813479572123630511003D-02,
     -      0.635268304031491589370853554861179226D-02,
     -      0.634647185357232077133865985209009272D-02,
     -      0.635181249249814242892031795292014390D-02,
     -      0.428165119949283277010087093781073727D-02,
     -      0.422202232490998637495728144620809755D-02,
     -      0.634778348531546779745432585508973398D-02,
     -      0.634951267513164769003766176048340372D-02,
     -      0.422642339558080654601589967624231413D-02,
     -      0.427805571097472327257310666979908530D-02,
     -      0.428333001472544586774899200922055888D-02,
     -      0.422628099746630357215619276748588928D-02,
     -      0.396229061169888000120684999242402148D-02,
     -      0.686463339110859910120088065362584293D-02,
     -      0.717663611447148279230459923807272151D-02,
     -      0.509015849891619466539716595367215603D-02,
     -      0.686223898740277065992493168831595826D-02,
     -      0.686273300910811802169336272306256605D-02,
     -      0.718303033372546953697766007687574368D-02,
     -      0.508760182929440736667384093669800088D-02,
     -      0.686472446391627854427555516302513847D-02,
     -      0.327475927063665676223969703436310710D-02,
     -      0.326835077302642517638686690480074000D-02,
     -      0.396123261679317611599680404263463936D-02,
     -      0.508806477119148878570286611696198677D-02,
     -      0.395325556879133573679459879801594935D-02,
     -      0.509102183232654039277118799799541890D-02,
     -      0.718509027657585643063467730130883603D-02,
     -      0.686517322792251751904309341076769165D-02,
     -      0.686170217205470684316759257706795371D-02,
     -      0.508996395191803560970095269701199089D-02,
     -      0.508873287805668622899586021094981072D-02,
     -      0.327319666953132751206977596216411102D-02,
     -      0.531453867049914881342555710164802843D-02,
     -      0.680170722521946089093751148280223119D-02,
     -      0.346960632301519628634803327869107779D-02,
     -      0.680333361291858039242917720250734287D-02,
     -      0.326899820613170288025725976566851806D-02,
     -      0.531300031095839902557365575813037813D-02,
     -      0.680540039454236605793833948988597230D-02,
     -      0.345285586209528503198810806438781584D-02,
     -      0.680142764083496140840849008043703527D-02,
     -      0.130577703307857034670747597673899822D-02,
     -      0.531678862437533747388155628073366155D-02,
     -      0.125163646759674717521859181700352195D-02,
     -      0.326729022672028303230887755073656770D-02,
     -      0.126222059976477077010586559438559644D-02,
     -      0.347158088200464352220120811993034854D-02,
     -      0.531376187324818505323497155634793245D-02,
     -      0.327447242122590558056495235589383616D-02,
     -      0.345577389843253203194014150235912284D-02,
     -      0.287235346310599641255117148679403505D-02,
     -      0.531510805954690544110098001449818049D-02,
     -      0.462070022010985670410865596874632753D-02,
     -      0.680370056738623135069735449636356628D-02,
     -      0.680390396382448614455763251765908005D-02,
     -      0.130749587780054861449525894098189773D-02,
     -      0.531452909747575770888054438314872082D-02,
     -      0.345665583326944610488246671203906535D-02,
     -      0.130462486759994760973621046008011366D-02,
     -      0.126181288963053095133907105614160226D-02,
     -      0.544920131510599526030878557998158200D-02,
     -      0.462599690228259155289523313923118017D-02,
     -      0.288789285448278073330396283878024070D-02,
     -      0.125156700490553983126110891039815467D-02,
     -      0.130548303629390221193892504071146948D-02,
     -      0.346791643377580626344780671427811438D-02,
     -      0.130574493979119991150091680776535378D-02,
     -      0.288367511822773909650732993030050463D-02,
     -      0.130391454039544023971083642156913296D-02,
     -      0.462283116938616980158380296619551978D-02,
     -      0.544939477661376274162825076066612329D-02,
     -      0.125274899101441211532673883171680537D-02,
     -      0.462153119726970214524545774716525373D-02,
     -      0.126158603376555510050731656158950978D-02,
     -      0.287378125411675186635252358416306760D-02,
     -      0.130741911688867713587289517642419239D-02,
     -      0.287750473990907676145949350543849288D-02,
     -      0.462397861948027216791662469251848722D-02,
     -      0.130880430766876434863321600631261302D-02,
     -      0.544941952223871253890909282428769782D-02,
     -      0.462316590644250306758600517788527285D-02,
     -      0.288390139264970407002077209974125458D-02,
     -      0.130410120798403522648127989030946839D-02,
     -      0.226343068798326831796074185280855672D-02,
     -      0.119342717677821346498637157605538147D-02,
     -      0.357928081580724932708627877996240838D-02,
     -      0.220600676631921821524096443953072252D-02,
     -      0.358132870825846711940260268357300738D-02,
     -      0.118089163539464511026984290209958781D-02,
     -      0.118247686678762973822946129448926548D-02,
     -      0.221411559083155538408166719322490906D-02,
     -      0.357924254285960197032977310718325785D-02,
     -      0.225400064117617033768471444300703778D-02,
     -      0.357911267676336751223647297689124717D-02,
     -      0.119114636516610104965261822932288583D-02,
     -      0.119045857245308578169541530516664555D-02,
     -      0.224983715210105549482451085396824842D-02,
     -      0.358134424759141061886964572111711169D-02,
     -      0.222029814161669649524491553051072711D-02,
     -      0.358071525107892680509314520670055525D-02,
     -      0.118404950206102163732911715573532533D-02,
     -      0.122928151077421916823471295648241243D-02,
     -      0.123128672190215494982597564383948895D-02,
     -      0.205019056938324444720954138117396134D-02,
     -      0.123130987700785108053277526570667544D-02,
     -      0.208849648544066768046760322265337673D-02,
     -      0.204003407968779081185046134374210218D-02,
     -      0.204197949827597453492282506653954262D-02,
     -      0.123029911102204980664611447263441122D-02,
     -      0.208798880280743345219283747150131167D-02,
     -      0.205215475053813024853952664628748986D-02,
     -      0.843108319826394253803392674949080238D-03,
     -      0.123039962091063677676489560524557216D-02,
     -      0.858359862272573165917749116729682889D-03,
     -      0.204892283198645386965658114235965143D-02,
     -      0.156255369766936823528508703609842680D-02,
     -      0.158663408014799327118791029351373953D-02,
     -      0.203858818624715016093160090144135829D-02,
     -      0.123057780670162411014458500543014896D-02,
     -      0.208868739039593229330447510556347046D-02,
     -      0.856795718696459034876874689766058275D-03,
     -      0.846720118117551053724768413666574640D-03,
     -      0.158239692935673180491933152128507506D-02,
     -      0.156756231801493427604768446669834708D-02,
     -      0.617353631599068808146037193560715795D-03,
     -      0.847080877013730944837416464695472415D-03,
     -      0.118787678878234344887869569260825173D-02,
     -      0.635388333474683536434533053858652919D-03,
     -      0.156871919276858974826051859134328140D-02,
     -      0.854448082604863543762843295361272144D-03,
     -      0.157999273558868639780585839083763155D-02,
     -      0.633458558380327694661465035567304836D-03,
     -      0.620050578858709493961265128232825681D-03,
     -      0.118775643052277432238221484134673260D-02,
     -      0.329989903996933139639348607001511169D-03,
     -      0.308967349789859722967613560318850055D-03,
     -      0.621557914798468051320691686320547847D-03,
     -      0.310741821158802300498984635254799492D-03,
     -      0.326544041932225886544605194587578351D-03,
     -      0.632594187932047180411001277034866090D-03,
     -      0.118830140001896349723409279015990652D-02,
     -      0.209005147070237509388146123600904928D-03,
     -      0.325269839938683099078898923496166297D-03,
     -      0.312427828812517784325925863181397236D-03,
     -      0.208918222836621349025701046319822640D-03,
     -      0.208846915881281042691343595084920636D-03/

        ier = 0
c
c       Return the 4th order quadrature.
c
        if (norder .eq. 4) then
        do 1004 i = 1,nquad4
        xsdisc(i)   = xs4(i)
        ysdisc(i)   = ys4(i)
        whtsdisc(i) = whts4(i)
 1004   continue
        return
        endif
c
c       Return the 8th order quadrature.
c
        if (norder .eq. 8) then
        do 1008 i = 1,nquad8
        xsdisc(i)   = xs8(i)
        ysdisc(i)   = ys8(i)
        whtsdisc(i) = whts8(i)
 1008  continue
        return
        endif
c
c       Return the 12th order quadrature.
c
        if (norder .eq. 12) then
        do 1012 i = 1,nquad12
        xsdisc(i)   = xs12(i)
        ysdisc(i)   = ys12(i)
        whtsdisc(i) = whts12(i)
 1012   continue
        return
        endif
c
c
c       Return the 16th order quadrature.
c
        if (norder .eq. 16) then
        do 1016 i = 1,nquad16
        xsdisc(i)   = xs16(i)
        ysdisc(i)   = ys16(i)
        whtsdisc(i) = whts16(i)
 1016 continue
        return
        endif

        ier = 1024
        end
