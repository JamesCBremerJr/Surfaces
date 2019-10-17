c
c       !!!!!111111FIX KERNEL DOCUMENTATION FOR WEAKLY SINGULAR ONLY!!!!!!!!!!!
c
c       REWRITE KERNELS WITH STARS USING PRIMES FOR TRANSPOSE!!!!!!!!!!!!!
c
c       KERNELS ARE NOT ADJOINTS BUT RATHER TRANSPOSES!!!!!!!!
c
c
        implicit double precision (a-h,o-z)
c
        dimension errl2s(100),ns(100),times(100),entries(100)
c
        double precision, allocatable :: disc(:),eval(:,:),zs(:,:)
c
        double complex, allocatable   :: x(:),y(:),amatr(:,:),vals(:,:)
        double complex, allocatable   :: amatr2(:,:)
        double precision, allocatable :: tris(:,:),w(:)
        integer, allocatable          :: itars(:),isrcs(:)
c
        double complex val,sum,zk,val0
        external surface
c
        pi = acos(-1.0d0)
c
        norder  = 6
        nself   = 0
        dnear   = 4.0d0
        alpha   = 0.5d0
        call surface_init(alpha)
c
        ldisc   = 200 000 000
        maxtris = 1 000 000
c
        lw = 400 000 000
        allocate(w(lw))
c
        allocate(disc(ldisc),tris(6,maxtris))
c
        niters = 4
        do 1000 iter=1,niters
        call elapsed(t1)
c
        call disc_init(ier,norder,nself,dnear,disc,ldisc)
        if (ier .ne. 0) then
        call prinf("after disc_init, ier = *",ier,1)
        stop
        endif
c
c       Add a decomposition of the torus 
c
        x1 = 0
        x2 = 2*pi
        y1 = 0
        y2 = 2*pi
c
        iparam  = 1
        nlevels = iter-1
c
        call triangulate(ier,iparam,x1,y1,x2,y2,nlevels,
     -    maxtris,ntris,tris)
        call disc_add(ier,disc,iparam,ntris,tris)
c
c       Fetch and display information about the decomposition stored
c       in the disc structure.
c
        call disc_info(disc,ntris,nquad,n,dnear,lkeep)
c
        call prinf("total # of triangles ntris = *",ntris,1)
        call prinf("total # of nodes n = *",n,1)
c
        dd = lkeep*8
        dd = dd/(1024d0**2)
c
        call prin2("disc structure size (in MB) = *",dd,1)
c
c       Fetch evaluation data and plot the discretization nodes.
c
        if (allocated(eval)) deallocate(eval)
        allocate(eval(13,n))
        call disc_data(disc,eval)
c
c       Plot the discretization nodes.
c
        if (allocated(zs)) deallocate(zs)
        allocate(zs(3,n))
        call disc_points(disc,zs)
c
        iplot=1
        call plot3d_points("discretization nodes *",iplot,n,zs)
c
        if( allocated(itars)) deallocate(itars,isrcs)
        allocate(itars(n),isrcs(n))
c
        ntars    = 10
        nsrcs    = n
c
        do 1100 i=1,n
        isrcs(i) = i
 1100 continue
c
        do 1200 i=1,n
        itars(i) = i
 1200 continue
        call corrand_integer_knuth(n,itars)
        call prinf("ntars = *",itars,ntars)
c
        if (allocated(amatr)) deallocate(amatr)
        allocate(amatr(ntars,nsrcs))
c
c       Double layer potential for Laplace's equation.
c
        ikernel  = 2
        zk       = 0.0d0
c
        call elapsed(t1)
        call disc_eval(ier,disc,ikernel,zk,ntars,itars,
     -    nsrcs,isrcs,amatr,w,lw)
        if (ier .ne. 0 ) then
        call prinf("after disc_eval, ier = *",ier,1)
        stop
        endif
c
        call elapsed(t2)
c
        tavg = (t2-t1)/(ntars*nsrcs)
        call prin2("average evaluation time = *",tavg,1)
c
        ents = 1/tavg
        call prin2("matrix entries per second = *",ents,1)
c
        entries(iter) = ents
c$$$c
c$$$c       Form a random submatrix and compare.
c$$$c
c$$$        nsrcs = min(nsrcs,100)
c$$$        ntars = min(ntars,100)
c$$$c
c$$$        if (allocated(amatr2)) deallocate(amatr2)
c$$$        allocate(amatr2(nsrcs,ntars))
c$$$        call corrand_integer_knuth(nsrcs,isrcs)
c$$$        call corrand_integer_knuth(ntars,itars)
c$$$c
c$$$        call disc_eval(ier,disc,ikernel,zk,ntars,itars,
c$$$     -    nsrcs,isrcs,amatr2,w,lw)
c$$$        if (ier .ne. 0 ) then
c$$$        call prinf("after disc_eval, ier = *",ier,1)
c$$$        stop
c$$$        endif
c$$$c
c$$$        errl2 = 0
c$$$        do 1300 ii=1,ntars
c$$$        do 1400 jj=1,nsrcs
c$$$        i = itars(ii)
c$$$        j = isrcs(jj)
c$$$        errabs = abs (amatr(i,j)-amatr2(ii,jj) )
c$$$        errl2  = errl2 + errabs**2
c$$$ 1400 continue
c$$$ 1300 continue
c$$$        errl2 = sqrt(errl2)
c$$$        call prin2("submatrix errl2 = *",errl2,1)
c
c        if (errl2 .gt. 0) stop
c
c       Form the vector to which the matrix will be applied and choose a
c       set of random targets.
c
        errmax = 0
        errl2  = 0
        dnorm  = 0
        do 3000 i=1,ntars
        sum = 0
        do 3100 j=1,n
        sum = sum + amatr(i,j)*eval(1,j)
 3100 continue
        errabs = abs(sum - eval(1,itars(i))/2)
        errmax = max(errmax,errabs)
        errl2  = errl2 + errabs**2
 3000 continue
        errl2 = sqrt(errl2)
        call prin2("errl2  = *",errl2,1)
        call prin2("errmax = *",errmax,1)
c
        call elapsed(t2)
        call prin2("time = *",t2-t1,1)
c
        call prina("*")
        call prina("*")
c
        errl2s(iter) = errl2
        ns(iter)     = n
        times(iter)  = t2-t1 
 1000 continue
c      
 5000 format(I6.6," ",D10.3,"   ",D10.3,"  (",F12.2,")    ",F12.2)
c
c        call prina("*")
c        call prina("*")
c        call prina("*")
c
        call prina("-------------------------------------------------*")
c
        do 5100 i=1,niters
        dratio = 0
        if (i .gt. 1) dratio = errl2s(i-1)/errl2s(i)
c
        write(*,5000)  ns(i),times(i),errl2s(i),dratio,entries(i)
        write(13,5000) ns(i),times(i),errl2s(i),dratio,entries(i)
 5100 continue
c
        call prina("*")
        call prina("*")

        end


        subroutine surface(iparam,s,t,r,dr)
        implicit double precision (a-h,o-z)
        dimension r(3),dr(3,2)        
        save
c
c       Supply a parameterization of a torus.
c
        a = 2.0d0
        b = alpha
c
        r(1)    = (a+b*cos(s))*sin(t)
        r(2)    = (a+b*cos(s))*cos(t)
        r(3)    = b*sin(s)
c
        dr(1,1) = -b*sin(s)*sin(t)
        dr(2,1) = -b*sin(s)*cos(t)
        dr(3,1) = b*cos(s)
c
        dr(1,2) = (a+b*cos(s))*cos(t)
        dr(2,2) = -(a+b*cos(s))*sin(t)
        dr(3,2) = 0
c
        return
c
        entry surface_init(alpha0)
        alpha=alpha0
        end



        subroutine sphere_data(n,x0,y0,z0,r,nn,eval)
        implicit double precision (a-h,o-z)
        dimension eval(13,1)
        data pi    / 3.14159265358979323846264338327950288d0 /
c
        dimension xslege(n),whtslege(n)
c
c       Return evaluation data for a user-specified sphere.
c
c
c       Fetch the n-point Legendre quadrature rule.
c
        call legequad(n,xslege,whtslege)
c
c       Construct the quadrature for spherical harmonics.
c
        nn = 0
c
        dd = n
        dd = 1.0d0/dd
c
        do 1000 i=1,n
        do 1100 j=1,n
        theta  = (j-1)*2*pi*dd
        phi    = acos(xslege(i))
c
        nn  = nn + 1
        eval(1,nn)  = sqrt(2*pi*r**2*dd*whtslege(i))
c
        eval(2,nn)  = x0+r*cos(theta)*sin(phi)
        eval(3,nn)  = y0+r*sin(theta)*sin(phi)
        eval(4,nn)  = z0+r*cos(phi)
c
        eval(5,nn)  = -r*sin(theta)*sin(phi)
        eval(6,nn)  =  r*cos(theta)*sin(phi)
        eval(7,nn)  =  0
c
        eval(8,nn)  = r*cos(theta)*cos(phi)
        eval(9,nn)  = r*sin(theta)*cos(phi)
        eval(10,nn) = -r*sin(phi)
c
        eval(11,nn) = cos(theta)*sin(phi)
        eval(12,nn) = sin(theta)*sin(phi)
        eval(13,nn) = cos(phi)
c
 1100 continue
 1000 continue
c
        end



        subroutine triangulate(ier,iparam,x1,y1,x2,y2,nlevels,
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
        call ttsplit_tri(tris(1,itri),children)
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





        subroutine ttsplit_tri(tri,children)
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



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the 
c       discretization code proper.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains code for discretizing the L^2 action of
c       certain singular integral operators given on surfaces.  The
c       discretization procedure is described in some detail in the
c       paper 
c
c         ``A Nystr\"om method for the discretization of weakly singular
c           integral operators on surfaces''
c          Journal of Computational Physics, 231 (2012), pg. 4885 - 4903
c
c       In brief, the user selects a discretization quadrature and 
c       provides a decomposition of a surface \Sigma.  Together, these 
c       give rise to a scheme for discretizing integral operators
c
c         Tf(x) = \int   K(x,y) \sigma(y) dS(y)
c                  \Sigma \ B_eps(x)
c
c       whose kernels K(x,y) are linear combinations of the classical 
c       layer potentials:
c
c                      exp(i k |x-y|)
c         S_k(x,y)   = --------------
c                           |x-y|
c  
c         D_k(x,y)   = \nabla_y S(x,y) \cdot \eta_y
c                                                                         (1)
c         D_k'(x,y)  = \nabla_x S(x,y) \cdot \eta_x
c
c       Here:
c
c         \nabla_p denotes the gradient with respect to the variable p;
c
c         \eta_p denotes the outward-pointing normal unit vector to 
c         the surface at the point p.
c
c       A decomposition of a surface \Sigma is a collection 
c
c         {\rho_j: T_j \to \mathbb{R}^3}
c
c       of smooth parameterizations each given on a triangle T_j in the
c       plane and such that \rho_j(T_j) form a disjoint union of \Sigma.
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
c       Such a structure is constructed by the user by first calling the
c       disc_init routine (below).  Information about the T_j is
c       added to the disc structure through calls to the routine
c       disc_add (below).
c
c       A discretization quadrature of order N is a quadrature rule on 
c       the simplex which integrates polynomials of order 2*N.  The
c       subroutine discquad in this file makes available a collection 
c       of discretization quadratures of various orders.

c       Associated with a discretization quadrature 
c
c          x_1,...,x_m,w_1,...,w_m                                        (2)
c
c       and a decomposition
c
c         D = {\rho_j: T_j \to \mathbb{R}^3}_{j=1}^k                      (3)
c
c       for a surface \Sigma is a scheme for discretizing integral 
c       operators whose kernels are linear combinations of those 
c       appearing in (1).
c
c       More specifically, we define for each j the mapping
c
c                       ( f(\rho_j(x_1)) \sqrt{|dp(x_1)|} \sqrt{w_1} )
c                       ( f(\rho_j(x_2)) \sqrt{|dp(x_2)|} \sqrt{w_2} )
c         \Phi_j(f) =   (                    .                       )
c                       (                    .                       )
c                       ( f(\rho_j(x_m)) \sqrt{|dp(x_m)|} \sqrt{w_m} )
c
c       which takes functions defined on \rho_j(T_j) to \mathbb{C}^m.
c
c       We now let \Phi be the amalgamated mapping which takes functions
c       on \Sigma to \mathbb{C}^{mk} and which maps the restriction f_j 
c       of f: \Sigma \to \mathbb{R}^3 to \rho(T_j) to \Phi_j(f_j).
c
c       The discretization A of the operator T associated with (2) and
c       (3) is the matrix A which makes the diagram
c
c                                   T
c                L^2(\Sigma) --------------->  L^2(\Sigma)
c                     |                             |
c                \Phi |                             | \Phi                (4)
c                     |                             |
c                    \_/            A              \_/
c                  C^{mk}     --------------->    C^{mk}
c       
c       commute.  
c
c       Once a full description of a surface has been input into a disc
c       structure, the user can call the routine disc_eval (see below)
c       to construct any desired submatrix of the matrix A discretizing
c       the operator T.
c
c       Throughout this file, we refer to the notion of "evaluation
c       data."  In order to evaluate many of the integral kernels (1)
c       K(x,y), first order derivative information at the points
c       x and y is required.  See disc_init for a discussion of how
c       evaluation data is organized.
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
c   disc_eval - evaluate a user-specified submatrix of the matrix
c        discretizing an integral operator with one of the kernels (1)
c
c   disc_eval_far - evaluate a specified entry of the matrix discretizing
c        an integral operator with one of the kernels (1); this routine
c        is only capable of handling "far" interactions, but it is 
c        considerably faster than disc_eval in that case
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine disc_init(ier,norder,nself,dnear,disc,ldisc)
        implicit double precision (a-h,o-z)
        dimension disc(1)
c
c       This subroutine intializes a discretization structure which
c       contains the necessary quadrature data and which will contain
c       a description of the decomposition.
c
c                            Input Parameters:
c
c   norder - order of the discretization quadrature used to represent 
c       solutions; at the present time, the possible orders are 4, 8,
c       12 and 16
c
c   nself - an integer parameter indicating which quadrature rules
c       should be used for self-interactions
c
c       nself = 0  means use a precomputed table of quadrature rules
c       nself > 0  means use piecewise double exponential quadrature
c                  rules of order nself in lieu of the precomputed
c                  tables
c           
c       NOTE: it is highly recommend that nself be set to zero
c
c   dnear - this parameter is the "near" interaction radius which
c       controls when a point is considered to be near a surface
c       region; typically, this parameter should be set to the value
c       2.0d0
c
c       EXPLANATION:  Let \Gamma be a surface region and p a point
c       outside of \Gamma.  To determine if p is near to \Gamma,
c       an approximate bounding ball B_r(x) for \Gamma is found 
c       and p is said to be near to \Gamma if it is contained in
c       B_{dnear *r} (x).
c
c   ldisc - the length of the user-supplied array disc
c
c                            Output Parameters:
c
c   ier - an error return code;
c       ier = 0     indicates successful execution
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
c   disc(7)  - ratio which defines near interactions                      (dnear)
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
c
c   disc(40) - a pointer to the rad structure which contains the tables   (irad)
c              used to construct self interaction quadratures
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
c   stored in a (13,nquad) array whose jth entry which contains the following
c    
c   eval(1,j)  - the *SQUARE ROOT* of the quadrature weight at this point (w)
c   eval(2,j)  - the x-coordinate of the point                            (x)
c   eval(3,j)  - the y-coordinate of the point                            (y)
c   eval(4,j)  - the z-coordinate of the point                            (z)
c   eval(5,j)  - the derivative of the x coordinate w.r.t. the first      (dx/ds)
c                parameterization variable
c   eval(6,j)  - the derivative of the y coordinate w.r.t. the first      (dy/ds)
c                parameterization variable
c   eval(7,j)  - the derivative of the z coordinate w.r.t. the first      (dz/ds)
c                parameterization variable
c   eval(8,j)  - the derivative of the x coordinate w.r.t. the second     (dx/dt)
c                parameterization variable
c   eval(9,j)  - the derivative of the y coordinate w.r.t. the second     (dy/dt)
c                parameterization variable
c   eval(10,j) - the derivative of the z coordinate w.r.t. the second     (dz/dt)
c                parameterization variable
c   eval(11,j) - the x component of the outward pointing unit normal      (dnx)
c   eval(12,j) - the y component of the outward pointing unit normal      (dny)
c   eval(13,j) - the z component of the outward pointing unit normal      (dnz)
c
        ier = 0
c
        if (ldisc .lt. 1 000 000) then
        ier = 4
        return
        endif
c
c       Allocate memory for the discretization quadrature and fetch it.
c
        ixs       = 100
        lxs       = 5000
c
        iys       = ixs+lxs
        lys       = 5000
c
        iwhts     = iys+lys
        lwhts     = 5000
c
c       Fetch the quadrature formula.
c
        npoly = (norder+1)*(norder+2)/2
        call discquad(ier,norder,nquad,disc(ixs),disc(iys),disc(iwhts))
        if (ier .ne. 0) return
c
        call prinf("in disc_init, norder = *",norder,1)
        call prinf("in disc_init, npoly = *",npoly,1)
        call prinf("in disc_init, nquad = *",nquad,1)
c
        iu        = iwhts+lwhts
        lu        = npoly*nquad
c
        iv        = iu+lu
        lv        = npoly*nquad
c
        lsofar    = iv+lv
        if (lsofar .gt. ldisc) then
        ier = 4
        return
        endif
c
        call disc_init0(norder,npoly,nquad,disc(ixs),disc(iys),
     -   disc(iwhts),disc(iu),disc(iv))
c
        iamatrin  = iv+lv
        lamatrin  = 4*nquad*nquad
c
        irad      = iamatrin+lamatrin
        lrad      = ldisc-irad
c
        if (lrad .le. 1 000 000) then 
        ier = 4
        return
        endif
c
c       Build the interpolation matrix for near evaluations.
c
        call disc_amatrin(norder,npoly,nquad,disc(ixs),disc(iys),
     -   disc(iwhts),disc(iu),disc(iamatrin))
c
c       Initialize the self-interaction quadratures.
c
        if (nself .eq. 0) then
        norder0 = norder
        if (norder .eq. 6)  norder0 = 8
        if (norder .eq. 10) norder0 = 12
        call radinit(ier,norder0,disc(irad),lrad,lkeep)
        if (ier .ne. 0) return
        else
        lkeep = 0
        endif
c
        lrad    = lkeep
c
        itris   = irad+lrad
        ltris   = ldisc-itris
c
        len     = 11+13*nquad
        maxtris = ltris/len
        ntris   = 0
c
        call prinf("in disc_init, maxtris = *",maxtris,1)
c
        disc(1)  = norder       
        disc(2)  = nself
        disc(3)  = maxtris
        disc(4)  = ntris
        disc(5)  = len
        disc(6)  = itris
        disc(7)  = dnear
c
        disc(20) = nquad
        disc(21) = ixs
        disc(22) = iys
        disc(23) = iwhts
        disc(24) = iu
        disc(25) = iv
c
        disc(34) = iamatrin
        disc(40) = irad
c
        end


        subroutine disc_init0(norder,npoly,nquad,xs,ys,whts,u,v)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1),u(npoly,nquad),v(nquad,npoly)
c
        dimension vals0(npoly)
c
c       Explicity construct v.
c
        do 1000 i=1,nquad
c
        x   = xs(i)
        y   = ys(i)
        wht = whts(i)
c
        call koorn(norder,x,y,vals0)
c
        do 1100 j=1,npoly
        v(i,j) = vals0(j)*sqrt(wht)
 1100 continue
 1000 continue
c
c       Construct u.
c
        do 2000 i=1,npoly
        do 2100 j=1,nquad
        u(i,j) = v(j,i)
 2100 continue       
 2000 continue
c
        end



        subroutine disc_amatrin(norder,npoly,nquad,xs,ys,whts,u,amatrin)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1)
        dimension amatrin(4*nquad,nquad),u(npoly,nquad)
        dimension simp(6),tris(6,4),amap(2,3)
        dimension vals0(npoly),vals(4*nquad,npoly)
c
c       Split the standard simplex.
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
        call koorn(norder,x,y,vals0)
        idx = (itri-1)*nquad+j
        do 1200 i=1,npoly
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
        do 2200 l=1,npoly
        sum = sum + vals(i,l)*u(l,j)
 2200 continue
        amatrin(i,j) = sum
 2100 continue
 2000 continue
c
        end



        subroutine disc_eval(ier,disc,ikernel,zk,ntars,itars,
     -    nsrcs,isrcs,amatr,w,lw)
        implicit double precision (a-h,o-z)
        dimension disc(1),itars(1),isrcs(1),w(1)
        double complex val,zk,amatr(ntars,nsrcs)
c
        external ksingle0,kdouble0,ksingle0prime
        external ksingle,kdouble,ksingleprime
c
c       Evaluate a user-specified submatrix of the matrix A in the
c       diagram (4).
c
c                            Input Parameters:
c
c   disc - the discretization structure describing a decomposition
c       of the surface and containing a discretization quadrature
c
c   ikernel - integer parameter indicating which kernel is to be
c        evaluated:
c
c      ikernel = 1  indicates the kernel 
c
c                         exp(i zk |x-y|)
c            K(x,y) =   ---------------
c                          4\pi |x-y|
c
c        of the single layer operator S_zk for the Helmholtz equation
c
c      ikernel = 2  means the kernel 
c
c                        (y-x) \cdot n_y
c            K(x,y) =   -----------------  exp( i zk |x-y| ) (1- i zk |x-y|)
c                           4\pi |x-y|^3
c
c        of the double layer operator D_zk for the Helmholtz equation
c
c      ikernel = 3 means the kernel 
c
c                        (x-y) \cdot n_x
c            K(x,y) =   -----------------  exp( i zk |x-y| ) (1- i zk |x-y|)
c                           4\pi |x-y|^3
c
c        of the transpose D_zk^' of the double layer operator for
c        the Helmholtz equation
c
c   zk - the complex-valued wavenumber for the problem;
c
c        NOTE: when this value is set to 0, the kernels above 
c        simplify to those for Laplace's equation
c
c   ntars - the number of target discretization nodes
c   itars - an integer array specifying the indices of the target 
c       discretization nodes
c
c   nsrcs - the number of source discretization nodes
c   isrcs - an integer array specifying the indices of the source
c       discretization nodes
c
c   w - a user-supplied work array
c   lw - the length of the array w
c
c                             Output Parameters:
c
c   ier - an error return code;
c
c       ier = 0     indicates successful execution
c       ier = 4     means that the work array was of insufficient length
c       ier = 16    means that the preset limits on the number of 
c                   subdivisions allowed in the near interaction 
c                   computations was exceeded
c       ier = 1024  means that a problem was encounterd
c
c   amatr - the (ntars,nsrcs) matrix whose (i,j)th entry is the
c       ( itars(i), isrcs(j) )th entry of the matrix A
c
        ier = 0
c     
c       Fetch data from the disc array.
c
        norder   = disc(1)
        nself    = disc(2)
        maxtris  = disc(3)
        ntris    = disc(4)
        len      = disc(5)
        itris    = disc(6)
        dnear    = disc(7)
c
        nquad    = disc(20)
        ixs      = disc(21)
        iys      = disc(22)
        iwhts    = disc(23)
        iu       = disc(24)
        iv       = disc(25)
c
        irad     = disc(40)
        npoly    = (norder+1)*(norder+2)/2
c
c       Allocate memory from the work array.
c
        ijtris  = 1
        ljtris  = nsrcs
c
        iiorder = ijtris+ljtris
        liorder = nsrcs
c
        ivals   = iiorder+liorder
        lvals   = nquad*ntars*2
c
        iw2     = ivals+lvals
        lw2     = lw-iw2
c
        if (lw2 .le. 0) then
        ier = 4
        return
        endif
c
        call disc_eval00(ier,disc,nquad,ikernel,zk,ntars,itars,
     -    nsrcs,isrcs,amatr,w(ijtris),w(ivals),w(iiorder),w(iw2),
     -    lw2)
c
c
        end


        subroutine disc_eval00(ier,disc,nquad,ikernel,zk,ntars,itars,
     -    nsrcs,isrcs,amatr,jtris,vals,iorder,w,lw)
        implicit double precision (a-h,o-z)
c
        dimension itars(1),isrcs(1)
        dimension jtris(nsrcs),iorder(nsrcs),w(1)
        double complex vals(nquad,ntars),amatr(ntars,nsrcs)
c
c       This is an auxilliary subroutine for disc_eval.  Its purpose
c       is to group evaluations by source triangle, call another
c       auxilliary routine to perform evaluation, and then return only 
c       the requested values.
c
c
c                           Input Parameters:
c
c  disc - the discretization structure
c  nquad - the number of quadrature nodes per source triangle
c  ikernel - the index of the kernel of interest
c  zk - the wavenumber for the kernel
c  
c
c
c                           Output Parameters:
c
c
c       Record the index of the triangle containing each source
c       point, reorder them by index triangle, and keep track of the
c       original ordering.
c
        do 1000 i=1,nsrcs
        j         = isrcs(i)
        jtri      = (j-1)/nquad
        jpt       = j-(jtri*nquad)
        jtri      = jtri+1
        jtris(i)  = jtri
        iorder(i) = i
 1000 continue
c
        call insorti2(nsrcs,jtris,iorder)
c
c       Process each source triangle, evaluating the interactions
c       with all targets in each pass.
c
        idx1 = 0
        idx2 = 0
        jtri = 0
c
 2000 continue
c
c       Find the extent of the points containing the first triangle
c       idx1:idx2
c
 2100 continue
        idx1 = idx2 + 1
        if (idx1 .gt. nsrcs) goto 3000
        jtri = jtris(idx1)
c
        idx2 = idx1
 2200 continue
        if (idx2 .ge. nsrcs) goto 2300
c
        if (jtris(idx2+1) .eq. jtri) then
        idx2 = idx2+1
        goto 2200
        endif
c
 2300 continue
c
c       Evaluate the interactions for jtri.
c
        call disc_eval2(ier,disc,ikernel,zk,jtri,ntars,itars,
     -    vals,w,lw)
        if (ier .ne. 0) return
c
 2400 continue
c
c       Copy out the values.
c       
        do 1500 jj=idx1,idx2
        idx         = iorder(jj)
        j           = isrcs(idx)-(jtri-1)*nquad
        do 1400 i=1,ntars
        amatr(i,idx) = vals(j,i)
 1400 continue
 1500 continue
c
        goto 2000
c
c       We are done processing triangles.
c      
 3000 continue
c
        end


        subroutine disc_eval_self(ier,disc,i,jtri,vals,ikernel,zk,w,lw)
        implicit double precision (a-h,o-z)
        dimension disc(1),w(lw)
c
        double complex vals(1),zk
c
c       This is an auxilliary subroutine for disc_eval which handles
c       "self" evaluations -- that is, the interaction of a triangle
c       with itself.
c     
c                             Input Parameters:
c
c   disc - the discretization structure describing
c   i - the index of the target node
c   jtri - the index of the triangle
c   
c   w - a user-supplied work array 
c   lw - the length of the user-specified work array
c
c                           Output Parameters:
c
c   vals - the 
c
        ier = 0
c
c       Fetch data from the discretization structure.
c
        norder   = disc(1)
        nself    = disc(2)
        maxtris  = disc(3)
        ntris    = disc(4)
        len      = disc(5)
        itris    = disc(6)
        dnear    = disc(7)
c
        nquad    = disc(20)
        ixs      = disc(21)
        iys      = disc(22)
        iwhts    = disc(23)
        iu       = disc(24)
        iv       = disc(25)
c
        irad     = disc(40)
        npoly    = (norder+1)*(norder+2)/2
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
c       Allocate memory from the work array.
c
        ivals0    = 1
        lvals0    = nquad
c
        ivals2    = ivals0+lvals0
        lvals2    = 2*(nquad+npoly)
c
        iw2       = ivals2+lvals2
        lw2       = lw-iw2
c
        maxquad   = (lw2-10)/5
c
        if (maxquad .lt. 2000) then
        ier = 4
        return
        endif
c
        ixsself   = iw2
        lxsself   = maxquad
c
        iysself   = ixsself+lxsself
        lysself   = maxquad
c
        iwhtsself = iysself+lysself
        lwhtsself = maxquad
c
        ivals1    = iwhtsself+lwhtsself
        lvals1    = 2*maxquad
c
c       Call an auxilliary routine to shape arrays.
c
        call disc_eval_self1(ier,norder,npoly,nself,nquad,disc(ixs),
     -    disc(iys),disc(iwhts),disc(iu),disc(iptr),len,disc(itris),
     -    jtri,ipt,maxquad,w(ixsself),w(iysself),w(iwhtsself),
     -    w(ivals0),w(ivals1),w(ivals2),disc(irad),vals,
     -    ikernel,zk)
c
        end



        subroutine disc_eval_self1(ier,norder,npoly,nself0,nquad,xs,ys,
     -    whts,u,x,len,tris,jtri,ipt,maxquad,xsself,ysself,whtsself,
     -    vals0,vals1,vals2,rad,vals,ikernel,zk)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1),u(npoly,nquad),x(13),y(13)
        dimension tris(len,1),xsself(maxquad),ysself(maxquad)
        dimension whtsself(maxquad)
        dimension amap(2,3),ainv(2,3),r(3),dr(3,2),dn(3)
        double complex vals1(maxquad),vals2(npoly+nquad)
        dimension vals0(maxquad)
        double complex val,sum,zk,vals(1)
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
        call selfquad(ier,rad,tris(1,jtri),s0,t0,x(5),
     -    maxquad,nself,xsself,ysself,whtsself)
        else
        call selfquad0(ier,nself0,norder,tris(1,jtri),s0,t0,
     -    maxquad,nself,xsself,ysself,whtsself)
        endif
c
        if (ier .ne. 0) return
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
        call eval_kernel(ikernel,zk,x,y,vals1(i))
 1000 continue
c
        do 1100 i=1,nquad
        vals2(i) = 0
 1100 continue
c
c       Form the product of the values vector with the matrix of 
c       polynomials values.
c
        do 1200 i=1,nself
        s   = ainv(1,1)*xsself(i)+ainv(1,2)*ysself(i)+ainv(1,3)
        t   = ainv(2,1)*xsself(i)+ainv(2,2)*ysself(i)+ainv(2,3)
        wht = whtsself(i)/det
c
        call koorn(norder,s,t,vals0)
c
        do 1300 j=1,npoly
        vals2(j) = vals2(j) + vals0(j)*sqrt(wht)*vals1(i)
 1300 continue
 1200 continue
c
c       Take the product vals2 = vals1 * amatr
c
        do 1400 j=1,nquad
        sum = 0
        do 1500 i=1,npoly
        sum = sum + vals2(i)*u(i,j)
 1500 continue
        vals(j) = sum
 1400 continue
        return
        end



        subroutine disc_eval2(ier,disc,ikernel,zk,jtri,ntars,itars,vals,
     -   w,lw)
        implicit double precision (a-h,o-z)
        dimension disc(1),itars(1),w(1)
        double complex vals(1),zk
c
c       This subroutine returns a submatrix of the matrix discretizing an
c       integral operator of the form (1).  The block is that associated
c       with a specified source triangle and a user-specified list of 
c       target nodes.
c
c                              Input Parameters:
c
c                              Output Parameters:
c
c   ier - an error return code
c       ier = 0    indicates successful execution
c       ier = 4    means that the work array w was of insufficient 
c                  length
c       ier = 16   means that one of the preset limits was exceeded 
c                  during the near interaction computation; 
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
        dnear      = disc(7)
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
        npoly      = (norder+1)*(norder+2)/2
c
c       Call an auxilliary routine to construct the tree structure.
c
        maxlevels  = 30
        maxtris0   = 500 000
c
        ilevels    = 1
        llevels    = 2*maxlevels
c
        itree      = ilevels+llevels
        ltree      = lw - itree
c
        maxtris    = ltree/(16+13*nquad)
c
        if (ltree .le. 0) then
        ier = 4
        return
        endif
c
        call disc_eval20(ier,disc,norder,npoly,nquad,disc(ixs),
     -   disc(iys),
     -   disc(iwhts),disc(iamatrin),ntars,itars,len,
     -   disc(itris),jtri,vals,ikernel,zk,dnear,maxlevels,maxtris,
     -   w(ilevels),w(itree),nlevels,nnodes,maxtris0)
c
        if (ier .ne. 0) return
c
        ltree = nnodes*(16+13*nquad)
c
        iistack      = itree+ltree
        listack      = nnodes
c
        ivals0       = iistack+listack
        lvals0       = 2*nquad*nnodes
c
        ivals1       = ivals0+lvals0
        lvals1       = 2*4*nquad
c
        iifprocessed = ivals1+lvals1
        lifprocessed = nnodes
c
        iw2          = iifprocessed+lifprocessed
        lw2          = lw-iw2
c
        if (lw2 .le. 0) then
        ier = 4
        return
        endif
c
        call disc_eval21(ier,disc,norder,npoly,nquad,disc(ixs),
     -   disc(iys),disc(iwhts),disc(iamatrin),ntars,itars,len,
     -   disc(itris),jtri,vals,ikernel,zk,dnear,maxlevels,maxtris,
     -   w(ilevels),w(ivals1),w(itree),nlevels,nnodes,
     -   w(iistack),w(ivals0),w(iifprocessed),w(iw2),lw2)
c
        if (ier .ne. 0) then
        print *,"after disc_eval21, ier = ",ier
        return
        endif
c
        end

c
        subroutine disc_eval21(ier,disc,norder,npoly,nquad,xs,ys,whts,
     -    amatrin,ntars,itars,len,tris,jtri,vals,ikernel,zk,dnear,
     -    maxlevels,maxtris,levels,vals1,tree,nlevels,nnodes,
     -    istack,vals0,ifprocessed,w,lw)
        implicit double precision (a-h,o-z)
        dimension disc(1),w(1)
        dimension xs(1),ys(1),whts(1),itars(1)
        dimension amatrin(4*nquad,nquad),tris(len,1)
c
        double complex vals(nquad,ntars)
c
        dimension children(6,4),xeval(13),verts(6),yeval(13)
        double complex val,zk,sum
c
        dimension tree(16+13*nquad,1),ifprocessed(nnodes)
        dimension levels(2,maxlevels),istack(nnodes)
        double complex vals1(4*nquad),vals0(nquad,nnodes)
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
        call disc_eval_self(ier,disc,i,jtri,vals(1,ii),ikernel,zk,w,lw)
        if (ier .ne. 0) return
        goto 2000
        endif
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
        dd = sqrt((bx-x)**2+(by-y)**2+(bz-z)**2)
c     
c       If the node is distant, evaluate the kernel at the
c       appropriate nodes and mark the node as processed.
c
        if (itri .eq. jtri .OR.
     -    (dd .gt. br*dnear .AND. level .ge. minlevels) ) then
        iptr = 16
c
        do 2400 i=1,nquad

        call eval_kernel(ikernel,zk,xeval,tree(iptr,inode),val)

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
c
 2000 continue
c
c
        end


        subroutine disc_eval20(ier,disc,norder,npoly,nquad,xs,ys,whts,
     -    amatrin,ntars,itars,len,tris,jtri,vals,ikernel,zk,dnear,
     -    maxlevels,maxtris,levels,tree,nlevels,nnodes,maxtris0)
        implicit double precision (a-h,o-z)
        dimension disc(1)
        dimension xs(1),ys(1),whts(1),itars(1)
        dimension amatrin(4*nquad,nquad),tris(len,1)
c
        double complex vals(nquad,ntars)
c
        dimension children(6,4),xeval(13),verts(6),yeval(13)
        double complex val,zk
c
        dimension tree(16+13*nquad,maxtris)
        dimension levels(2,maxlevels)
c
        ier =0 
c
        iparam = tris(7,jtri)
c
c       The tree structure is organized as follows:
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
        dd   = sqrt((x-bx)**2+(y-by)**2+(z-bz)**2)
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
           print *,nnodes,maxtris
        ier = 4 
        return
        endif
c
        if (nnodes+4 .ge. maxtris0) then
        ier = 16
        return
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
c
 1400 continue
 1100 continue
c
        newnodes = nnodes - levels(1,level+1)+1
        levels(2,level+1) = newnodes
        if (newnodes .eq. 0) goto 1999
 1000 continue
c
c       We get here if the maximum number of levels is exceeded.
c
        ier = 16
        return
 1999 continue
c
        nlevels = level
c
        end




        subroutine disc_eval1(nquad,xs,ys,whts,verts,tree,iparam)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1),verts(6)
        dimension tree(1),amap(2,3),dn(3),r(3),dr(3,2)
 
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
        bx = 0
        by = 0
        bz = 0
c
        do 1000 i=1,nquad
        s   = amap(1,1)*xs(i)+amap(1,2)*ys(i)+amap(1,3)
        t   = amap(2,1)*xs(i)+amap(2,2)*ys(i)+amap(2,3)
        wht = whts(i)*det
c
        call surface(iparam,s,t,r,dr)
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
c
 1000 continue
c
        call disc_bounding_ball(nquad,tree(16),bx,by,bz,br)
c
        tree(7)  = bx
        tree(8)  = by
        tree(9)  = bz
        tree(10) = br
c
        end



        subroutine disc_eval_far(disc,ikernel,zk,i,j,val)
        implicit double precision (a-h,o-z)
        dimension disc(1)
        double complex val
        external kernel
c
c       Evaluate an entry of a discretization matrix corresponding to
c       a single source point and a single target in the far regime.
c
c       WARNING: it is incumbent on the user to ensure that the requested
c       value of the matrix corresponds to a ``far-field'' interaction.
c       Erroneous results will be returned if this is not the case.
c
c
c                          Input Parameters:
c
c   disc - the discretization structure
c   (i,j) - entry of the matrix to evaluate
c
c   ikernel - integer parameter indicating which kernel is to be
c        evaluated; see the disc_eval routine for the meanings of
c        the possible values
c   zk - the complex-valued wavenumber for the problem;
c
c        NOTE: when this value is set to 0, the kernels above 
c        simplify to those for Laplace's equation
c
c                         Output Parameters:
c
c   val - the requested complex-valued entry of the matrix
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
        call eval_kernel(ikernel,zk,disc(iptr),disc(jptr),val)
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
c   ntris0 - the number of triangles to add the the decomposition
c   tris0 - a (2,3,ntris0) array which gives the coordinates of the
c        vertices of the triangles in the parameterization domain
c        corresponding to the surface elements to add to the 
c        decomposition
c
c                          Output Parameters:
c 
c   ier - an error return code;
c
c       ier = 0   indicates successful execution
c       ier = 4   means that the disc array is of insufficient length
c   
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
        dnear      = disc(7)
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
     -    len,disc(itris),iparam,ntris0,tris0)
c
c       Update the header.
c
        disc(4) = ntris
c
        end

       
        subroutine disc_add0(nquad,xs,ys,whts,ntris,len,tris,
     -    iparam,ntris0,tris0)
        implicit double precision (a-h,o-z)
        dimension xs(1),ys(1),whts(1),tris(len,1),dn(3)
        dimension tris0(6,1),r(3),dr(3,2),amap(2,3)
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
c       Construct the evaluation data for the surface element and
c       find a bounding ball.
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
c       Record the parameters for the bounding ball.
c
        call disc_bounding_ball(nquad,tris(12,itri),bx,by,bz,br)
c
        tris(8, itri) = bx
        tris(9, itri) = by
        tris(10,itri) = bz
        tris(11,itri) = br
c
 1000 continue
        ntris = ntris+ntris0
c
        end


        subroutine disc_info(disc,ntris,nquad,n,dnear,lkeep)
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
c   n - the total number of discretization nodes
c   dnear - the near interaction radius
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
        dnear      = disc(7)       
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
c       decomposition.  See the documentation at the top of this
c       file for a discussion of ``evaluation data.''
c
c                       Input Parameters:
c
c   disc - the discretization structure
c
c                       Output Parameters:
c
c   eval - a (13,n) array each column of which gives the evaluation
c       data for a single point
c
c       Fetch data from the disc structure.
c
        norder     = disc(1)
        nself      = disc(2)
        maxtris    = disc(3)
        ntris      = disc(4)
        len        = disc(5)
        itris      = disc(6)
        dnear      = disc(7)
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
c       Return the coordinates of the discretization nodes on the 
c       surface described in the disc structure.
c
c                          Input Parameters:
c
c   disc - the discretization structure
c 
c                         Output Parameters:
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
        dnear    = disc(7)
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


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Subroutines for finding approximate minimum bounding balls for
c       collections of points in R^3.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine disc_bounding_ball(n,eval,bx,by,bz,br)
        implicit double precision (a-h,o-z)
        dimension eval(13,1)
c
c       Given evaluation data for a collection of points, find an 
c       approximate minimum bounding ball for them.  The resulting ball
c       is guaranteed to enclose the points, but not guaranteed to be
c       an enclosing ball of minimum radius.
c
c
c                           Input parameters:
c
c   n - the number of points
c   eval - an (13,n) array containing evaluation data for the points
c
c
c                         Output parameters:
c
c   (bx,by,bz) - the center of the bounding ball
c   br - the radius of the bounding ball
c
c
c       Pick an arbitrary point.
c
        i1 = min(n,max(1,n/2))
c
c       Find the point y most distant from it.
c
        i2   = 0
        dmax = 0
        do 1000 j=1,n
        dd = (eval(2,i1) - eval(2,j))**2
        dd = (eval(3,i1) - eval(3,j))**2 + dd
        dd = (eval(4,i1) - eval(4,j))**2 + dd
c
        if (dd .gt. dmax) then
        i2   = j
        dmax = dd
        endif
 1000 continue
c
c       Find the point z most distant from y.
c
        i3   = 0
        dmax = 0
        do 1100 j=1,n
        dd = (eval(2,i2) - eval(2,j))**2
        dd = (eval(3,i2) - eval(3,j))**2 + dd
        dd = (eval(4,i2) - eval(4,j))**2 + dd
        if (dd .gt. dmax) then
        i3   = j
        dmax = dd
        endif
 1100 continue
c
c       Take the initial bounding ball to be the ball centered at the
c       midpoint of y and z with radius equal to half their distance.
c
        bx = (eval(2,i2) + eval(2,i3))/2
        by = (eval(3,i2) + eval(3,i3))/2
        bz = (eval(4,i2) + eval(4,i3))/2
c
        br = (eval(2,i2) - eval(2,i3))**2
        br = (eval(3,i2) - eval(3,i3))**2 + br
        br = (eval(4,i2) - eval(4,i3))**2 + br
c
c       Iteratively construct a larger bounding ball.  Check to see
c       if all points are enclosed; if not, extend the ball and move
c       it closer to the point which is not enclosed.
c
        iters = 0
 2000 continue
        iters = iters+1
c
        if (iters .gt. 100) then
        print *,"disc_bounding_ball bombed !!!!!"
        stop
        endif
c
        do 2100 i=1,n
        dd = (bx-eval(2,i))**2 + (by-eval(3,i))**2 + (bz-eval(4,i))**2
        if (dd .gt. br) goto 2200        
 2100 continue
        goto 3000
 2200 continue
c
        bx = (bx + eval(2,i))/2
        by = (by + eval(3,i))/2
        bz = (bz + eval(4,i))/2
        br = (eval(2,i)-bx)**2 + (eval(3,i)-by)**2 + (eval(4,i)-bz)**2
        goto 2000

 3000 continue
        br = sqrt(br)
c
        end



        subroutine disc_bounding_ball0(n,eval,bx,by,bz,br)
        implicit double precision (a-h,o-z)
        dimension eval(13,1)
c
c       Given evaluation data for a collection of points, find an 
c       approximate minimum bounding ball for them.  The resulting ball
c       is guaranteed to enclose the points, but not guaranteed to be
c       an enclosing ball of minimum radius.
c
c
c                           Input parameters:
c
c   n - the number of points
c   eval - an (13,n) array containing evaluation data for the points
c
c
c                         Output parameters:
c
c   (bx,by,bz) - the center of the bounding ball
c   br - the radius of the bounding ball
c
        bx = 0
        by = 0
        bz = 0
        dd = 1.0d0/n
        do 1000 j=1,n
c
c        dd = eval(1,j)**2
c        bx = bx + eval(2,j)*dd
c        by = by + eval(3,j)*dd
c        bz = bz + eval(4,j)*dd
c
        bx = bx + eval(2,j)*dd
        by = by + eval(3,j)*dd
        bz = bz + eval(4,j)*dd
c
 1000 continue
c
        br = 0
        do 1100 i=1,n
        dd = (bx-eval(2,i))**2 + (by-eval(3,i))**2 + (bz-eval(4,i))**2
        br = max(br,dd)
 1100 continue
c
        br = sqrt(br)
c
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Main kernel evaluation routine.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine eval_kernel(ikernel,zk,x,y,val)
        implicit double precision (a-h,o-z)
        dimension x(13),y(13)
        double complex zk,val
c
c       Evaluate one of the kernels associated with the Helmholtz
c       equation at wavenumber zk.         
c
        if (ikernel .eq. 1) then
        call ksingle(x,y,zk,val)
        return
        endif
c
        if (ikernel .eq. 2) then
        call kdouble(x,y,zk,val)
        return
        endif
c
        if (ikernel .eq. 3) then
        call ksingleprime(x,y,zk,val)
        return
        endif
c
        end


        subroutine ksingle(x,y,zk,val)
        implicit double precision (a-h,o-z)
        double complex zk,val,ima
c
        data ima     / (0.0d0,1.0d0) /
        data pi      / 3.14159265358979323846264338327950288d0 /
        data over4pi / 0.0795774715459476678844418816862571810d0 /
c
        dimension x(13),y(13)
c
        xwht  =  x(1)
        x1    =  x(2)
        x2    =  x(3)
        x3    =  x(4)
c        dx11  =  x(5)
c        dx21  =  x(6)
c        dx31  =  x(7)
c        dx12  =  x(8)
c        dx22  =  x(9)
c        dx32  =  x(10)
c        dnx1  =  x(11)
c        dnx2  =  x(12)
c        dnx3  =  x(13)
c
        ywht  =  y(1)
        y1    =  y(2)
        y2    =  y(3)
        y3    =  y(4)
c        dy11  =  y(5)
c        dy21  =  y(6)
c        dy31  =  y(7)
c        dy12  =  y(8)
c        dy22  =  y(9)
c        dy32  =  y(10)
        dny1  =  y(11)
        dny2  =  y(12)
        dny3  =  y(13)
c
        dd    = sqrt((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)
c
        val   = exp(ima*zk*dd)/dd
        val   = -val*xwht*ywht*over4pi
c                
        end


        subroutine kdouble(x,y,zk,val)
        implicit double precision (a-h,o-z)
        dimension x(13),y(13)
        double complex val,ima,zk
c
        data ima     / (0.0d0,1.0d0) /
        data pi      / 3.14159265358979323846264338327950288d0 /
        data over4pi / 0.0795774715459476678844418816862571810d0 /
c
        xwht  =  x(1)
        x1    =  x(2)
        x2    =  x(3)
        x3    =  x(4)
c        dx11  =  x(5)
c        dx21  =  x(6)
c        dx31  =  x(7)
c        dx12  =  x(8)
c        dx22  =  x(9)
c        dx32  =  x(10)
c        dnx1  =  x(11)
c        dnx2  =  x(12)
c        dnx3  =  x(13)
c
        ywht  =  y(1)
        y1    =  y(2)
        y2    =  y(3)
        y3    =  y(4)
c        dy11  =  y(5)
c        dy21  =  y(6)
c        dy31  =  y(7)
c        dy12  =  y(8)
c        dy22  =  y(9)
c        dy32  =  y(10)
        dny1  =  y(11)
        dny2  =  y(12)
        dny3  =  y(13)
c
        dd    = sqrt((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)
        val   = ((x1-y1)*dny1 + (x2-y2)*dny2 + (x3-y3)*dny3)/dd**3
        val   = val * exp(ima*zk*dd)*(1-ima*zk*dd)
        val   = -val*xwht*ywht*over4pi
c
        end




        subroutine ksingleprime(x,y,zk,val)
        implicit double precision (a-h,o-z)
        dimension x(13),y(13)
        double complex val,ima,zk
c
        data ima     / (0.0d0,1.0d0) /
        data pi      / 3.14159265358979323846264338327950288d0 /
        data over4pi / 0.0795774715459476678844418816862571810d0 /
c
        xwht  =  x(1)
        x1    =  x(2)
        x2    =  x(3)
        x3    =  x(4)
c        dx11  =  x(5)
c        dx21  =  x(6)
c        dx31  =  x(7)
c        dx12  =  x(8)
c        dx22  =  x(9)
c        dx32  =  x(10)
        dnx1  =  x(11)
        dnx2  =  x(12)
        dnx3  =  x(13)
c
        ywht  =  y(1)
        y1    =  y(2)
        y2    =  y(3)
        y3    =  y(4)
c        dy11  =  y(5)
c        dy21  =  y(6)
c        dy31  =  y(7)
c        dy12  =  y(8)
c        dy22  =  y(9)
c        dy32  =  y(10)
c        dny1  =  y(11)
c        dny2  =  y(12)
c        dny3  =  y(13)
c
        dd    = sqrt((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)
        val   = ((y1-x1)*dnx1 + (y2-x2)*dnx2 + (y3-x3)*dnx3)/dd**3
        val   = val * exp(ima*zk*dd)*(1-ima*zk*dd)
        val   = -val*xwht*ywht*over4pi
c
        end

        subroutine kdouble0(x,y,val)
        implicit double precision (a-h,o-z)
        dimension x(13),y(13)
c
        data pi      / 3.14159265358979323846264338327950288d0 /
        data over4pi / 0.0795774715459476678844418816862571810d0 /
c
        double complex val
c
        xwht  =  x(1)
        x1    =  x(2)
        x2    =  x(3)
        x3    =  x(4)
c        dx11  =  x(5)
c        dx21  =  x(6)
c        dx31  =  x(7)
c        dx12  =  x(8)
c        dx22  =  x(9)
c        dx32  =  x(10)
c        dnx1  =  x(11)
c        dnx2  =  x(12)
c        dnx3  =  x(13)
c
        ywht  =  y(1)
        y1    =  y(2)
        y2    =  y(3)
        y3    =  y(4)
c        dy11  =  y(5)
c        dy21  =  y(6)
c        dy31  =  y(7)
c        dy12  =  y(8)
c        dy22  =  y(9)
c        dy32  =  y(10)
        dny1  =  y(11)
        dny2  =  y(12)
        dny3  =  y(13)
c
        dd    = sqrt((x1-y1)**2+(x2-y2)**2+(x3-y3)**2)
        val   = ((x1-y1)*dny1 + (x2-y2)*dny2 + (x3-y3)*dny3)/dd**3
        val   = -val*xwht*ywht*over4pi
c
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Koornwinder polynomials evaluation routines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine koorn(n,x,y,pols)
        implicit double precision (a-h,o-z)
        dimension pols(1)
c
        dimension plege(n+1),pjac(n+1)
c
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



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Discretization quadrature table
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine discquad(ier,norder,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension xs(*),ys(*),whts(*)
c
c       Return one of the discretization quadratures for integrating
c       polynomials over the standard simplex.  The returned quadrature
c       will integrate products of the polynomials of the specified
c       order.
c
c              discretization    quadrature
c                 order            order        nodes     Gaussian
c              --------------    ----------     -----     --------
c                   4                8            17         15
c                   6               12            32         28
c                   8               16            52         45
c                  12               24           112         91
c                  16               32           192         153
c            
c       All quadrature formulas returned by this routine are accurate
c       to around 30 digits, have positive weights and interior nodes.
c
c                            Input Parameters:
c
c   norder - the order of the quadrature to return
c
c                            Output Parameters:
c
c   ier - an error return code;
c      ier = 0     indicates successful execution
c      ier = 1024  means the specified parameter norder was invalud
c
c   (nquad,xs,ys,whts) - the requested quadrature formula will be 
c      returned in these user-supplied parameters
c
c       4th order
c
        parameter ( nquad4  =  18 )
        dimension xs4( 18)
        dimension ys4( 18)
        dimension whts4( 18)
        data xs4 / 
     -    0.111510689055420320758371932973043818D+00,
     -    0.205785376558420858705704866257697464D+00,
     -    0.679821402152797400322428843504781994D+00,
     -    0.809370810339840651983923586728957624D-01,
     -    0.913158580969700885327005941830180654D+00,
     -    0.257306011352595793532499998036283573D-01,
     -    0.106259825952878440913038723553464688D+00,
     -    0.832993122198875709269739460322035329D-02,
     -    0.268944523777391920136175712909689334D+00,
     -    0.873800804180638866412103742559710143D+00,
     -    0.405395602394820033507763708136243590D+00,
     -    0.673699314472323647289115087022576234D+00,
     -    0.412257117298964921429989954791549002D+00,
     -    0.176198777999669354675412722900667663D+00,
     -    0.555922628161725058090311212523691165D-01,
     -    0.884525178986688614088777268604863572D-02,
     -    0.417036739497255038597412945682500780D+00,
     -    0.574936088667004367325608906906876131D+00 /
        data ys4 / 
     -    0.682990142202969793092685844706760447D+00,
     -    0.113992628707034030323071388254327075D+00,
     -    0.209352371212608321842695475061743206D+00,
     -    0.913649518493695627772833722987951447D+00,
     -    0.265867392329577964963305429059887938D-02,
     -    0.872488688034925654539505035869161689D+00,
     -    0.240776928055565033960604368483818363D-01,
     -    0.797776771967703610671107477750121448D-01,
     -    0.675611395682081854894144617450485668D+00,
     -    0.991667708880568818600282966018722061D-01,
     -    0.184356213133506468406111067480624818D+00,
     -    0.575293015528198627436039107997520030D-01,
     -    0.411386095302821876062050225309057698D+00,
     -    0.417447411019152874705623622712643179D+00,
     -    0.270983993096316974992450524844352587D+00,
     -    0.576320687299224518539508433403681129D+00,
     -    0.128466493679543560514759440479375799D-01,
     -    0.416364441545200713343822338877151376D+00 /
        data whts4 /
     -    0.378582318218716976642047432125302763D-01,
     -    0.384561846415558453341262110015790847D-01,
     -    0.381545581496914418755590250750621446D-01,
     -    0.680865335378137294866146645745330185D-02,
     -    0.643248615292015290164033173467939966D-02,
     -    0.126042690923014720779779765113741890D-01,
     -    0.123278785956719274268900889011554758D-01,
     -    0.747626094357094436108643392498068890D-02,
     -    0.338327485393527309153674830610553581D-01,
     -    0.131876003476571734867848464944241576D-01,
     -    0.585539384963104933422048556552705365D-01,
     -    0.348689216044728480884515191054490889D-01,
     -    0.578727064016056482491941071142408602D-01,
     -    0.579968247065765105103569343925627687D-01,
     -    0.340203625729134656323532741777056864D-01,
     -    0.160099331691533763415395912577018785D-01,
     -    0.175575850052951555652317094957090150D-01,
     -    0.159808564052977432783694024270671019D-01 /
c
c       6th order
c
        parameter ( nquad6  =  32 , norder6 =  12)
        dimension xs6( 32)
        dimension ys6( 32)
        dimension whts6( 32)
        data xs6 / 
     -    0.183513326653622185111446165515263088D+00,
     -    0.103586042548336847857731791011550210D+00,
     -    0.736233314262211881314007548241130836D+00,
     -    0.347166116435620294236330997799043892D+00,
     -    0.240518094255921132157666909369109644D-01,
     -    0.881610884728015720207995361887434604D+00,
     -    0.549927013527146686476676871247294031D+00,
     -    0.989939382387382249730279834844157818D+00,
     -    0.275916649046619493624396578044038177D+00,
     -    0.498403582189943103872958494460632159D+00,
     -    0.111920755631506668506637215833217146D+00,
     -    0.488155780103293880906919139215303254D+00,
     -    0.708023777480999534673920751518388510D+00,
     -    0.197488192620153143581475832504511794D-01,
     -    0.125679160353626808138944898616833843D+00,
     -    0.261719938191235412434297734195244654D-01,
     -    0.245399353379381620824743858396482351D-01,
     -    0.262465226331651614935134489813952854D-01,
     -    0.728598494409173534850529149711395789D+00,
     -    0.686140089682223671607992650554843666D+00,
     -    0.518743992530373426152407768389661877D+00,
     -    0.294311220857264016670774897960434226D+00,
     -    0.288648589715406979448454849532880741D+00,
     -    0.886618507714876873456674395101099712D+00,
     -    0.301746948809257881741338370069750250D+00,
     -    0.308569705049188308647222997683396171D+00,
     -    0.248831874041978507784532403672498834D-01,
     -    0.489994383563345996135582164621987630D+00,
     -    0.127446516821654592409604002171720246D+00,
     -    0.254621751302731936951266534540104189D-01,
     -    0.133586062223426582426937820161510058D+00,
     -    0.131281077553662432681040175809219943D+00 /
        data ys6 / 
     -    0.681862126824319831730088521396117929D+00,
     -    0.747068194142398202367800358727193960D+00,
     -    0.162814213983355766067207639667682689D+00,
     -    0.531477463584557554761349811767572313D+00,
     -    0.869031292334338966869665275841664206D+00,
     -    0.977048268426027595285544289912362268D-01,
     -    0.345004606775824624996379656587207622D+00,
     -    0.706854454348007784655002180868635406D-02,
     -    0.698440759151665627172258742398610989D+00,
     -    0.120639181826794350322361310529093784D+00,
     -    0.859876678198155631209517545745753400D+00,
     -    0.490064993688693502546215405473968903D+00,
     -    0.271927113450600784782161327383899715D+00,
     -    0.962962361709883022459675790687922104D+00,
     -    0.547672718857204310443398148594106097D+00,
     -    0.986400896754916178570132011540117469D-01,
     -    0.698808210898627912870953540363081732D+00,
     -    0.138612222073217563736055024717633078D-01,
     -    0.207277020765661146051575550200051935D-01,
     -    0.103747620648624502983247990848385487D+00,
     -    0.241334135762768027438066303292102930D-01,
     -    0.253080202801288760044391688781869514D+00,
     -    0.431157254088421078920543942972886995D+00,
     -    0.272838942192273874375794510160565550D-01,
     -    0.103237601805295370643464521156013011D+00,
     -    0.187803250540523372265626432527017722D-01,
     -    0.483249249181520380950754730447593357D+00,
     -    0.264508002701460293245945973499403079D+00,
     -    0.339509895505890367643601623230058081D+00,
     -    0.266602002442582265003443146794883358D+00,
     -    0.302799789561626639503248932524401863D-01,
     -    0.152013372350200824824418320157985569D+00 /
        data whts6 /
     -    0.166902288033607042632363953365499226D-01,
     -    0.122192985047995315374895554170623125D-01,
     -    0.174082476103751060212369518402162945D-01,
     -    0.249636024838736593236024015271984328D-01,
     -    0.794590787870612493792626740844970872D-02,
     -    0.711187951547014455918236432885992081D-02,
     -    0.226474840765953036376004010874490991D-01,
     -    0.890853247581781975112830906447458664D-03,
     -    0.125600305850699454762408949485926400D-01,
     -    0.252415325183105663543102672486808132D-01,
     -    0.935670884289935580676961744983307021D-02,
     -    0.123999384567669915033652235186072119D-01,
     -    0.104604734083625576260267917743580605D-01,
     -    0.248146410876727165547678593853795862D-02,
     -    0.275149049757459092577733312676451273D-01,
     -    0.863760099706868559214645234441559822D-02,
     -    0.123708707452694026883259072966325871D-01,
     -    0.273954898414442435097065684090104260D-02,
     -    0.102082869217206711767785679640168731D-01,
     -    0.173469709761512261552501680843574107D-01,
     -    0.132354383798579877584335826249388877D-01,
     -    0.321629165754492887287969971924782343D-01,
     -    0.319907606311028652014256960836987553D-01,
     -    0.879885906362217375681468697640489268D-02,
     -    0.229412755836343370488601736201431823D-01,
     -    0.975627049378412635768736649536614807D-02,
     -    0.142142367111125514545530392382839065D-01,
     -    0.306027260020876203216186965267089354D-01,
     -    0.282391287815149437116087129637136030D-01,
     -    0.129409701513860186937965165441897836D-01,
     -    0.111329103646817682069814640065333291D-01,
     -    0.227886736207269548606012351987280270D-01 /
c
c       8th order
c
        parameter ( nquad8  =  52 )
        dimension xs8( 52)
        dimension ys8( 52)
        dimension whts8( 52)
        data xs8 / 
     -    0.519767719034468700213084323369720028D-01,
     -    0.211258744427687119322488961610729198D+00,
     -    0.122449382192729489361928862684160417D+00,
     -    0.717253120688961360343722482416616123D-01,
     -    0.901105923295648842430323311731358268D-02,
     -    0.696622357473493654323582659757981007D-01,
     -    0.787901113737806778269328663722378344D-01,
     -    0.338258237218030053917737448607493836D+00,
     -    0.165543163238460326185605625210038364D+00,
     -    0.500431753616941046908512361434360785D+00,
     -    0.152529901125185637200723371275502276D-01,
     -    0.167600841291828906597132379009694039D-01,
     -    0.631397846320871230338539874861617474D+00,
     -    0.176292874820775809992775935717534327D-01,
     -    0.301132311718762033773240298786681881D+00,
     -    0.199156956991012573363723368885381586D+00,
     -    0.880957343495858479405632043072009160D-01,
     -    0.164559201775063868958487075886138893D-01,
     -    0.311758943102047526080787765099928905D+00,
     -    0.144708019550043019122070579906036093D-01,
     -    0.200315829598290329064572410607878181D+00,
     -    0.213595011643334540691063135544537433D+00,
     -    0.709506076854025941720863033591359619D+00,
     -    0.477714035273900016464466481472710276D+00,
     -    0.675453626354009817384805512644925858D+00,
     -    0.187120823764174936787329814901944719D+00,
     -    0.194012066711161934202702372356502229D+00,
     -    0.461018993386758488268950923469021166D+00,
     -    0.148165955897092972208649903945106039D-01,
     -    0.826278781347896016697543939908396320D+00,
     -    0.372095890924438919283480147897196900D+00,
     -    0.799425078667219774397928012484754674D+00,
     -    0.543797598593706606352030245442951201D+00,
     -    0.810587932046778344077922580271860159D-01,
     -    0.901075509252472185397176995088630884D-01,
     -    0.339347906715904739216511798625088770D+00,
     -    0.926820087115153932730499165042136312D+00,
     -    0.170675773533975512484546272289632021D+00,
     -    0.509146287205555863912881021973778281D+00,
     -    0.141376412025942196019301231911754247D-01,
     -    0.857758097597551106604210345794347173D+00,
     -    0.342921232003291773062785224272638283D+00,
     -    0.990530381544746459926942817891645093D+00,
     -    0.349832648205172564504569681068118833D+00,
     -    0.648572178061644397558533411758126167D+00,
     -    0.677816629159252916802803389766245827D+00,
     -    0.733546298514833225922504759913852056D-01,
     -    0.867713282157954102822869438126893218D-01,
     -    0.803836602340422634781789381620287079D+00,
     -    0.925335102751862491856218650717776017D+00,
     -    0.514622693272293861716759121846551483D+00,
     -    0.158619856235392864202338577665451322D-01 /
        data ys8 / 
     -    0.192149790930970547102978292326053978D+00,
     -    0.193687229807571176370382490044191347D+00,
     -    0.203348676205388475369176844704154530D+00,
     -    0.615127315284195311969320897498241294D+00,
     -    0.162187678217058001551416992196650746D+00,
     -    0.855350540644524887736232090202859574D-01,
     -    0.345526845862314954435757841173348135D+00,
     -    0.194828806734926685287077815847156210D+00,
     -    0.842182233782586404294775262262580899D-01,
     -    0.183139351904140555263346539536821417D+00,
     -    0.624973735382955992758953261976076849D-01,
     -    0.806366572427141394133828193012694531D+00,
     -    0.730131919470696119184716450618648503D-01,
     -    0.111766102357329829303814775872259453D-01,
     -    0.822101516496926381263130811976032527D-01,
     -    0.504867884843073384488973064451216316D+00,
     -    0.690135727303891249465845577482517864D+00,
     -    0.482902981742844411676684094122907773D+00,
     -    0.674863951302845104190298122516693178D+00,
     -    0.313006923145335091231543218826138267D+00,
     -    0.646098574370965201354517508309403903D+00,
     -    0.161744635663579101972062123217274739D-01,
     -    0.131494829863217564871343079157695086D-01,
     -    0.507002212972828444271945977565420263D+00,
     -    0.150485752974188123553076701508976899D+00,
     -    0.756991658662128924490678776077041641D+00,
     -    0.343363151950373007820887093443229191D+00,
     -    0.816642774089177825350034761794449326D-01,
     -    0.920360941555398236436919472326586322D+00,
     -    0.107081338819927095798302438275380132D+00,
     -    0.159448745884964983959152015096681057D-01,
     -    0.478426665718762392255509473475882570D-01,
     -    0.156669785681175642759539393237372337D-01,
     -    0.820140266028819565436832650794091159D+00,
     -    0.167702482909426330236570066195601922D-01,
     -    0.587164170611884575199575148796415881D+00,
     -    0.211680195756237837818963868914339594D-01,
     -    0.821071518554994572373394217606949502D+00,
     -    0.410627338218868891146934420270033092D+00,
     -    0.977146484499110220791357272310888007D+00,
     -    0.267759590516013136914996315446611409D-02,
     -    0.332978420232052109409181330435959229D+00,
     -    0.367222643319820750628715370889156918D-02,
     -    0.468259706469813914713835514885783305D+00,
     -    0.335782411160853982454335858246249064D+00,
     -    0.244042792899463513525466200224446435D+00,
     -    0.906523493846325446698675709637525680D+00,
     -    0.504679547915022381967593716786125820D+00,
     -    0.181932069588624646913684124887713876D+00,
     -    0.643219209918314451968190137848612826D-01,
     -    0.297373997605183206068709932984292614D+00,
     -    0.655270932891651324566605566822520517D+00 /
        data whts8 /
     -    0.775969325534846352478245692121494492D-02,
     -    0.137198985605302663019750601259205823D-01,
     -    0.105620762786817358871660855756871620D-01,
     -    0.438874565320160964397521680626513792D-02,
     -    0.312404464756714747253831704461383289D-02,
     -    0.693425556376469421386269662869965453D-02,
     -    0.142398136881103122901457943568942678D-01,
     -    0.193069092957917908495963915501886757D-01,
     -    0.108615436204903318623443015318373769D-01,
     -    0.196584026964320480067877524127210073D-01,
     -    0.286954932381207330952238706734519242D-02,
     -    0.577107949136726646680283001492357171D-02,
     -    0.130126277566500750045396100704250008D-01,
     -    0.132143848391148730531753774119163216D-02,
     -    0.136437783476604058813089770341687027D-01,
     -    0.211475244359060152275876075528208106D-01,
     -    0.126211514903489776519172662490114145D-01,
     -    0.734453143083447672266067710872383647D-02,
     -    0.539593596277463533296620325260501999D-02,
     -    0.601981237262828773312692890842863820D-02,
     -    0.163700479831639305117379259407248791D-01,
     -    0.596095270251309016473199387074266969D-02,
     -    0.520117299387137695164533091268039942D-02,
     -    0.673052313646667729160593711125050864D-02,
     -    0.156509808335783272727689331629697388D-01,
     -    0.953834972727612515879869164543940664D-02,
     -    0.213504443965654993968547618770178303D-01,
     -    0.148992541540567379655486582903802266D-01,
     -    0.328067057693242124774117037952649542D-02,
     -    0.890649743244088218138214166119878443D-02,
     -    0.686786361362015796781614481631863126D-02,
     -    0.883352756498781560283609235386784552D-02,
     -    0.683558937145799700133349243586997356D-02,
     -    0.953125392623532488410900229114902469D-02,
     -    0.426694414709939765651219640389099246D-02,
     -    0.137614270158639175555395337589694929D-01,
     -    0.372883956132380278469443318218521299D-02,
     -    0.281482980507945316454534063266262277D-02,
     -    0.153094867816189920344129154463877190D-01,
     -    0.952761106659411008795318716124918639D-03,
     -    0.174292188523195709887220103578834878D-02,
     -    0.232497938575005274848179755284289899D-01,
     -    0.414907115489048361550471173644081503D-03,
     -    0.206582053651551312139319234175387680D-01,
     -    0.662377406202887703343322582935592872D-02,
     -    0.134322675480357718990450996415175325D-01,
     -    0.410442681747783878597947796123077576D-02,
     -    0.145828593965065003824020616176967131D-01,
     -    0.505898275678972493875407184589202050D-02,
     -    0.257068676907789846563417164613851532D-02,
     -    0.203868459244809664781885677546534471D-01,
     -    0.668009931560231936505663970507130367D-02 /
c
c       12th order
c
        parameter ( nquad12  = 112)
        dimension xs12(112)
        dimension ys12(112)
        dimension whts12(112)
        data xs12 / 
     -    0.212007984771972848950259429081967053D-01,
     -    0.334518912987752716191069463986953654D-01,
     -    0.836096397616232326266075142986624256D+00,
     -    0.751637629819907463890568484669755920D-01,
     -    0.826620172019117014880452564440682633D+00,
     -    0.897622250564933297580794021396885131D+00,
     -    0.683252089200133155216296356353614741D+00,
     -    0.893717820633065829711745143830922470D+00,
     -    0.466994557996387790065439988779282169D+00,
     -    0.602352940951354264657129062659814173D+00,
     -    0.740619571703650229755537661598976251D+00,
     -    0.352566734924697496849778422288604515D+00,
     -    0.387414352434552600267074497694992955D-01,
     -    0.644409216537094862163061749637319837D+00,
     -    0.115955698909983826656405389689046946D+00,
     -    0.568171586556750175583961484176157028D+00,
     -    0.802301396113893457399941306437122445D-01,
     -    0.257702527324955559447136882517471015D+00,
     -    0.591494612727228496328852152049547030D+00,
     -    0.313663973186239823437922647840056984D-01,
     -    0.336982371947612532937918819038471568D-02,
     -    0.694173482165086503124682625592507144D+00,
     -    0.742162617421569270132798462301910086D+00,
     -    0.624559824560955414129464354798179452D+00,
     -    0.804207229831245880012054285264899019D+00,
     -    0.608888379140937635638691211867905691D-02,
     -    0.113855942603054370472947933046767158D+00,
     -    0.227664539878483925826601243522330216D-01,
     -    0.369028037974573517673330360179443337D+00,
     -    0.865263890896172247762031453464161677D-01,
     -    0.866112490098715501684267902688798158D-01,
     -    0.794167791008842367284095954418163570D+00,
     -    0.476391520952460657987245666840279072D-01,
     -    0.173587501861882536666913420075868889D+00,
     -    0.154112352251690303719368190627375432D+00,
     -    0.200443343280572849673541902390281634D+00,
     -    0.241407519633375089039598392674971248D+00,
     -    0.275234031992592729652731969299223708D-01,
     -    0.969700768731799763002314625896007685D-01,
     -    0.352100324973462508942497050966815023D+00,
     -    0.202223341788925361778773704454064407D+00,
     -    0.781354202505159318950259515365253544D+00,
     -    0.481867186234012534063997132943560348D+00,
     -    0.383570792747980778862042977927606579D-01,
     -    0.135121734432126021420753081067964129D+00,
     -    0.235215993437063641010456022311318381D-01,
     -    0.464541206694505201287950013893138187D+00,
     -    0.412474511295396353287012172539870357D+00,
     -    0.266776953032586646871408687579658499D+00,
     -    0.519811155466883610653498507203627507D+00,
     -    0.476367203326364325540779924900917130D-01,
     -    0.518878791437239454909467310590528033D+00,
     -    0.170136776639897086657344270898367582D+00,
     -    0.773168330990482422591621140962346010D-01,
     -    0.597890027612330935498449953124808458D-02,
     -    0.340688085125802802870160323657312948D+00,
     -    0.222039833695936679499234204062375540D+00,
     -    0.104605551400143849195142729432326131D+00,
     -    0.146830460267175816058501557697732587D+00,
     -    0.711516207646938601380239860104316143D+00,
     -    0.328287570914417034416626193044741166D+00,
     -    0.942521907219667722111023873332305594D+00,
     -    0.873551944932777565764508848175474983D+00,
     -    0.658143426076844057238967242416248332D+00,
     -    0.250776679904002413938491440982057904D+00,
     -    0.689467898024509684537474418911259583D-01,
     -    0.375356140152911097265703249681803738D-01,
     -    0.930786138995389289754329912941140905D-01,
     -    0.539801434741713223205593576321352208D+00,
     -    0.433999700763260698237979943423809824D+00,
     -    0.141811503330297824641158462067058224D+00,
     -    0.877551008681931861956596611518861668D+00,
     -    0.689628956374518731403217986600498007D+00,
     -    0.366554097539169970280141947004217299D+00,
     -    0.261865592062846383105012774701550708D+00,
     -    0.673673458203306699174433494621652120D+00,
     -    0.237197853523820430303698108024890749D+00,
     -    0.427948435339640990686212035050676309D+00,
     -    0.938765427233566897966155110964725075D-01,
     -    0.297873168666635326000686218495431161D+00,
     -    0.819980582833908210849697375705317255D-02,
     -    0.406263082875665133585440349727664937D-01,
     -    0.423815192295754920308539528833603803D-01,
     -    0.233714580115707238602852153224948795D+00,
     -    0.554306719034746673158377814234885624D+00,
     -    0.950044876799911079397448265050040963D+00,
     -    0.739088513149581679445334613603807473D-01,
     -    0.445059990169966074669737751161364455D+00,
     -    0.569762094758820984853429269045589778D+00,
     -    0.595877024585143181373339994353206419D+00,
     -    0.789958621069065904500101630082982174D-02,
     -    0.187565143936474041798160565840313276D+00,
     -    0.315564519710598263021957408410952300D+00,
     -    0.174564980904398380104869334722150790D+00,
     -    0.386374023997383268499830766804175622D-01,
     -    0.156888693175283457580548127956663726D+00,
     -    0.787791593465877471202345078785322468D+00,
     -    0.273525573904712633776423023115470743D+00,
     -    0.312617243740232826209116572942469876D+00,
     -    0.398427408996808790933862447852584198D+00,
     -    0.148794656074623026750450035785659202D+00,
     -    0.423500072935007177327595445466109247D-02,
     -    0.636925144148089210437233348614361874D-03,
     -    0.471878863699833819385974885139208808D+00,
     -    0.946808916924728797590763414087000654D+00,
     -    0.353780747780002081582798178392718755D+00,
     -    0.616627361484751313932121999527658952D-02,
     -    0.745484012236257341391678454366085812D-02,
     -    0.775304812943010855021902236105325588D-02,
     -    0.983477544716303033675107523004640106D+00,
     -    0.851582400112790086670199367658855773D-02,
     -    0.851603488386221933904411932648994301D-02 /
        data ys12 / 
     -    0.921310657538613517792326107150123387D+00,
     -    0.940865136710725688383728523523223035D+00,
     -    0.702805676786681798504454732088700445D-01,
     -    0.894778283439697258022676869083536855D+00,
     -    0.115085504824427211765834381657915421D+00,
     -    0.290973360489551893184005560821144877D-01,
     -    0.711598190161278995020047889878112615D-01,
     -    0.752220063094018365043178856721686556D-01,
     -    0.486957042440043445246297618478156545D+00,
     -    0.136871697225858278400149739396671008D+00,
     -    0.982514480313270555801020035854112563D-01,
     -    0.611456460180077059601783463819253060D+00,
     -    0.866539711072106515519873762770887135D+00,
     -    0.172993414350105311631460695437798988D+00,
     -    0.571843856281274550038111274940783991D-02,
     -    0.881165092510220025315395228198393282D-01,
     -    0.280029047457005651985584659097983372D-01,
     -    0.716447789303908586573413837980056545D+00,
     -    0.379008130166744199997177555571288667D+00,
     -    0.281116444742220583144007594928825904D-01,
     -    0.951092223809558294893285609034848616D+00,
     -    0.285559854889598896787548722644076623D+00,
     -    0.158190782730182155707425948631361762D+00,
     -    0.245317106553968346806811191223646218D+00,
     -    0.359767199922500977247034565091608897D-01,
     -    0.658807292814385981912332145257529650D+00,
     -    0.878431228839042923291144713818926398D+00,
     -    0.712175383701717769704461711956364383D+00,
     -    0.150780916602064253030440112155927437D+00,
     -    0.841003698574732327647339374420887566D+00,
     -    0.554558906677782215866196786031754130D+00,
     -    0.173397056005448719776632932157613170D+00,
     -    0.949285115547783199486900389680419459D+00,
     -    0.152399066884997434815404038478962302D+00,
     -    0.347638476314395627794360839457185850D-01,
     -    0.793101208966036497316997438481697120D+00,
     -    0.748347409373782604322410662846186919D-01,
     -    0.726933601057982951073900696766689108D-01,
     -    0.242898826052194426348160894430555258D+00,
     -    0.299912885312231440679008866799602014D-01,
     -    0.627915690252356729291472212162322834D-02,
     -    0.213287451653352188822348115600175188D+00,
     -    0.173706434451914887958570041266484181D+00,
     -    0.590696496939078439731750341061875920D+00,
     -    0.569675531571458629351411024287691344D+00,
     -    0.800593001672240891665084945033168532D+00,
     -    0.973317271766847749691987610708116695D-01,
     -    0.378055004010849992039285805936697475D+00,
     -    0.142077499797059716029764786580832198D+00,
     -    0.241185324883598680454325145526153667D+00,
     -    0.405366585859484144888818468458374928D-02,
     -    0.323207370356106396933529082054493571D+00,
     -    0.238835270973650539557129978119060859D+00,
     -    0.769102503636360750293214042443230885D-01,
     -    0.112002999684139074892627438284797215D+00,
     -    0.567842244398767029617611738964353613D+00,
     -    0.534187658441230397218076172038965374D+00,
     -    0.450664799128007028927064617665134468D+00,
     -    0.746454319597525301085351556495133288D+00,
     -    0.290580306542870877135306999043099796D-01,
     -    0.495420155835724174576233670987295599D+00,
     -    0.283187331680849521821957135877109832D-01,
     -    0.119161007982616189722373609458419230D+00,
     -    0.340236992030594905021038310439339365D+00,
     -    0.284347143232910137060247130959207417D-01,
     -    0.677954552208482096199247557850374343D+00,
     -    0.143215860468262268742384390895513426D+00,
     -    0.153306388243057682086569810732043866D+00,
     -    0.452346236929728670944985430752150292D+00,
     -    0.724460081258809586052051381215353659D-02,
     -    0.663555243763182654288015388315370720D+00,
     -    0.724140242035501197170033000520163189D-02,
     -    0.247932692578902374384149549399218083D+00,
     -    0.230193914400912711507729202169232140D+00,
     -    0.230258205021103567976749490305634330D+00,
     -    0.624997366024122532147715131271281605D-02,
     -    0.691099697358728234249019189279836839D+00,
     -    0.563494277746820231670626687439678722D+00,
     -    0.345065948968649731180844599809409210D+00,
     -    0.413747254067906410245538102276259852D+00,
     -    0.541935394603519813854442356010899825D+00,
     -    0.238867794194192442447617995224142095D+00,
     -    0.464738184451583023904017627254843056D+00,
     -    0.620526771464334320687888938874877069D+00,
     -    0.765762492771085843648300880572453593D-02,
     -    0.361963112318588443171649443118483471D-02,
     -    0.771872007247715751534846818922182661D+00,
     -    0.442381981042948379648457675312123968D+00,
     -    0.351976600465541592208410245764543788D+00,
     -    0.371500521902699254332354679807347447D-01,
     -    0.198413089733604904016045380175844495D+00,
     -    0.446393510774999060614089033487347348D+00,
     -    0.679785492260373217019768616404092347D+00,
     -    0.338429771703118990314403157637586492D+00,
     -    0.347334225104676397563699952859027390D+00,
     -    0.803924326020873204856871813751380320D+00,
     -    0.598007509757986293640399842913793008D-02,
     -    0.327340971817736000786679227429053349D+00,
     -    0.437578185983409923029938771504765309D-02,
     -    0.297388114626245537073506621983368356D+00,
     -    0.858435391848871472337882947447562818D-01,
     -    0.440864151648885364612893493467658341D-01,
     -    0.782623843633260319798772560762777950D+00,
     -    0.400419152677640393624487579579841231D-01,
     -    0.486394631343715321359969356226748125D-01,
     -    0.781710698519761043072132677663555804D-01,
     -    0.878526596395423554518556561809349787D+00,
     -    0.303065953227123597629454668170874103D+00,
     -    0.421626979967290019683966140728105260D+00,
     -    0.850855503172922851627730984186269802D-02,
     -    0.982899521605245386787339540600396953D+00,
     -    0.790241486122170751525098614815975064D-02 /
        data whts12 /
     -    0.115770259962358461498826714303361630D-02,
     -    0.117832437830862475199853530412886648D-02,
     -    0.363062791621321464402893147315672693D-02,
     -    0.197287052432453849572824753030262926D-02,
     -    0.265816052583096653624342054735779444D-02,
     -    0.208001143727411939549404365280311160D-02,
     -    0.526688375751911717156016078094878849D-02,
     -    0.227165593694590219069961899513833487D-02,
     -    0.625958552597740633956407856969452864D-02,
     -    0.686016632509918049282437566761572267D-02,
     -    0.580655499310519524691375832821840627D-02,
     -    0.459282850123998205021953182517354616D-02,
     -    0.261404146001352398348617538704662441D-02,
     -    0.671214131749702989262693822536830399D-02,
     -    0.115975438185289135166681453653643213D-02,
     -    0.603019341919750060205803462561723067D-02,
     -    0.215439304385399193229824807724376639D-02,
     -    0.309949940242650413850669453032527551D-02,
     -    0.391730753606741829888378735154982514D-02,
     -    0.131947384150053795956941066699474889D-02,
     -    0.534588767373359696005156568161995249D-03,
     -    0.313641057874468429258656350733732483D-02,
     -    0.585950195072932714867168588109558016D-02,
     -    0.694403030604171906854371522233536122D-02,
     -    0.389267539156616579104167086160929269D-02,
     -    0.176900523919831519926844890367823345D-02,
     -    0.147826648962994265946479455519358082D-02,
     -    0.314950479078666916287770830922168016D-02,
     -    0.823552323227079651777746097979862568D-02,
     -    0.352732692452334870740036578673862182D-02,
     -    0.521346217905350800853254994287738848D-02,
     -    0.356425993107281387226931033003769280D-02,
     -    0.560371047745585916229743366450974078D-03,
     -    0.689568050155775108829542155574338783D-02,
     -    0.351177937489407811546018474489023086D-02,
     -    0.163404635101788973694406646521235094D-02,
     -    0.604597787538863107160421636082823396D-02,
     -    0.199742357775635336758254996077771012D-02,
     -    0.631163417425175446768611151994367350D-02,
     -    0.402920575224166454615728570936579440D-02,
     -    0.154990004753523716319528840876777247D-02,
     -    0.147033921546298989466462187154710721D-02,
     -    0.988424053892361278615524058063056711D-02,
     -    0.506813816097167964622170453156786147D-02,
     -    0.727183178952195435757698247138159396D-02,
     -    0.311792548807300277218462637477118058D-02,
     -    0.717777479199477824160535998249422357D-02,
     -    0.949374538590281791374248719405747360D-02,
     -    0.777499120226100569184221065176268381D-02,
     -    0.101329735520882671800982044703107032D-01,
     -    0.640428232908581323199583693763215143D-03,
     -    0.797531796165930928280451788365105646D-02,
     -    0.771431193075750031187679531603764019D-02,
     -    0.387536642150138671444819483396083836D-02,
     -    0.121245239979236951723392437389476375D-02,
     -    0.769860120554101416875288728081892190D-02,
     -    0.955363532225328225702159954581827598D-02,
     -    0.784144667771844838624190872335174074D-02,
     -    0.620321677603631942399449928766816514D-02,
     -    0.333234041012865785974791825430598195D-02,
     -    0.972145678521743411173038946201328212D-02,
     -    0.122327720561874647247402216667372311D-02,
     -    0.151779823197093763554098845304702295D-02,
     -    0.948955705762713724134174458740413040D-03,
     -    0.335942826295028033878115802418643661D-02,
     -    0.625615040187627289201870819427959332D-02,
     -    0.361956057847485965648306500749237170D-02,
     -    0.565329034384501033636859126705519372D-02,
     -    0.222571119353413515926558975012315586D-02,
     -    0.223333922625983796544204059041993346D-02,
     -    0.772600455885893644896115356307656451D-02,
     -    0.144989762776022536268598721953384188D-02,
     -    0.562991025288700245896760974308568589D-02,
     -    0.876344815033089344363826311866367054D-02,
     -    0.944791414817708718942135494513786176D-02,
     -    0.187491413909015712245950247254475056D-02,
     -    0.593879527123554959347197140315983743D-02,
     -    0.245367008790029044322207567061227093D-02,
     -    0.724768128608164195780480494054310709D-02,
     -    0.108462496767359860407612930072777627D-01,
     -    0.248706084866858695932081120141785679D-02,
     -    0.467972826942943578866060295311895417D-02,
     -    0.584072958537811953558446470530054462D-02,
     -    0.834035465526379358586313287325845570D-02,
     -    0.235075219020532194193977438924710291D-02,
     -    0.616578062226889110505865984519958669D-03,
     -    0.540671776114902850929605749550351236D-02,
     -    0.865437068867206865536312354105979001D-02,
     -    0.706601924083463892661926984539762095D-02,
     -    0.490351699659677463809901194377775750D-02,
     -    0.192839759929153481728307549877482110D-02,
     -    0.104800477442472066557575257756571186D-01,
     -    0.150357164595949091577444088486067032D-02,
     -    0.959326226501388433965366756047401232D-02,
     -    0.491085940821317887918907325058598774D-02,
     -    0.414173751766085558204220454263752688D-02,
     -    0.160622922533836156456558201921888627D-02,
     -    0.106696353848284641403077979767635177D-01,
     -    0.146309413684628812738335669187115707D-02,
     -    0.104232786648009638922265022157993420D-01,
     -    0.490898137604637403190038989438496788D-02,
     -    0.638095842145505562140941559636776316D-03,
     -    0.767336876653061430851143719852549921D-03,
     -    0.557669188374741186299813150046274993D-02,
     -    0.723630824520942598859508927087774983D-03,
     -    0.691143210948163115895365494059041740D-02,
     -    0.129938939786654432650381718867944389D-02,
     -    0.216338409354786631541051561167954944D-02,
     -    0.240949938459188612601955819122635587D-02,
     -    0.452573596872072442591479876267600373D-03,
     -    0.478878762601534950276798718934314067D-03,
     -    0.444882055984480858059335584171554970D-03 /
c
c       16th order
c
        parameter ( nquad16  = 192 )
        dimension xs16(192)
        dimension ys16(192)
        dimension whts16(192)
        data xs16 / 
     -    0.605912786488662267270246341906806959D+00,
     -    0.731644135567609711107739688244948012D+00,
     -    0.183843477281472319981041896880769009D+00,
     -    0.687172392410141350310170568829470048D+00,
     -    0.688919951282826843643167044105365958D+00,
     -    0.553375373691711939570307596548756011D+00,
     -    0.580687557066870729258903080943575016D+00,
     -    0.547246866267507394113210902944955001D+00,
     -    0.416154160127948357871225578275645012D+00,
     -    0.111403188782795548047269582171098005D+00,
     -    0.535318860605049441867647072711736976D+00,
     -    0.632684312240306584702487482226216993D+00,
     -    0.841087037949741274858611530058367583D-01,
     -    0.106164610868507645721580286186153001D+00,
     -    0.140929097757753467212799730165840009D+00,
     -    0.405568815963952563804601063312441978D+00,
     -    0.604746142463519226462945556169910023D+00,
     -    0.461027767204503615195702339625426013D+00,
     -    0.487605173433380523475521013667839990D+00,
     -    0.323088618749324676790746885991835126D-01,
     -    0.542602859995790382080404915858406003D+00,
     -    0.690199413548752575161860898840081012D+00,
     -    0.480160952140214825136890517459004008D+00,
     -    0.425766569117567507462646038519799973D-01,
     -    0.739266638142273306954829689388879128D-01,
     -    0.240191834262341052490807017296627997D+00,
     -    0.388213961386427801921354358670596096D-01,
     -    0.222071324081617355231322243514046996D+00,
     -    0.628510743443898478962192407039094041D+00,
     -    0.285943664489726862070903416746572985D+00,
     -    0.867875927228837516076161362741195042D+00,
     -    0.632333411073021626722806391862876047D+00,
     -    0.156542851268511285667123753376824314D-01,
     -    0.737107261980060871973977056770216594D-01,
     -    0.880348693840611892532813466350488028D+00,
     -    0.740374997087746396980671413727872975D+00,
     -    0.711798038113870445971129105332038752D-01,
     -    0.818284607995712545322890376135588976D+00,
     -    0.755368027576413568712622418197106189D-01,
     -    0.765708196402146059096991291876298961D+00,
     -    0.712148990502203068359025768156805020D+00,
     -    0.381466670984379775357048186952104020D-01,
     -    0.206442237023029627770063919656938087D-01,
     -    0.175428476480742704445726154712855612D-01,
     -    0.670197974417353823610475440136084572D-01,
     -    0.873596167454914717186556254652092974D+00,
     -    0.787813500731596206949092747049493603D-01,
     -    0.329100670513266811275579172108039986D+00,
     -    0.464226735883768300300490274339508016D+00,
     -    0.486007595528132716656658828511890001D+00,
     -    0.649884433447743045375863540928391978D+00,
     -    0.757966843059074240772336112898394032D+00,
     -    0.912271651144903888012326287010473429D-01,
     -    0.792125291251883639864383050022401097D-02,
     -    0.247727024997153464134212163948177009D+00,
     -    0.538614834666816940013340278087955027D+00,
     -    0.658295615903097296359809945478110498D-03,
     -    0.398610901756426100537823715458008327D-01,
     -    0.569225542724627097186079297929736963D+00,
     -    0.485304515172750895695909014676280288D-01,
     -    0.375580496633109838076737696587919001D+00,
     -    0.158581953150788363711132048073407002D+00,
     -    0.319193293363308887210767265001485006D+00,
     -    0.279002224945076572401811529864888988D+00,
     -    0.135515013235447280030672325511529203D-01,
     -    0.213664948058237511670442146529442989D-02,
     -    0.212857355310781953702006845888376005D+00,
     -    0.506661386975151137332584125391475375D-01,
     -    0.360837014789707714102890899245826986D+00,
     -    0.937748440954206178228020017228262066D-01,
     -    0.380415383345867725850154367270684986D+00,
     -    0.819371955952570507161791160050433961D+00,
     -    0.145598681026674127481052356935756001D+00,
     -    0.120768433165047830839053472884073602D-01,
     -    0.212818663073603123915416391996941298D-01,
     -    0.639215630323290008182474429014806028D+00,
     -    0.475690348959208729216365529848655514D-01,
     -    0.297765443227622082630303506323754980D+00,
     -    0.884627741853965219075796024862288975D-01,
     -    0.864568727226756720859042451617123004D+00,
     -    0.279412353076146831222213980420310006D+00,
     -    0.168749509237686325520737867472848715D-01,
     -    0.209748584816818487161761308346682191D-02,
     -    0.502670407737302874199422935391178305D-01,
     -    0.780489008513990372041837604067253011D+00,
     -    0.142492286440007663905930135434243000D+00,
     -    0.901031527723605120613757294845585027D-01,
     -    0.410480066638944284137934714664981992D+00,
     -    0.856202769144277889223548421087161026D+00,
     -    0.204342394717842768398142943127910995D+00,
     -    0.465137085304281541729823578694959012D+00,
     -    0.151683628452550059432585485092295696D-01,
     -    0.789230987429405382722731826775979595D-01,
     -    0.454389956735023004611331408689076993D+00,
     -    0.435618725732091186169102643170516785D-01,
     -    0.515021842007392172236641876624314020D-01,
     -    0.284989597875659677043541929284721993D+00,
     -    0.336331819346883487199167173182319982D+00,
     -    0.470404395595510514028781234817318006D+00,
     -    0.246418943701211268700795899139862692D-01,
     -    0.144712047300640119662479961332532997D+00,
     -    0.402179114155203041065658429394425768D-02,
     -    0.739457082638371491967319874261742952D+00,
     -    0.455759269307039919812080865174016415D-01,
     -    0.129646343682152546593667150408271997D+00,
     -    0.638448167293646443443452480492845016D+00,
     -    0.286345445487152095196167039696979994D+00,
     -    0.137392954647170860379899906783722011D+00,
     -    0.117419925668302343520869342018377001D-01,
     -    0.206854988071210976230202521645211008D+00,
     -    0.205584366522247567326830588293763999D+00,
     -    0.823420877006223468060838948384157766D-01,
     -    0.811692312247131560898653796888907999D+00,
     -    0.350385286136594968902774263426271987D+00,
     -    0.127858591059562080564618754660417999D+00,
     -    0.161989059160254356316452401979306003D+00,
     -    0.346400688648167273449325807880035993D-02,
     -    0.394226797413142774043907754948695418D-01,
     -    0.806024691703380701042686444059717002D+00,
     -    0.253126663126863050746377020777658004D+00,
     -    0.477628922940507030739413767057487866D-02,
     -    0.181478072778990474088980211357202395D-01,
     -    0.107912238924124534476585905518748005D+00,
     -    0.225229642297073545096227820575845007D+00,
     -    0.281612911790340402103483656254168008D+00,
     -    0.127050101599883966285515157566574994D+00,
     -    0.227213961479294520644883756210449009D+00,
     -    0.791475326633070078624049495854001005D+00,
     -    0.112889328008704715704619864663721999D+00,
     -    0.229576423747032144808547500410285011D+00,
     -    0.195174677073452388141187134468539989D+00,
     -    0.745609684703029417339626942081401017D+00,
     -    0.368490651724789295700933162946037977D+00,
     -    0.158799587118674438121347442580034003D+00,
     -    0.494198487011383741111979259438551924D-01,
     -    0.128325665686733791188966009858396010D+00,
     -    0.399757453718107406259445127902838980D+00,
     -    0.989396711580881718701348717505816985D+00,
     -    0.872550650999380735991271436863438962D+00,
     -    0.199051244499529433788411746586506406D-01,
     -    0.289898865268452903955429605130509002D+00,
     -    0.372671627101357607320333381483813989D+00,
     -    0.988838892891072409943494547524988045D+00,
     -    0.809852574910397715721930274269688033D+00,
     -    0.959517564033955715510294412060827997D+00,
     -    0.256613227059912011542622078883787789D-02,
     -    0.156628122225589180347509632651789988D+00,
     -    0.104983989371973945846442876773593003D+00,
     -    0.304851279695973001310071164040394999D+00,
     -    0.677305794971507074240104467108506966D+00,
     -    0.328982686146986699272936935639195593D-02,
     -    0.220194723474798340284852310756359006D+00,
     -    0.510101006488170321743491591105754212D-01,
     -    0.511667348103994715623611073437519413D-02,
     -    0.369547951043731542844402414401074977D+00,
     -    0.508217375138171759789056595153359967D+00,
     -    0.959700604223827012997807483591734012D+00,
     -    0.238894238878598439050628797722080206D-02,
     -    0.201643983687079837791788719630966010D+00,
     -    0.251632230458502938176278099441554020D+00,
     -    0.309942394629587121153135754269524985D+00,
     -    0.918288125128679140915853459284461980D+00,
     -    0.173538348281990733200622344973130996D+00,
     -    0.168225707507016182948269798148049004D+00,
     -    0.915300125392637737527204468430463004D+00,
     -    0.377593875744033754692494758299553817D-02,
     -    0.546529853190871884976520295244047044D+00,
     -    0.192608966291757621233849426039386693D-01,
     -    0.540089567543424776347570639028433964D+00,
     -    0.603804002975863155815787144896986999D+00,
     -    0.170862157994026910966300280099599002D+00,
     -    0.924450552586783959471812111387757037D+00,
     -    0.701069537628429861415970244516379962D+00,
     -    0.411719250654044037813284307231778007D+00,
     -    0.442936682857231227751128799770395663D-02,
     -    0.416784762040788528201993650762057817D-02,
     -    0.961085893342229736737290909289620598D-01,
     -    0.442754176550653772103180535788750014D+00,
     -    0.966850122173295600197469386081433955D+00,
     -    0.210940789207804314380706521776354212D-01,
     -    0.357787682422142190832880061989011706D-02,
     -    0.692583603374638815201123650150461007D+00,
     -    0.927936136633538343088642956024104959D+00,
     -    0.229809006858059047891690428097647998D-01,
     -    0.439400177203575625604094976960839002D+00,
     -    0.344880529980326873016151421878899986D+00,
     -    0.307345906994825981952940701087617005D+00,
     -    0.205531317017316746544614307519434008D+00,
     -    0.596971218406940879903882578057223971D+00,
     -    0.173385270626270955130422365396236595D-02,
     -    0.401575648433094865001469193621449984D+00,
     -    0.499289974080756797611296679126066982D+00 /
        data ys16 / 
     -    0.347746335068449506089255973161153023D+00,
     -    0.249124305813420746407557802620655990D+00,
     -    0.678638119666369956308515371369069955D+00,
     -    0.228616213048335400890280079095279999D+00,
     -    0.268138238938005000262019573849479000D+00,
     -    0.166102254279560855396814759546007991D+00,
     -    0.210722990539668442453564735283762997D+00,
     -    0.117404350919289294850301027383527004D+00,
     -    0.354950336994960122172979081813428984D+00,
     -    0.327314634303517110340632120665311016D+00,
     -    0.293354867973821514521611578995137976D+00,
     -    0.294363542131984561768210662377892016D+00,
     -    0.289642608025081220975631902451983994D+00,
     -    0.464598952754336140977973880353465000D+00,
     -    0.691342096217473902649924389427414977D+00,
     -    0.250822222259440292174171306889972000D+00,
     -    0.442317805729638871565175554464427875D-01,
     -    0.373833919079848555411990969945724979D+00,
     -    0.213226754238766311864133461974102000D+00,
     -    0.398055960837352505471914622006683532D-02,
     -    0.394782317414308175303386940205512006D+00,
     -    0.424432943255081191721109668107646611D-01,
     -    0.274033144309442900928997236770828014D+00,
     -    0.605008844638221507806520891507629990D+00,
     -    0.882673527801923548944486461187660990D+00,
     -    0.624231815894449387076770155627111967D+00,
     -    0.772807739871492273429139725643584013D+00,
     -    0.285388843922390442435179961779420976D+00,
     -    0.815027220741467960629115673909690688D-01,
     -    0.628957728014298718551601815421699006D+00,
     -    0.565763079573138219594309909613167808D-01,
     -    0.233739220662737140342161458691360008D+00,
     -    0.164146167767711896390557775961394487D-01,
     -    0.213739265914298826838158785470373001D+00,
     -    0.741750691841393089960468072059364157D-01,
     -    0.157965527137765076995865143701553985D-01,
     -    0.484118266501870851110022068944007872D-02,
     -    0.179368442769058993270160141829324113D-01,
     -    0.635475533017380240435553993959058964D+00,
     -    0.417140624445412467324951527803375974D-01,
     -    0.795807997416920182500932500937356819D-01,
     -    0.698686959780470062934074161211014040D+00,
     -    0.876358252899108217627073365278373997D+00,
     -    0.966293480073970737048447101388475998D+00,
     -    0.538405448640924366212114874668856030D+00,
     -    0.242688584325602043254990421903600414D-01,
     -    0.715469099194186913555480303574128007D+00,
     -    0.557364452267135450023467104126040953D+00,
     -    0.106694789672772147828146328258997004D+00,
     -    0.405925844019251088870885970341289979D+00,
     -    0.177650752551755457046958481578201801D-01,
     -    0.192737384633270154735408668766472998D+00,
     -    0.397556433561737351014040645659457999D+00,
     -    0.101162007712719576657013197429574801D-02,
     -    0.105191085736651501499353901185174698D-01,
     -    0.684878193555037014052379674945796033D-01,
     -    0.795137373629545109153199509000285973D+00,
     -    0.191276378394229519892171340625675994D+00,
     -    0.321806745669773464563196957669252011D+00,
     -    0.877808489216210324002788844310920016D+00,
     -    0.321412555350567643097124146746119019D+00,
     -    0.278945537178516519192015456014389984D+00,
     -    0.100317959665767875302583919263129998D+00,
     -    0.672319448611382349828053333374931810D-01,
     -    0.821128889931103073203033831286029029D+00,
     -    0.136202368027122018958677154755662998D+00,
     -    0.702654623238985886106201812180646000D+00,
     -    0.127774961087675812896256959236118002D+00,
     -    0.212979666128456003397805829488784002D+00,
     -    0.766672541547368687449292819941549983D+00,
     -    0.458682039911005966569852923835870996D+00,
     -    0.127595106414520600835555993569128009D+00,
     -    0.757274719179733110402911010398487031D+00,
     -    0.745896352467815712397601879304311006D+00,
     -    0.481367803506762179438230591449971122D-01,
     -    0.130365745099968945029737341976139993D+00,
     -    0.238864244111282443792447053312081590D-01,
     -    0.383939274598983491210827998438399008D+00,
     -    0.825957014286031614691455913842362019D+00,
     -    0.369530005889348688690671697329945402D-02,
     -    0.672005280900258520261643771786105953D+00,
     -    0.649955825420479076716513449479155988D+00,
     -    0.994801627858906451651341412435819440D-02,
     -    0.928129224182144231218202261707839984D+00,
     -    0.215845553931234533380831316594072997D+00,
     -    0.854515925533726723476678864960919960D+00,
     -    0.869962701166936711438189364699215997D-01,
     -    0.490712064434827373771873534247530993D+00,
     -    0.140221969786639329250400936254534000D+00,
     -    0.770752674010485986544707431495776268D-01,
     -    0.479863426582325044578141738620963006D+00,
     -    0.164317784998813362914417864907109989D+00,
     -    0.442445783981742972957366225672710921D-01,
     -    0.167109796107968575328775853843436002D+00,
     -    0.280607817410325944923682004603726984D+00,
     -    0.467524553761069149954871063468822017D+00,
     -    0.316476042849100199990125419514409714D-01,
     -    0.429081651464320009335725452317011996D+00,
     -    0.499731268688493454954006615110019573D-01,
     -    0.927562123487102412308123547676268004D+00,
     -    0.954878483872474167258555962745391833D-01,
     -    0.963915774373009881587700336567385999D+00,
     -    0.161306354500723310305726457114112005D+00,
     -    0.708729776488155722941998932431874221D-01,
     -    0.819530107022585998752153198498119973D+00,
     -    0.340959352072116873236545241391477995D+00,
     -    0.159598789639420208342981485513765995D+00,
     -    0.469203488027839599638858616211488116D-01,
     -    0.986234660194573673978066262772936958D+00,
     -    0.373415537240117074053433238233375493D-01,
     -    0.590855355533391335964700957351950398D-03,
     -    0.912906413572819026435975906207960961D+00,
     -    0.558218623538691409542834960307886875D-01,
     -    0.146409642568681391688209827795249998D-01,
     -    0.179760274199994701072707586551643805D-02,
     -    0.456203907964812279556106538572447020D+00,
     -    0.370596422682144621087407973469047829D-01,
     -    0.956867117325088188060877061785325024D+00,
     -    0.998322510250815431256113587321544947D-01,
     -    0.493881893777173379309254493092821984D+00,
     -    0.786749252028865572373111789249006704D-01,
     -    0.248645381988795363498235199043381997D+00,
     -    0.871775651123070336638988648661011965D+00,
     -    0.443556443410185910290491894261036021D+00,
     -    0.224424458657791268773962390053186004D+00,
     -    0.219623183346770069754488928136823004D+00,
     -    0.363687541516096044235990705309129987D+00,
     -    0.295712340178232784884256150720443308D-02,
     -    0.552049242206578061348693455243761979D+00,
     -    0.124794297548958954095769561271618006D+00,
     -    0.758168754519413133621649767416000037D+00,
     -    0.104353533711567300293722247063482004D+00,
     -    0.459645855256725349649650420327827627D-01,
     -    0.151205439215303773948879218674651008D+00,
     -    0.821307187774205591431752918011378037D+00,
     -    0.631863374503895563223181162936136016D+00,
     -    0.846639197699236995886118458115955508D-01,
     -    0.662013396719534948773189996980627730D-03,
     -    0.106191298984453163121618475670908997D+00,
     -    0.104676204099281778094939910769930995D+00,
     -    0.533017987412201400524694575680893969D+00,
     -    0.146268476767384004231261930230397010D+00,
     -    0.892167966883408302633603257553958703D-02,
     -    0.170458167949039339872990969087545992D+00,
     -    0.429068914424591745773916942821609672D-02,
     -    0.871601652117496329948883070969238021D+00,
     -    0.369928511302431348107891342152302997D+00,
     -    0.178378941945487395416815700454291296D-01,
     -    0.296700232023351134788889809562344987D+00,
     -    0.161434560444058581699187549178914989D+00,
     -    0.215283809743802390517027081445249992D+00,
     -    0.776124827768380177559034935614923984D+00,
     -    0.369551797270255538810333098562716019D+00,
     -    0.923640857706019157239373440242610969D+00,
     -    0.573512164619628283261040436633048037D+00,
     -    0.453900621017250703294335606629930067D-02,
     -    0.365606065171750277481149559941028175D-01,
     -    0.695137705830528736868927257932220016D+00,
     -    0.592904563960382032720546763124526035D+00,
     -    0.727807453040706382914478502278920953D+00,
     -    0.225469647390484112735205839249892298D-02,
     -    0.526091671355519141405741020195970568D-02,
     -    0.537359219168477558796471069556037979D+00,
     -    0.135374390360700112135551445295180404D-01,
     -    0.800289106988620684858851039903480752D-01,
     -    0.303839748848598984420426752838278985D+00,
     -    0.231609014318146814131091079704200915D-01,
     -    0.344871615572919060476823911711653007D+00,
     -    0.436217418398070522676739316587647976D+00,
     -    0.411288133198589760101163719764492672D-02,
     -    0.811167808592809292081596146779662030D+00,
     -    0.258225648583376070394983595785942814D-01,
     -    0.288659891945216145069729948793236104D-02,
     -    0.286687823178524043223703277129742482D-02,
     -    0.499932487536325282030027796654960006D+00,
     -    0.596929696728276650489953822771095024D+00,
     -    0.148944418544028000565638680742865012D+00,
     -    0.200593820914494138127913601947469786D-01,
     -    0.163680314022836674793996920679753491D-01,
     -    0.439704598525130881547551863564951990D+00,
     -    0.401069437953194061277925315617758984D+00,
     -    0.303404931139574879600817759825886977D+00,
     -    0.500590411112857001591331378009790173D-01,
     -    0.542626371187898068752394478213878975D+00,
     -    0.537300864806688438880626597572697994D+00,
     -    0.632478321750155611381765500551299989D+00,
     -    0.688275979958745506937044459866086997D+00,
     -    0.202505344062487821911565867347754999D+00,
     -    0.398770376436230603535097640554604989D+00,
     -    0.990983889795537350467779830827030962D+00,
     -    0.594124819895014456738908905636889994D+00,
     -    0.496068650580086066143071596824139019D+00 /
        data whts16 /
     -    0.202545007112460065595710866636814582D-02,
     -    0.177099417788321245447414749584142492D-02,
     -    0.253413642673988526523671077241491697D-02,
     -    0.292963963689327017495398082906118596D-02,
     -    0.216472548032116133811346350437904088D-02,
     -    0.426236723884360160808762380054224378D-02,
     -    0.556109537121606702422531779146068126D-02,
     -    0.458902953756652699803725604533704094D-02,
     -    0.527607435728611742665364707772547279D-02,
     -    0.309724518799810518852834134505567492D-02,
     -    0.527709813617109311027453306821431380D-02,
     -    0.258467957034259660149898855434217782D-02,
     -    0.357643723641519102899931182741675204D-02,
     -    0.430755020293215222046422025626081634D-02,
     -    0.376979851492682841757305845827738395D-02,
     -    0.542712738547393623970200129371155088D-02,
     -    0.258876881955502729088897366657263418D-02,
     -    0.489546745879779160240341813564056620D-02,
     -    0.436877021631995825472141623366586929D-02,
     -    0.325571141927113107930222587211997000D-03,
     -    0.323300826626098929451957346703502616D-02,
     -    0.264308704804184842578591459275888984D-02,
     -    0.576753900315109311092176860981875633D-02,
     -    0.260523349301709159908930496728313198D-02,
     -    0.143488320908003195416526322474821291D-02,
     -    0.409880031568426605654916517743309632D-02,
     -    0.225385619279428779857437918250494304D-02,
     -    0.571222438876046708102248113204676271D-02,
     -    0.386999765552240432236879878078923912D-02,
     -    0.296138920596027620957868528419429909D-02,
     -    0.172908322557166476504783871352658808D-02,
     -    0.464135758622042978383126923720172417D-02,
     -    0.479278868951271219347906426522794008D-03,
     -    0.314067773225574431811590055673310402D-02,
     -    0.138083029342344362203813624736730204D-02,
     -    0.163847658094329721639968854290953206D-02,
     -    0.568874523997443954343363777187827569D-03,
     -    0.149950818789440480461671506395946591D-02,
     -    0.376633458885686416588761943052455114D-02,
     -    0.215724319298799857615944056287187494D-02,
     -    0.329358247934874903088620842143352595D-02,
     -    0.266014013752526305801956594086428097D-02,
     -    0.140682304329325420440457859650633708D-02,
     -    0.506770827738736533756405417616780907D-03,
     -    0.343196189060053660434465103313830313D-02,
     -    0.156762107162307932725812521230721292D-02,
     -    0.361598976879869877268081403069375415D-02,
     -    0.419398791689941700068999887894310168D-02,
     -    0.447846625797415353840440161074327105D-02,
     -    0.442715255712207482689075239846568633D-02,
     -    0.190344012903233112028369138357785399D-02,
     -    0.262071246158760364346171011063905083D-02,
     -    0.387841273554196430212468422435351717D-02,
     -    0.754547955107086019625321601603675931D-04,
     -    0.142401683235890641457800869636277697D-02,
     -    0.351590787407807924636513866042394604D-02,
     -    0.355344809610963288707847814527449806D-03,
     -    0.229641537190335610614394309070452911D-02,
     -    0.441704607486966779641376319198363434D-02,
     -    0.150971671892481519513895493889084396D-02,
     -    0.557163367348741550479375073089227385D-02,
     -    0.482010696160380961631397938170665263D-02,
     -    0.412291261435361904846269135827001667D-02,
     -    0.332144817805565906007221750665721406D-02,
     -    0.140211634346426562368028188208684309D-02,
     -    0.449638169801363449120872417638460483D-03,
     -    0.317716088223456244797089060312207413D-02,
     -    0.244686422988696583850769082992660588D-02,
     -    0.549111623804690636895512309857608106D-02,
     -    0.304423464111767519270077365708694309D-02,
     -    0.561102910251051104312725701267686833D-02,
     -    0.237516164757093157812194468290131496D-02,
     -    0.334789604559289821938050321931522118D-02,
     -    0.154094104085841980141114702019774209D-02,
     -    0.989480872013132832579061548033418589D-03,
     -    0.463094341578931427926966942241502468D-02,
     -    0.111004253685202986546190127845585806D-02,
     -    0.613289731592851291550267885693975064D-02,
     -    0.227149253819346354653767183917864002D-02,
     -    0.618222591564695180323419481541261215D-03,
     -    0.281783703879632740444159769489521001D-02,
     -    0.189941878951863560544250391592024002D-02,
     -    0.110882296512258251943790712453862795D-03,
     -    0.100850643840161503423608690404735393D-02,
     -    0.772730129022300545751720979412349479D-03,
     -    0.568949021505929079305523575105430843D-03,
     -    0.259357102011102771239198665945488809D-02,
     -    0.419978035075237638959049979915000691D-02,
     -    0.642875262206786599957913819440150886D-03,
     -    0.315131459666799893603912052586520390D-02,
     -    0.349964647555664519376228996921663612D-02,
     -    0.140395128130126063656281691770107900D-02,
     -    0.159411373462121598768998395025228904D-02,
     -    0.571446445413174770045571181803388387D-02,
     -    0.280068999388987305289466904458196417D-02,
     -    0.352057298304286878694587572664259310D-02,
     -    0.229687869043391938367781046842040484D-02,
     -    0.605680543930789477208818860282633667D-02,
     -    0.348100630858888432615161221449678418D-02,
     -    0.112584361077856789355139379399452003D-02,
     -    0.316422134512586168259100727668129209D-02,
     -    0.336214552755681398264028501204709423D-03,
     -    0.365348427418058086213013745797192391D-02,
     -    0.162732578994347600103328208883679407D-02,
     -    0.234414576027200372597711976284189414D-02,
     -    0.207774683482242655558372456149241187D-02,
     -    0.486240755965413770820754905512172475D-02,
     -    0.256637723891001120915820713564265300D-02,
     -    0.115606004230772096764713112725127102D-03,
     -    0.243590090776891193921829088809759694D-02,
     -    0.333562880226322711784149330237906886D-03,
     -    0.605628863273078566048285100838318898D-03,
     -    0.255361520886864583040619605820232301D-02,
     -    0.180395084389973684475452565734150705D-02,
     -    0.401244010798833019287448063943147808D-03,
     -    0.528665885694523431848639480675989731D-02,
     -    0.317385486080104918699406131976678683D-03,
     -    0.342440026576420793593040564319500108D-03,
     -    0.262179948293656702450871510338150186D-02,
     -    0.569748846604771011330257332638357664D-02,
     -    0.583029510406660207470936036566126546D-03,
     -    0.179673508608818062650297779947105105D-02,
     -    0.133034271263040074938301672682482792D-02,
     -    0.491550980193299967796682753558662832D-02,
     -    0.537446846408153574839880064271375226D-02,
     -    0.460767531158990009543095528383727368D-02,
     -    0.597024030826267275132501488109038535D-02,
     -    0.640651342737739072815997881299394839D-03,
     -    0.471757783416240420519629080243214634D-02,
     -    0.462965242971584582920670789398694210D-02,
     -    0.257088308778015098598038456686648618D-02,
     -    0.330209335991852312353702206936142013D-02,
     -    0.351777588930488405508392580768707497D-02,
     -    0.424874110865248681192514395942359267D-02,
     -    0.225274281674937606953028622993258518D-02,
     -    0.471459066423673027365759281633612573D-02,
     -    0.385524029353912897461202019421593882D-02,
     -    0.791662164533434581599924597903805969D-04,
     -    0.139302651859959615990699603705353498D-02,
     -    0.133739926707291519373975305422829503D-02,
     -    0.554698876762924674767523452133026692D-02,
     -    0.542506974278358676847665240340786278D-02,
     -    0.119183190751704915804454894008195195D-03,
     -    0.168271938290142743273015286435710401D-02,
     -    0.368662693066499105359174092283769595D-03,
     -    0.469942176107648846738318225624092384D-03,
     -    0.573955734994784132385254591779341371D-02,
     -    0.128863493381297705844680132414480503D-02,
     -    0.681433625358497196889934553563217884D-02,
     -    0.428718208076615037407965394013743976D-02,
     -    0.707164249125183880423456464149468307D-03,
     -    0.786383558625025438567929159183709265D-03,
     -    0.358522921014349073760788031137666291D-02,
     -    0.599821680160495894449653037167069013D-03,
     -    0.386041040829531203699369233567746406D-02,
     -    0.110530191922246650516646683876510108D-02,
     -    0.347684817044839119257848284592062899D-03,
     -    0.636578486132780531922271001534334674D-03,
     -    0.548208505389258030602698551474999926D-02,
     -    0.197322430073675956249841773150773109D-02,
     -    0.609394726076270754924285668147694093D-03,
     -    0.635711157680036100236957199236704073D-03,
     -    0.536068994648287324253780305986663218D-02,
     -    0.149186553515444989102812296078134303D-02,
     -    0.606789883414910355464388102905976981D-03,
     -    0.892837288015536252068846515567821390D-03,
     -    0.260959106844814447594577338749982191D-02,
     -    0.211932359994639416121555887898328795D-02,
     -    0.266249699062065499044958171217061904D-02,
     -    0.100497752175056801547933731334091694D-02,
     -    0.158744657847805129203058477359027508D-02,
     -    0.123164149824977872022555550569261795D-02,
     -    0.725340320235501511079278912975118146D-03,
     -    0.776292579968383798968621273263023184D-03,
     -    0.110702217783615527684812493535078798D-02,
     -    0.102357709092582919717865380031207296D-02,
     -    0.357891371503332517530072353351544988D-02,
     -    0.233553405468887576107173903230405291D-02,
     -    0.526310251509651340140425616767042539D-03,
     -    0.239176928161419348425127965545995911D-02,
     -    0.931139833742874878561985580178056053D-03,
     -    0.937772838802478577177074924812912043D-03,
     -    0.106132304576097004259561046472006205D-02,
     -    0.263101985820015101913886778746157003D-02,
     -    0.255449620492014307287417101600081881D-02,
     -    0.238435628443987163439343659251816506D-02,
     -    0.101259013780165945475815345984906606D-02,
     -    0.578385345410589365092509896762021531D-02,
     -    0.106094973494820244113667430555976498D-02,
     -    0.928409812480481694422135322458798110D-04,
     -    0.107256741780964829365322539297673399D-02,
     -    0.116080783360129741951412195941370700D-02 /
c
        ier = 0
c
c       Return the 4th order quadrature.
c
        if (norder .eq. 4) then
        nquad = nquad4
        do 1004 i = 1,nquad4
        xs(i)   = xs4(i)
        ys(i)   = ys4(i)
        whts(i) = whts4(i)
 1004 continue
        return
        endif
c
c       Return the 6th order quadrature.
c
        if (norder .eq. 6) then
        nquad = nquad6
        do 1006 i = 1,nquad6
        xs(i)   = xs6(i)
        ys(i)   = ys6(i)
        whts(i) = whts6(i)
 1006 continue
        return
        endif
c
c       Return the 8th order quadrature.
c
        if (norder .eq. 8) then
        nquad = nquad8
        do 1008 i = 1,nquad8
        xs(i)   = xs8(i)
        ys(i)   = ys8(i)
        whts(i) = whts8(i)
 1008  continue
        return
        endif
c
c       Return the 12th order quadrature.
c
        if (norder .eq. 12) then
        nquad = nquad12
        do 1012 i = 1,nquad12
        xs(i)   = xs12(i)
        ys(i)   = ys12(i)
        whts(i) = whts12(i)
 1012 continue
        return
        endif
c
c       Return the 16th order quadrature.
c
        if (norder .eq. 16) then
        nquad = nquad16
        do 1016 i = 1,nquad16
        xs(i)   = xs16(i)
        ys(i)   = ys16(i)
        whts(i) = whts16(i)
 1016 continue
        return
        endif

c
        ier = 1024
        end
