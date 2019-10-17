        implicit double precision (a-h,o-z)
c
        double precision, allocatable :: w(:),tris(:,:)
        double precision, allocatable :: evalin(:,:),evalout(:,:)
c
        double complex zk
c
        pi  = acos(-1.0d0)
c
c       Set discretization and solve parameters.
c
        norder  = 4
        nself   = 0
        dnear   = 2.0d0
        zk      = 1.0d0
        nlevels = 1
c
c       Allocate a big work array.
c
        lw = 50 000 000
c
 0001 continue
        if (allocated(w)) deallocate(w)
        allocate(w(lw),STAT=istat)
        if (istat. eq. 0) then
        lw = lw + 50 000 000
        goto 0001
        else
        lw = lw - 50 000 000
        endif
c
        if (allocated(w)) deallocate(w)
        allocate(w(lw))
c
        call prinf("in scattering, lw = *",lw,1)
c
c       Construct a discretization structure.
c
        idisc = 1
        ldisc = lw-idisc
c
        pi    = acos(-1.0d0)
        nself = norder
c
        call disc_init(ier,norder,nself,dnear,w(idisc),ldisc)
        if (ier .ne. 0) then
        call prinf("after disc_init, ier = *",ier,1)
        stop
        endif
c
c       Add a decomposition of the torus.
c
        maxtris = 100 000
        allocate(tris(6,maxtris))
c
        x1 = 0
        x2 = 2*pi
        y1 = 0
        y2 = 2*pi
c
        iparam  = 1
c
        call triangulate(ier,iparam,x1,y1,x2,y2,nlevels,
     -    maxtris,ntris,tris)
c
c$$$        M = nlevels
c$$$        N = 4*nlevels
c$$$        call torus_triangles(M,N,ntris,tris)
c
        call disc_add(ier,w(idisc),iparam,ntris,tris)
        deallocate(tris)
c
c       Fetch and display information about the decomposition stored
c       in the disc structure.
c
        call disc_info(w(idisc),ntris,nquad,n,dnear,lkeep)
        ldisc = lkeep
c
        call prinf("total # of triangles ntris = *",ntris,1)
        call prinf("total # of nodes n = *",n,1)
c
        dd = lkeep*8
        dd = dd/(1024d0**2)
c
        call prin2("disc structure size (in MB) = *",dd,1)
c
c       Fetch the evauation data and coordinates of the discretization
c       nodes.
c
        ieval = idisc+ldisc
        leval = 13*n
c
        izs   = ieval+leval
        lzs   = 3*n
c
        lsofar = izs+lzs
        if (lsofar .gt. lw) then
        call prina("out of memory*")
        stop
        endif
c
        call disc_points(w(idisc),w(izs))
        call disc_data(w(idisc),w(ieval))
c
        iplot=1
        call plot3d_points("discretization nodes *",iplot,n,w(izs))  
c
c       Build the octree and near interaction lists.
c
        ioc = izs+lzs
        loc = lw-ioc
c
        k   = 10 000
        call elapsed(t1)
        call octree(ier,n,w(izs),k,w(ioc),loc,lkeep)
        call elapsed(t2)
        if (ier .ne. 0) then
        call prinf("after octree, ier = *",ier,1)
        stop
        endif
        loc = lkeep
c
        iinear = ioc+loc
        linear = lw-iinear
c
        call near_lists(ier,w(idisc),w(ioc),w(iinear),linear,lkeep)
        if (ier .ne. 0) then
        call prinf("after near_lists, ier = *",ier,1)
        stop
        endif
c
        linear = lkeep
c
c       Construct incoming and outcoming evaluation data.
c
        ievalin  = iinear+linear
        levalin  = 13 * 10 000
c
        ievalout = ievalin+levalin
        levalout = 13 * 10 000
c
        iw2      = ievalout+levalout
        lw2      = lw-iw2
c
        if (lw2 .le. 0) then
        call prina("out of memory*")
        stop
        endif
c
        n0 = 30
        x0 = 2.0d0
        y0 = 0.0d0
        z0 = 0.0d0
        r0 = 0.1d0
c
        call sphere_data(n0,x0,y0,z0,r0,nin,w(ievalin))
c
        n0 = 30
        x0 = 10.0d0
        y0 = 0.0d0
        z0 = 0.0d0
        r0 = 1.0d0
c
        call sphere_data(n0,x0,y0,z0,r0,nout,w(ievalout))
c
c       Build the scattering matrix.
c
        naux       = 30
        eps        = 1.0d-13
c
        call scattering_matrix(ier,eps,zk,w(idisc),w(ioc),w(iinear),
     -    nin,w(ievalin),nout,w(ievalout),naux,w(iw2),lw2)
c
        end


        subroutine surface(iparam,s,t,r,dr)
        implicit double precision (a-h,o-z)
        dimension r(3),dr(3,2)        
c
c       Supply a parameterization of the torus formed by revolving a
c       circle of radius b in the zy plane about the circle of radius
c       a in the xy plane.
c
        a = 2.0d0
        b = 1.0d0
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
        end


        subroutine near_lists(ier,disc,oc,inear,linear,lkeep)
        implicit double precision (a-h,o-z)
        dimension disc(1),oc(1),inear(1)
c
        integer, allocatable :: iwork(:)
c      
        ier = 0
c
c       Fetch some data from the disc array.
c
        ntris = disc(4)
        len   = disc(5)
        itris = disc(6)
        nquad = disc(20)
        ixs   = disc(21)
        iys   = disc(22)
        iwhts = disc(23)
        dnear = disc(99)
c
        n = nquad*ntris
c
c       Allocate work arrays from the heap.
c
        allocate(iwork(n))
c
c       Allocate some memory from inear array.
c
        iinearptrs  = 100
        linearptrs  = 4*n
c
        iinearidxs  = iinearptrs+linearptrs
        linearidxs  = linear-iinearidxs
c
        if (linearidxs .le. 0) then
        ier = 4
        return
        endif
c
c       Construct near target lists.
c
        call near_lists0(ier,oc,n,nquad,ntris,len,disc(itris),iwork,
     -    inear(iinearptrs),inear(iinearidxs),linearidxs,dnear)
        if (ier .ne. 0) return
c
        lkeep = iinearidxs+linearidxs
c
        inear(1) = n
        inear(2) = iinearptrs
        inear(3) = iinearidxs
c
        end


        subroutine near_lists0(ier,oc,n,nquad,ntris,len,tris,iwork,
     -    inearptrs,inearidxs,linearidxs,dnear)
        implicit double precision (a-h,o-z)
        dimension tris(len,1),inearptrs(4,n),inearidxs(1),iwork(1)
c
        dimension idxs(nquad)
c
        iptr = 1
c
c       First, construct the lists of near targets for each point.
c
        do 1000 jtri=1,ntris
c
        x0 = tris(8,jtri)
        y0 = tris(9,jtri)
        z0 = tris(10,jtri)
c
        r0 = sqrt(tris(11,jtri))*dnear
c
        call ocrange(oc,x0,y0,z0,r0,nn,iwork)
        if (iptr + nn .gt. linearidxs) then
        ier = 4
        return
        endif
c
        call quicksorti(nn,iwork)
c
        do 1111 i=1,nquad
        idxs(i)=(jtri-1)+i
 1111 continue
c
        call iremove(nn,iwork,nquad,idxs)
c
        do 1100 j=1,nn
        inearidxs(iptr+j-1) = iwork(j)
 1100 continue
c    
        do 1200 j=1,nquad
        jj = (jtri-1)*nquad+j
        inearptrs(3,jj) = nn
        inearptrs(4,jj) = iptr
 1200 continue
c
        iptr=iptr+nn
 1000 continue
c
c       Count the number of near sources for each target point.
c
        do 2000 i=1,n
        inearptrs(1,i) = 0
 2000 continue
c
        do 2100 i=1,n
        nn = inearptrs(3,i)
        i1 = inearptrs(4,i)
c
        do 2200 l=1,nn
        idx = inearidxs(i1+l-1)
        inearptrs(1,idx) = inearptrs(1,idx)+1
 2200 continue
 2100 continue
c
c       Allocate memory for the lists.
c
        do 2300 i=1,n
        nn = inearptrs(1,i)
        inearptrs(2,i) = iptr
        iptr = iptr + nn
 2300 continue
        if (iptr .gt. linearidxs) then
        ier = 4
        return
        endif
c
c       Form the lists of near sources.
c
        do 3000 i=1,n
        inearptrs(1,i) = 0
 3000 continue
c
        do 3100 i=1,n
        nn = inearptrs(3,i)
        i1 = inearptrs(4,i)
c
        do 3200 l=1,nn
        idx = inearidxs(i1+l-1)
        mm = inearptrs(1,idx)
        j1 = inearptrs(2,idx)
        inearidxs(j1+mm) = i
        mm = mm + 1
        inearptrs(1,idx) = mm
 3200 continue
 3100 continue
c
        linearidxs = iptr
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
c       The number of points generated is 
        
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


        subroutine torus_triangles(M,N,ntris,tris)
        implicit double precision (a-h,o-z)
        dimension tris(6,1)
        pi = acos(-1.0d0)
c
        ntris = 0
c
        dd1 = 2*pi/M
        dd2 = 2*pi/N
c
c     M
c     N       
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
        tris(5,ntris) = x2
        tris(6,ntris) = y2
c
        ntris = ntris + 1
        tris(1,ntris) = x1
        tris(2,ntris) = y1
        tris(3,ntris) = x1
        tris(4,ntris) = y2
        tris(5,ntris) = x2
        tris(6,ntris) = y2
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
c       scattering matrix code proper.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains code for compressing the integral
c       operator associated with the exterior Neumann problem
c
c
c             \Delta u + k^2 u = 0         in  \Omega
c                         dN u = f         on  \partial\Omega
c 
c       where \Omega a sufficiently, \partial\Omega denotes the boundary
c       of \Omega and dN denotes the derivative w.r.t. the outward
c       pointing normal vector to
c
c       This solver ignores
c
c       This file contains a code for compressing integral operators
c       of scattering theory on surfaces.
c
c
c
c          T\sigma(x) =  \sigma(x)  + \int K(x,y) \sigma(y) ds(y)         (1)
c                                        \Sigma
c
c       The following routines are user-callable:
c      
c   scattering_matrix - construct a scattering matrix
c
c   scattering_near_lists - construct the near interaction lists for
c       the discretization nodes on the 
c
c       The following subroutines are used internally by the scattering
c       matrix code.
c
c   scattering_garbage_collection - perform garbage collection on
c       the heap used by the scattering matrix code
c
c   scattering_create - construct
c
c   scattering_combine - combine scattering matrices on the regions
c       \Gamma_1 and \Gamma_2 in order to form a new scattering
c       matrix for their union
c
c   scattering_bounding_ball0 - construct an approximate minimum 
c       bounding ball for a collection of points in R^3 using
c       a crude method
c
c   scattering_bounding_ball - construct an approximate minimum
c       bounding ball for a collection of points in R^3 using a 
c       somewhat less crude but surprisingly less effective
c       mechanism
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine scattering_matrix(ier,eps,zk,disc,oc,inear,
     -    nin,evalin,nout,evalout,naux,w,lw)
        implicit double precision (a-h,o-z)
        dimension disc(1),oc(1),inear(1),w(1),evalin(13,1)
        dimension evalout(13,1)
c
        integer, allocatable :: iboxes(:)
c
c       This subroutine is the user-callable interface to the
c       accelerated direct solver.  It constructs a scattering matrix
c       for the problem (1).
c
c                            Input Parameters:
c
c  zk - the wavenumber for the problem
c  
c  ikernel - the index of the kernel for the integral operator which
c       is being compressed
c
c  disc - structure describing the discretization of the scatterer
c
c  oc - an octree containing the discretization nodes on the scatterer
c       (see octree.f)
c
c        This structure is used to 
c
c  inear  a structure containing the near interaction lists (see 
c       nearlists.f)
c
c  nin - the number of quadrature nodes on the surface which specifies
c       incoming potentials
c  evalin - the evaluation data for the incoming surface
c
c  nout - the number of quadrature nodes on the surface which specifies
c       outgoing potentials
c  evalout - the evaluation data for the outgoing surface
c
c  naux - the square root of the number of points per wavelength^2 
c      which should be sampled on auxilliary contours --- in other words,
c      naux^2 points are used to discretize an auxilliary contour of
c      area 1 wavelength^2
c
c
c                             Output Parameters:
c
c   ier - an error return code;
c
c       ier = 0     indiciates successful execution
c       ier = 4     means that the work array w was of insufficient length
c       ier = 1024  means that a matrix evaluation failed
c
c
c   w - upon successful execution, the beginning of the work array will
c       contain a scattering matrix for the exterior Neumann problem
c       (1)
c
c   w(1) - the number of points in the source skeleton                    (nin)
c   w(2) - the number of points in the target skeleton                    (nout)
c   w(3) - a pointer to the list of the incoming skeleton                 (iin)
c   w(4) - a pointer to the list of the outgoing skeleton                 (iiout)
c   w(5) - a pointer into the work array to the nout x nin scattering     (iskel)
c          matrix
c
c
c       Fetch some data from the near interaction structure.
c
        iinearptrs = inear(2)
        iinearidxs = inear(3)
c
c       Fetch some data from the octree.
c
        call ocinfo(oc,nlevels,nboxes,n)
c
c       Allocate a few variables.
c
        allocate(iboxes(nboxes+n))
c
c       Allocate global variables from the work array.
c
        iiflags     = 1
        liflags     = 2*n
c
        ieval       = iiflags+liflags
        leval       = 13*n
c
        iismatrptrs = ieval+leval
        lismatrptrs = nboxes
c
        iheap       = iismatrptrs+lismatrptrs
        lheap       = lw-iheap
c
        iheapptr    = 1
c
c       Zero all of the pointers and the flags.
c
        call izero(nboxes,w(iiptrs))
        call izero(2*n,w(iiflags))
c
        call disc_data(disc,w(ieval))
c
c       Construct scattering matrices on the lowest level.
c
        call oclevel(oc,nlevels,nboxes,iboxes)
c
        call prinf("in scattering_matrix, lowest level boxes = *",
     -    iboxes,nboxes)
c
        do 1000 ii=1,nboxes
        ibox = iboxes(ii)
        call prinf("in scattering_matrix, ibox = *",ibox,1)
c
        call create_smatr(eps,n,zk,w(ieval),w(iiflags),
     -    nin,evalin,
     -    nout,evalout,
     -    naux,disc,oc,inear(iinearidxs),
     -    inear(iinearptrs),ibox,w(iheap),
     -    w(iismatrptrs),iheapptr,lheap)
c
 1000 continue
c
c       Recursively combine scattering matrices level-by-level.
c
c        do 1100 ilevel=nlevels-1,1
c        call oclevel(oc,ilevel,nboxes,iboxes)
c        call prinf("in scattering_matrix, ilevel = *",ilevel,1)
c 1100 continue

        end


        subroutine create_smatr(eps,n,zk,eval,iflags,
     -    nin,evalin,nout,evalout,
     -    naux0,disc,oc,
     -    inearidxs,inearptrs,ibox,heap,
     -    ismatrptrs,iheapptr,lheap)
        implicit double precision (a-h,o-z)
        dimension disc(1),oc(1),heap(1)
        dimension eval(13,n),evalin(13,nin),evalout(13,nout)
        dimension ismatrptrs(1)
        dimension inearidxs(1),inearptrs(4,n)
c
        dimension iflags(2,n)
c
        double complex val,zk
c
        double precision, allocatable :: evalaux(:,:)
c
        double complex, allocatable   :: smatr(:,:),sinv(:,:)
        double complex, allocatable   :: xx(:,:),yy(:,:)
c
        double complex, allocatable   :: aout(:,:),rout(:)
        integer, allocatable          :: icolsout(:)
c
        double complex, allocatable   :: aint(:,:),rint(:)
        integer, allocatable          :: icolsin(:)
c
c
        double precision, allocatable :: w(:)
c
        integer, allocatable :: idxs(:),idxssrc(:),idxstar(:)
        integer, allocatable :: iproximate(:)
c
        integer, allocatable :: ilocal(:),inearsrcs(:),ineartars(:)
c
c       Construct a scattering matrix for a user-specified portion
c       of the surface.  This subroutine constructs such a matrix
c       at the lowest level of the hiearchical decomposition of the
c       input surface.
c
c       This is done by first constructing matrices Ain and Aout which
c       account for incoming and outgoing potentials and factoring
c       them as 
c
c         Ain   = R_in B_in
c 
c         Aout  = Bout R_out
c
c       where Bin and Bout are submatrices of Ain and Aout, 
c       respectively.
c
c       Next the restriction T of the integral operator (1) to the 
c       specified region is constructed and inverted.  The scattering
c       matrix is 
c
c         S     = Rout T^{-1} R_in.
c
c
c                               Input Parameters:
c
c   n - the number of total discretization nodes on the boundary surface
c   zk - the complex-valued wavenumber for the problem
c   eval - a (13,n) array providing evaluation data for all 
c       discretization nodes on the surface
c
c   iflags - a (2,n) integer array keeps track of which target and 
c       source nodes are still extant
c
c   inearidxs - an integer array storing the indices
c   inearptrs - a (2,n) array which specifies
c
c
c                              Output parameters:
c
c
c   w - upon successful execution, the beginning of the work array
c       will contain the scattering matrix; 
c
c   w(1) - the number of points in the target skeleton                    (nout)
c   w(2) - the number of points in the source skeleton                    (nin)
c   w(3) - a pointer to the indices of the target skeleton                (iiout)
c   w(4) - a pointer to the indices of the source skeleton                (iin)
c   w(5) - a pointer to the (nout,nin) scattering matrix
c

c
c       Allocate temporary variables from the beginning of the work
c       array.
c


c
c
c       Fetch some parameters from the structures.
c
        dnear = disc(99)
c
c       Find the indices of the points in the box.
c
        allocate(idxs(n),iproximate(n),ilocal(n))
        allocate(inearsrcs(n),ineartars(n))
        call ocpoints(oc,ibox,nidxs,idxs)
        call prinf("in create_smatr, idxs = *",idxs,nidxs)
c
c       Mark all nodes in the box as being unavailable as both
c       targets and sources.
c
        do 1000 i=1,nidxs
        iflags(1,idxs(i)) = 1
        iflags(2,idxs(i)) = 1
 1000 continue
c
c       Construct a bounding ball for near interactions.
c
        call scattering_bounding_ball0(nidxs,idxs,eval,bx,by,bz,br)
c
c       Find the set of "proximate" near interactions.
c
        br  = br*dnear
        call ocrange(oc,bx,by,bz,br,nproximate,iproximate)
c
c       Construct the list of "discretization" near interactions.
c
        nneartars = 0
        do 2000 i=1,n
        ilocal(i) = iflags(1,i)
 2000 continue
c
        do 2100 ii=1,nproximate
        i = iproximate(ii)
        if (ilocal(i) .eq. 0) then
        nneartars = nneartars+1
        ineartars(nneartars) = i
        ilocal(i) = 1
        endif
 2100 continue
c
        do 2200 i=1,nidxs
        ii = idxs(i)
c
        nn   = inearptrs(3,i)
        iptr = inearptrs(4,i)
        do 2300 jj=1,nn
        j = inearidxs(jj)
        if (ilocal(j) .eq. 0) then
        ilocal(j) = 1
        nneartars = nneartars+1
        ineartars(nneartars) = j
        endif
 2300 continue
 2200 continue
c
        call quicksorti(nneartars,ineartars)
        call prinf("nneartars = *",nneartars,1)
        call prinf("ineartars = *",ineartars,nneartars)
c
c       Construct the auxilliary target contour.
c
        nauxeval = 0
c
c       Form the matrix Aout and factor it as Aout = Bout * Rin
c
        nnout = nout
        allocate(rout(nidxs*nnout),icolsout(nidxs))
        allocate(aout(nnout,nidxs))
c
        do 3000 jj=1,nidxs
        j = idxs(jj)
        do 3200 i=1,nout
        call ksingle(evalout(1,i),eval(1,j),zk,val)
        aout(i,jj) = val
 3200 continue
 3000 continue
c
        allocate(w(2*(3*nidxs+1)*nnout+1))
        call factor_right(eps,nnout,nidxs,aout,krankout,icolsout,rout,w)
        deallocate(aout,w)
c
c       Form the matrix Ain^t.
c
        nnin = nin
        allocate(icolsin(nnin),rint(nnin*nidxs))
        allocate(aint(nnin,nidxs))
        do 3300 jj=1,nidxs
        j = idxs(jj)
        do 3400 i=1,nin
        call ksingle(eval(1,j),evalin(1,i),zk,val)
        aint(i,jj) = val
 3400 continue
 3300 continue
c
        allocate(w(2*(3*nidxs+1)*nnin+1))
        call factor_right(eps,nnin,nidxs,aint,krankin,icolsin,rint,w)
        deallocate(w,aint)
c
c       Form and invert the diagonal operator.
c
        allocate(smatr(nidxs,nidxs),sinv(nidxs,nidxs))
        ikernel = 3
        call disc_eval(ier,disc,ikernel,zk,nidxs,idxs,
     -    nidxs,idxs,smatr,w,lw)
c
        stop
c
        stop
        do 4000 i=1,nidxs
        smatr(i,i) = smatr(i,i) + 0.5d0
 4000   continue
c
        do 4100 i=1,nidxs
        do 4200 j=1,nidxs
        sinv(i,j) = smatr(i,j)
 4200 continue
 4100 continue
        call lapack_invert(nidxs,sinv)
c
c       Form the scattering matrix proper by taking the product
c
c          R_out (krank) S^{-1} R_in  = R_out S^{-1} ((R_in)^t)^t
c
c
        allocate(xx(krankout,nidxs),yy(krankout,krankin))
        call apply_left(nidxs,nidxs,krankout,icolsout,rout,smatr,xx)
        call apply_rightt(nidxs,krankout,krankin,icolsin,rint,xx,yy)
c
c       Form the lists of retained targets and sources.
c
        allocate(idxssrc(krankout),idxstar(krankin))
c
        stop
c
c       Form the lists of near sources and targes.
c

c       Copy the scattering matrix and its associated data into the heap,
c       store the length of the data.
c
        ismatr     = iheap+20
        lsmatr     = 2*krankout*krankin
c
        iisrcs     = ismatr+lsmatr
        lisrcs     = krankout
c
        iitars     = iisrcs+lisrcs
        litars     = krankin
c
c        iineartars = iitars+litars
c        lineartars = 0
c
c        iinearsrcs = iineartars+lineartars
c        linearsrcs = 0
c
c        idneartars = iinearsrcs+linearsrcs
c        ldneartars = 2*
c 
c        idnearsrcs = idneartars+ldneartars
c        ldnearsrcs = 2*
c
c
c        heap(iheapptr)   = krankout
c        heap(iheapptr+1) = krankin
c        heap(iheapptr+2) = ismatr
c        heap(iheapptr+3) = iisrcs
c        heap(iheapptr+4) = iitars
c        heap(iheapptr+5) = iineartars
c        heap(iheapptr+6) = iinearsrcs
c        heap(iheapptr+7) = 
c
c       Mark the retained sources and targets in the flags array.
c

c
c       Move the scattering matrix to the beginning of the work array.
c

        end


c
c      
c

        subroutine combine_smatrs()
        implicit double precision (a-h,o-z)
c
c       This routine combines too scattering matrices.
c
        end



        subroutine find_near_targets()
        implicit double precision (a-h,o-z)
c
c       Compile a list of the target points which are "near," either
c       because of proximaty or because of quadrature evaluation,
c       to a specified collection of source points.
c
        end


        subroutine find_near_sourcs()
        implicit double precision (a-h,o-z)
        end


c$$$        subroutine mark_points(iiflags,idxs)
c$$$        implicit double precision (a-h,o-z)
c$$$c
c$$$c       Mark a collection of
c$$$c
c$$$        end
c$$$
c$$$
c$$$        subroutine unmark_targets()
c$$$        implicit double precision (a-h,o-z)
c$$$c
c$$$c       Mark a collection of points as missing.
c$$$c
c$$$        end



        subroutine lapack_invert(n,zmatr)
        implicit double precision (a-h,o-z)
        double complex zmatr(n,n)
        integer*4 lw4,n4,info4
        integer*4, allocatable      :: ipiv4(:)
        double complex, allocatable :: w(:)
c
c       Invert an (n,n) double complex matrix via LAPACK.
c
        allocate(ipiv4(n+10 000))
        n4  = n
        lw4 = -1
c
c        call zgetri(n4,zmatr,n4,ipiv4,ww,lw4,info4)
c        lw4 = ww
c        allocate(w(lw4))
c
c        call zgetri(n4,zmatr,n4,ipiv4,w,lw4,info4)
        info = info4
c
        call prinf("after dgetri, info = *",info,1)
c        
        end



        subroutine scattering_bounding_ball(n,idxs,eval,bx,by,bz,br)
        implicit double precision (a-h,o-z)
        dimension eval(13,1),idxs(1)
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
c   idxs - an integer array specifying the indices of the input points
c   eval - an (13,*) array containing evaluation data for the points
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
        i1 = min(n,max(1,idxs(n/2)))
c
c       Find the point y most distant from it.
c
        i2   = 0
        dmax = 0
        do 1000 jj=1,n
        j  = idxs(jj)
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
        do 1100 jj=1,n
        j  = idxs(jj)
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
        do 2100 ii=1,n
        i  = idxs(ii)
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
        end



        subroutine scattering_bounding_ball0(n,idxs,eval,bx,by,bz,br)
        implicit double precision (a-h,o-z)
        dimension eval(13,1),idxs(1)
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
c   idxs - an integer array specifying the indices of the input points
c   eval - an (13,*) array containing evaluation data for the points
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
        do 1000 jj=1,n
c
        j = idxs(jj)
c
        bx = bx + eval(2,j)*dd
        by = by + eval(3,j)*dd
        bz = bz + eval(4,j)*dd
c
 1000 continue
c
        br = 0
        do 1100 ii=1,n
        i  = idxs(ii)
        dd = (bx-eval(2,i))**2 + (by-eval(3,i))**2 + (bz-eval(4,i))**2
        br = max(br,dd)
 1100 continue
c
        br = sqrt(br)
c
        end

