
c       octree_decomp
c
c       memory management
c
c       CHANGE NIN TO NEVALIN AND NOUT TO NEVALOUT?
c
c      -2. Make flags (2,n) instead of (4,n)
c
c      -1.  BLAS AND LAPACK EVERYTHING??  This might make parallelization more
c           difficult??
c
c       0.  Add auxilliary routines to form Ain and Aout to reuse for both???
c
c       1.  I should probably store near interactions (this will trade memory for speed)
c
c       2.  scattering_combine should combine more than one matrix
c
c       3.  the scattering_matrix code will construct the near lists eventually
c
        implicit double precision (a-h,o-z)
c
        dimension xout(13)
c
        double precision, allocatable :: w(:)
c
        double complex, allocatable   :: xx(:),yy(:),smatr(:,:)
        integer, allocatable          :: iin(:),iout(:)
c
        double complex zk,val,der,sum
c
        pi  = acos(-1.0d0)
c
c       Set discretization and solve parameters.
c
        norder  = 6
        nself   = 0
        dnear   = 2.00d0
        zk      = 0.00d0
        nlevels = 2
        alpha   = 1.00d0
        call surface_init(alpha)
c
c       Allocate a big work array.
c
        ld = 100 000 000
        lw = 100 000 000
c
 0001 continue
        if (allocated(w)) deallocate(w)
        allocate(w(lw),STAT=istat)
        if (istat. eq. 0) then
        lw = lw + ld
        goto 0001
        else
        lw = lw - ld
        endif
c
        if (allocated(w)) deallocate(w)
        allocate(w(lw))
c
        call prinf("in scattering, lw = *",lw,1)
c
        maxtris = 100 000
        itris   = 1
        ltris   = 6*maxtris
c
c       Construct a discretization structure.
c
        idisc = itris+ltris
        ldisc = lw-idisc
c
        call disc_init(ier,norder,nself,dnear,w(idisc),ldisc)
        if (ier .ne. 0) then
        call prinf("after disc_init, ier = *",ier,1)
        stop
        endif
c
c       Add a decomposition of the torus.
c
        x1 = 0
        x2 = 2*pi
        y1 = 0
        y2 = 2*pi
c
        iparam  = 1
c
        call triangulate(ier,iparam,x1,y1,x2,y2,nlevels,
     -    maxtris,ntris,w(itris))
c
        call disc_add(ier,w(idisc),iparam,ntris,w(itris))
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
c       Build a trivial decomposition structure.
c
        iidecomp = izs+lzs
        lidecomp = lw-iidecomp
c
        nlevels_decomp  = 3
        nbranch         = 2
        call decomp_trivial(ier,n,w(izs),nlevels_decomp,nbranch,
     -    w(iidecomp),lidecomp,lkeep)
        if (ier .ne. 0) then
        call prinf("after decomp_trivial, ier = *",ier,1)
        stop
        endif
        lidecomp = lkeep
c
c       Build the octree and near interaction lists.
c
        ioc = iidecomp+lidecomp
        loc = lw-ioc
c
        k   = 1000
        call elapsed(t1)
        call octree(ier,n,w(izs),k,w(ioc),loc,lkeep)
        call elapsed(t2)
c
        if (ier .ne. 0) then
        call prinf("after octree, ier = *",ier,1)
        stop
        endif
        call prin2("octree time = *",t2-t2,1)
c
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
        n0 = 10
        x0 = 2.00d0
        y0 = 0.00d0
        z0 = 0.00d0
        r0 = 0.1d0
c
        call sphere_data(n0,x0,y0,z0,r0,nin,w(ievalin))
c
        n0 = 10
        x0 = 100.0d0
        y0 = 0.0d0
        z0 = 0.0d0
        r0 = 1.0d0
c
        call sphere_data(n0,x0,y0,z0,r0,nout,w(ievalout))
c
c       Build the scattering matrix.
c
        naux       = 0
        eps        = 1.0d-12
c
        call scattering_matrix(ier,eps,zk,w(idisc),w(iidecomp),
     -    w(ioc),w(iinear),nin,w(ievalin),nout,w(ievalout),naux,w(iw2),
     -    lw2,lkeep)
        if (ier .ne. 0) then
        call prinf("after scattering_matrix, ier = *",ier,1)
        stop
        endif
c
c       Extract the scattering matrix from the work array.
c
        len    = w(iw2)
        nout   = w(iw2+1)
        nin    = w(iw2+2)
        iiout  = w(iw2+3)+iw2-1
        iiin   = w(iw2+4)+iw2-1
        ismatr = w(iw2+5)+iw2-1
c
        allocate(iout(nout),iin(nin),smatr(nout,nin))
c
        call scattering_imove(nout,w(iiout),iout)
        call scattering_imove(nin,w(iiin),iin)
c
        call scattering_move(2*nin*nout,w(ismatr),smatr)
c
        call prinf("iout = *",iout,nout)
        call prinf("iin = *",iin,nin)
c
c       Evaluate the incoming potential at the nodes of the incoming
c       skeleton.
c
        allocate(xx(nin),yy(nout))
        do 2000 ii=1,nin
        i = iin(ii)
        call incoming(w(ieval+(i-1)*13),zk,val,der)
        xx(ii) = der
 2000 continue
c
c       Apply the scattering matix.
c        
        do 2100 i=1,nout
        sum = 0
        do 2200 j=1,nin
        sum = sum + smatr(i,j)*xx(j)
 2200 continue
        yy(i) = sum
 2100 continue
c
c       Test the outgoing potential at a point distant from the
c       domain.
c
        sum = 0
        do 2300 ii=1,nout
        i = iout(ii)
        call ksingle(w(ievalout),w(ieval+(i-1)*13),zk,val)
        sum = sum + val*yy(ii)
 2300 continue
c
        call incoming(w(ievalout),zk,val,der)
c
        print *,""
        print *,sum
        print *,val
        print *,abs(sum-val)
c
        end


        subroutine incoming(x,zk,val,der)
        implicit double precision (a-h,o-z)        
        double complex zk,val,der
        dimension y(13)
c
c       Evaluate the incoming potential at a collection of evaluation
c       points.
c

c
c       Unit point souce at (2,0,0)
c
        y(1) = 1.0d0
        y(2) = 2.0d0
        y(3) = 0.0d0
        y(4) = 0.0d0
        call ksingle(x,y,zk,val)
        call ksingleprime(x,y,zk,der)
c
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

        entry surface_init(alpha0)
        alpha = alpha0
        end



        subroutine near_lists(ier,disc,oc,inear,linear,lkeep)
        implicit double precision (a-h,o-z)
        dimension disc(1),oc(1),inear(1)
c
        integer, allocatable :: iwork(:)
c
c       Construct the near interaction lists for the discretization
c       nodes.  For each point, the 
c
c                               Input Parameters:
c
c                               Output Parameters:
c 
c     
        ier = 0
c
c       Fetch some data from the disc array.
c
        ntris = disc(4)
        len   = disc(5)
        itris = disc(6)
        dnear = disc(7)
        nquad = disc(20)
        ixs   = disc(21)
        iys   = disc(22)
        iwhts = disc(23)
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

        double precision, allocatable :: xslege(:),whtslege(:)
c
c       Return evaluation data for a user-specified sphere.
c
c       The number of points generated is 
c
c
c

        
c
c       Fetch the n-point Legendre quadrature rule.
c
        allocate(xslege(n),whtslege(n))
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


        subroutine decomp_octree()
        implicit double precision (a-h,o-z)
        end


        subroutine decomp_trivial(ier,n,zs,nlevels,nbranch,idecomp,
     -   lidecomp,lkeep)
        implicit double precision (a-h,o-z)
        dimension zs(3,n),idecomp(1)
c
        double precision, allocatable :: xs(:)
c
c       This subroutine contains code for constructing a trivial
c       decomposition of the 
c
c                             Input Parameters:
c
c   n - the number of discretization nodes on the surface
c   zs - a (3,n) array specifying the coordinates of the discretization
c       nodes
c
c   nlevels - the number of levels for the decomposition
c   nbranch - 
c
c   lidecomp - the length of the user-supplied output array idecomp
c       
c                             Output Parameters:
c
c   ier - an error return code;
c       ier = 0   indicates successful execution
c       ier = 4   means that the idecomp array was of insufficient length
c
c   idecomp - upon successful return, this array will contain a structure
c       which describes a hiearchical decomposition of the discretization
c       nodes on the surface.  
c
c       It is organized as follows:
c
c   idecomp(1) - the number of levels in the decomposition                (nlevels)
c   idecomp(2) - the number of regions in the decomposition               (nregions)
c   idecomp(3) - a pointer to the iregions array                          (iiregions)
c   idecomp(4) - a pointer to the ilevels array                           (iilevels)
c
c   ilevels(1,j) - the number of region on level j of the decomposition   (nn)
c   ilevels(2,j) - the index of the first region on level j               (i1)
c
c   iregions(1,j) - the number of children of this region                 (nchildren)
c   iregions(2,j) - a pointer to the list of children                     (iichildren)
c   iregions(3,j) - the number of points in the region; this is set only  (npoints)
c                   if nchildren = 0
c   iregions(4,j) - a pointer to the list of points                       (iipts)
c 
c   lkeep - the length of the idecomp array produced by this subroutine
c
        ier = 0
c
        nregions   = (nbranch**(nlevels)-1)/(nbranch-1)
c
c       Allocate memory iregion sand ilevels arrays.
c
        iidxs      = 1000
        lidxs      = n
c
        iilevels   = iidxs+lidxs
        lilevels   = 2*nlevels
c
        iiregions  = iilevels+lilevels
        liregions  = 4*nregions
c
        iichildren = iiregions+liregions
        lichildren = nregions
c
        lkeep      = iichildren+lichildren
c
        if (lkeep .gt. lidecomp) then
        ier = 4
        return
        endif
c
c       Sort the points by x-coordinate in order to divide them into
c       two regions.
c
        allocate(xs(n))
        do 1000 i=1,n
        xs(i)              = zs(1,i)
        idecomp(iidxs+i-1) = i
 1000 continue
        call quicksort(n,xs,idecomp(iidxs+i-1))
c
c       Construct the ilevels and iregions arrays.
c
        call decomp_trivial0(n,nlevels,nregions,idecomp(iilevels),
     -   idecomp(iiregions),idecomp,iidxs,idecomp(iichildren),
     -   iichildren,nbranch)
c
        idecomp(1) = nlevels
        idecomp(2) = nregions
        idecomp(3) = iiregions
        idecomp(4) = iilevels
c
        end


        subroutine decomp_trivial0(n,nlevels,nregions,ilevels,iregions,
     -   idecomp,iidxs,ichildren,iichildren,nbranch)
        implicit double precision (a-h,o-z)
        dimension ilevels(2,nlevels),iregions(4,nregions),idecomp(1)
c
c       Build the levels array.
c
        i1      = 1
        do 1000 i=1,nlevels
        nn = nbranch**(i-1)
        ilevels(1,i) = nn
        ilevels(2,i) = i1
c
        i1 = i1 + nn
c
 1000 continue
c
c       Construct the regions array for all but the final level.  Make
c       the list of children.
c
        iregion = 1
        iptr    = iichildren
c
        do 1100 i=1,nlevels-1
        ichild     = ilevels(2,i+1)
        nn         = ilevels(1,i)
        nchildren  = nbranch
c
        do 1200 j=1,nn
        
        iregions(1,iregion) = nchildren
        iregions(2,iregion) = iptr
        do 1300 l=1,nchildren
        idecomp(iptr+l-1) = ichild
        ichild=ichild+1
 1300 continue
c
c        call prinf("iregion = *",iregion,1)
c        call prinf("ichildren = *",idecomp(iptr),nchildren)
c
        iptr    = iptr+nchildren
        iregion = iregion+1
 1200 continue
 1100 continue
c
        nn = n / (nbranch**(nlevels-1))
        i1 = iidxs
c
        do 2000 i=1,nbranch**(nlevels-1)
c
        iregions(1,iregion) = 0
        iregions(2,iregion) = 0
        iregions(3,iregion) = nn
        iregions(4,iregion) = i1
c
c        call prinf("iregion = *",iregion,1)
c        call prinf("idxs = *",idecomp(i1),nn)
c
        i1 =i1 + nn
        iregion=iregion+1
 2000 continue
c
        end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the
c       scattering matrix code proper.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       The exerior Neumann problem
c
c             \Delta u + k^2 u = 0         in  \Omega
c                         dN u = f         on  \partial\Omega,             (1)
c                      condition at infinity 
c 
c       where \Omega a suitable domain, \partial\Omega denotes the
c       boundary of \Omega and dN denotes the derivative w.r.t. the
c       outward pointing normal vector to \partial\Omega, can be
c       reformulated as the integral equation
c
c          T\sigma(x) =  1/2 \sigma(x)  + \int K(x,y) \sigma(y) ds(y)      (2)
c                                          \partial\Omega
c
c       with
c
c
c
c       Here, we are using a formulation which is not resistant to
c       resonances.  In the case
c
c       The compression scheme is based on recursively constructing
c       and combining scattering matrices.
c
c       EXPLAIN MEMORY MODEL
c
c       This file contains a code for compressing integral operators
c       of scattering theory on surfaces.
c
c
c          T\sigma(x) =  \sigma(x)  + \int K(x,y) \sigma(y) ds(y)          (1)
c                                        \Sigma
c
c
c
c       The following routines are user-callable:
c      
c   scattering_matrix - construct a scattering matrix; this subroutine
c
c   scattering_decomp_octree - construct a decomposition structure for
c       the surface using an octree
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       The following are internal subroutines which are important
c       enough to be documented carefully:
c
c   scattering_create - create a scattering matrix for a lowest-level
c       region of the surface decomposition
c
c   scattering_combine - combine scattering matrices on a set of 
c       regions in order to form a scattering matrix for their union
c
c   scattering_near_lists - build a structure which stores lists of
c       near targets and sources
c
c   scattering_targets- construct the list of near targets for a 
c       given set of discretization nodes
c
c   scattering_sources- construct the list of near sources for a 
c       given set of discretization nodes
c
c   scattering_garbage_collection - perform garbage collection on
c       the heap used by the scattering matrix code; that is to say,
c       the code moves the collection of extant scattering matrices
c       to the beginning of the work array thus hopefully making room
c       at the end for new scattering matrices.  It is called by the
c       scattering_matrix code if insufficient
c
c   scattering_bounding_ball0 - construct an approximate minimum 
c       bounding ball for a collection of points in R^3 using
c       a crude method
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



        subroutine scattering_matrix(ier,eps,zk,disc,idecomp,oc,inear,
     -    nin,evalin,nout,evalout,naux,w,lw,lkeep)
        implicit double precision (a-h,o-z)
        double complex zk,val,sum
        dimension disc(1),idecomp(1),oc(1),inear(1),w(1)
        dimension evalin(13,1),evalout(13,1)
c
c       This subroutine is the user-callable interface to the
c       accelerated direct solver.  It constructs a scattering matrix
c       for the problem (1).  That is to say, given user-specified
c       contours R_in and R_out such that the right-hand side f in (1)
c       is of the form
c
c         f(x) = \int D_y G(x,y) \sigma(y) ds(y)
c
c       and the 
c
c         u(x) = \int D_
c
c       constructs a matrix which takes the values of f(x) at a 
c       collection of points
c
c
c
c        
c 
c
c
c       The user must specify a contour \Gamma_in which generates 
c       the space of possible incoming potentials and a contour
c       \Gamma_out which generates the 
c       
c
c       This subroutine is 
c
c                             Input Parameters:
c
c  zk - the wavenumber for the problem
c  disc - structure describing the discretization of the scatter's
c       boundary; see discretize.f for a description of this structure
c  idecomp - structure describing a hiearchical decomposition of the 
c       discretization nodes on the boundary of the scatter; see
c       below for a description of this structure
c
cccccccccccccccc   to be eliminated and constructed here?? cccccccccc
c
c  oc - a structure describing an octree containing the discretization
c       nodes on the boundary of the scatterer (see octree.f).  This
c       octree is used to perform range searches
c
c        This structure is used to 
c
c
c  inear  a structure containing the near interaction lists 
c        DEFINITELY ELIMINATE THIS
c
cccccccccccccccc
c
c  nin - the number of quadrature nodes on the surface which specifies
c       incoming potentials (the "incoming surface")
c  evalin - the evaluation data for the incoming surface
c
c  nout - the number of quadrature nodes on the surface which specifies
c       outgoing potentials (the "outgoing surface")
c  evalout - a (13,nout) array specifying the evaluation data for the  
c      outgoing surface
c
c  naux  - this parameter controls the number of points per wavelength
c      which are used to discretize auxilliary contours
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
c       ier = 128   means that an auxilliary contour required too many
c                   discretization nodes
c       ier = 400   means that the number of children of a region in
c                   the decomposition exceeded 100
c       ier = 1024  means that a matrix evaluation failed
c
c
c   w - upon successful execution, the beginning of the work array will
c       contain a smatrix structure describing a scattering matrix for the
c       exterior Neumann problem (1)
c
c       The structure used to describe scattering matrices is organized
c       as a header containing pointers which are relative to the beginning
c       of the structure.  The header has the following information:
c
c   w(1) - the length of all smatr data including this element            (len)
c   w(2) - the number of points in the outgoing skeleton                  (nout)
c   w(3) - the number of points in the incoming skeleton                  (nin)
c   w(4) - a pointer to the list of the outgoing skeleton nodes           (iout)
c   w(5) - a pointer to the list of the incoming skeleton nodes           (iiin)
c   w(6) - a pointer into the work array to the nout x nin scattering     (ismatr)
c          matrix
c
c      THESE ARE POTENTIAL FUTURE ADDITIONS:
c
c   w(7)  - the number of near targets                                    (nneartars)
c   w(8)  - a pointer to the list of near targets                         
c   w(9)  - a pointer to the (nneartars,nout) array containing the
c           near target interactions
c
c   w(10) - the number of near sources                                    (nnearsrcs)
c   w(11) - a pointer to the list of near sources
c   w(12) - a pointer to the (nnearsrcs,nin) array containg the near
c           source interactions
c
c
c   ???Store the indices????
c
c   lkeep - the size of the structure describing the resulting scattering matrix

c
c       Fetch some data from the discretization structure.
c
        call disc_info(disc,ntris,nquad,n,dnear,lkeep)
c
c       And from the decomposition structure.
c
        nlevels   = idecomp(1)
        nregions  = idecomp(2)
        iiregions = idecomp(3)
        iilevels  = idecomp(4)
c
        call prinf("in scattering_matrix, nlevels = *",nlevels,1)
        call prinf("in scattering_matrix, nregions = *",nregions,1)
c
c       Construct an octree for range searches.
c

c
c       Construct the near interaction lists using the octree.
c


c
c       Fetch some data from the near interaction structure.
c
        iinearptrs     = inear(2)
        iinearidxs     = inear(3)
c
c       Allocate some global variables from the work array.
c
        iiflags        = 1
        liflags        = 2*n
c
        iinearsrcs     = iiflags+liflags
        linearsrcs     = 2*n
c
        iineartars     = iinearsrcs+linearsrcs
        lineartars     = 2*n
c
        iidxs          = iineartars+lineartars
        lidxs          = n
c
        iidxsin        = iidxs+lidxs
        lidxsin        = n
c
        iidxsout       = iidxsin+lidxsin
        lidxsout       = n
c
        iiwork         = iidxsout+lidxsout
        liwork         = 2*n
c
        ievalaux       = iiwork+liwork
        levalaux       = 13*100 000
c
        ieval          = ievalaux+levalaux
        leval          = 13*n
c
        iisptrs        = ieval+leval
        lisptrs        = 2*nregions
c
        iiout          = iisptrs+lisptrs
        liout          = n
c
        iiin           = iiout+liout
        liin           = n
c
        iw2            = iiin+liin
        lw2            = lw-iw2
c
        if (lw2 .le. 0) then
        ier = 4
        return
        endif
c
c       Call an auxilliary subroutine to form arrays.
c
        call scattering_matrix0(ier,n,eps,zk,naux,nlevels,nregions,
     -   idecomp,idecomp(iiregions),idecomp(iilevels),dnear,disc,oc,
     -   inear(iinearptrs),inear(iinearidxs),w(iiflags),
     -   w(iinearsrcs),w(iineartars),w(iidxs),w(iidxsin),w(iidxsout),
     -   w(iiwork),w(ievalaux),w(ieval),w(iisptrs),w(iiout),w(iiin),nin,
     -   evalin,nout,evalout,w(iw2),lw2,idecomp(iiregions),
     -   idecomp(iilevels))
c
c       Copy the top-level scattering matrix to the beginning of the
c       work array.
c
c
c       ------------- WARNING ------------------------------------------
c
c       BECAUSE THE TARGET AND SOURCE BLOCKS MAY OVERLAP, THIS OPERATION
c       CANNOT BE PERFORMED USING A CALL TO AN EXTERNAL SUBROUTINE.  
c
c       DUE TO FORTRAN'S ALIASING RULES USING SUCH AN APPROACH CAN 
c       RESULT IN INCORRECT RESULTS.
c
c       ----------------------------------------------------------------
c
        len   = w(iw2)
        lkeep = len
        do 3000 i=1,len
        w(i) = w(iw2+i-1)
 3000 continue
c
        end


        subroutine scattering_matrix0(ier,n,eps,zk,naux0,nlevels,
     -   nregions,idecomp,iregions,ilevels,dnear,disc,oc,inearptrs,
     -   inearidxs,iflags,inearsrcs,ineartars,idxs,idxsin,idxsout,
     -   iwork,evalaux,eval,isptrs,iout,iin,nin,evalin,nout,evalout,
     -   w,lw)
        implicit double precision (a-h,o-z)
        double complex zk,val,sum
        dimension disc(1),oc(1),inearptrs(4,1),inearidxs(1),iflags(2,1)
        dimension idecomp(1),ilevels(2,1),iregions(4,1)
        dimension inearsrcs(1),ineartars(1),idxs(1),iwork(1)
        dimension idxsout(1),idxsin(1)
        dimension evalaux(13,1),eval(13,1),isptrs(2,1)
        dimension evalin(13,1),evalout(13,1)
        dimension iin(1),iout(1),w(1)
c
        dimension ichildptrs(1000)
c
        ier = 0
        pi  = acos(-1.0d0)
c
c       Fetch some data from the disc array.
c
        call disc_data(disc,eval)
c
c       Initialize the scattering matrix pointers.
c
        do 0100 i=1,n
        isptrs(1,i) = 0
        isptrs(2,i) = 0
 0100 continue
c
c       Initialize the pointer used to keep track of the first free
c       entry in the work array.
c
        ifree = 1
c
c       Mark all source and target nodes as present.
c
c       flags(1,i) - source 
c       flags(2.i) - target
c     
        do 1000 i=1,n
        iflags(1,i) = 1
        iflags(2,i) = 1
 1000 continue
c
c       Construct scattering matrices for the boxes at the lowest
c       level of the decomposition.
c       
        nn        = ilevels(1,nlevels)
        iregion0  = ilevels(2,nlevels)
c
        call prinf("in scattering_matrix, level = *",nlevels,1)
        call prinf("in scattering_matrix, nn = *",nn,1)
c
        do 1100 iregion=iregion0,iregion0+nn-1
        call prinf("in scattering_matrix, iregion = *",iregion,1)
c
c       Find the indices in the region and construct a bounding ball for
c       the corresponding points.
c
        nidxs = iregions(3,iregion)
        iptr  = iregions(4,iregion)
c
        call scattering_move(nidxs,idecomp(iptr),idxs)
c
c       Mark the nodes as absent.
c
        do 1200 i0=1,nidxs
        i = idxs(i0)
        iflags(1,i) = 0
        iflags(2,i) = 0
 1200 continue
c
        call scattering_targets(n,oc,eval,dnear,nidxs,idxs,iwork,
     -    inearptrs,inearidxs,nneartars,ineartars,iflags)
c
        call scattering_sources(n,oc,eval,dnear,nidxs,idxs,iwork,
     -    inearptrs,inearidxs,nnearsrcs,inearsrcs,iflags)
c
        call prinf("in scattering_matrix, nidxs = *",nidxs,1)
        call prinf("in scattering_matrix, idxs = *",idxs,nidxs)
c        call prinf("nnearsrcs = *",nnearsrcs,1)
c        call prinf("inearsrcs = *",inearsrcs,nnearsrcs)
c        call prinf("nneartars = *",nneartars,1)
c        call prinf("ineartars = *",ineartars,nneartars)
c
c       Build the evaluation data for the auxilliary contour, but only
c       if we are not on the top level.
c
c       Compute the wavelength.

        wl   = 2*pi/abs(zk)
        dd   = sqrt(4*pi*br0**2 / wl**2)
c
        naux = 0
        nn0  = naux0*max(1.0d0,dd)
c
        if (nn0 .gt. 100 000) then
        ier = 128
        return
        endif
c
c       Construct the auxilliary contour.
c
        call scattering_bounding_ball0(nidxs,idxs,eval,bx,by,bz,br)
c     
        call scattering_sphere_data(nn0,bx,by,bz,br0,naux,evalaux)
        call prinf("naux = *",naux,1)
c
c       Call a subroutine to construct the scattering matrix structure 
c       in the work array, record its position and advance the heap.
c
        lfree = lw-ifree
c
        call scattering_create(ier,disc,eps,zk,eval,nidxs,idxs,nin,
     -    evalin,nout,evalout,naux,evalaux,nneartars,ineartars,
     -    nnearsrcs,inearsrcs,krankout,iout,krankin,iin,w(ifree),lfree,
     -    lkeep)
        if (ier .ne. 0) return
c
c       Set flags; outgoing nodes are potential sources which incoming
c       nodes are potential targets.
c
        do 1800 i=1,krankin
        iflags(2,iin(i)) = 1
 1800 continue
c
        do 1900 i=1,krankout
        iflags(1,iout(i)) = 1
 1900 continue
c
        if (ier .ne. 0) return
        isptrs(1,iregion) = ifree
        isptrs(2,iregion) = lkeep
c
        ifree = ifree + lkeep
 1100 continue
c
c       Combine scattering matrices level-by-level.
c
        do 2000 ilevel=nlevels-1,1,-1
        call prina("*")
        call prina("--------------------------*")
        call prinf("in scattering_matrix, ilevel = *",ilevel,1)
c
        nn       = ilevels(1,ilevel)
        iregion0 = ilevels(2,ilevel)
c
        do 2100 iregion=iregion0,iregion0+nn-1
        nchildren = iregions(1,iregion)
        iptr      = iregions(2,iregion)
        call prinf("in scattering_matrix, iregion = *",iregion,1)
        call prinf("in scattering_matrix, children = *",idecomp(iptr),
     -   nchildren)
c
c       Collect pointers to the children's scattering matrices and
c       construct the list of incoming and outgoing indices.
c
        if (nchildren .gt. 100) then
        ier = 400
        return
        endif
c
        nidxsout = 0
        nidxsin  = 0
        do 2200 i=1,nchildren        
        ichild        = idecomp(iptr+i-1)
        iptr0          = isptrs(1,ichild)
        ichildptrs(i) = iptr0
c
c
        krankout      = w(iptr0+1)
        krankin       = w(iptr0+2)
        iiout         = w(iptr0+3)+iptr0-1
        iiin          = w(iptr0+4)+iptr0-1
c
        call scattering_imove(krankout,w(iiout),idxsout(nidxsout+1))
        call scattering_imove(krankin,w(iiin),idxsin(nidxsin+1))
c
        nidxsout = nidxsout+krankout
        nidxsin  = nidxsin+krankin
 2200 continue
c
c
        call prinf("idxsin = *",idxsin,nidxsin)
        call prinf("idxsout = *",idxsout,nidxsout)
c
c       Incoming nodes are potential targets while outgoing nodes
c       are potential sources.
c
        do 2300 i0=1,nidxsout
        i           = idxsout(i0)
        iflags(1,i) = 0
 2300 continue
c
        do 2400 i0=1,nidxsin
        i           = idxsin(i0)
        iflags(2,i) = 0
 2400 continue
c      
c       Construct the near lists; we need the list of near targets
c       for the source points and near sources for the target points.
c
        call scattering_targets(n,oc,eval,dnear,nidxsout,idxsout,iwork,
     -    inearptrs,inearidxs,nneartars,ineartars,iflags)
c
        call scattering_sources(n,oc,eval,dnear,nidxsin,idxsin,iwork,
     -    inearptrs,inearidxs,nnearsrcs,inearsrcs,iflags)
c
        call prinf("nnearsrcs = *",nnearsrcs,1)
        call prinf("inearsrcs = *",inearsrcs,nnearsrcs)
c
        call prinf("nneartars = *",nneartars,1)
        call prinf("ineartars = *",ineartars,nneartars)
c
c
c       Build the evaluation data for the auxilliary contour, but only
c       if we are not on the top level.
c
c       Compute the wavelength.

        wl   = 2*pi/abs(zk)
        dd   = sqrt(4*pi*br0**2 / wl**2)
c
        naux = 0
        nn0  = naux0*max(1.0d0,dd)
c
        if (ilevel .eq. 1) nn0 = 0
c
        if (nn0 .gt. 100 000) then
        ier = 128
        return
        endif
c
c       Construct the auxilliary contour.
c
        call scattering_bounding_ball0(nidxs,idxs,eval,bx,by,bz,br)
        call scattering_sphere_data(nn0,bx,by,bz,br0,naux,evalaux)
        call prinf("naux = *",naux,1)

        lfree = lw-ifree
c
        call scattering_combine(ier,disc,eps,zk,n,eval,nchildren,
     -    ichildptrs,nidxsin,idxsin,nidxsout,idxsout,nneartars,
     -    ineartars,nnearsrcs,inearsrcs,
     -    nin,evalin,nout,evalout,naux,evalaux,
     -    krankout,iout,krankin,iin,w,w(ifree),lfree,
     -    lkeep)
c
        if (ier .ne. 0) return
c
        isptrs(1,iregion) = ifree
        isptrs(2,iregion) = lkeep
c
        ifree = ifree+lkeep
c
c       Mark the children as free.
c
        do 2500 i=1,nchildren
        ichild = idecomp(iptr+i-1)
        isptrs(1,ichild) = 0
 2500 continue
c
c       Set flags; outgoing nodes are potential sources which incoming
c       nodes are potential targets.
c
        do 2600 i=1,krankin
        iflags(2,iin(i)) = 1
 2600 continue
c
        do 2700 i=1,krankout
        iflags(1,iout(i)) = 1
 2700 continue
        
 2100 continue
 2000 continue
c
c       Copy out the top-level scattering matrix.
c
c       ------------- WARNING ------------------------------------------
c
c       BECAUSE THE TARGET AND SOURCE BLOCKS MAY OVERLAP, THIS OPERATION
c       CANNOT BE PERFORMED USING A CALL TO AN EXTERNAL SUBROUTINE.  
c
c       DUE TO FORTRAN'S ALIASING RULES USING SUCH AN APPROACH CAN 
c       RESULT IN INCORRECT RESULTS.
c
c       ----------------------------------------------------------------
c
        ismatr = isptrs(1,1)
        len    = isptrs(2,1)
c
        do 3000 i=1,len
        w(i) = w(ismatr+i-1)
 3000 continue
        end



        subroutine scattering_garbage_collection(heap,isptrs,iwork,
     -    ifree)
        implicit double precision (a-h,o-z)
c
c       This routine moves all extant scattering matrices to the
c       beginning of the work array, thus (hopefully) freeing space
c       for
c
c
c                            Input Parameters:
c
c
c                            Output Parameters:
c
        end


        subroutine scattering_combine(ier,disc,eps,zk,n,eval,k,
     -   iptrs,nidxsin,idxsin,nidxsout,idxsout,nneartars,ineartars,
     -   nnearsrcs,inearsrcs,nin,evalin,nout,evalout,naux,evalaux,
     -   krankout,iout,krankin,iin,heap,w,lw,lkeep)
        implicit double precision (a-h,o-z)
        dimension disc(1),eval(13,1),iptrs(1),w(1),heap(1)
        dimension evalaux(13,1),inearsrcs(1),ineartars(1)
        dimension evalin(13,nin),evalout(13,nout)
        dimension iout(1),iin(1),idxsin(1),idxsout(1)
        double complex zk,val,sum
c
        dimension kranksin(100),kranksout(100)
c
c       Combine a collection of scattering matrices S1,...,Sk given
c       on regions R1,...,Rk to form a scattering matrix for the union
c       of the R_j.
c
c       This is down by first constructing the matrix 
c
c                    ( S1           ) ( 0    T12  ... T1k )
c          A = I +   (    S2        ) ( T21  0    ... T2k )
c                    (       ...    ) (           ...     )
c                    (           Sk ) ( Tk1  Tk2  ...  0  )
c
c       where T_ij is the interaction matrix with targets on the
c       region Ri and sources on the region Rj.  The matrix A is
c       then inverted and the uncompressed scattering matrix
c
c                     (S1           )
c         U = A^{-1}  (   S2        )                                       (3)
c                     (      ...    )
c                     (          Sk )
c
c       is computed.  
c
c       In the next step, matrices Ain and Aout which account for incoming 
c       and outgoing potentials on the combined region are formed and 
c       factored as 
c
c         Ain   = R_in B_in
c 
c         Aout  = Bout R_out
c
c       where Bin and Bout are submatrices of Ain and Aout, 
c       respectively.  Then the compressed scattering matrix
c
c         S = R_out U R_in
c
c       is formed.
c
c       It is possible that the inversion of the matrix A could be
c       sped up using various Sherman-Morrison-Woodbury type tricks,
c       but for the time-being this routine inverts it via dense
c       LAPACK routines.
c
c
c                            Input Parameters:
c
c   disc - a structure describing the discretization of the boundary
c       surface
c   eps - precision for the scattering matrix computations
c   zk - wavenumber for the problem
c   n - the number of discretization nodes on the boundary
c   eval - a (13,n) array with the evaluation data
c   
c   k - the number of scattering matrices which are to be combined
c   iptrs - a integer array of length k whose entries are pointers
c       into the heap to the scattering matrices which are to be
c       combined
c
c   heap - the pointer to the beginning of the work array containing
c       the scattering matrices
c
c   lw - the length of the user-supplied work array w
c
c                            Output Parameters:
c
c   ier - an error return code;
c       ier = 0  indicates successful execution
c       ier = 4  means that the work array w was of insufficient 
c                length
c
c   w - upon successful execution, 
c
c   krankout - the number of 
c   iout - this user-supplied integer array will contain the indices
c       of the outgoing nodes for the combined region
c
c   krankin - 
c
c
c   lkeep - the length of the newly created scattering matrix structure
c
c     
c
c   w(1) - the length of all smatr data including this element            (len)
c   w(2) - the number of points in the outgoing skeleton                  (nout)
c   w(3) - the number of points in the incoming skeleton                  (nin)
c   w(4) - a pointer to the list of the outgoing skeleton nodes           (iout)
c   w(5) - a pointer to the list of the incoming skeleton nodes           (iiin)
c   w(6) - a pointer into the work array to the nout x nin scattering     (ismatr)
c
        ier   = 0
        lkeep = 0
c
        call prinf("in scattering_combine, nidxsin = *",nidxsin,1)
        call prinf("idxsin = *",idxsin,nidxsin)
        call prinf("in scattering_combine, nidxsout = *",nidxsout,1)
        call prinf("idxsout = *",idxsout,nidxsout)
c
c       Allocate memory.
c
        nin0      = nin  + naux + nnearsrcs
        nout0     = nout + naux + nneartars
c
        ismatr    = 1
        lsmatr    = 2*nidxsout*nidxsin
c
        iamatr    = ismatr+lsmatr
        lamatr    = 2*nidxsout*nidxsout
c
        iaint     = iamatr+lamatr
        laint     = 2*nidxsin*nin0
c
        irint     = iaint+laint
        lrint     = 2*nidxsin*nin0
c
        iatars    = irint+lrint
        latars    = 2*nneartars*nidxsin
c
        iasrcs    = iatars+latars
        lasrcs    = 2*nnearsrcs*nidxsout
c
        iaout     = iasrcs+lasrcs
        laout     = 2*nidxsout*nout0
c
        irout     = iaout+laout
        lrout     = 2*nidxsin*nout0
c
        iicolsin  = irout+lrout
        licolsin  = nidxsin
c
        iicolsout = iicolsin+licolsin
        licolsout = nidxsout
c
        iw2       = iicolsout+licolsout
        lw2       = lw-iw2
c
        if (lw2 .lt. 0) then
        ier = 4
        return
        endif
c
c       Call an auxilliary routine to shape arrays.
c
        call scattering_combine0(ier,eps,zk,disc,eval,
     -   k,nidxsout,idxsout,nidxsin,idxsin,w(iamatr),
     -   w(ismatr),heap,iptrs,w(iw2),lw2,krankout,iout,krankin,iin,
     -   naux,evalaux,nout0,nneartars,ineartars,nout,evalout,
     -   w(iaout),w(irout),w(iatars),w(iicolsout),
     -   nin0,nnearsrcs,inearsrcs,nin,evalin,w(iaint),w(irint),
     -   w(iasrcs),w(iicolsin))
c
c       Copy out the scattering matrix.
c
        ismatr = 7
        lsmatr = 2*krankout*krankin
c
        iiout  = ismatr+lsmatr
        liout  = krankout
c
        iiin   = iiout+liout
        liin   = krankin
c
        lkeep = iiin+liin
c
        if (lkeep .gt. lw) then
        ier = 4
        return
        endif
c
        do 3000 i=1,lsmatr
        w(ismatr+i-1) = w(iw2+i-1)
 3000 continue
c
        call scattering_imove(krankout,iout,w(iiout))
        call scattering_imove(krankin,iin,w(iiin))
c
        w(1) = lkeep
        w(2) = krankout
        w(3) = krankin
        w(4) = iiout
        w(5) = iiin
        w(6) = ismatr
c
        end


        subroutine scattering_combine0(ier,eps,zk,disc,eval,
     -   k,nidxsout,idxsout,nidxsin,idxsin,amatr,smatr,
     -   heap,iptrs,w,lw,krankout,iout,krankin,iin,
     -   naux,evalaux,nout0,nneartars,ineartars,nout,evalout,
     -   aout,rout,atars,icolsout,
     -   nin0,nnearsrcs,inearsrcs,nin,evalin,aint,rint,asrcs,
     -   icolsin)
        implicit double precision (a-h,o-z)
        dimension disc(1),heap(1),iptrs(1),w(1),iout(1),iin(1)
        double complex amatr(nidxsout,nidxsout),smatr(nidxsout,nidxsin)
        double complex aout(nout0,nidxsout),rout(nout0,nidxsout)
        double complex atars(nneartars,nidxsout)
        double complex aint(nin0,nidxsin),rint(nidxsin,1)
        double complex asrcs(nidxsin,nnearsrcs)
        dimension ineartars(1),evalout(13,1),evalaux(13,1),icolsout(1)
        dimension idxsin(1),idxsout(1),eval(13,1)
c
        dimension inearsrcs(1),evalin(13,1),icolsin(1)

        double complex zk,sum,val
c
        ier = 0
c
c       Form the matrix A.
c
        amatr = 0
        i1    = 1
        i2    = 1
        do 1000 i=1,k
        iptri     =  iptrs(i)
        krankouti =  heap(iptri+1)
        krankini  =  heap(iptri+2)
        iiini     =  heap(iptri+4)+iptri-1
        ismatri   =  heap(iptri+5)+iptri-1
c
        j1        = 1
        j2        = 1
        do 1100 j=1,k
c
        iptrj     =  iptrs(j)
        krankoutj =  heap(iptrj+1)
        iioutj    =  heap(iptrj+3)+iptrj-1
c
        if (i .ne. j) then
c
c       Evaluate T_ij which captures the outgoing potential from region 
c       j to region i.
c
        iw2       = 2*krankini*krankoutj+1
        lw2       = lw-iw2
c
        ikernel = 3
        call disc_eval(ier,disc,ikernel,zk,krankini,heap(iiini),
     -    krankoutj,heap(iioutj),w,w(iw2),lw2)
        if (ier .ne. 0) return
c
c       Multipole by S_i on the left.
c
        call scattering_mult(krankouti,krankoutj,krankini,heap(ismatri),
     -    w,amatr(i2,j2),nidxsout)
        endif
        j2 = j2+krankoutj
 1100 continue
        i2 = i2+krankouti
 1000 continue
c
        do 1200 i=1,nidxsout
        amatr(i,i) = amatr(i,i)+1.0d0
 1200 continue
c
c       Invert A.
c
        call lapack_invert(ier,nidxsout,amatr,w,lw)
        if (ier .ne. 0) return
c
        call prinf("after lapack_invert, ier = *",ier,1)
c
c       Form the uncompressed scattering matrix u.
c
        i1=1
        j1=1
        do 1300 j=1,k
        iptr     =  iptrs(j)
        krankout =  heap(iptr+1)
        krankin  =  heap(iptr+2)
        ismatr   =  heap(iptr+5)+iptr-1
c
        call scattering_mult(nidxsout,krankin,krankout,amatr(1,i1),
     -    heap(ismatr),smatr(1,j1),nidxsout)
c
        i1=i1+krankout
        j1=j1+krankin
 1300 continue
c
c       Form the matrix Aout.
c
c
c       Form the matrix Aout(nidxsout,nout0) and factor it as 
c
c         Aout(nout0,nidxsout) = Bout(nout0,krankout) * Rout(krankout,nidxsout)
c
        ikernel = 3
        call disc_eval(ier,disc,ikernel,zk,nneartars,ineartars,
     -    nidxsout,idxsout,atars,w,lw)
        if (ier .ne. 0) return
c
        do 2000 j=1,nidxsout
        do 2100 i=1,nout
        call ksingle(evalout(1,i),eval(1,idxsout(j)),zk,val)
        aout(i,j) = val
 2100 continue
        do 2200 i=1,naux
        call ksingleprime(evalaux(1,i),eval(1,idxsout(j)),zk,val)
        aout(nout+i,j) = val
 2200 continue
        do 2300 i=1,nneartars
        aout(nout+naux+i,j) = atars(i,j)
 2300 continue
 2000 continue
c
        if (lw .le. 2*(3*nidxsout+1)*nout+1) then
        ier = 4
        return
        endif
c
        call factor_right(eps,nout0,nidxsout,aout,krankout,icolsout,
     -   rout,w)
        do 2500 i=1,krankout
        iout(i) = idxsout(icolsout(i))
 2500 continue
c
c       Form the matrix Aint (nin0,nidxs) and factor it as 
c
c            Aint (nin0,nidxs) = Aint (nin0,krankin) Rint (krankin, nidxs)
c
        ikernel = 3
        call disc_eval(ier,disc,ikernel,zk,nidxsin,idxsin,nnearsrcs,
     -    inearsrcs,asrcs,w,lw)
        if (ier .ne. 0 ) return
c
        do 3000 j=1,nidxsin
        do 3100 i=1,nin
        call ksingleprime(eval(1,idxsin(j)),evalin(1,i),zk,val)
        aint(i,j) = val
 3100 continue
        do 3200 i=1,naux
        call ksingle(eval(1,idxsin(j)),evalaux(1,i),zk,val)
        aint(nin+i,j) = val
 3200 continue
        do 3300 i=1,nnearsrcs
        aint(nin+naux+i,j) = asrcs(j,i)
 3300 continue
 3000 continue
c
        if (lw .le. 2*(3*nidxsin+1)*nin+11) then
        ier = 4
        return
        endif
c
        call factor_right(eps,nin0,nidxsin,aint,krankin,icolsin,rint,w)
c
        do 3600 i=1,krankin
        iin(i) = idxsin(icolsin(i))
 3600 continue
        call prinf("in scattering_combine, iout = *",iout,krankout)
        call prinf("in scattering_combine, iin = *",iin,krankin)
c
c       Multiply smatr on the left by Rout and on the right by Rin.
c
        ib = 1
        lb = 2*krankout*krankin+1 000 000
c
        ia = ib+lb
        la = 2*krankout*nidxsin
c
        lneeded = ia+la
c
        if (lneeded  .gt. lw) then
        ier = 4
        return
        endif
c
        call apply_left(nidxsout,nidxsin,krankout,icolsout,
     -    rout,smatr,w(ia))
c
        call apply_rightt(nidxsin,krankout,krankin,icolsin,rint,w(ia),
     -   w(ib))
c
c       Form the list of incoming columns.
c
c        krankin = nidxsin
c        call scattering_imove(krankin,idxsin,iin)
c        call prinf("in scattering_combine, iin = *",iin,krankin)
c
c$$$c
c$$$        krankin  = nidxsin
c$$$c
c$$$        idx2     = 1
c$$$c
c$$$        do 4000 i=1,k
c$$$        iptr      = iptrs(i)
c$$$        krankin0  = heap(iptr+2)
c$$$        iiin0     = heap(iptr+4)+iptr-1
c$$$c
c$$$        call scattering_imove(krankin0,heap(iiin0),iin(idx2))
c$$$c
c$$$        idx2 = idx2+krankin0
c$$$ 4000 continue
c$$$c
c$$$        call prinf("in scattering_combine, iin = *",iin,krankin)
c$$$        stop
c
        end



        subroutine scattering_targets(n,oc,eval,dnear,nidxs,idxs,iwork,
     -    inearptrs,inearidxs,nneartars,ineartars,iflags)
        implicit double precision (a-h,o-z)
        dimension oc(1),eval(13,1),idxs(1),iwork(1),inearptrs(4,1)
        dimension inearidxs(1),ineartars(1),iflags(2,n)
c        
c       Construct the list of near targets for a user-specified 
c       collection of points.
c
c                              Input Parameters:
c
c   n - the total number of discretization nodes
c   oc - a structure describing an octree containing the discretization
c       nodes
c   eval - the evaluation data for the discretization nodes
c   dnear - the near interaction radius
c   nidxs - the number of input points
c   idxs - the indices of the input points
c   inearptrs - the (4,1) pointers array from the inear structure
c   inearidxs - the list of indices from the inear structure
c   iflags - the (2,n) array specifying which incoming and outgoing
c       nodes are still extant
c
c                             Output Parameters:
c
c   nneartars - the number of near targets
c   ineartars - this user-supplied array will contain the list of
c       near points 
c
c       IMPORTANT: the ineartars array must be of length at least
c       2*n!!!
c
        call scattering_bounding_ball0(nidxs,idxs,eval,bx,by,bz,br)
c
c       Add the nodes which are near by virtue of geometry to the
c       lists of near target and source nodes.
c
        nneartars = 0
c
        br0  = br*dnear
        call ocrange(oc,bx,by,bz,br0,n1,iwork)
c
        do 1000 i0=1,n1
        i = iwork(i0)
        if (iflags(2,i) .eq. 1) then
        nneartars            = nneartars+1
        ineartars(nneartars) = i
        endif
 1000 continue
c
c       Add the target nodes which are near by virtue of discretization
c       to the list.
c
        do 1100 i0=1,nidxs
        i    = idxs(i0)
        nn   = inearptrs(3,i)
        iptr = inearptrs(4,i)
        do 1200 j0=iptr,iptr+nn-1
        j = inearidxs(j0)           
        if (iflags(2,j) .eq. 1) then
        nneartars = nneartars+1
        ineartars(nneartars) = j
        endif
 1200 continue
        call iduplicates(nneartars,ineartars)
 1100 continue
c
        end



        subroutine scattering_sources(n,oc,eval,dnear,nidxs,idxs,iwork,
     -    inearptrs,inearidxs,nnearsrcs,inearsrcs,iflags)
        implicit double precision (a-h,o-z)
        dimension oc(1),eval(13,1),idxs(1),iwork(1),inearptrs(4,1)
        dimension inearidxs(1),inearsrcs(1),iflags(2,n)
c        
c       Construct the list of near sources for a user-specified 
c       collection of points.
c
c                              Input Parameters:
c
c   n - the total number of discretization nodes
c   oc - a structure describing an octree containing the discretization
c       nodes
c   eval - the evaluation data for the discretization nodes
c   dnear - the near interaction radius
c   nidxs - the number of input points
c   idxs - the indices of the input points
c   inearptrs - the (4,1) pointers array from the inear structure
c   inearidxs - the list of indices from the inear structure
c   iflags - the (2,n) array specifying which incoming and outgoing
c       nodes are still extant
c
c                             Output Parameters:
c
c   nnearsrcs - the number of near sources
c   inearsrcs - this user-supplied array will contain the list of
c       near sources
c
c       IMPORTANT: the inearsrcs array must be of length at least
c       2*n!!!
c
        call scattering_bounding_ball0(nidxs,idxs,eval,bx,by,bz,br)
c
c       Add the nodes which are near by virtue of geometry to the
c       lists of near target and source nodes.
c
        nnearsrcs = 0
c
        br0  = br*dnear
        call ocrange(oc,bx,by,bz,br0,n1,iwork)
c
        do 1000 i0=1,n1
        i = iwork(i0)
        if (iflags(1,i) .eq. 1) then
        nnearsrcs            = nnearsrcs+1
        inearsrcs(nnearsrcs) = i
        endif
 1000 continue
c
c       Add the target nodes which are near by virtue of discretization
c       to the list.
c
        do 1100 i0=1,nidxs
        i    = idxs(i0)
        nn   = inearptrs(1,i)
        iptr = inearptrs(2,i)
        do 1200 j0=iptr,iptr+nn-1
        j = inearidxs(j0)           
        if (iflags(1,j) .eq. 1) then
        nnearsrcs = nnearsrcs+1
        inearsrcs(nnearsrcs) = j
        endif
 1200 continue
        call iduplicates(nnearsrcs,inearsrcs)
 1100 continue
c
        end




        subroutine scattering_mult(n,m,k,a,b,c,ldc)
        implicit double precision (a-h,o-z)
        double complex a(n,k),b(k,m),c(ldc,1)
        double complex sum
        do 1000 i=1,n
        do 1100 j=1,m
        sum = 0
        do 1200 l=1,k
        sum = sum + a(i,l)*b(l,j)
 1200 continue
        c(i,j) = sum
 1100 continue
 1000 continue
        end



        subroutine scattering_copy(n,m,a,lda,b)
        implicit double precision (a-h,o-z)
        double complex b(n,m),a(lda,1)
c
c       Copy the matrix B into a block of the matrix A.
c
        do 1000 j=1,m
        do 1100 i=1,n
        a(i,j) = b(i,j)
 1100 continue
 1000 continue
c
        end

      
c$$$        subroutine scattering_create(ier,disc,eps,zk,eval,nidxs,idxs,
c$$$     -    nin,evalin,nout,evalout,naux,evalaux,nneartars,ineartars,
c$$$     -    nnearsrcs,inearsrcs,krankout,iout,krankin,iin,
c$$$     -    w,lw,lkeep)
c$$$        implicit double precision (a-h,o-z)
c$$$        double complex zk
c$$$        dimension eval(13,1),evalin(13,1),evalout(13,1),evalaux(13,1)
c$$$        dimension idxs(1),ineartars(1),inearsrcs(1),w(1)
c$$$        dimension iin(1),iout(1)
c$$$c
c$$$c       This subroutine constructs a scattering matrix for a region on
c$$$c       the lowest level of the decomposition of the surface.
c$$$c
c$$$c       This is done by first constructing matrices Ain and Aout which
c$$$c       account for incoming and outgoing potentials and factoring
c$$$c       them as 
c$$$c
c$$$c         Ain   = R_in B_in
c$$$c 
c$$$c         Aout  = Bout R_out
c$$$c
c$$$c       where Bin and Bout are submatrices of Ain and Aout, 
c$$$c       respectively.
c$$$c
c$$$c       Next the restriction T of the integral operator (2) to the 
c$$$c       specified region is constructed and inverted.  The scattering
c$$$c       matrix is 
c$$$c
c$$$c         S     = Rout T^{-1} R_in.
c$$$c
c$$$c       This routine returns a structure describing the scattering
c$$$c       matrix and its associated data at the beginning of the work
c$$$c       array.
c$$$c
c$$$c                               Input Parameters:
c$$$c
c$$$c   n - the number of total discretization nodes on the boundary surface
c$$$c   zk - the complex-valued wavenumber for the problem
c$$$c   eval - a (13,n) array providing evaluation data for all 
c$$$c       discretization nodes on the surface
c$$$c
c$$$c   nin - the number of evaluation nodes on the contour which specifies
c$$$c       incoming potentials
c$$$c   evalin - the (13,n) array which specifies the evaluation data for
c$$$c       incoming potentials
c$$$c
c$$$c   nout - the number of evaluation nodes on the contout which specifies
c$$$c        outgoing potentials
c$$$c   evalout - 
c$$$c
c$$$c   iflags - a (2,n) integer array keeps track of which target and 
c$$$c       source nodes are still extant
c$$$c
c$$$c   inearidxs - an integer array storing the indices
c$$$c   inearptrs - a (2,n) array which specifies
c$$$c
c$$$c   w - a user-supplied work array
c$$$c   lw - the length of the user-supplied work array
c$$$c
c$$$c                              Output parameters:
c$$$c
c$$$c
c$$$c   ier - an error return code;
c$$$c       ier = 0      indicates successful execution
c$$$c       ier = 4      means that the work array was of insufficient length
c$$$c       ier = 1000
c$$$c
c$$$c   w - upon successful execution, the beginning of the work array
c$$$c       will contain the scattering matrix; 
c$$$c
c$$$c   w(1) - the length of the scattering matrix structure
c$$$c   w(2) - the number of points in the target skeleton                    (nout)
c$$$c   w(3) - the number of points in the source skeleton                    (nin)
c$$$c   w(4) - a pointer to the indices of the target skeleton                (iiout)
c$$$c   w(5) - a pointer to the indices of the source skeleton                (iin)
c$$$c   w(6) - a pointer to the (nout,nin) scattering matrix                  (ismatr)
c$$$c
c$$$c   lkeep - the length of the scattering matrix structure
c$$$c
c$$$c       NOTE: all pointers are relative to the beginning of the structure
c$$$c
c$$$        ier = 0
c$$$c
c$$$c       Allocate some temporary variables from the beginning of the work
c$$$c       array.
c$$$c
c$$$        ier   = 0
c$$$        lkeep = 0
c$$$c
c$$$        nin0  = nin  + nnearsrcs + naux
c$$$        nout0 = nout + nneartars + naux
c$$$c
c$$$c       Allocate memory for the routine.
c$$$c        
c$$$        iaint     = 1
c$$$        laint     = 2*nidxs*nin0
c$$$c
c$$$        irint     = iaint+laint
c$$$        lrint     = 2*nidxs*nin0
c$$$c
c$$$        iatars    = irint+lrint
c$$$        latars    = 2*nneartars*nidxs
c$$$c
c$$$        iasrcs    = iatars+latars
c$$$        lasrcs    = 2*nnearsrcs*nidxs
c$$$c
c$$$        iaout     = iasrcs+lasrcs
c$$$        laout     = 2*nidxs*nout0
c$$$c
c$$$        irout     = iaout+laout
c$$$        lrout     = 2*nidxs*nout0
c$$$c
c$$$        it        = irout+lrout
c$$$        lt        = 2*nidxs*nidxs
c$$$c
c$$$        ismatr    = it+lt
c$$$        lsmatr    = 2*nidxs*nidxs
c$$$c
c$$$        iicolsin  = ismatr+lsmatr
c$$$        licolsin  = nidxs
c$$$c
c$$$        iicolsout = iicolsin+licolsin
c$$$        licolsout = nidxs
c$$$c
c$$$        iw2       = iicolsout+licolsout
c$$$        lw2       = lw-iw2
c$$$c
c$$$        if (lw2 .le. 0) then
c$$$        ier = 4
c$$$        return
c$$$        endif
c$$$c
c$$$c       Call an auxilliary subroutine to shape arrays.
c$$$c
c$$$        call scattering_create0(ier,eps,disc,zk,nidxs,idxs,eval,naux,
c$$$     -   evalaux,nin0,nin,evalin,nnearsrcs,inearsrcs,nout0,nout,
c$$$     -   evalout,nneartars,ineartars,w(it),w(iaint),w(irint),w(iaout),
c$$$     -   w(irout),w(iicolsout),krankout,iout,w(iicolsin),krankin,iin,
c$$$     -   w(iw2),lw2,w(iasrcs),w(iatars))
c$$$        if (ier .ne. 0) return
c$$$c
c$$$c       Form the scattering matrix structure.
c$$$c
c$$$        ismatr0  = 7
c$$$        lsmatr0  = 2*krankin*krankout
c$$$c
c$$$        iiout0   = ismatr0+lsmatr0
c$$$        liout0   = krankout
c$$$c
c$$$        iiin0    = iiout0+liout0
c$$$        liin0    = krankin
c$$$c
c$$$        lkeep    = iiin0+liin0
c$$$c
c$$$        if (lkeep .gt. lw) then
c$$$        ier = 4
c$$$        return
c$$$        endif
c$$$c
c$$$c       Move the scattering matrix to the beginning of the
c$$$c       work array for output.
c$$$c
c$$$        do 1000 i=1,lsmatr0
c$$$        w(ismatr0+i-1) = w(iw2+i-1)
c$$$ 1000 continue
c$$$c
c$$$c       Copy over the lists of chosen source and target nodes.
c$$$c
c$$$        call scattering_imove(krankout,iout,w(iiout0))
c$$$        call scattering_imove(krankin,iin,w(iiin0))
c$$$c
c$$$        w(1)     = lkeep
c$$$        w(2)     = krankout
c$$$        w(3)     = krankin
c$$$        w(4)     = iiout0
c$$$        w(5)     = iiin0
c$$$        w(6)     = ismatr0
c$$$c
c$$$        end
c$$$
c$$$
c$$$        subroutine scattering_create0(ier,eps,disc,zk,nidxs,idxs,eval,
c$$$     -   naux,evalaux,nin0,nin,evalin,nnearsrcs,inearsrcs,nout0,nout,
c$$$     -   evalout,nneartars,ineartars,t,aint,rint,aout,rout,icolsout,
c$$$     -   krankout,iout,icolsin,krankin,iin,w,lw,asrcs,atars)
c$$$        implicit double precision (a-h,o-z)
c$$$        double complex zk,val,sum
c$$$c
c$$$        double complex t(nidxs,nidxs)
c$$$        double complex aint(nin0,nidxs),rint(nidxs,nidxs)
c$$$        double complex aout(nout0,nidxs),rout(nidxs,nidxs)
c$$$        double complex asrcs(nidxs,nnearsrcs),atars(nneartars,nidxs)
c$$$c
c$$$        dimension inearsrcs(1),ineartars(1)
c$$$c
c$$$        dimension eval(13,1),evalaux(13,1),evalin(13,1),evalout(13,1)
c$$$        dimension w(1),idxs(1),disc(1),icolsout(1),iout(1),iin(1)
c$$$        dimension icolsin(1)
c$$$c
c$$$c       Form and invert the restricted operator T.
c$$$c
c$$$        ikernel = 3
c$$$        call disc_eval(ier,disc,ikernel,zk,nidxs,idxs,
c$$$     -    nidxs,idxs,t,w,lw)
c$$$        if (ier .ne. 0) return
c$$$        do 1000 i=1,nidxs
c$$$        t(i,i) = t(i,i) + 0.5d0 
c$$$ 1000 continue
c$$$        call lapack_invert(ier,nidxs,t,w,lw)
c$$$        if (ier .ne. 0) return
c$$$c
c$$$c       Form the matrix Aout(nidxs,nout0) and factor it as 
c$$$c
c$$$c           Aout(nout0,nidxs) = Bout(nout0,krankout) * Rout(krankout,nidxs)
c$$$c
c$$$        ikernel = 3
c$$$        call disc_eval(ier,disc,ikernel,zk,nneartars,ineartars,
c$$$     -    nidxs,idxs,atars,w,lw)
c$$$        if (ier .ne. 0) return
c$$$c
c$$$        do 2000 j=1,nidxs
c$$$        do 2100 i=1,nout
c$$$        call ksingle(evalout(1,i),eval(1,idxs(j)),zk,val)
c$$$        aout(i,j) = val
c$$$ 2100 continue
c$$$        do 2200 i=1,naux
c$$$        call ksingleprime(evalaux(1,i),eval(1,idxs(j)),zk,val)
c$$$        aout(nout+i,j) = val
c$$$ 2200 continue
c$$$        do 2300 i=1,nneartars
c$$$        aout(nout+naux+i,j) = atars(i,j)
c$$$ 2300 continue
c$$$ 2000 continue
c$$$c
c$$$        if (lw .le. 2*(3*nidxs+1)*nout+1) then
c$$$        ier = 4
c$$$        return
c$$$        endif
c$$$c
c$$$        call factor_right(eps,nout0,nidxs,aout,krankout,icolsout,rout,w)
c$$$        do 2500 i=1,krankout
c$$$        iout(i) = idxs(icolsout(i))
c$$$ 2500 continue
c$$$c
c$$$c       Form the matrix Aint (nin0,nidxs) and factor it as 
c$$$c
c$$$c            Aint (nin0,nidxs) = Aint (nin0,krankin) Rint (krankin, nidxs)
c$$$c
c$$$        ikernel = 3
c$$$        call disc_eval(ier,disc,ikernel,zk,nidxs,idxs,nnearsrcs,
c$$$     -    inearsrcs,asrcs,w,lw)
c$$$        if (ier .ne. 0 ) return
c$$$c
c$$$        do 3000 j=1,nidxs
c$$$        do 3100 i=1,nin
c$$$        call ksingleprime(eval(1,idxs(j)),evalin(1,i),zk,val)
c$$$        aint(i,j) = val
c$$$ 3100 continue
c$$$        do 3200 i=1,naux
c$$$        call ksingleprime(eval(1,idxs(j)),evalaux(1,i),zk,val)
c$$$        aint(nin+i,j) = val
c$$$ 3200 continue
c$$$        do 3300 i=1,nnearsrcs
c$$$        aint(nin+naux+i,j) = asrcs(j,i)
c$$$ 3300 continue
c$$$ 3000 continue
c$$$c
c$$$        if (lw .le. 2*(3*nidxs+1)*nin+11) then
c$$$        ier = 4
c$$$        return
c$$$        endif
c$$$c
c$$$        call factor_right(eps,nin0,nidxs,aint,krankin,icolsin,rint,w)
c$$$c
c$$$        do 3600 i=1,krankin
c$$$        iin(i) = idxs(icolsin(i))
c$$$ 3600 continue
c$$$c
c$$$c       Form the matrix Rout(krankout,nidxs) * T^{-1} * Rint(nidxs,krankin)
c$$$c
c$$$        ib = 1
c$$$        lb = 2*krankout*krankin
c$$$c
c$$$        ia = ib+lb
c$$$        la = 2*krankout*nidxs
c$$$c
c$$$        call apply_left(nidxs,nidxs,krankout,icolsout,rout,t,w(ia))
c$$$        call apply_rightt(nidxs,krankout,krankin,icolsin,rint,w(ia),
c$$$     -   w(ib))
c$$$c
c$$$        call prinf("in create_smatr, krankout = *",krankout,1)
c$$$        call prinf("in create_smatr, iout = *",iout,krankout)
c$$$        call prinf("in create_smatr, krankin = *",krankin,1)
c$$$        call prinf("in create_smatr, iin = *",iin,krankin)
c$$$c
c$$$        end


        subroutine scattering_create(ier,disc,eps,zk,eval,nidxs,
     -    idxs,nin,evalin,nout,evalout,naux,evalaux,nneartars,ineartars,
     -    nnearsrcs,inearsrcs,krankout,iout,krankin,iin,
     -    w,lw,lkeep)
        implicit double precision (a-h,o-z)
        double complex zk
        dimension eval(13,1),evalin(13,1),evalout(13,1),evalaux(13,1)
        dimension idxs(1),ineartars(1),inearsrcs(1),w(1)
        dimension iin(1),iout(1)
c
c       This subroutine constructs a scattering matrix for a region on
c       the lowest level of the decomposition of the surface.
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
c       Next the restriction T of the integral operator (2) to the 
c       specified region is constructed and inverted.  The scattering
c       matrix is 
c
c         S     = Rout T^{-1} R_in.
c
c       This routine returns a structure describing the scattering
c       matrix and its associated data at the beginning of the work
c       array.
c
c                               Input Parameters:
c
c   n - the number of total discretization nodes on the boundary surface
c   zk - the complex-valued wavenumber for the problem
c   eval - a (13,n) array providing evaluation data for all 
c       discretization nodes on the surface
c
c   nin - the number of evaluation nodes on the contour which specifies
c       incoming potentials
c   evalin - the (13,n) array which specifies the evaluation data for
c       incoming potentials
c
c   nout - the number of evaluation nodes on the contout which specifies
c        outgoing potentials
c   evalout - 
c
c   iflags - a (2,n) integer array keeps track of which target and 
c       source nodes are still extant
c
c   inearidxs - an integer array storing the indices
c   inearptrs - a (2,n) array which specifies
c
c   w - a user-supplied work array
c   lw - the length of the user-supplied work array
c
c                              Output parameters:
c
c
c   ier - an error return code;
c       ier = 0     indicates successful execution
c       ier = 4     means that the work array was of insufficient length
c       ier = 1000
c
c   w - upon successful execution, the beginning of the work array
c       will contain the scattering matrix; 
c
c   w(1) - the length of the scattering matrix structure
c   w(2) - the number of points in the target skeleton                    (nout)
c   w(3) - the number of points in the source skeleton                    (nin)
c   w(4) - a pointer to the indices of the target skeleton                (iiout)
c   w(5) - a pointer to the indices of the source skeleton                (iin)
c   w(6) - a pointer to the (nout,nin) scattering matrix                  (ismatr)
c
c   lkeep - the length of the scattering matrix structure
c
c       NOTE: all pointers are relative to the beginning of the structure
c
        ier = 0
c
c       Allocate some temporary variables from the beginning of the work
c       array.
c
        ier   = 0
        lkeep = 0
c
        ismatr = 7
        lsmatr = 2*nidxs*nidxs
c
        iw2    = ismatr+lsmatr
        lw2    = lw-iw2
c
        ikernel = 3
        call disc_eval(ier,disc,ikernel,zk,nidxs,idxs,
     -    nidxs,idxs,w(ismatr),w(iw2),lw2)
c
        call scattering_addi(nidxs,w(ismatr))
c
        call lapack_invert(ier,nidxs,w(ismatr),w(iw2),lw2)
        if (ier .ne. 0) return
c
        krankin  = nidxs
        krankout = nidxs
c
        do 1000 i=1,nidxs
        iout(i) = idxs(i)
        iin(i)  = idxs(i)
 1000 continue
c
        iiin    = ismatr+lsmatr
        liin    = nidxs
c
        iiout   = iiin+liin
        liout   = nidxs
c
        lkeep   = iiout+liout
c
        call scattering_move(krankin,iin,w(iiin))
        call scattering_move(krankout,iout,w(iiout))
c
        w(1)    = lkeep
        w(2)    = krankout
        w(3)    = krankin
        w(4)    = iiout
        w(5)    = iiin
        w(6)    = ismatr
c
        end

        subroutine scattering_addi(n,a)
        implicit double precision (a-h,o-z)
        double complex a(n,n)
        do 1000 i=1,n
        a(i,i) = a(i,i)+0.5d0
 1000 continue
        end


        subroutine lapack_invert(ier,n,zmatr,w,lw)
        implicit double precision (a-h,o-z)
        double complex zmatr(n,n)
        dimension w(1)
        integer*4 lw4,n4,ier4
c
c       Invert an (n,n) double complex matrix via LAPACK.
c
c                          Input Parameters:  
c
c                         Output parameters:
c
        ier   = 0
        n4    = n
c
        iipiv = 1
        lipiv = n+10
c
        iw2   = iipiv+lipiv
        lw2   = lw-iw2
c
        lw4    = min(1 000 000 000,lw2/4)
c
        call lapack_invert0(ier4,n4,zmatr,w(iipiv),w(iw2),lw4)
c
        ier = ier4
        end


        subroutine lapack_invert0(ier,n,zmatr,ipiv,w,lw)
        implicit none
        double complex zmatr(n,n),w(lw)
        integer*4 ipiv(n+10)
        integer*4 n,ier,lw,lw0,lneeded,info,ipiv0(1)
c
        ier = 0
c
        call zgetrf(n,n,zmatr,n,ipiv,info)
        if (info .ne. 0) then
        ier = 1000
        return
        endif
c
        lw0 = -1
        call zgetri(n,zmatr,n,ipiv,w,lw0,info)
        if (info .ne. 0) then
        ier = 1000
        return
        endif
        lneeded = w(1)
        if (lw .le. lneeded) then
        ier = 4
        return
        endif
c
        call zgetri(n,zmatr,n,ipiv,w,lw,info)
c
        if (info .ne. 0) then
        ier = 1000
        return
        endif
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


        subroutine scattering_sphere_data(n,x0,y0,z0,r,nn,eval)
        implicit double precision (a-h,o-z)
        dimension eval(13,1)
        data pi    / 3.14159265358979323846264338327950288d0 /
c
c        dimension xslege(n),whtslege(n)
c
c       Return evaluation data for a user-specified sphere.
c       
c
c       WAVENUMBER NEEDS TO BE INCLUDED HERE???
c
c       USE SOMETHING OTHER THAN LEGENDRE TO MAKE THIS FAST????
c
c
c       Fetch the n-point Legendre quadrature rule.
c
c        call legequad(n,xslege,whtslege)
c
c       Construct the quadrature for spherical harmonics.
c
        nn = 0
c
        dd = n
        dd = 1.0d0/(dd)
c
        do 1000 i=1,n
        do 1100 j=1,n
        theta  = (j-1)*2*pi*dd
        phi    = (i-1)*pi*dd
c
        nn          = nn + 1
        eval(1,nn)  = sqrt(2*pi*r**2*dd*dd)
c        eval(1,nn)  = sqrt(2*pi*r**2*dd*whtslege(i))
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



        subroutine scattering_imove(k,ia,ib)
        implicit double precision (a-h,o-z)
        dimension ia(1),ib(1)
        do 1000 i=1,k
        ib(i) = ia(i)
 1000 continue
        end
      

        subroutine scattering_move(k,a,b)
        implicit double precision (a-h,o-z)
        dimension a(1),b(1)
        do 1000 i=1,k
        b(i) = a(i)
 1000 continue
        end


