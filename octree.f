        implicit double precision (a-h,o-z)
        double precision, allocatable :: zs(:,:),oc(:),zs0(:,:)
        integer, allocatable          :: ipoints1(:),ipoints2(:)
        dimension r(3),dr(3,2)
c
c       Sample a collection of points on a torus
c
        m = 100
        n = m*m
c
        allocate(ipoints1(n),ipoints2(n))
        allocate(zs(3,n),zs0(3,n+10 000))
c
        call make_points2(m,n,zs)
        call prinf("before octree, n =*",n,1)
c
c       Plot the points.
c     
c        iplot = 1
c        call plot3d_points("input points*",iplot,n,zs)
c
c       Construct the octree.
c
        loc = 400 000 000
        allocate(oc(loc))
c
        k = 100
c
        call elapsed(t1)
        call octree(ier,n,zs,k,oc,loc,lkeep)
        call elapsed(t2)
        if (ier .ne. 0) then
        call prinf("after octree, ier = *",ier,1)
        stop
        endif
c
        call prinf("after octree, lkeep = *",lkeep,1)
        call prin2("octree time = *",t2-t1,1)
c
c       Perform a range search with the structure.
c
        bx = 3.0d0
        by = 0.0d0
        bz = 1.0d0
        br = 1.0d0
c
        call ocrange(oc,bx,by,bz,br,npoints1,ipoints1)
        call quicksorti(npoints1,ipoints1)
        call prinf("after ocrange, npoints1 = *",npoints1,1)
        call prinf("after ocrange, ipoints1 = *",ipoints1,npoints1)
c
c       Do the same range search with brute force.
c
        npoints2=0
        br0 = br**2
c
        do 2000 i=1,n
        x = zs(1,i)
        y = zs(2,i)
        z = zs(3,i)
c
        dd = (x-bx)**2+(y-by)**2+(z-bz)**2
        if (dd .lt. br0) then
        npoints2=npoints2+1
        ipoints2(npoints2)=i
        endif
 2000 continue
c
        call prinf("after ocrange, npoints2 = *",npoints2,1)
        call prinf("after ocrange, ipoints2 = *",ipoints2,npoints2)
c
        call prina("*")
        call prina("*")
        call prina("*")
c
        call ocinfo(oc,nlevels,nboxes,n)
c
        call ocleaves(oc,nleaves,ipoints1)
        call prinf("ileaves = *",ipoints1,nleaves) 
c
        nn = 0
        do 3100 i=1,nleaves
        ibox = ipoints1(i)
        call ocpoints(oc,ibox,npoints,ipoints2)
        nn = nn + npoints
 3100 continue
        call prinf("# points in leaves = *",nn,1)
c
        call prinf("nboxes=*",nboxes,1)

        do 3000 ibox=1,nboxes

        call ocbox(oc,ibox,x1,y1,z1,x2,y2,z2,nchildren,level,iparent)
        call ocpoints(oc,ibox,npoints,ipoints1)
c
        call prinf("ibox = *",ibox,1)
        call prinf("ipoints = *",ipoints1,npoints)
c
        call oczs(oc,npoints,ipoints1,zs0)
c
        do 4000 i=1,npoints
        x = zs0(1,i)
        y = zs0(2,i)
        z = zs0(3,i)

        eps = 1.0d-12
        if (x .gt. x1-eps .AND. x .lt. x2+eps .AND.
     -      y .gt. y1-eps .AND. y .lt. y2+eps .AND.
     -      z .gt. z1-eps .AND. z .lt. z2+eps) goto 4000
c
        print *,ibox,i
        print *,x1,x,x2
        print *,y1,y,y2
        print *,z1,z,z2
        print *,""
        stop
 4000 continue
c
        call occhildren(oc,ibox,nchildren,ipoints1)
        call prina("*")
        call prinf("ibox = *",ibox,1)
        call prinf("npoints=*",npoints,1)
        call prinf("level = *",level,1)
        call prinf("iparent = *",iparent,1)
        call prinf("ichildren = *",ipoints1,nchildren)
c
        call ocneighbors(oc,ibox,nneighbors,ipoints1)
        call prinf("ineigbhors = *",ipoints1,nneighbors)
c
        call ocpatrons(oc,ibox,npatrons,ipoints1)
        call prinf("ipatrons = *",ipoints1,npatrons)
c
        call occlose(oc,ibox,nclose,ipoints1)
        call prinf("iclose = *",ipoints1,nclose)
c        
 3000 continue
c$$$
c$$$c
c$$$        ibox = 1
c$$$        call occhildren(oc,ibox,npoints1,ipoints1)
c$$$        call prinf("ichildren = *",ipoints1,npoints1)
c$$$c
c$$$c       Check some of the support routines.
c$$$c
c$$$        ibox = 3
c$$$        call ocpoints(oc,ibox,npoints1,ipoints1)
c$$$c
c$$$        call prinf("ibox = *",ibox,1)
c$$$        call prinf("after ocrange, npoints1 = *",npoints1,1)
c$$$        call prinf("after ocrange, ipoints1 = *",ipoints1,npoints1)
c$$$c
c$$$        ilevel = 3
c$$$        call oclevel(oc,ilevel,nboxes,ipoints1)
c$$$        call prinf("nboxes = *",nboxes,1)
c$$$        call prinf("iboxes1 = *",ipoints1,nboxes)
c$$$c
c$$$        call ocleaves(oc,npoints1,ipoints1)
c$$$        call prinf("ileaves = *",ipoints1,npoints1)
        end



        subroutine make_points2(m,n,zs)
        implicit double precision (a-h,o-z)
        dimension zs(3,1)
c
c       Sample points on a torus.
c
        a = 2.000d0
        b = 1.000d0

        pi    = acos(-1.0d0)
        delta = 2*pi/m
c        
        n=0
        do 1000 i=1,m
        t = delta*(i-1)
        do 1100 j=1,m
        s = delta*(j-1)
        n=n+1
        zs(1,n)    = (a+b*cos(s))*sin(t)
        zs(2,n)    = (a+b*cos(s))*cos(t)
        zs(3,n)    = b*sin(s)
 1100 continue
 1000 continue
c
        end

        subroutine make_points(n,zs)
        implicit double precision (a-h,o-z)
        dimension zs(3,1)
c
        pi    = acos(-1.0d0)
        delta = 2*pi/n
c        
        do 1000 i=1,n
c
        theta = delta * (i-1)
c
        x = cos(theta)
        y = sin(theta)
        z = 0
ccos(30*theta)
c
        zs(1,i) = x
        zs(2,i) = y
        zs(3,i) = z
c
 1000 continue
c
        end

        subroutine sphere(s,t,r,dr,iface,par2,par3,par4)
        implicit double precision (a-h,o-z)
        dimension r(3),dr(3,2)
c
c       This subroutine supplies a discretization of the unit sphere
c       over a collection of rectangles in the plane.
c
        dd  = 1.0d0/sqrt(1.0d0+s**2+t**2)
        dds = -dd**3 * s
        ddt = -dd**3 * t
c
        if (iface .eq. 1) then
c
        r(1) = s * dd
        r(2) = t * dd
        r(3) = 1 * dd
c
        dr(1,1) = dd + s*dds
        dr(2,1) = t * dds
        dr(3,1) = dds
c
        dr(1,2) = s * ddt
        dr(2,2) = dd + t * ddt
        dr(3,2) = ddt
c
        return
        endif
c
        if (iface .eq. 2) then
c        
        r(1) = s * dd
        r(2) = t * dd
        r(3) =-1 * dd
c
        dr(1,1) = dd + s*dds
        dr(2,1) = t * dds
        dr(3,1) = -dds
c
        dr(1,2) = s * ddt
        dr(2,2) = dd + t * ddt
        dr(3,2) = -ddt
        return
        endif
c
        if (iface .eq. 3) then
        r(3) = s * dd
        r(1) = t * dd
        r(2) = 1 * dd
c
        dr(3,1) = dd + s*dds
        dr(1,1) = t * dds
        dr(2,1) = dds
c
        dr(3,2) = s * ddt
        dr(1,2) = dd + t * ddt
        dr(2,2) = ddt
        return
        endif
c
        if (iface .eq. 4) then
        r(3) = s * dd
        r(1) = t * dd
        r(2) = -1 * dd
c
        dr(3,1) = dd + s*dds
        dr(1,1) = t * dds
        dr(2,1) = -dds
c
        dr(3,2) = s * ddt
        dr(1,2) = dd + t * ddt
        dr(2,2) = -ddt
        return
        endif
c
        if (iface .eq. 5) then
        r(2) = s * dd
        r(3) = t * dd
        r(1) = 1 * dd
c
        dr(2,1) = dd + s*dds
        dr(3,1) = t * dds
        dr(1,1) = dds
c
        dr(2,2) = s * ddt
        dr(3,2) = dd + t * ddt
        dr(1,2) = ddt
        return
        endif
c
        if (iface .eq. 6) then
        r(2) = s * dd
        r(3) = t * dd
        r(1) = -1 * dd
c
        dr(2,1) = dd + s*dds
        dr(3,1) = t * dds
        dr(1,1) = -dds
c
        dr(2,2) = s * ddt
        dr(3,2) = dd + t * ddt
        dr(1,2) = -ddt
        return
        endif
c
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c       This file houses an "octree" code.  That is, it  contains 
c       subroutines for constructing a data structure which describes 
c       a simple spatial decomposition of a collection of points in 
c       R^3.
c
c       The principal purpose of this code is to allow for "fast" range
c       searches on collections of points in space.  Once an octree
c       containing n points has been constructed (this takes O(n log n)
c       operations) a range search can be conducted in O(log n) 
c       operations.
c
c       It also constructs the interaction lists necessary for a
c       standard multipole code (that is, lists of neighbors,
c       lists of patrons, and lists of close boxes).
c
c       The following subroutines are user-callable:
c
c   octree - construct an octree data structure containing a user-
c       specified collection of points
c
c   ocrange - perform a range search; in this case, return the list of
c       all points in the octree which are inside of a user-specified
c       ball
c
c   ocrangebox - perform a range search; in this case, return the
c       list of all points in the octree which are inside of a
c       user-specified box
c
c   ocinfo - return information about the octree structure
c
c   oclevel - return the list of boxes on a user-specified level of
c       the octree
c
c   ocbox - return information about a specified box
c
c   ocpoints - return the list of points in a box
c
c   occhildren - return the list of children of a box
c
c   oczs - return the coordinates of a collection of points in the
c       octree
c
c   ocneighbors - return the list of neighbors of a specified box
c       
c
c       The neighbors of a box ibox are the set of all boxes
c       which are children of neighbors of its parent which
c       intersect with it.  The root box has no neighbors by
c       definition.
c
c       NOTE: the box ibox is a neighbor of itself by definition
c
c   ocpatrons - return the list of patrons of a box
c
c       Any box which is a child of a neighbor of the parent of
c       ibox but which does not intersect ibox is a patron of
c       ibox.
c
c
c   occlose - for a specified *leaf* box ibox, return the list of 
c       leaves which are "close" to ibox
c
c       A box jbox is close to ibox if either
c
c         ibox is not well-separated from jbox 
c
c       or
c
c         jbox is contained in a neighbor of ibox
c  
c       
c   ocleaves - return the indices of all leaf (childless) boxes in the 
c       tree
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



        subroutine octree(ier,n,zs,k,oc,loc,lkeep)
        implicit double precision (a-h,o-z)
        dimension zs(3,1),oc(1)
c
c       Construct an octree structure for a collection of user-
c       specified points in R^3.  Also, construct all of the lists
c       necessary for standard multipole codes.
c
c       This code operates by constructing a bounding cube for the
c       input points and recursively subdividing it until each
c       subcube contains fewer than k points.
c  
c
c                           Input Parameters:
c
c   n - the number of points
c   zs - a (3,n) array each column of which specified the coordinates 
c       of one point
c   k - the maximum number of points which a leaf node may contain
c   loc - the length of the user-supplied oc array
c
c                           Output Parameters:
c
c   ier - an error return code;
c       ier = 0   indicates successful execution
c       ier = 4   means that the oc structure was of insufficient length
c       ier = 16  means that the preset maximum number of levels was
c                 exceeded (this usually occurs because of a problem
c                 with the user-supplied input points)
c       ier = 32  means that the number of boxes exceeded the number of
c                 points
c
c   oc - upon return, this user-supplied array will contain a structure
c       describing the octree
c   lkeep - length of the data stored in the oc array
c
c                             Octree Structure:
c
c   The header of the octree structure is formatted as follows:
c
c   oc(1)  - the number of levels in the tree                             (nlevels)
c   oc(2)  - the number of boxes in the tree                              (nboxes)
c   oc(3)  - the number of points contained in the tree                   (n)
c   oc(4)  - a pointer to the ilevels array                               (iilevels)
c   oc(5)  - a pointer to the boxes array                                 (iboxes)
c   oc(6)  - a pointer to the zs array                                    (izs)
c   oc(7)  - a pointer to the list of indices                             (iidxs)
c   oc(8)  - a pointer to the ilists array                                (iilists)
c   oc(9)  - a pointer to a work array of length 3*n                      (iwork)
c   oc(10) - the number of leaves (childless nodes) in the tree           (nleaves)
c   oc(11) - a pointer into the ilists array to the list of leaf nodes    (ileaves)
c
c   The ilevels array is dimensioned (2,nlevels).  Each record 
c   describes one level of the octree:
c
c   ilevels(1,j) - the number of boxes on the jth level
c   ilevels(2,j) - the index of the first box on the jth level
c
c   NOTE: the records for the boxes on each level are stored in 
c   consecutive positions in the boxes array
c
c   The boxes array is dimensioned (30,nboxes) with the jth entry containing
c   the following information about the jth box:
c
c   boxes(1,j)  - the x-coord of the lower left point of the box          (x1)
c   boxes(2,j)  - the y-coord of the lower left point of the box          (y1)
c   boxes(3,j)  - the z-coord of the lower left point of the box          (z1)
c   boxes(4,j)  - the x-coord of the upper right point of the box         (x2)
c   boxes(5,j)  - the y-coord of the upper right point of the box         (y2)
c   boxes(6,j)  - the z-coord of the upper right point of the box         (z2)
c   boxes(7,j)  - the level of the box                                    (level) 
c   boxes(8,j)  - the index of the boxes parent                           (iparent)
c   boxes(9,j)  - the number of points in the box                         (nn)
c   boxes(10,j) - a pointer into the idxs array of the list of points     (i1)
c                 in the box
c
c   boxes(11,j) - the number of neighboring boxes                         (nneighbors)
c   boxes(12,j) - a pointer into the ilists array to the list of          (ineighbor1)
c                 neighboring boxes
c   boxes(13,j) - the number of patron boxes boxes                        (npatrons)
c   boxes(14,j) - a pointer into the ilists array to the first entry      (ipatron1)
c                 of the list of patrons
c   boxes(15,j) - the number of leaf nodes which are close                (nclose)
c   boxes(16,j) - a pointer into the ilists array to the list of          (iclose1)
c                 close nodes
c           
c   boxes(20,j)    - the number of children of the box                    (nchildren)
c   boxes(21-28,j) - indices of the child boxes
c     
c   The variable ilists contains various integer lists which are 
c   associated with the octree and with its boxes.
c
c   NOTE: the root box is always box 1 and it is the sole box on level 1
c   of the tree.
c
        ier = 0
c
        maxlevels = 500
c
c       Allocate memory from the oc array.
c
        iwork    = 100
        lwork    = 3*n
c
        izs      = iwork+lwork
        lzs      = 3*n
c
        iilevels = izs+lzs
        lilevels = 2*maxlevels
c
        iidxs    = iilevels+lilevels
        lidxs    = n
c
        iboxes   = iidxs+lidxs
        lboxes   = loc-iboxes
c
        maxboxes = lboxes/30
c
        if (maxboxes .le. 1) then
        ier = 4
        return
        endif
c
        call prinf("in octree, maxboxes = *",maxboxes,1)
c
c       Make a copy of the coordinates of the points.
c
        call ocmove(3*n,zs,oc(izs))
c
c       Form the list of boxes by recursively subdividing
c
        call octree0(ier,k,maxlevels,maxboxes,n,nboxes,nlevels,
     -   oc(iilevels),oc(iidxs),oc(iboxes),oc(izs),oc(iwork))
        if (ier .ne. 0) return
c
        call prinf("in octree, nlevels = *",nlevels,1)
        call prinf("in octree, nboxes = *",nboxes,1)
c
c       Construct the "lists" of children, neighbors and patrons.
c
        lboxes  = nboxes*30
c
        iilists = iboxes+lboxes
        lilists = loc-iilists
c
        call octree_lists(ier,nlevels,oc(iilevels),nboxes,oc(iboxes),
     -    oc(iilists),lilists,lused,oc(iwork),oc(iwork+n),oc(iwork+2*n),
     -    nleaves,ileaves)
        if (ier .ne. 0) return
c
        lkeep = iilists+lused
c
c       Construct the header of the structure.
c
        oc(1)  = nlevels
        oc(2)  = nboxes
        oc(3)  = n
        oc(4)  = iilevels
        oc(5)  = iboxes
        oc(6)  = izs
        oc(7)  = iidxs
        oc(8)  = iilists
        oc(9)  = iwork
        oc(10) = nleaves
        oc(11) = ileaves
c
        end



        subroutine octree0(ier,k,maxlevels,maxboxes,n,nboxes,nlevels,
     -    ilevels,idxs,boxes,zs,work)
        implicit double precision (a-h,o-z)
        dimension ilevels(2,maxlevels),boxes(30,1),idxs(1),zs(3,1)
        dimension work(1)
c
        ier = 0
c
        nboxes  = 0
        nlevels = 0
c
c       Initialize the idxs array.
c
        do 1000 i=1,n
        idxs(i) = i
 1000 continue
c
c       Initialize the ilevels array.
c
        do 1100 i=1,nlevels
        ilevels(1,i) = 0
        ilevels(2,i) = 0
 1100 continue
c
c       Find the extents of the root node.
c
        xmin = 1d50
        xmax =-1d50
        ymin = 1d50
        ymax =-1d50
        zmin = 1d50
        zmax =-1d50
c
        do 1200 i=1,n
        x = zs(1,i)
        y = zs(2,i)
        z = zs(3,i)
        xmin = min(xmin,x)
        ymin = min(ymin,y)
        zmin = min(zmin,z)
        xmax = max(xmax,x)
        ymax = max(ymax,y)
        zmax = max(zmax,z)
 1200 continue
c
c       Make sure it is a cube.
c
        dd =max(xmax-xmin,ymax-ymin,zmax-zmin)
c
        xx   = (xmin+xmax)/2
        xmin = xx-dd/2
        xmax = xx+dd/2
        yy   = (ymin+ymax)/2
        ymin = yy - dd/2
        ymax = yy + dd/2
        zz   = (zmin+zmax)/2
        zmin = zz - dd/2
        zmax = zz + dd/2
c
c       Initialize the root box
c
        nboxes=1
        boxes(1,nboxes)  = xmin
        boxes(2,nboxes)  = ymin
        boxes(3,nboxes)  = zmin
        boxes(4,nboxes)  = xmax
        boxes(5,nboxes)  = ymax
        boxes(6,nboxes)  = zmax
        boxes(7,nboxes)  = 1
        boxes(8,nboxes)  = 0
        boxes(9,nboxes)  = n
        boxes(10,nboxes) = 1
        boxes(20,nboxes) = 0
c
        ilevels(1,1) = 1
        ilevels(2,1) = 1
c
c       Construct the boxes level-by-level.
c
        do 2000 ilevel=1,maxlevels-1
c
c       Split each box on level ilevel.
c
        nn    = ilevels(1,ilevel)
        ibox1 = ilevels(2,ilevel)
c
        ilevels(1,ilevel+1) = 0
        ilevels(2,ilevel+1) = nboxes+1
c
        if (nn .eq. 0) goto 2200
c
        do 2100 ii=1,nn
        ibox    = ibox1+ii-1
        xmin    = boxes(1,ibox)
        ymin    = boxes(2,ibox)
        zmin    = boxes(3,ibox)
        xmax    = boxes(4,ibox)
        ymax    = boxes(5,ibox)
        zmax    = boxes(6,ibox)
        level   = boxes(7,ibox)
        iparent = boxes(8,ibox)
        nn      = boxes(9,ibox)
        i1      = boxes(10,ibox)
c
c       Add a record of the box to its parent's record.
c
        if (iparent .ne. 0) then
        nchildren = boxes(20,iparent)
        nchildren = nchildren+1
        boxes(20+nchildren,iparent) = ibox
        boxes(20,iparent)           = nchildren
        endif
c
c       Determine if the box should be split or not.
c
        if (nn .lt. k) goto 2100
c
        if (nboxes .gt. n) then
        ier = 32
        return
        endif
c
        if (level .eq. maxlevels-1) then
        ier = 16
        return
        endif
c
        if (nboxes+8 .gt. maxboxes) then
        ier  = 4
        return
        endif
c
        xx = (xmin+xmax)/2
        yy = (ymin+ymax)/2
        zz = (zmin+zmax)/2
c
        call octree_split(xmin,xmax,ymin,ymax,zmin,zmax,nn,i1,
     -    idxs,work(1),work(n+1),zs,nn1,nn2,nn3,nn4,nn5,nn6,nn7,nn8,
     -    xx,yy,zz)
c
        if (nn1 .gt. 0) then
        nboxes=nboxes+1
        ilevels(1,ilevel+1) = ilevels(1,ilevel+1)+1
        boxes(1,nboxes)  = xmin
        boxes(2,nboxes)  = ymin
        boxes(3,nboxes)  = zmin
        boxes(4,nboxes)  = xx
        boxes(5,nboxes)  = yy
        boxes(6,nboxes)  = zz
        boxes(7,nboxes)  = level+1
        boxes(8,nboxes)  = ibox
        boxes(9,nboxes)  = nn1
        boxes(10,nboxes) = i1
        boxes(20,nboxes) = 0
        endif
c
        if (nn2 .gt. 0) then
        nboxes=nboxes+1
        ilevels(1,ilevel+1) = ilevels(1,ilevel+1)+1
        boxes(1,nboxes)  = xmin
        boxes(2,nboxes)  = ymin
        boxes(3,nboxes)  = zz
        boxes(4,nboxes)  = xx
        boxes(5,nboxes)  = yy
        boxes(6,nboxes)  = zmax
        boxes(7,nboxes)  = level+1
        boxes(8,nboxes)  = ibox
        boxes(9,nboxes)  = nn2
        boxes(10,nboxes) = i1+nn1
        boxes(20,nboxes) = 0
        endif
c
        if (nn3 .gt. 0) then
        nboxes=nboxes+1
        ilevels(1,ilevel+1) = ilevels(1,ilevel+1)+1
        boxes(1,nboxes)  = xmin
        boxes(2,nboxes)  = yy
        boxes(3,nboxes)  = zmin
        boxes(4,nboxes)  = xx
        boxes(5,nboxes)  = ymax
        boxes(6,nboxes)  = zz
        boxes(7,nboxes)  = level+1
        boxes(8,nboxes)  = ibox
        boxes(9,nboxes)  = nn3
        boxes(10,nboxes) = i1+nn1+nn2
        boxes(20,nboxes) = 0
        endif
c
        if (nn4 .gt. 0) then
        nboxes=nboxes+1
        ilevels(1,ilevel+1) = ilevels(1,ilevel+1)+1
        boxes(1,nboxes)  = xmin
        boxes(2,nboxes)  = yy
        boxes(3,nboxes)  = zz
        boxes(4,nboxes)  = xx
        boxes(5,nboxes)  = ymax
        boxes(6,nboxes)  = zmax
        boxes(7,nboxes)  = level+1
        boxes(8,nboxes)  = ibox
        boxes(9,nboxes)  = nn4
        boxes(10,nboxes) = i1+nn1+nn2+nn3
        boxes(20,nboxes) = 0
        endif
c
        if (nn5 .gt. 0) then
        nboxes=nboxes+1
        ilevels(1,ilevel+1) = ilevels(1,ilevel+1)+1
        boxes(1,nboxes)  = xx
        boxes(2,nboxes)  = ymin
        boxes(3,nboxes)  = zmin
        boxes(4,nboxes)  = xmax
        boxes(5,nboxes)  = yy
        boxes(6,nboxes)  = zz
        boxes(7,nboxes)  = level+1
        boxes(8,nboxes)  = ibox
        boxes(9,nboxes)  = nn5
        boxes(10,nboxes) = i1+nn1+nn2+nn3+nn4
        boxes(20,nboxes) = 0
        endif
c
        if (nn6 .gt. 0) then
        nboxes=nboxes+1
        ilevels(1,ilevel+1) = ilevels(1,ilevel+1)+1
        boxes(1,nboxes)  = xx
        boxes(2,nboxes)  = ymin
        boxes(3,nboxes)  = zz
        boxes(4,nboxes)  = xmax
        boxes(5,nboxes)  = yy
        boxes(6,nboxes)  = zmax
        boxes(7,nboxes)  = level+1
        boxes(8,nboxes)  = ibox
        boxes(9,nboxes)  = nn6
        boxes(10,nboxes) = i1+nn1+nn2+nn3+nn4+nn5
        boxes(20,nboxes) = 0
        endif
c
        if (nn7 .gt. 0) then
        nboxes=nboxes+1
        ilevels(1,ilevel+1) = ilevels(1,ilevel+1)+1
        boxes(1,nboxes)  = xx
        boxes(2,nboxes)  = yy
        boxes(3,nboxes)  = zmin
        boxes(4,nboxes)  = xmax
        boxes(5,nboxes)  = ymax
        boxes(6,nboxes)  = zz
        boxes(7,nboxes)  = level+1
        boxes(8,nboxes)  = ibox
        boxes(9,nboxes)  = nn7
        boxes(10,nboxes) = i1+nn1+nn2+nn3+nn4+nn5+nn6
        boxes(20,nboxes) = 0
        endif
c
        if (nn8 .gt. 0) then
        nboxes=nboxes+1
        ilevels(1,ilevel+1) = ilevels(1,ilevel+1)+1
        boxes(1,nboxes)  = xx
        boxes(2,nboxes)  = yy
        boxes(3,nboxes)  = zz
        boxes(4,nboxes)  = xmax
        boxes(5,nboxes)  = ymax
        boxes(6,nboxes)  = zmax
        boxes(7,nboxes)  = level+1
        boxes(8,nboxes)  = ibox
        boxes(9,nboxes)  = nn8
        boxes(10,nboxes) = i1+nn1+nn2+nn3+nn4+nn5+nn6+nn7
        boxes(20,nboxes) = 0
        endif        
c
 2100 continue
 2000 continue
 2200 continue
c
        nlevels = ilevel-1
c
        end


        subroutine octree_split(xmin,xmax,ymin,ymax,zmin,zmax,nn,i1,
     -    idxs,iwork1,iwork2,zs,nn1,nn2,nn3,nn4,nn5,nn6,nn7,nn8,
     -    xx,yy,zz)
        implicit double precision (a-h,o-z)
        dimension idxs(1),iwork1(1),iwork2(1),zs(3,1),nns(8),mms(8)
c
        do 0100 i=1,8
        nns(i) = 0
 0100 continue
c
c       Decide which box each point in the parent belongs too.
c    
        do 1000 ii=1,nn
        i = idxs(i1+ii-1)
        x = zs(1,i)
        y = zs(2,i)
        z = zs(3,i)
c
        ichild = 1
c
        if (x .gt. xx) ichild=ichild+4
        if (y .gt. yy) ichild=ichild+2
        if (z .gt. zz) ichild=ichild+1
c
        nns(ichild) = nns(ichild)+1
        iwork1(ii)  = ichild
c
 1000 continue
c
        mms(1) = 0
        do 1100 i=2,8
        mms(i) = mms(i-1)+nns(i-1)
 1100 continue
c
        do 1200 ii=1,nn
        i = idxs(i1+ii-1)
        ichild = iwork1(ii)
        j = mms(ichild)
        j = j+1
        iwork2(j) = i
        mms(ichild)=j
 1200 continue
c
        do 1300 ii=1,nn
        idxs(i1+ii-1) = iwork2(ii)
 1300 continue
c
        nn1 = nns(1)
        nn2 = nns(2)
        nn3 = nns(3)
        nn4 = nns(4)
        nn5 = nns(5)
        nn6 = nns(6)
        nn7 = nns(7)
        nn8 = nns(8)
c
        end


        subroutine oclevel(oc,ilevel,nboxes,iboxes)
        implicit double precision (a-h,o-z)
        dimension oc(1),iboxes(1)
c
        nlevels  = oc(1)
        nboxes   = oc(2)
        n        = oc(3)
        iilevels = oc(4)
        iboxes   = oc(5)
        izs      = oc(6)
        iidxs    = oc(7)
        iilists  = oc(8)
        iwork    = oc(9)
c
        call oclevel0(ilevel,oc(iilevels),nboxes,iboxes)
c        
        end


        subroutine oclevel0(ilevel,ilevels,nboxes,iboxes)
        implicit double precision (a-h,o-z)
        dimension ilevels(2,1),iboxes(1)
        nboxes = ilevels(1,ilevel)
        i1     = ilevels(2,ilevel)
c
        do 1000 i=1,nboxes
        iboxes(i) = i1+i-1
 1000 continue
        end


        subroutine ocpoints(oc,ibox,npoints,ipoints)
        implicit double precision (a-h,o-z)
        dimension oc(1),ipoints(1)
c
        nlevels  = oc(1)
        nboxes   = oc(2)
        n        = oc(3)
        iilevels = oc(4)
        iboxes   = oc(5)
        izs      = oc(6)
        iidxs    = oc(7)
        iilists  = oc(8)
        iwork    = oc(9)
c
        call ocpoints0(ibox,oc(iboxes),oc(iidxs),npoints,ipoints)
c
        end


        subroutine ocpoints0(ibox,boxes,idxs,npoints,ipoints)
        implicit double precision (a-h,o-z)
        dimension boxes(30,1),idxs(1),ipoints(1)
c
        npoints = 0
c
        nn = boxes(9,ibox)
        i1 = boxes(10,ibox)
c
        do 1000 ii=i1,i1+nn-1
        npoints = npoints+1
        ipoints(npoints) = idxs(ii)
 1000 continue
c
        end


        subroutine oczs(oc,npoints,ipoints,zs)
        implicit double precision (a-h,o-z)
        dimension oc(1),ipoints(1),zs(3,1)
c
c       Fetch data from the octree structure.
c        
        nlevels  = oc(1)
        nboxes   = oc(2)
        n        = oc(3)
        iilevels = oc(4)
        iboxes   = oc(5)
        izs      = oc(6)
        iidxs    = oc(7)
        iilists  = oc(8)
        iwork    = oc(9)
c
        call oczs0(npoints,ipoints,oc(izs),zs)
c
        end


        subroutine oczs0(npoints,ipoints,zs0,zs)
        implicit double precision (a-h,o-z)
        dimension ipoints(1),idxs(1),zs0(3,1),zs(3,1)
        do 1000 ii=1,npoints
        i = ipoints(ii)
        zs(1,ii) = zs0(1,i)
        zs(2,ii) = zs0(2,i)
        zs(3,ii) = zs0(3,i)
 1000 continue
        end


        subroutine ocinfo(oc,nlevels,nboxes,n)
        implicit double precision (a-h,o-z)
        dimension oc(1)
        nlevels = oc(1)
        nboxes  = oc(2)
        n       = oc(3)
        end


        subroutine occhildren(oc,ibox,nchildren,ichildren)
        implicit double precision (a-h,o-z)
        dimension oc(1),ichildren(1)
c
c       Fetch data from the octree structure.
c        
        nlevels  = oc(1)
        nboxes   = oc(2)
        n        = oc(3)
        iilevels = oc(4)
        iboxes   = oc(5)
        izs      = oc(6)
        iidxs    = oc(7)
        iilists  = oc(8)
        iwork    = oc(9)
c
c       Call an auxilliary routine to shape arrays.
c
        call occhildren0(ibox,oc(iboxes),nchildren,ichildren)
c
c
        end


        subroutine occhildren0(ibox,boxes,nchildren,ichildren)
        implicit double precision (a-h,o-z)
        dimension boxes(30,1),ichildren(1)
c
        nchildren = boxes(20,ibox)
        do 1000 i=1,nchildren
        ichildren(i) = boxes(20+i,ibox)
 1000 continue
        end


        subroutine ocleaves(oc,nleaves,ileaves)
        implicit double precision (a-h,o-z)
        dimension oc(1),ileaves(1)
c
c       Fetch data from the octree structure.
c        
        nlevels  = oc(1)
        nboxes   = oc(2)
        n        = oc(3)
        iilevels = oc(4)
        iboxes   = oc(5)
        izs      = oc(6)
        iidxs    = oc(7)
        iilists  = oc(8)
        iwork    = oc(9)
        nleaves  = oc(10)
        i1       = oc(11)
c
c       Call an auxilliary routine to shape arrays.
c
        call ocleaves0(nleaves,ileaves,oc(iilists),i1)
c       
        end
 
        subroutine ocleaves0(nleaves,ileaves,ilists,i1)
        implicit double precision (a-h,o-z)
        dimension ileaves(1),ilists(1)
        do 1000 i=1,nleaves
        ileaves(i) = ilists(i1+i-1)
 1000 continue
        end


        subroutine ocneighbors(oc,ibox,nneighbors,ineighbors)
        implicit double precision (a-h,o-z)
        dimension oc(1),ineighbors(1)
c
c
c       Fetch data from the octree structure.
c        
        nlevels  = oc(1)
        nboxes   = oc(2)
        n        = oc(3)
        iilevels = oc(4)
        iboxes   = oc(5)
        izs      = oc(6)
        iidxs    = oc(7)
        iilists  = oc(8)
        iwork    = oc(9)
c
c       Call an auxilliary routine to shape arrays.
c
       call ocneighbors0(ibox,oc(iboxes),oc(iilists),nneighbors,
     -    ineighbors)
c
        end



        subroutine ocneighbors0(ibox,boxes,ilists,nneighbors,ineighbors)
        implicit double precision (a-h,o-z)
        dimension boxes(30,1),ilists(1),ineighbors(1)
c
        nneighbors = boxes(11,ibox)
        i1         = boxes(12,ibox)
c
        do 1000 i=1,nneighbors
        ineighbors(i) = ilists(i1+i-1)
 1000 continue
        end


        subroutine ocpatrons(oc,ibox,npatrons,ipatrons)
        implicit double precision (a-h,o-z)
        dimension oc(1),ipatrons(1)
c
c
c       Fetch data from the octree structure.
c        
        nlevels  = oc(1)
        nboxes   = oc(2)
        n        = oc(3)
        iilevels = oc(4)
        iboxes   = oc(5)
        izs      = oc(6)
        iidxs    = oc(7)
        iilists  = oc(8)
        iwork    = oc(9)
c
c       Call an auxilliary routine to shape arrays.
c
        call ocpatrons0(ibox,oc(iboxes),oc(iilists),npatrons,
     -    ipatrons)
c
        end


        subroutine ocpatrons0(ibox,boxes,ilists,npatrons,ipatrons)
        implicit double precision (a-h,o-z)
        dimension boxes(30,1),ilists(1),ipatrons(1)
c
        npatrons = boxes(13,ibox)
        i1       = boxes(14,ibox)
        do 1000 i=1,npatrons
        ipatrons(i) = ilists(i1+i-1)
 1000 continue
        end


        subroutine occlose(oc,ibox,nclose,iclose)
        implicit double precision (a-h,o-z)
        dimension oc(1),iclose(1)
c
c
c       Fetch data from the octree structure.
c        
        nlevels  = oc(1)
        nboxes   = oc(2)
        n        = oc(3)
        iilevels = oc(4)
        iboxes   = oc(5)
        izs      = oc(6)
        iidxs    = oc(7)
        iilists  = oc(8)
        iwork    = oc(9)
c
c       Call an auxilliary routine to shape arrays.
c
        call occlose0(ibox,oc(iboxes),oc(iilists),nclose,
     -    iclose)
c
        end


        subroutine occlose0(ibox,boxes,ilists,nclose,iclose)
        implicit double precision (a-h,o-z)
        dimension boxes(30,1),ilists(1),iclose(1)
c
        nclose   = boxes(15,ibox)
        i1       = boxes(16,ibox)
        do 1000 i=1,nclose
        iclose(i) = ilists(i1+i-1)
 1000 continue
        end


        subroutine ocrange(oc,bx,by,bz,br,npoints,ipoints)
        implicit double precision (a-h,o-z)
        dimension oc(1),ipoints(1)
c
c       Perform a "fast" range search using the octree structure; 
c       in particular, use the octree structure to find all points 
c       inside of a user-specified ball in O(log(n)) time.
c 
c                           Input Parameters:
c
c  oc - the array containing the octree structure
c  (bx,by,bz) - the center of the ball which is the subject of the
c       query
c  br - the radius of the ball which is the subject of the query
c
c                          Output Parameters:
c
c   npoints - the number of points in the ball
c   ipoints - the first npoints entries of this user-specified
c       array will contain the indices of the points inside of the 
c       ball
c
c       Fetch data from the octree structure.
c        
        nlevels  = oc(1)
        nboxes   = oc(2)
        n        = oc(3)
        iilevels = oc(4)
        iboxes   = oc(5)
        izs      = oc(6)
        iidxs    = oc(7)
        iilists  = oc(8)
        iwork    = oc(9)
c
c       Call an auxilliary routine to do the dirty work.
c
        call ocrange0(oc(izs),oc(iboxes),oc(iilists),oc(iidxs),
     -   bx,by,bz,br,npoints,ipoints,oc(iwork))
c
        end


       subroutine ocrange0(zs,boxes,ilists,idxs,bx,by,bz,br,
     -    npoints,ipoints,istack)
        implicit double precision (a-h,o-z)
        dimension boxes(30,1),zs(3,1)
        dimension ipoints(1),idxs(1),istack(1)
c
        npoints = 0
c
c
c       Found a bounding box for the ball.
c     
        br0 = br**2
c
        x3 = bx-br
        x4 = bx+br
        y3 = by-br
        y4 = by+br
        z3 = bz-br
        z4 = bz+br
c
c       Start with the top box.
c
        nstack = 1
        istack(1) = 1
c
 1000 continue
        if (nstack .eq. 0) goto 1100
        ibox = istack(nstack)
        nstack = nstack-1
        x1 = boxes(1,ibox)
        y1 = boxes(2,ibox)
        z1 = boxes(3,ibox)
        x2 = boxes(4,ibox)
        y2 = boxes(5,ibox)
        z2 = boxes(6,ibox)
c
c       Check for an intersection of the two boxes.
c
        call octree_ifinters(x1,x2,y1,y2,z1,z2,x3,x4,y3,y4,z3,z4,
     -    ifinter)
c
        if (ifinter .eq. 0) goto 1000
c
        nchildren = boxes(20,ibox)
c
c       Process childless boxes.
c
        if (nchildren .eq. 0) then
        nn = boxes(9,ibox)
        i1 = boxes(10,ibox)
        do 1300 ii=1,nn
        i = idxs(i1+ii-1)
        x = zs(1,i)
        y = zs(2,i)
        z = zs(3,i)
c
        dd = (x-bx)**2+(y-by)**2+(z-bz)**2
        if (dd .lt. br0) then
        npoints=npoints+1
        ipoints(npoints) = i
        endif
 1300 continue
c
        goto 1000
        endif
c
        do 1200 i=1,nchildren
        ichild = boxes(20+i,ibox)
        nstack=nstack+1
        istack(nstack)=ichild
 1200 continue
c
        goto 1000
 1100 continue
c        
        end


        subroutine ocrangebox(oc,x1,y1,z1,x2,y2,z2,npoints,ipoints)
        implicit double precision (a-h,o-z)
        dimension oc(1),ipoints(1)
c
c       Perform a "fast" range search using the octree structure; that
c       is, use the octree structure to find all points inside of 
c       a user-specified box.
c 
c                           Input Parameters:
c
c  oc - the array containing the octree structure
c  (x1,y1,z1) - the coordinates of the lower left corner of the
c       box
c  (x2,y2,z2) - the coordinates of the upper  right corner of the
c       box
c  br - the radius of the ball which is the subject of the query
c
c                          Output Parameters:
c
c   npoints - the number of points in the box
c   ipoints - the first npoints entries of this user-specified
c       array will contain the indices of the points inside of the 
c       ball
c
c       Fetch data from the octree structure.
c        
        nlevels  = oc(1)
        nboxes   = oc(2)
        n        = oc(3)
        iilevels = oc(4)
        iboxes   = oc(5)
        izs      = oc(6)
        iidxs    = oc(7)
        iilists  = oc(8)
        iwork    = oc(9)
c
c       Call an auxilliary routine to do the dirty work.
c
        call ocrangebox0(oc(izs),oc(iboxes),oc(iilists),oc(iidxs),
     -   x1,y1,z1,x2,y2,z2,npoints,ipoints,oc(iwork))
c
        end


       subroutine ocrangebox0(zs,boxes,ilists,idxs,x1,y1,z1,x2,y2,z2,
     -    npoints,ipoints,istack)
        implicit double precision (a-h,o-z)
        dimension boxes(30,1),zs(3,1)
        dimension ipoints(1),idxs(1),istack(1)
c
        eps0 = 1.0d-12
        npoints = 0
c
c       Start with the top box.
c
        nstack = 1
        istack(1) = 1
c
 1000 continue
        if (nstack .eq. 0) goto 1100
        ibox = istack(nstack)
        nstack = nstack-1
        x3 = boxes(1,ibox)
        y3 = boxes(2,ibox)
        z3 = boxes(3,ibox)
        x4 = boxes(4,ibox)
        y4 = boxes(5,ibox)
        z4 = boxes(6,ibox)
c
c       Check for an intersection of the two boxes.
c
        call octree_ifinters(x1,x2,y1,y2,z1,z2,x3,x4,y3,y4,z3,z4,
     -    ifinter)
c
        if (ifinter .eq. 0) goto 1000
c
        nchildren = boxes(20,ibox)
c
c       Process childless boxes.
c
        if (nchildren .eq. 0) then
        nn = boxes(9,ibox)
        i1 = boxes(10,ibox)
        do 1300 ii=1,nn
        i = idxs(i1+ii-1)
        x = zs(1,i)
        y = zs(2,i)
        z = zs(3,i)
c
        if ( (x .gt. x1-eps0 .AND. x .lt. x2+eps0) .AND.
     -       (y .gt. y1-eps0 .AND. y .lt. y2+eps0) .AND.
     -       (z .gt. z1-eps0 .AND. z .lt. z2+eps0)) then
c
        npoints=npoints+1
        ipoints(npoints) = i
        endif
c
 1300 continue
c
        goto 1000
        endif
c
        do 1200 i=1,nchildren
        ichild = boxes(20+i,ibox)
        nstack=nstack+1
        istack(nstack)=ichild
 1200 continue
c
        goto 1000
 1100 continue
c        
        end



        subroutine octree_lists(ier,nlevels,ilevels,nboxes,boxes,ilists,
     -   lilists,lused,iwork1,iwork2,ileaf,nleaves,ileaves)
        implicit double precision (a-h,o-z)
        dimension boxes(30,nboxes),ilists(1),ilevels(2,1)
        dimension iwork1(1),iwork2(1),ileaf(1)
c
c       Zero the counts of patrons, neighbors and close boxes.
c
        nleaves = 0
        do 1000 ibox=1,nboxes
        boxes(11,ibox) = 0
        boxes(12,ibox) = 0
        boxes(13,ibox) = 0
        boxes(14,ibox) = 0
        boxes(15,ibox) = 0
        boxes(16,ibox) = 0
 1000 continue
c
        iptr = 1
c
c       Make lists of neighbors and patrons; this must proceed in
c       a top down fashion.
c
        do 1200 ibox=1,nboxes
        x1         = boxes(1,ibox)
        y1         = boxes(2,ibox)
        z1         = boxes(3,ibox)
        x2         = boxes(4,ibox)
        y2         = boxes(5,ibox)
        z2         = boxes(6,ibox)
        iparent    = boxes(8,ibox)
        nchildren  = boxes(20,ibox)
c
        if (nchildren .eq. 0) then
        nleaves = nleaves+1
        ileaf(nleaves) = ibox
        endif
c
c       The root box is its own neighbor; it has no patrons.
c
        if (ibox .eq. 1) then
        boxes(11,1)  = 1
        boxes(12,1)  = iptr
        ilists(iptr) = 1
        iptr=iptr+1
        goto 1200
        endif
c
        nneighbors = 0
        npatrons   = 0
c
c       Process the children of the neighbors of the parent box.
c
        nn         = boxes(11,iparent)
        i1         = boxes(12,iparent)
c
        do 1400 ii = 1,nn
        in         = ilists(i1+ii-1)
        nchildren  = boxes(20,in)
        do 1500 jj = 1,nchildren
        ichild     = boxes(20+jj,in)
c
        x3         = boxes(1,ichild)
        y3         = boxes(2,ichild)
        z3         = boxes(3,ichild)
        x4         = boxes(4,ichild)
        y4         = boxes(5,ichild)
        z4         = boxes(6,ichild)
c
        call octree_ifinters(x1,x2,y1,y2,z1,z2,x3,x4,y3,y4,z3,z4,
     -    ifinter)
c
        if (ifinter .eq. 1) then
        nneighbors = nneighbors+1
        iwork1(nneighbors) = ichild
        else
        npatrons = npatrons+1
        iwork2(npatrons) = ichild
        endif
 1500 continue           
 1400 continue
c
c       Copy the lists of neigbhors and patrons into ilists
c
        if ( iptr + nneighbors + npatrons .gt. lilists) then
        ier = 4
        return
        endif
c
        boxes(11,ibox) = nneighbors
        boxes(12,ibox) = iptr
c
        do 1600 ii=1,nneighbors
        ilists(iptr) = iwork1(ii)
        iptr=iptr+1
 1600 continue
c
        boxes(13,ibox) = npatrons
        boxes(14,ibox) = iptr
        do 1700 ii=1,npatrons
        ilists(iptr) = iwork2(ii)
        iptr=iptr+1
 1700 continue
c
 1200 continue
c
c
c       Copy the list of leaves into ilists array.
c
        if (iptr+nleaves .gt. lilists) then
        ier = 4
        return
        endif

        ileaves = iptr
        do 2000 i=1,nleaves
        ilists(iptr+i-1) = ileaf(i)
 2000 continue
        iptr=iptr+nleaves
c
c       Form a list of "close" leaves for each leaf box.
c
        do 3000 ii=1,nleaves
        ibox = ileaf(ii)
c
        nneighbors = boxes(11,ibox)
        i1         = boxes(12,ibox)
c
        x3         = boxes(1,ibox)
        y3         = boxes(2,ibox)
        z3         = boxes(3,ibox)
        x4         = boxes(4,ibox)
        y4         = boxes(5,ibox)
        z4         = boxes(6,ibox)
c
        nclose = 0
c
c       First, find all leaf boxes jbox such that ibox is not
c       well-separated from jbox.
c
        nstack = 1
        iwork1(1) = nstack
c
 3100 continue
        if (nstack .eq. 0) goto 3200
        jbox = iwork1(nstack)
        nstack = nstack-1
c
        x1         = boxes(1,jbox)
        y1         = boxes(2,jbox)
        z1         = boxes(3,jbox)
        x2         = boxes(4,jbox)
        y2         = boxes(5,jbox)
        z2         = boxes(6,jbox)
        nchildren  = boxes(20,jbox)
c
        call octree_ifseparated(x1,x2,y1,y2,z1,z2,
     -    x3,x4,y3,y4,z3,z4,ifsep)
c
        if (ifsep .eq. 1) goto 3100
c
        if (nchildren .eq. 0) then
        nclose = nclose+1
        iwork2(nclose) = jbox
        goto 3100
        endif
c
        do 3300 i=1,nchildren
        ichild = boxes(20+i,jbox)
        nstack=nstack+1
        iwork1(nstack)=ichild
 3300 continue
        goto 3100
c
 3200 continue
c
c       Now add to the list all leaf boxes which are contained
c       in neighbors of ibox.
c
        nstack     = 0
        nneigbhors = boxes(11,ibox)
        i1         = boxes(12,ibox)
c
        do 3400 i=1,nneighbors
c
        jbox = ilists(i1+i-1)
c
        nstack=nstack+1
        iwork1(nstack)=jbox
 3400 continue
c
 3500 continue
        if (nstack .eq. 0) goto 3600
        jbox   = iwork1(nstack)
        nstack = nstack-1
c
        nchildren0 = boxes(20,jbox)
c
        if (nchildren0 .eq. 0) then
        nclose         = nclose+1
        iwork2(nclose) = jbox
        else
c
        do 3700 i=1,nchildren0
        ichild0 = boxes(20+i,jbox)
        nstack=nstack+1
        iwork1(nstack) = ichild0
 3700 continue
c
        endif
        goto 3500
 3600 continue
c
        call quicksorti(nclose,iwork2)
        call iduplicates(nclose,iwork2)
c
c       Copy the list of close leaf boxes into ilists.
c
        boxes(15,ibox) = nclose
        boxes(16,ibox) = iptr
c
        if (iptr + nclose .gt. lilists) then
        ier = 4
        return
        endif
c
        do 3900 i=1,nclose
        ilists(iptr) = iwork2(i)
        iptr=iptr+1
 3900 continue
c
 3000 continue
c
        lused = iptr
c
        end




        subroutine octree_ifseparated(x1,x2,y1,y2,z1,z2,
     -    x3,x4,y3,y4,z3,z4,ifsep)
        implicit double precision (a-h,o-z)
c
c       Determine if the box B with corners (x3,y3,z3) and (x4,y4,z4) is
c       separated from the box C with corners (x1,y1,z1) and
c       (x2,y2,z2) 
c
        eps0 = 1.d-12
c
        ifsep = 0
c
        delta = x2-x1
        delta = delta-eps0
c
        if (x3 .gt. x2+delta) ifsep=1
        if (x4 .lt. x1-delta) ifsep=1
        if (y3 .gt. y2+delta) ifsep=1
        if (y4 .lt. y1-delta) ifsep=1
        if (z3 .gt. z2+delta) ifsep=1
        if (z4 .lt. z1-delta) ifsep=1
c
        end



        subroutine octree_ifinters(x1,x2,y1,y2,z1,z2,x3,x4,y3,y4,z3,z4,
     -    ifinter)
        implicit double precision (a-h,o-z)
c
        eps0 = 1.0d-12
c
c       Determine if two boxes intersect.
c
        eps = eps0
        ifinter = 1
c
        if (x1 .gt. x4+eps) ifinter=0
        if (x3 .gt. x2+eps) ifinter=0
        if (y1 .gt. y4+eps) ifinter=0
        if (y3 .gt. y2+eps) ifinter=0
        if (z1 .gt. z4+eps) ifinter=0
        if (z3 .gt. z2+eps) ifinter=0
c
        end





        subroutine ocbox(oc,ibox,x1,y1,z1,x2,y2,z2,nchildren,level,
     -    iparent)
        implicit double precision (a-h,o-z)        
        dimension oc(1)
c
c       Fetch data from the octree structure.
c        
        nlevels  = oc(1)
        nboxes   = oc(2)
        n        = oc(3)
        iilevels = oc(4)
        iboxes   = oc(5)
        izs      = oc(6)
        iidxs    = oc(7)
        iilists  = oc(8)
        iwork    = oc(9)
c
c       Call an auxilliary routine to shape arrays.
c
        call ocbox0(ibox,oc(iboxes),x1,y1,z1,x2,y2,z2,nchildren,level,
     -    iparent)
c
        end


        subroutine ocbox0(ibox,boxes,x1,y1,z1,x2,y2,z2,nchildren,level,
     -   iparent)
        implicit double precision (a-h,o-z)
        dimension boxes(30,1)
        x1        = boxes(1,ibox)
        y1        = boxes(2,ibox)
        z1        = boxes(3,ibox)
        x2        = boxes(4,ibox)
        y2        = boxes(5,ibox)
        z2        = boxes(6,ibox)        
        level     = boxes(7,ibox)
        iparent   = boxes(8,ibox)
        nchildren = boxes(20,ibox)
c
        end


        subroutine ocmove(n,a,b)
        implicit double precision (a-h,o-z)
        dimension a(1),b(1)
        do 1000 i=1,n
        b(i) = A(i)
 1000 continue
        end


