


        subroutine triangulate_split(ier,maxtris,ntris,tris)
        implicit double precision (a-h,o-z)
        dimension tris(6,1),children(6,4)
c
        ier = 0
        pi  = acos(-1.0d0)
c
        ntris0 = ntris
        do 1000 itri=1,ntris0
        if (ntris+4 .gt. maxtris) then
        ier = 4
        return
        endif
c
        call ttsplit_tri(tris(1,itri),children)
c
        do 1100 i=1,6
        tris(i,itri) = children(i,1)
 1100 continue
c
        do 1200 l=2,4
        ntris=ntris+1
        do 1300 i=1,6
        tris(i,ntris) = children(i,l)
 1300 continue
 1200 continue
 1000 continue
c
        end



        subroutine triangulate3(ier,iparam,x1,y1,x2,y2,eps,
     -    maxtris,ntris,tris,diammax)
        implicit double precision (a-h,o-z)
        dimension tris(6,1),simp(3,2),r(3),dr(3,2),dn(3),dr2(2,2)
        dimension children(6,4),amap(2,3)
c
        double precision, allocatable :: istack(:,:)
        double precision, allocatable :: xs(:),ys(:),whts(:)
        double precision, allocatable :: zs(:,:)
c
        ier = 0
        pi  = acos(-1.0d0)
        diammax = 0
c
        allocate(xs(10 000),ys(10 000),whts(10 000))
        allocate(zs(3,10 000))
c        
        norder = 12
        call discquad(ier,norder,ndisc,xs,ys,whts)
        if (ier .ne. 0) then
        call prinf("in triangulate1, ier = *",ier,1)
        stop
        endif
c
        allocate(istack(2,maxtris))
c
        if (maxtris .lt. 10) then
        ier = 4
        return
        endif
c
c       Construct the two top-level triangles.
c 
        nstack        = 0
        ntris         = 0

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
 1000 continue
c
        if (nstack .eq. 0) goto 2000
c
        itri   = istack(1,nstack)
        level  = istack(2,nstack)
        nstack = nstack-1
c
        call disc_trimap(tris(1,itri),amap,det)
c
c       Approximate the diameter of the region and compute
c       the largest mean curvature on the triangle.
c
        dmaxmean = 0
        area     = 0
        do 3000 i=1,ndisc
        s   = amap(1,1)*xs(i)+amap(1,2)*ys(i)+amap(1,3)
        t   = amap(2,1)*xs(i)+amap(2,2)*ys(i)+amap(2,3)
        wht = whts(i)*det
c
        call surface(iparam,s,t,r,dr) 
c
        dn(1) = dr(2,1)*dr(3,2) - dr(3,1)*dr(2,2)
        dn(2) = dr(3,1)*dr(1,2) - dr(1,1)*dr(3,2)
        dn(3) = dr(1,1)*dr(2,2) - dr(2,1)*dr(1,2)
        da    = sqrt(dn(1)**2+dn(2)**2+dn(3)**2)
        area  = area + da*wht
c
        zs(1,i) = r(1)
        zs(2,i) = r(2)
        zs(3,i) = r(3)
c
 3000 continue
c
c        call disc_bounding_ball(ndisc,zs,bx,by,bz,br)
        diam    = br*2
        if (diam .lt. eps) then
        diammax = max(diam,diammax)
        goto 1000
        endif
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



        subroutine triangulate1(ier,iparam,x1,y1,x2,y2,epsarea,
     -    maxtris,ntris,tris)
        implicit double precision (a-h,o-z)
        dimension tris(6,1),simp(3,2),r(3),dr(3,2),dn(3)
        dimension children(6,4),amap(2,3)
c
        double precision, allocatable :: istack(:,:)
        double precision, allocatable :: xs(:),ys(:),whts(:)
c
        ier = 0
        pi  = acos(-1.0d0)
c
        allocate(xs(10 000),ys(10 000),whts(10 000))
c        
        norder = 12
        call discquad(ier,norder,ndisc,xs,ys,whts)
        if (ier .ne. 0) then
        call prinf("in triangulate1, ier = *",ier,1)
        stop
        endif
c        ndisc = (norder+1)*(norder+2)/2
c
        allocate(istack(2,maxtris))
c
        if (maxtris .lt. 10) then
        ier = 4
        return
        endif
c
c       Construct the two top-level triangles.
c 
        nstack        = 0
        ntris         = 0

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
c
        if (nstack .eq. 0) goto 2000
c
        itri   = istack(1,nstack)
        level  = istack(2,nstack)
        nstack = nstack-1
c
        call disc_trimap(tris(1,itri),amap,det)
c
        area = 0
        do 3000 i=1,ndisc
        s   = amap(1,1)*xs(i)+amap(1,2)*ys(i)+amap(1,3)
        t   = amap(2,1)*xs(i)+amap(2,2)*ys(i)+amap(2,3)
        wht = whts(i)*det
        call surface(iparam,s,t,r,dr)
c
        dn(1) = dr(2,1)*dr(3,2) - dr(3,1)*dr(2,2)
        dn(2) = dr(3,1)*dr(1,2) - dr(1,1)*dr(3,2)
        dn(3) = dr(1,1)*dr(2,2) - dr(2,1)*dr(1,2)
        da    = sqrt(dn(1)**2+dn(2)**2+dn(3)**2)
        area  = area + da*wht
 3000 continue
c
        if (area .lt. epsarea) goto 1000
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


        subroutine triangulate2(ier,iparam,x1,y1,x2,y2,nlevels,
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


