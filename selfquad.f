        implicit double precision (a-h,o-z)
        dimension verts(2,3),dr(3,2),r(3),dn(3)
        double precision, allocatable :: rad(:)
        double precision, allocatable :: xs(:),ys(:),whts(:)
        dimension xeval(13),yeval(13)
        double complex sum1,sum2,val
        external kernel
c
        pi = acos(-1.0d0)
c
        lrad   = 5000 000
        allocate(rad(lrad))
        norder = 8
c
        maxquad = 5 000 000
        allocate(xs(maxquad),ys(maxquad),whts(maxquad))

        call radinit(ier,norder,rad,lrad,lkeep)
c
        if (ier .ne. 0) then
        call prinf("after radinit, ier = *",ier,1)
        stop
        endif
c
        iparam = 1
        x0     = 0.90d0
        y0     = 0.01d0
c
        verts(1,1) = cos(0.0d0)
        verts(2,1) = sin(0.0d0)
c
        verts(1,2) = cos(2*pi/3)
        verts(2,2) = sin(2*pi/3)
c
        verts(1,3) = cos(4*pi/3)
        verts(2,3) = sin(4*pi/3)
c
        scale = .15d0
c
        do 0100 i=1,3
        verts(1,i) = verts(1,i)*scale
        verts(2,i) = verts(2,i)*scale
 0100 continue
        x0 = x0*scale
        y0 = y0*scale
c
c       Build evaluation data for the target point.
c
        call surface(iparam,x0,y0,r,dr)
c
        call prin2("r=*",r,3)
        call prin2("dr=*",dr,6)
c
        dn(1) = dr(2,1)*dr(3,2) - dr(3,1)*dr(2,2)
        dn(2) = dr(3,1)*dr(1,2) - dr(1,1)*dr(3,2)
        dn(3) = dr(1,1)*dr(2,2) - dr(2,1)*dr(1,2)
        da    = sqrt(dn(1)**2+dn(2)**2+dn(3)**2)
        dn(1) = dn(1)/da
        dn(2) = dn(2)/da
        dn(3) = dn(3)/da
c
        xeval(1)  = 1
        xeval(2)  = r(1)
        xeval(3)  = r(2)
        xeval(4)  = r(3)
        xeval(5)  = dr(1,1)
        xeval(6)  = dr(2,1)
        xeval(7)  = dr(3,1)
        xeval(8)  = dr(1,2)
        xeval(9)  = dr(2,2)
        xeval(10) = dr(3,2)
        xeval(11) = dn(1)
        xeval(12) = dn(2)
        xeval(13) = dn(3)
c
c       Adaptively construct the quadrature.
c
        call mach_zero(eps)
        eps = eps*100
c
        n      = 70 000
        call elapsed(t1)
        call selfquad0(ier,n,norder,verts,x0,y0,
     -    maxquad,nquad,xs,ys,whts)
        call elapsed(t2)
        call prin2("selfquad0 time = *",t2-t1,1)
c
        if (ier .ne. 0 ) then
        call prinf("after selfquad0, ier = *",ier,1)
        stop
        endif
c
        call prinf("after selfadap0, nquad = *",nquad,1)
c
        sum1 = 0
c
        do 1000 i=1,nquad
        x   = xs(i)
        y   = ys(i)
        wht = whts(i)
c
        call surface(iparam,x,y,r,dr)
c
        dn(1) = dr(2,1)*dr(3,2) - dr(3,1)*dr(2,2)
        dn(2) = dr(3,1)*dr(1,2) - dr(1,1)*dr(3,2)
        dn(3) = dr(1,1)*dr(2,2) - dr(2,1)*dr(1,2)
        da    = sqrt(dn(1)**2+dn(2)**2+dn(3)**2)
        dn(1) = dn(1)/da
        dn(2) = dn(2)/da
        dn(3) = dn(3)/da
c
        yeval(1)  = 1
        yeval(2)  = r(1)
        yeval(3)  = r(2)
        yeval(4)  = r(3)
        yeval(5)  = dr(1,1)
        yeval(6)  = dr(2,1)
        yeval(7)  = dr(3,1)
        yeval(8)  = dr(1,2)
        yeval(9)  = dr(2,2)
        yeval(10) = dr(3,2)
        yeval(11) = dn(1)
        yeval(12) = dn(2)
        yeval(13) = dn(3)
c
        call kernel(xeval,yeval,val,par1,par2,par3,par4) 
        sum1 = sum1 + val*wht
 1000 continue
c      
c       Now use the precomputed quadrature.
c
        call surface(iparam,x0,y0,r,dr)     
        call elapsed(t1)
        call selfquad(ier,rad,verts,x0,y0,dr,maxquad,
     -   nquad,xs,ys,whts)
        call elapsed(t2)
        call prin2("selfquad time = *",t2-t1,1)
c
        call prinf("after selfquad, nquad = *",nquad,1)
c
        sum2 = 0
        do 2000 i=1,nquad
        x   = xs(i)
        y   = ys(i)
        wht = whts(i)
c
        call surface(iparam,x,y,r,dr)
c
        dn(1) = dr(2,1)*dr(3,2) - dr(3,1)*dr(2,2)
        dn(2) = dr(3,1)*dr(1,2) - dr(1,1)*dr(3,2)
        dn(3) = dr(1,1)*dr(2,2) - dr(2,1)*dr(1,2)
        da    = sqrt(dn(1)**2+dn(2)**2+dn(3)**2)
        dn(1) = dn(1)/da
        dn(2) = dn(2)/da
        dn(3) = dn(3)/da
c
        yeval(1)  = 1
        yeval(2)  = r(1)
        yeval(3)  = r(2)
        yeval(4)  = r(3)
        yeval(5)  = dr(1,1)
        yeval(6)  = dr(2,1)
        yeval(7)  = dr(3,1)
        yeval(8)  = dr(1,2)
        yeval(9)  = dr(2,2)
        yeval(10) = dr(3,2)
        yeval(11) = dn(1)
        yeval(12) = dn(2)
        yeval(13) = dn(3)
c
        call kernel(xeval,yeval,val,par1,par2,par3,par4)
        sum2 = sum2 + val*wht
 2000 continue
c
c
        print *,""
        print *,sum1
        print *,sum2
        print *,""
c
        print *,abs(sum2-sum1)
        print *,abs(sum2-sum1)/abs(sum1)
c        
        end


        subroutine disk0(n,nlege,rcut,
     -    maxquad,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension xslege(nlege+10),whtslege(nlege+10)
        dimension xs(1),ys(1),whts(1)
c
        data pi      / 3.14159265358979323846264338327950288d0 /
        external funuser
c
c       Construct a quadrature for a circle of a given radius.

        ier   = 0
        nquad = 0 
c
        dd = n
        dd = 2*pi/dd
c
        call legequad(nlege,xslege,whtslege)
c
        r0 = 0
        r1 = rcut
c
        do 1000 i=1,nlege
        r    = (r1-r0)/2*xslege(i) + (r1+r0)/2
        rwht = (r1-r0)/2*whtslege(i)
c
        do 1100 j=1,n
        t    = dd*(j-1)
        twht = dd
c
        x   = r*cos(t)
        y   = r*sin(t)
        wht = twht*rwht*r
c
        nquad = nquad+1
        xs(nquad)   = x
        ys(nquad)   = y
        whts(nquad) = wht
 1100 continue
 1000 continue
        end




        subroutine surface(iparam,s,t,r,dr)
        implicit double precision (a-h,o-z)
        dimension r(3),dr(3,2)        
c
c       Supply a parameterization of a torus.
c
        a = 2.000d0
        b = 0.001d0
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
        end


        subroutine kernel(x,y,val,par1,par2,par3,par4)
        implicit double precision (a-h,o-z)
        double complex val
        dimension x(1),y(1)
c        call kdouble0(x,y,val)
        call ksingle0(x,y,val)
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
        val   = ((x1-y1)*dny1 + (x2-y2)*dny2 + (x3-y3)*dny3)
        val   = val /dd**3
        val   = -val*xwht*ywht*over4pi
c
        end



        subroutine ksingle0(x,y,val)
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
        val   = -xwht*ywht*over4pi/dd
c
        end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the
c       code proper.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains code for generating quadratures for the
c       evaluation of integrals of the form
c
c           \int    K(p(x),p(y)) \sigma(p(y)) |dp(y)|                     (1)
c               T
c       where
c
c       1. T is a triangle in the plane containing the point x;
c
c       2. p: T --> R^3 is a smooth one-to-one mapping;
c
c       3. \sigma(x) is a smooth function;
c
c       4. K is a linear combination of the kernels
c
c                    exp(i k |x-y|)
c         S(x,y)   = --------------
c                         |x-y|
c  
c         D(x,y)   = \nabla_y S(x,y) \cdot \eta_y
c
c         D*(x,y)  = \nabla_x S(x,y) \cdot \eta_x
c
c         R(x,y)   = \nabla_x S(x,y) 
c
c       Here, \nabla_x denotes the gradient in the x-variable while
c       \eta_x denotes the outward-pointing normal unit vector to to the
c       surface under consideration at the point x.
c
c       The behavior of the Jacobian dp of the parameterization
c       p at the point x has a large influence on the form of the 
c       integrand of (1).  This code proceeds by first composing the 
c       given parameterization p with an appropriate linear mapping 
c
c                 A: \mathbb{R}^2 \to \mathbb{R}^2
c
c       in order to  form a new parameterization p' = p \ocirc A such 
c       that the Jacobian of p' at the point x is conformal.  
c
c       Then the quadrature rules generated by the code in radial.f 
c       are used to evaluate the resulting integral.
c
c       The following subroutines are user-callable:
c
c   selfquad - return a quadrature for evaluating an integral of
c       the form (1); a precomputed table of "generalized Legendre rules"
c       and a precomputed table of quadrature rules stored
c       in a "rad" structure (see radial.f) is used to generate the rule
c
c   selfquad0 - return a quadrature for evaluating an integral of
c       the form (1); the rule returned by this subroutine is 
c       constructed using double exponential quadrature rules instead
c       of a precomputed table
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        subroutine selfquad(ier,rad,verts0,x0,y0,dr,maxquad,
     -    nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension verts(2,3),dr(3,2),a(2,2),ainv(2,2),verts0(2,3)
        dimension xs(1),ys(1),whts(1),b(3,2)
c
c       Return a quadrature formula for evaluating integrals of the
c       form (1).
c
c                            Input Parameters:
c
c    rad - the rad structure generated by the radinit subroutine in
c       radial.f
c    verts - a (2,3) matrix which specifys the vertices of the triangle
c       over which the surface element is parameterized
c    (x0,y0) - the coordinates (in the parameterization variables)
c       of the target point x
c    dr - the (3,2) Jacobian of the parameterization at (x0,y0)
c    maxquad - the largest quadrature rule which can be accomodated
c       by the user-supplied output arrays xs, ys and whts
c
c                           Output Parameters:
c
c    ier - an error return code;
c       ier = 0     indicates succesful execution
c       ier = 16    means that the requested quadrature rule is of length
c                   larger than maxquad
c       ier = 64    means that one of the necessary quadrature rules is 
c                   missing from the table in the rad structure; typically, 
c                   this  indicates that the triangle is very close to 
c                   degenerate
c       ier = 128   means that one of more of the vertices is *very*
c                   close to the origin
c       ier = 256   means the origin is *very* close to the boundary
c                   of the triangle
c       ier = 512   means that the vertices are very nearly colinear
c
c    nquad - the number of nodes in the resulting quadrature formula
c    xs - this user-supplied array will contain the x coordinates of 
c       the quadrature nodes upon return
c    ys - this user-supplied array will contain the y coordinates of 
c       the quadrature nodes upone return
c    whts - this user-supplied array will contain the quadrature weights
c       upon return
c
        ier = 0
c
c       Find a linear map A which makes the Jacobian at the target node
c       conformal.  Also, compute the inverse of A.
c
        ainv=0
        call self_findmap(dr,a,ainv)
c
c       Apply the inverse of A to the vertices of the triangle T in
c       order to find the vertices of the triangle T_0.
c
        do 2000 i=1,3
        x = verts0(1,i)-x0
        y = verts0(2,i)-y0
c
        xx = ainv(1,1)*x+ainv(1,2)*y
        yy = ainv(2,1)*x+ainv(2,2)*y
c
        verts(1,i) = xx
        verts(2,i) = yy
 2000 continue
c
c       Fetch a quadrature on T_0 for radially singular functions.
c
        call radquad(ier,rad,verts,maxquad,nquad,xs,ys,whts)
        if (ier .ne. 0) return
c
c       Apply the mapping A to the quadrature formula.
c
        det = abs(a(1,1)*a(2,2)-a(1,2)*a(2,1))
c
        sum1=0
        sum2=0
c
        do 3000 i=1,nquad
        s   = xs(i)
        t   = ys(i)
        wht = whts(i)
c
        sum1=sum1+wht
c       
        u = a(1,1)*s + a(1,2)*t
        v = a(2,1)*s + a(2,2)*t
        wht = wht*det
c
        sum2=sum2+wht
c
        xs(i)   = u+x0
        ys(i)   = v+y0
        whts(i) = wht
 3000 continue
c

        return
        end


        subroutine selfquad0(ier,n,norder,verts0,x0,y0,
     -    maxquad,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension verts0(2,3),verts(2,3),xs(1),ys(1),whts(1)
c
c       This routine returns a quadrature formula for the evaluation
c       of an integral of the form (1).  The rule returned by this
c       routine is constructed using piecewise double exponential
c       quadrature rules and a precomputed "generalized Legendre"
c       quadrature.
c
c                              Input Parameters:
c
c   n - an integer parameter controlling the order of the double 
c       exponential  rule used
c   norder - the order of the generalized Legendre rule; the possible
c       values for norder are 4, 8, 12 and 16.
c   verts0 - a (2,3) array each column of which specifies a 
c   maxquad - the largest quadrature rule which can be accomodated by the
c       user-supplied output arrays xs, ys and whts
c
c                             Output Parameters:
c
c   ier - an error return code;
c       ier = 0     indicates successful execution
c       ier = 16    means that the requested quadrature is of length
c                   greater than maxquad
c       ier = 512   means that the specified triangle was too degenerate to
c                   be accomodated by the precomputed quadratures in the
c                   subroutine pvquad
c       ier = 1024  means that an invalid value of norder was specified
c
c   (nquad,xs,ys,whts) - upon successful execution, the requested 
c       quadrature rule
c
        do 1000 i=1,3
        verts(1,i) = verts0(1,i)-x0
        verts(2,i) = verts0(2,i)-y0
 1000 continue
c
        call radquad0(ier,verts,n,norder,maxquad,nquad,xs,ys,whts)
        if (ier .ne. 0) return
c
        do 2000 i=1,nquad
        xs(i) = xs(i)+x0
        ys(i) = ys(i)+y0
 2000 continue
c
        end


        subroutine self_findmap(a,b,binv)
        implicit double precision (a-h,o-z)
        dimension a(3,2),b(2,2),binv(2,2),q(3,2)
c
c       Given a (3,2) matrix A, return a (2,2) matrix B such that
c       A*B has orthonormal columns.  Also, return the inverse
c       BINV of the mapping B.
c
c
c                             --- WARNING ---
c       THIS ROUTINE IS A BIT ROUGH --- IT IS EXPECTED TO FAIL
c       IN THE EVENT THAT THE MATRIX A IS NEARLY SINGULAR.  BUT
c       WHAT ARE YOU DOING USING PARAMETERIZATION WITH JACOBIANS
c       WHICH ARE NEARLY SINGULAR ANYWAY?
c                              --------------
c
c
c                          Input Parameters:
c
c   a - the (3,2) input matrix 
c
c                         Output Parameters:
c
c   b - a (2,2) matrix such that a*b has orthonormal columns
c   binv - the (2,2) inverse of b
c
c
c       Orthonormalize the columns of the matrix a.
c
        sum1 = sqrt(a(1,1)**2 + a(2,1)** 2 + a(3,1)**2)
        q(1,1) = a(1,1) / sum1
        q(2,1) = a(2,1) / sum1
        q(3,1) = a(3,1) / sum1
c
        dip = a(1,2)*q(1,1)+a(2,2)*q(2,1)+a(3,2)*q(3,1)
c
        q(1,2) = a(1,2) - dip * q(1,1)
        q(2,2) = a(2,2) - dip * q(2,1)
        q(3,2) = a(3,2) - dip * q(3,1)
c
        sum2 = sqrt(q(1,2)**2 + q(2,2)** 2 + q(3,2)**2)
        q(1,2) = q(1,2) / sum2
        q(2,2) = q(2,2) / sum2
        q(3,2) = q(3,2) / sum2
c
c       Compute BINV = Q'*A.
c     
        binv(1,1) = q(1,1)*a(1,1) + q(2,1)*a(2,1) + q(3,1)*a(3,1)
        binv(1,2) = q(1,1)*a(1,2) + q(2,1)*a(2,2) + q(3,1)*a(3,2)
        binv(2,1) = q(1,2)*a(1,1) + q(2,2)*a(2,1) + q(3,2)*a(3,1)
        binv(2,2) = q(1,2)*a(1,2) + q(2,2)*a(2,2) + q(3,2)*a(3,2)
c
c       Compute B = (BINV)^(-1).
c
        det = binv(1,1)*binv(2,2) - binv(1,2)*binv(2,1)
        b(1,1) = binv(2,2)/det
        b(2,2) = binv(1,1)/det
        b(1,2) = -binv(1,2)/det
        b(2,1) = -binv(2,1)/det
c
        end

