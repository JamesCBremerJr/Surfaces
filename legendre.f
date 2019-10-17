        implicit double precision (a-h,o-z)
        parameter (n = 220, m = 200)
c
        dimension xs(n),whts(n),u(n,n),v(n,n)
        dimension eye(n,n)
        dimension pols(n),ders(n),whtslog(n)
        dimension vals(n),coefs(n)
c
        dimension ys(m),ywhts(m),amatrin(m,n)
        dimension vals2(m)
c        
        parameter (mm = 1 000 000)
        dimension zs(mm),vs(mm)
        external fun
c
        double complex values(10 000),z,ima,h1,h0,sumz
        data ima /(0,1)/
c
        call prini(6,14)
        pi = acos(-1.0d0)
c
c       Fetch the n-point quadrature and the transformation matrix.
c
        call legendre(n,xs,whts,u,v)
c
        call prin2("legendre nodes as constructed = *",xs,n)
        call prin2("legendre weights as constructed = *",whts,n)
c
c       Test that u*v = I.
c
        call matmul(n,n,n,u,v,eye)
        errl2=0
        do 1000 i=1,n
        do 1100 j=1,n
        d=0
        if (i .eq. j) d = 1        
        errl2=errl2+abs(eye(i,j)-d)**2
 1100 continue
 1000 continue
        errl2=sqrt(errl2)
        call prin2("u times v difference from identity = *",errl2,1)
c
c       Compute the integral of Sin(x)^2+Arctan[x] on [0,Pi].
c
        sum = 0
        do 1300 j=1,n
        x   = (xs(j)+1)*pi/2
        wht = whts(j)*pi/2
        sum = sum + wht*(sin(x)**2+atan(x))
 1300 continue
        val = 0.5d0*(pi+2*pi*atan(pi)-log(1.0d0+pi**2))
c
        call prina("*")
        call prin2("integral as computed = *",sum,1)
        call prin2("true value = *",val,1)
        call prin2("difference = *",abs(sum-val),1)
c
c       Test the evaluation of a Legendre expansion.
c
        do 2000 j=1,n
        x       = xs(j)
        wht     = whts(j)
        val     = sin(x)*exp(-x)
        vals(j) = val*sqrt(wht)
 2000 continue
        call apply(n,n,u,vals,coefs)
c
        x = .999999d0
        call legeeval(n,coefs,x,sum,sumder)
c
        call prina("*")
        call prin2("evaluation error = *",sum-sin(x)*exp(-x),1)
        call prin2("derivative error = *",
     1      sumder-(cos(x)*exp(-x)-sin(x)*exp(-x)),1)
c
c       Test the construction of interpolation matrices.
c
        call legequad(m,ys,ywhts)
        call legeinterp(n,xs,whts,u,v,m,ys,ywhts,amatrin)
c
        do 2100 j=1,n
        x   = xs(j)
        wht = whts(j)
        vals(j)=sin(x)*sqrt(wht)
 2100 continue
        call apply(m,n,amatrin,vals,vals2)
c
        errl2 = 0
        do 2200 j=1,m
        y     = ys(j)
        ywht  = ywhts(j)
        val   = sin(y)*sqrt(ywht)
        errl2 = errl2 + (vals2(j)-val)**2
c
 2200 continue
        errl2=sqrt(errl2)
        call prin2("interp error = *",errl2,1)
c
        call print_sing("singular values of amatrin = *",m,n,
     1    amatrin)
c
c       Test the log quadratures ... the values have been verified
c       via Mathematica.
c
        z = 1+1.0d-3
        call legelog(n,xs,whts,u,z,whtslog)
c
        call prina("*")
        call prin2("z = *",z,2)
        call prin2("log|z-x| weights = *",whtslog,n)
c
        sum  = 0
        do 3000 j=1,n
        x   = xs(j)
        wht = whtslog(j)
c    
        sum = sum+wht*cos(x)
 3000 continue
c
        rval = -0.408391983349019983309917662897360892d0
c
        call prin3("integral of cos(x)log(z-x) = *",sum,1)
        call prin3("value as computed by Mathematica = *",rval,1)
        call prin2("difference = *",rval-sum,1)
c
        z = 1+ima
        call legelog(n,xs,whts,u,z,whtslog)
c
        call prin2("z = *",z,2)
        call prin2("log|z-x| weights = *",whtslog,n)
c
        sum  = 0
        do 3100 j=1,n
        x   = xs(j)
        wht = whtslog(j)
c    
        sum = sum+wht*cos(x)
 3100 continue
c
        rval = 0.598928475673910538940755387021777203d0
c
        call prin3("integral of cos(x)log(z-x) = *",sum,1)
        call prin3("value as computed by Mathematica = *",rval,1)
        call prin2("difference = *",rval-sum,1)
c
c       Fetch a large quadrature formula and time the results.
c
        call prinf("fetching legendre nodes and weights, m = *",mm,1)

        call elapsed(t1)
        call legequad(mm,zs,vs)
        call elapsed(t2)
c
        call prin2("CPU time = (in seconds) *",t2-t1,1)
c
        end



        subroutine apply(n,m,a,x,y)
        implicit double precision (a-h,o-z)
        dimension a(n,m),x(1),y(1)
        do 1000 i=1,n
        sum=0
        do 1100 j=1,m
        sum=sum+a(i,j)*x(j)
 1100 continue
        y(i)=sum
 1000 continue
        end


        subroutine matmul(n,k,m,a,b,c)
        implicit double precision (a-h,o-z)
        dimension c(n,m),a(n,k),b(k,m)
c
        do 1000 i=1,n
        do 1100 j=1,m        
        sum=0
        do 1200 l=1,k
        sum=sum+a(i,l)*b(l,j)
 1200 continue
        c(i,j)=sum
 1100 continue
 1000 continue
        end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains code for generating Gauss-Legendre quadratures
c       and for manipulating Legendre expansions.
c
c       The following subroutines are user-callable:
c
c    legendre - return the weights and nodes of the n-point Gauss-
c       Legendre quadrature on the interval [-1,1].  Also, return the 
c       associated transformation matrices which take scaled functions
c       values to coefficient expansions and vice-versa.
c
c    legequad - return only the nodes and weights of the n-point Gauss-
c       Legendre quadrature on the interval [-1,1]
c
c    legeeval - evaluate a Legendre expansion and its derivative at a 
c       user-specified point in the interval [-1,1]
c
c    legeinterp - return the matrix interpolating polynomials of degree
c       n-1 from the nodes of the n-point Gauss-Legendre quadrature to a 
c       collection of user-specified nodes in the interval [-1,1]
c
c    lege - evaluate the Legendre polynomials of a given order at a
c       a user-specified point in the interval [-1,1]
c
c    legeders - evaluate the Legendre polynomials of a given order and
c       their derivatives at a user-specified point in the interval 
c       [-1,1]
c
c    legelog - return the weights of an n-point Chebyshev-Legendre
c       quadrature for the kernel log|z-t| where z is an arbitrary
c       complex point
c
c    legeq - evaluate the Legendre functions of the second kind of 
c       order less than or equal to n at a point in the complex plane
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



        subroutine legendre(n,xs,whts,u,v)
        implicit double precision (a-h,o-z)
        dimension xs(1),whts(1),u(n,n),v(n,n)
c
c       Return the nodes and weights of the n-point Gauss-Legendre
c       quadrature as well as the associated transformation matrices 
c       which takes the scaled values of functions at the nodes of the 
c       quadrature to coefficient expansions and vice-versa.
c 
c                              Input Parameters:
c
c   n - the number of nodes in the quadrature formula to be generated
c
c                             Output Parameters:
c
c   xs - the nodes of the n-point quadrature
c   whts - the weights of the n-point quadrature
c   u - the (n,n) matrix which maps scaled function values to 
c       coefficient expansions
c   v - the (n,n) matrix which is the inverse of u
c
c
        if (n .lt. 150) then
        call legequad1(n,xs,whts)
        else
        call legequad2(n,xs,whts)
        endif
c
        do 1000 j=1,n
        call lege(n-1,xs(j),u(1,j))
        do 1100 i=1,n
        d=sqrt(whts(j))
        u(i,j)=u(i,j)*d
 1100 continue
 1000 continue
c
        do 2000 i=1,n
        do 2100 j=1,n
        v(i,j)=u(j,i)
 2100 continue
 2000 continue
c
        do 3000 i=1,n
        do 3100 j=1,n
        u(i,j)=u(i,j)*(2.0d0*i-1.0d0)/2.0d0
 3100 continue
 3000 continue
        end


        subroutine legequad(n,xs,whts)
        implicit double precision (a-h,o-z)
        dimension xs(1),whts(1)
c
c       Return the nodes and weights of the n-point Gauss-Legendre 
c       quadrature.
c
c                              Input Parameters:
c
c   n - the number of nodes in the quadrature formula to be generated
c
c                             Output Parameters:
c
c   xs - the nodes of the n-point quadrature
c   whts - the weights of the n-point quadrature
c
c
        if (n .lt. 150) then
        call legequad1(n,xs,whts)
        else
        call legequad2(n,xs,whts)
        endif
c
        end


        subroutine legeeval(n,coefs,x,val,der)
        implicit double precision (a-h,o-z)
        dimension coefs(n)
c
c       Evaluate an n-term Legendre expansion and its derivative at a 
c       user-specified point in the interval [-1,1].
c
c                          Input Parameters:
c
c   n - the number of terms in the expansion
c   coefs - the n coefficients in the expansion
c   x - the point in the interval [-1,1] at which to evaluate the
c       expansion
c
c                         Output Parameters:
c
c   val - the value of the expansion at the point x
c   der - the derivative of the expansion at the point x
c
        p0 = 1
        p1 = x
c
        d0 = 0
        d1 = 1
c
        val = coefs(1)+coefs(2)*x
        der = coefs(2)
c
        do 1000 k=1,n-2
        p2 = (2*k+1.0d0)*x*p1-k*p0
        p2 = p2/(k+1.0d0)
c
        d2 = (2*k+1.0d0)*x*d1+(2*k+1)*p1-k*d0
        d2 = d2/(k+1.0d0)
c
        val = val+p2*coefs(k+2)
        der = der+d2*coefs(k+2)
c
        p0 = p1
        p1 = p2
c
        d0 = d1
        d1 = d2
 1000 continue
c
        end



        subroutine legeinterp(n,xs,whts,u,v,m,ys,ywhts,amatrin)
        implicit double precision (a-h,o-z)
        dimension xs(1),whts(1),u(n,n),v(n,n)
        dimension ys(1),ywhts(1),amatrin(m,n)
c
        dimension w(n,m)
c
c       Construct the matrix taking the scaled values of polynomials of 
c       degree n-1 from the nodes of the n-point Legendre quadrature 
c       to their scaled values at the nodes of a user-specified 
c       quadrature on [-1,1].
c
c       That is, construct the matrix mapping the vector with entries
c
c          ( f(x_i) \sqrt{w_i} ),
c       
c       where {x_i,w_i} are the nodes and weights of the n-point Gauss-
c       Legendre quadrature, to the vector with entries
c
c          ( f(y_i) \sqrt{v_i} ),
c
c       where {y_i,v_i} are the nodes and weights of a user-specified
c       quadrature formula.
c
c       In the event that the quadrature {y_i,v_i} integrates products
c       of polynomials of degree n-1, this mapping will be perfectly
c       conditioned.
c
c       This routine uses m*n words of stack space.
c
c        
c                            Input Parameters:
c
c   n - the dimension of the space of polynomials which are to be 
c       interpolated
c   xs - the nodes of the n-point Gauss-Legendre quadrature
c   whts - the weights of the n-point Gauss-Legendre quadrature
c   u - the transformation matrix u returned by the subroutine legendre
c   v - the inverse of u, also returned by the subroutine legendre
c
c   m - the number of points in the output 
c   ys - the nodes of the output quadrature
c   ywhts - the weights of the output quadrature
c
c                           Output Parameters:
c
c   amatrin - the (m,n) interpolation matrix requested by the user
c
c
c       Compute the values of the Legendre polynomials at the nodes of
c       the specified quadrature.
c
        do 1000 i=1,m
        y    = ys(i)
        ywht = ywhts(i)
        dd   = sqrt(ywht)
c
        call lege(n-1,y,w(1,i))
        do 1100 j=1,n
        w(j,i)=w(j,i)*dd
 1100 continue
 1000 continue
c
c       Compute the product (w^t )* u and store it in amatrin.
c
        do 2000 i=1,m 
        do 2100 j=1,n
        sum = 0
        do 2200 l=1,n
        sum=sum+w(l,i)*u(l,j)
 2200 continue
        amatrin(i,j)=sum
 2100 continue
 2000 continue
c
        end


        subroutine legelog(n,xs,whts,u,z,whtslog)
        implicit double precision (a-h,o-z)
        dimension xs(1),whts(1),u(n,n),whtslog(1)
        dimension dints(n+1)
        double complex z,qs(n+2)
c
c       Return quadrature weights w_1,\ldots,w_n such that the formula
c
c              1                   n
c         \int   log|x-z| f(x) =  \sum  f(x_j) w_j,                     (1)
c             -1                   j=1
c
c       where x_1, \ldots, x_n are the nodes of the n-point Gauss-
c       Legendre quadrature and z is a user-specified complex-valued
c       point, holds whenever f is a polynomial of degree n-1.
c
c       This subroutine uses 3*n+7 double words of stack space.
c
c                             Input Parameters:
c
c   n - the number of nodes in the quadrature formula to be generated
c   xs - the nodes of the n-point Gauss-Legendre formula 
c   whts - the weights of the n-point Gauss-Legendre formula
c   u - the (n,n) transformation matrix returned by the subroutine 
c       legendre
c   z - the double complex value of z in (1)
c
c                            Output Parameters:
c
c   whtslog - the n weights which make formula (1) exact for polynomials
c       of degree n-1
c
        call mach_zero(eps0)
        eps0=eps0*10
c
c        eps0 = 1.0d-15
c
c       First compute the integrals
c
c          1       1              1
c         ---  \int    P_k(x) * -------    k=0,...,n+1
c          2      -1            z - x
c               
        call legeq(n+1,z,qs)
c
c       Now take the appropriate differences.
c
        do 0100 j=1,n-1
        dints(j+1)=2*(qs(j+2)-qs(j))/(2.0d0*j+1.0d0)
 0100 continue
c
        x = real(z)
        y = imag(z)
c
c       -1<=z<=1 and Im(z) = 0
c
        if(abs(y) .lt. eps0 .AND. abs(x) .lt. 1) then
        dints(1)=-2-z*log(z)+z*log((z*(1 + z))/(-1+z))+log(1-z**2)
        goto 1000
        endif
c
        dints(1) = -2 + Log(-1 + z) + z*(log(1+z)-log(1-z) + Log(-z) 
     1    - Log(z)) + Log(1 + z)
c
 1000 continue
c
        do 2100 j=1,n
        sum = 0
        do 2200 i=1,n
        sum=sum+dints(i)*u(i,j)*sqrt(whts(j))
 2200 continue
        whtslog(j)=sum
 2100 continue
c
        end



        subroutine legequad1(n,xs,whts)
        implicit double precision (a-h,o-z)
        dimension xs(n),whts(n)
        data pi /3.14159265358979323846264338327950288d0/
c
c       Use Newton's method and the recurrence relation to compute the
c       nodes and weights of the n-point Gauss-Legendre quadrature.  The
c       running time of this procedure grows as O(n**2), but it is 
c       quite efficient for small n.  
c
c                            Input Parameters:
c
c   n - the size of the Legendre quadrature to compute
c
c                           Output Parameters:
c
c   xs - the n nodes of the n-point Gauss-Legendre quadrature on the 
c       interval [-1,1]
c   whts - the n weights of the n-point Gauss-Legendre quadrature on 
c       the interval [-1,1]
c
        maxiters = 10
c        pi       = acos(-1.0d0)
c
        call mach_zero(eps0)
        eps0=eps0*10
c
c        eps0     = 1.0d-15
c
c       Use Newton and the recurrence relation to find half of the
c       roots --- the others are obtained via symmetry.
c
c       Note that we also store the value of the derivative at the 
c       obtained root for use in computing the weights below.
c
        ifodd = 0
        nn = (n+1)/2
        if (nn .ne. n/2) then
        ifodd=1
        nn=nn-1
        endif
c
        do 2000 i=nn+1,n
c
c       Use Cheybshev node as initial approximation ....
c
c        x0 = cos(-pi+(2*i-1)*pi/(2*n))
c
c       ... or use a somewhat improved estimate.
c
        dk = i-nn-1
        dn = n
        theta = (4*(ceiling(dn/2)-dk)-1) / (4*dn+2) * pi
c
        x0 = 1.0d0 - (dn-1)/(8*dn**3)-1.0d0/(384.0d0*dn**4)*
     1      (39.0d0-28.0d0/sin(theta)**2)
c        
        x0=x0*cos(theta)
c
c       Conduct Newton iterations.
c
        do 2100 iter=1,maxiters
        call lege0(n,x0,pol,der)
c
        dd = pol/der
        x0 = x0-dd
c
        if (abs(dd) .lt. eps0) then
        xs(i)=x0
        whts(i)=der
        goto 2000
        endif
 2100 continue
c
c       This doesn't happen.
c
        print *,"legeroots bombed!"
        stop
 2000 continue
c
c       Compute the weights using the derivatives we stored above.
c        
        do 3000 j=nn+1,n
        x       = xs(j)
        dd      = 2.0d0/(1-x**2)
        whts(j) = dd/whts(j)**2
 3000 continue
c
c       Reflect the quadrature on [0,1].
c
        do 4000 j=1,nn
        xs(j)   = -xs(n-j+1)
        whts(j) = whts(n-j+1)
 4000 continue
c
c       Handle the root at 0 if n is odd.
c
        if (ifodd .eq. 1) then
        x0          = 0
        call lege0(n,x0,pol,der)
c
        xs(nn+1)   = x0
        whts(nn+1) = 2.0d0/der**2
        endif
c
        end


        subroutine legequad2(n,xs,whts)
        implicit double precision (a-h,o-z)
        dimension xs(n),whts(n)
        data pi /3.14159265358979323846264338327950288d0/
c
c       Use Newton's method and local Taylor expansions to compute the 
c       n roots of the Legendre polynomial of degree n.  This routine 
c       is O(n) in the number of nodes n.
c
c                            Input Parameters:
c
c   n - the size of the Legendre quadrature to compute
c
c                           Output Parameters:
c
c   xs - the n nodes of the n-point Gauss-Legendre quadrature on the 
c       interval [-1,1]
c   whts - the n weights of the n-point Gauss-Legendre quadrature on 
c       the interval [-1,1]
c
        maxiters = 12
c        eps0     = 1.0d-15
        call mach_zero(eps0)
        eps0=eps0*10
        k        = 60
c        if (n .gt. 300) k=50
c
c       Find all of the roots which lie in [0,1].
c
        ifodd = 0
        nn = (n+1)/2
        if (nn .ne. n/2) then
        ifodd=1
        nn=nn-1
        endif
c
c       Use a negative value of x0 to indicate that the procedure has
c       not been initialized.
c
        x0 = -1
c
        do 2000 i=nn+1,n
c
c       Use Chebyshev node as initial guess for roots ...
c
c        x1 = cos(-pi+(2*i-1)*pi/(2*n))
c
c       ... or use this somewhat better approximation.
c
        dk = i-nn-1
        dn = n
        theta = (4*(ceiling(dn/2)-dk)-1) / (4*dn+2) * pi
c
        x1 = 1.0d0 - (dn-1)/(8*dn**3)-1.0d0/(384.0d0*dn**4)*
     1      (39.0d0-28.0d0/sin(theta)**2)
c        
        x1=x1*cos(theta)
c
c       Conduct Newton iterations.
c
        do 2100 iter=1,maxiters
c
c       Evaluate the nth polynomial at x1 using a Taylor expansion
c       and the recurrence relation.  The first evaluation must be
c       handled specially.
c
        if (x0 .lt. 0) then
        call lege0(n,x1,pol,der)
        else
        dx = x1-x0       
        call legetayl(n,k,x0,dx,pol,der)
        endif
c
c       Newton iteration.
c
        x0 = x1
        dd = -pol/der
        x1 = x0+dd
c
        if (abs(dd) .lt. eps0) then
        xs(i)=x1
        whts(i)=der
        goto 2000
        endif
c
 2100 continue
c
c       This doesn't happen.
c
        print *,"legeroots bombed!"
        stop
 2000 continue
c
c       Compute the weights using the derivatives we stored above.
c        
        do 3000 j=nn+1,n
        x       = xs(j)
        dd      = 2.0d0/(1-x**2)
        whts(j) = dd/whts(j)**2
 3000 continue
c
c       Reflect the quadrature on [-1,0] to [0,1].
c
        do 4000 j=1,nn
        xs(j)   = -xs(n-j+1)
        whts(j) = whts(n-j+1)
 4000 continue
c
c       Handle the root at 0 if n is odd.
c
        if (ifodd .eq. 1) then
        x0          = 0
        call lege0(n,x0,pol,der)
c
        xs(nn+1)   = x0
        whts(nn+1) = 2.0d0/der**2
        endif
c
        end


        subroutine legetayl(n,k,x,dx,pol,der)
        implicit double precision (a-h,o-z)
c
c       Evaluate the Legendre polynomial of order n using a Taylor
c       expansion and the recurrence relation for the power series
c       coefficients.
c
        k0 = min(k,n)
c
        p1     = pol
        p2     = der*dx
c
        sum    = p1+p2
        sumder = p2/dx
c
        do 1000 j=0,k0-2
c
        d = 2*x*(j+1.0d0)**2/dx*p2-(n*(n+1.0d0)-j*(j+1.0d0))*p1
        d = d/(j+1.0d0)/(j+2.0d0)*dx**2/(1.0d0-x**2)
c
        sum    = sum+d
        sumder = sumder+d*(j+2)/dx
c
        p1 = p2
        p2 = d
 1000 continue
c
        pol = sum
        der = sumder
c
        end


        subroutine legeq(n,z,vals)
        implicit double precision (a-h,o-z)
        double complex vals(1)
        double complex z,ima,q00,q11
        double complex q0,q1,q2,alpha
        data ima / (0,1) /
c
c       Compute the values of the Legendre functions of the second kind
c       of orders 0 through n at the (complex-valued) point z.
c
c       Note that the branch cut is chosen depending on the value of
c       z so as to make the identity
c
c                 1    P_n(x)
c       P.V.  \int    --------   = 2*Q_n(z)
c                -1     z-x
c
c       hold for all z.
c
c                          Input Parameters:
c
c   n - the order of the Legendre functions to evaluate
c   z - the point in the complex plane at which to evaluate 
c
c                          Output Parameters:
c
c   vals - an array containing the (n+1) requested values
c
c        eps0 = 1.0d-15
        call mach_zero(eps0)
        eps0=eps0*10
c
c
c       Pick the branch cut we want based on the placement of z:
c

c
c       For z in [-1,1]
c
        x = real(z)
        y = imag(z)
        if (abs(y) .le. eps0) then
        if (x .gt. 1) then
        q00 = (log(1+x)-log(x-1))/2
        goto 1000
        endif
        if (x .lt. -1.0d0) then
        q00 = (log(-1-x)-log(1-x))/2
        goto 1000
        endif
        q00=(log(1+x)-log(1-x))/2
        goto 1000
        endif
c
c       For z in the complement of the real line.
c
        q00=((log(1+z)-log(-1+z)))/2
c
 1000 continue
c
        q11 = (-1.0d0+z*q00)
c
c       If |z| is small we can run the recurrence forward.
c
        if (abs(z) .gt. 1.1d0) goto 3000
c
        vals(1) = q00
        vals(2) = q11
c
        do 2000 j=2,n-1
        vals(j+1)=((2*j-1.0d0)*z*vals(j)-(j-1.0d0)*vals(j-1))/j
 2000 continue
        return
 3000 continue
c
c       The backward recurrence is called when |z|>>1.
c
        q2 = 1
        q1 = 0
c
c       Run the recurrence backward ... we might need more than n 
c       values.
c        
        k = 2*n
c
        do 3100 j=k,1,-1
c
        dd1 = j
        dd1 = (j+1)/dd1
c
        dd2 = j
        dd2 = (2*j+1)/dd2
c
        q0 = -dd1*q2+dd2*z*q1
        if (j .le. n) vals(j)= q0
c
        q2 = q1
        q1 = q0
 3100 continue
c
c       Normalize the computed values.
c
        alpha = q00/q0
c
        do 3300 j=1,n
        vals(j)=vals(j)*alpha
 3300 continue
c
        end



        subroutine lege(n,x,vals)
        implicit double precision (a-h,o-z)
        dimension vals(1)
c     
c       Evaluate the (n+1) Legendre polynomials of order 0 through n at
c       the point x.
c
c                          Input Parameters:
c
c   n - the order of the Legendre polynomials to evaluate
c   x - the point in the interval [-1,1] at which to evaluate them
c
c                         Output Parameters:
c
c   vals - the requested values
c
c   
c
        vals(1) = 1
        if (n .eq. 0) return
c
        vals(2) = x
        if (n .eq. 1) return
c
        do 1000 j=2,n
        vals(j+1)=((2.0d0*j-1.0d0)*x*vals(j)-(j-1.0d0)*vals(j-1))/j
 1000 continue
c
        end



        subroutine legeders(n,x,pols,ders)
        implicit double precision (a-h,o-z)
        dimension pols(1),ders(1)
c
c       Evaluate the (n+1) Legendre polynomials of order 0 through n
c       at the point x as well as their derivatives.
c
c                          Input Parameters:
c
c   n - the order of the Legendre polynomials to evaluate
c   x - the point in the interval [-1,1] at which to evaluate them
c
c                         Output Parameters:
c
c   pols - the values of the polynomials
c   ders - the values of the derivatives
     
c       Evaluate the n+1 normalized Legendre polynomials of degree 0 
c       through n and their derivatives at the point x.
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
        pols(j+1)=((2.0d0*j-1.0d0)*x*pols(j)-(j-1.0d0)*pols(j-1))/j
c
c        ders(j+1)=((2.0d0*j-1.0d0)*(x*ders(j)+pols(j))-
c     1    (j-1.0d0)*ders(j-1))/j
c
 1000 continue
 2000 continue
c
c       Compute the derivatives.
c
        d=x**2.0d0-1.0d0
        do 3000 j=3,n+1
        ders(j)=(j-1.0d0)*(x*pols(j)-pols(j-1.0d0))/d
 3000 continue
c
        end


        subroutine lege0(n,x,pol,der)
        implicit double precision (a-h,o-z)
c
c       Evaluate the Legendre polynomial of degree n and its derivative
c       at the point x using well-known formulae.
c
c       NOTE: this routine does not evaluate the normalized polynomial
c 
c                           Input Parameters:
c
c   n - degree of the Legendre polynomial to evaluate
c   x - point in the interval [-1,1] at which to evaluate the polynomial
c
c                          Output Parameters:
c
c   pol - the value of the polynomial at that point
c   der - the value of the derivative of the polynomial at that point
c
c
c       Handle orders 0 and 1.
c
        if (n .eq. 0) then
        pol = 1
        der = 0
        return
        endif
c
        if (n .eq. 1) then
        pol = x
        der = 1
        return
        endif
c
c       Use the well-known recurrence relation.
c
        p1 = 1
        p2 = x
c
        do 1000 j=2,n
        p  = ((2*j-1)*x*p2-(j-1)*p1)/j
        p1 = p2
        p2 = p
 1000 continue
c
        pol = p
c
c       Compute the derivative using another well-known formula.
c
        d=x**2-1
        der=n*(x*p2-p1)/d
c
        end
