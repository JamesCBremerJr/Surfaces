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


