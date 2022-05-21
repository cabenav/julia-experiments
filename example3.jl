using LinearAlgebra, Polynomials, SpecialPolynomials
using QuadGK, Calculus, CalculusWithJulia
using Plots


function MatrixM(N,v2,a,b)
   Matt = zeros(N,N)
   aux= [-100 v2 v2;v2 -100 v2;v2 v2 -100]
   for i = 1:Int(N/3)
      Matt[3*(i-1)+1:3*i,3*(i-1)+1:3*i] = aux
   end
   for i = 1:Int(N/3-1)
      Matt[3*i,3*i+1]=b
      Matt[3*i+1,3*i]=b
      Matt[3*i-1,3*i+1]=a
      Matt[3*i,3*i+2]=a
      Matt[3*i+1,3*i-1]=a
      Matt[3*i+2,3*i]=a
   end
   return Matt
end

function eigens(N)
   eigen = zeros(41,N)
   eigen2 = zeros(41,N)
   x = -20:20
   for a=-20:20
      eigen[a+21,:] = eigvals(MatrixM(N,4,a,1))
   end
   for a=-20:20
      eigen2[a+21,:] = eigvals(MatrixM(N,4,1,a))
   end
   display(plot(x,eigen, legend = false))
end

eigens(600)
readline()
