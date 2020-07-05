/* This gp file contains the gp commands which compute the upper bounds for the symmetric functions of degree 8 fields, 
according to Hunter-Pohst-Martinet Theorem and to the estimate of Theorem 24 */



bounds(S1,a8,d)=		/* This command returns the upper bound for the second symmetric function */
{ my(dk,n,T);
n=8;
dk= abs(d);
T= S1^2/n + (dk*n)^(1/(n-1));   /* The estimate is the suitable one for degree 8 fields */
N=abs(a8);
if(N>(T/n)^(n/2), return(oo));
return(T);
}




hpeq(n,N,T)=			/* This command returns the vector y which is mandatory for optimizing the symmetric functions */ 
{ y=vector(n-1);
for(i=1,n-1, 
	y[i]=solve(x=0,1,(n-i)*x^(2*n/i)-T*x^(2*n/i-2)+i*N^(2/i)););  /* Here we use the solve function between the parameters 0 and 1. In our cases of interest, this always worked whenever N=1 */
return(y);
}


hpeqU(n,N,T,x0,x1)=								/* This command is equal to the previous one, the only difference being the choice of the parameters x0 and x1 for the solve 											function. This command is needed whenever N>1: the right choices for the parameters are written in the appendix */ 
{ y=vector(n-1);
for(i=1,n-1, 
	y[i]=solve(x=x0,x1,(n-i)*x^(2*n/i)-T*x^(2*n/i-2)+i*N^(2/i)););
return(y);
}


upbounds(n,y,a,b) =								/* This command returns the upper bounds for the symmetric functions, given the degree, a vector y for the optimization, and two 											parameters a and b which tell the number of symmetric functions which must be estimated. The functions with positive exponents are 											estimated by values in the vector UP, while the ones with negative exponent are bounded by values in vector UN */
{ my(c,maxi,temp);
if(!length(y)==n-1 || a<=2 || b>=0, return(oo));
UP=vector(a-2);
c=abs(b);
UN=vector(c);

for(i=3,a, maxi=0; 
	for(j=1, n-1, temp=j*y[j]^((j-n)*(i/j))+(n-j)*y[j]^i; 
		if(temp>=maxi, maxi=temp);
	); 
	UP[i-2]=maxi; 
); 

for(i=1,c, maxi=0; 
	for(j=1, n-1, temp=j*y[j]^((j-n)*(i/j))+(n-j)*y[j]^i; 
		if(temp>=maxi, maxi=temp);
	);
	UN[i]=maxi;
);
 
}



