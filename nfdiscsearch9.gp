install(ZpX_primedec,GG);					/*It installs an ad hoc command needed for the programs to work */

\r bounds_n9.gp;						/* Here we read the commands needed for the upper bounds of symmetric functions */

{print("\n  
	Welcome to the nfdiscsearch for degree 9 fields. \n 
        To launch the search, write: \n
	nfdiscsearch9(r1,S1,S2,S3,a9,d,valmod2)\n 
	where the parameters are as follows: \n 
	* r1:= real signature of the fields, so an integer among 1,3,5,7, and 9;
	* S1:= value of the trace of the HPM-elements, which is an integer between 0 and 4;
	* S2:= value of the second symmetric function: it must be an integer congruent to S1 mod 2;
	* S3:= value of the third symmetric function: it must be an integer congruent to S1 mod 3; 
	* a9:= value of the last coefficient of the polynomials, i.e. the opposite of the norm of the HPM-elements;
	*  d:= upper bound for the discriminants of the number fields in absolute value;
	* valmod2: write 1 to generate polynomials p(x) such that p(1) is odd, write 0 for p(1) even. ");
}

nfdiscsearch9(r1,S1,S2,S3,a9,d,valmod2)={

L=List();				/* List where we put the polynomials surviving to the tests */



my(T,a1,a,b,n,a7bound,a6bound,P1A,P1B,e,c,k2,p2,a2,L2,k3,p3,a3,L3,t1,k4,S4,p4,a4,L4,k5,S5,p5,a5,L5,k6,S6,p6,a6,L6,sym,k7,S7,p7,a7,L7,sym2,k8,S8,p8,a8,L8,adj,s,v,s1,sm1,G,t3,t4,valmod8);

a1=-S1; a=8; b=-2; n=9;

T=bounds(S1,a9,d);


/* Here we use the commands of bounds_n8.gp in order to compute the symmetric functions */
if(N==1, 
	hpeq(n,N,T);
);

x0=0;
x1=1.36;
if(N!=1,
	hpeqU(n,N,T);
);


upbounds(n,y,a,b);

/* Here we give some upper bounds for the last two coefficients of the generated polynomials G(x) and two numbers P1A and P1B, which will bound G(1) and G(-1) respectively */

a8bound=floor(UN[1])*abs(a9);
a7bound=floor(UN[2])*abs(a9);

P1A=floor(((T-2*S1+n)/n)^(n/2));
P1B=floor(((T+2*S1+n)/n)^(n/2));


c=0;							/* Counter for the number of polynomials generated during the tests */
e=0;							/* Counter for the number of polynomials survived to the tests */


k2=(-a1*S1) % 2 ;					/* Here we use the congruence relation between the first and the second symmetric function to generate the range of values for the symmetric function S2, 								including the number p2 of steps and the lower bound L2 for S2, together with the coefficient a2 */	
k22=S2%2;
if(k2!=k22, break;); 	
p2=floor((T-k2)/2)+floor((T-(2-k2))/2)+1;
S22=floor((T-k2)/2)*2+k2;					
L2=S22-2*p2+(2*a1^2)/n;
if(S2<L2, break;);
a2=(-S2-a1*S1)/2;


	k3=(-a2*S1-S2*a1) % 3;					/* Here we use the congruence relations to generate the range of values for the symmetric 									function S3, including the number p3 of steps and the lower bound L3 for S3, together with the coefficient a3 */
	k33=S1%3;
	if(k3!=k33, break;);
	p3=floor((UP[1]-k3)/3)+floor((UP[1]-(3-k3))/3)+1;
	S33=floor((UP[1]-k3)/3)*3+k3;
	L3=S3-3*p3;	
	if(a1==0,
		L3=0;
	);
	if(S3 < L3, break;);
	a3=(-S3-a2*S1-S2*a1)/3;
	

		k4=(-a3*S1-a2*S2-a1*S3) % 4;				/* Here we use the congruence relations to generate the range of values for the symmetric 										function S4, including the number p4 of steps and the lower bound L4 for S4, together with the coefficient a4 */
		S4=floor((UP[2]-k4)/4)*4+k4;
		p4=floor((UP[2]-k4)/4)+floor((UP[2]-(4-k4))/4)+1;
		a4=(-S4-a3*S1-a2*S2-a1*S3)/4;
		L4=S4-4*p4;
		t1=getabstime();					/* It begins to register the time needed for the computation */

		while(S4>=-2*(T-S2)^2 && abs(S3)<=sqrt((S2+T)*(S4+2*(T-S2)^2)/2) && S4>=L4,			/* While the value for S4 is acceptable and other conditions are satisfied, generate a range for 															the fifth symmetric function S5, including the number of steps p5 and the lower bound L5, together 															with the coefficient a5 */

			k5=(-a4*S1-a3*S2-a2*S3-a1*S4) % 5;
			S5=floor((UP[3]-k5)/5)*5+k5;
			p5=floor((UP[3]-k5)/5)+floor((UP[3]-(5-k5))/5)+1;
			a5=(-S5-a4*S1-a3*S2-a2*S3-a1*S4)/5;		
			L5=S5-5*p5;
			while(S5>=L5,										/* While the value for S5 is acceptable, generate a range for 															the sixth symmetric function S6, including the number of steps p6 and the lower bound L6, together 															with the coefficient a6 */

				k6=(-a5*S1-a4*S2-a3*S3-a2*S4-a1*S5) % 6;
				S6=floor((UP[4]-k6)/6)*6+k6;
				p6=floor((UP[4]-k6)/6)+floor((UP[4]-(6-k6))/6)+1;
				a6=(-S6-a5*S1-a4*S2-a3*S3-a2*S4-a1*S5)/6;	
				L6=S6-6*p6;

					

				
				while(S6>=L6,									/* While the value for S6 is acceptable, generate a range for 															the seventth symmetric function S7, including the number of steps p7 and the lower bound L7, 															together with the coefficient a7. Morevor, it begins to verify additional conditions */

					sym=-a6*S1-a5*S2-a4*S3-a3*S4-a2*S5-a1*S6;
					
					k7= sym%7;
					S7=floor((UP[5]-k7)/7)*7+k7;
					p7=floor((UP[5]-k7)/7)+floor((UP[5]-(7-k7))/7)+1;
					a7=(-S7+sym)/7;		
					L7=S7-7*p7;

					while(S7>=L7,								/* While the value for S7 is acceptable, generate a range for 															the seventth symmetric function S8, including the number of steps p8 and the lower bound L8, 															together with the coefficient a8. Morevor, it begins to verify additional conditions */
						sym2=-a7*S1-a6*S2-a5*S3-a4*S4-a3*S5-a2*S6-a1*S7;
						k8= sym2%8;
						S8=floor((UP[6]-k8)/8)*8+k8;
						p8=floor((UP[6]-k8)/8)+floor((UP[6]-(8-k8))/8)+1;
						a8=(-S8+sym2)/8;		
						L8=S8-8*p8;

						if(a8<-a8bound*abs(a9),						/* If the coefficient a8 is below -a8bound*abs(a9), then increase it in such a way to prevent this, 															modifying the values of S8 as consequence */
							adj=-a8-a8bound*abs(a9);
							a8=a8+adj;
							S8=S8-8*adj;
						);

						v=[1,a1,a2,a3,a4,a5,a6,a7,a8,a9];				/* Vector with the coefficients of the generated polynomial G(x), which must be tested */
						s1=vecsum(v);							/* This is G(1) */
						
						if((s1<-P1A),							/* If G(1) <-P1A,increase its value by increasing the value of a8 */
							adj=-s1-P1A;
							a8=a8+adj;
							S8=S8-8*adj;
							v[#v-1]=a8;
							s1=vecsum(v);
						);

						if(s1%2!=valmod2,						/* If G(1) has not the parity condition expressed by valmod2, then increase a7 of one in order to 															reach the correct parity */
							S8=S8-8;
							a8=a8+1; 
							v[#v-1]=a8;
							s1=vecsum(v);
						);

						valmod8=s1%8;							/* If G(1) is even, put a bonus parameter which allows to decrease of 4 units the value of S7, 															after choosing the correct class mod 8 which prevents p(1) from being an exact multiple of 2 or 														4. */		
						bonus=1;

						if( valmod2==0,
							S8=S8-8*valmod8;
							a8=a8+valmod8;
							v[#v-1]=a8;
							s1=vecsum(v);
							bonus=4;
						);				
					
						sm1=vecsum([-1,a1,-a2,a3,-a4,a5,-a6,a7,-a8,a9]);		/* This is G(-1) */
						
						G=Pol(v);							/* Here we give the instance of the polynomial G(x) which must be tested */
					
						while(abs(a8)<= a8bound &&  s1<=P1A && S8>=L8,			/* While the coefficient a8 is small, G(1) is less than P1A and S8 attains acceptable values, we 															test the polynomial G we have just created */

							c++;							/* Increase of one the number of created polynomials */
			
						/* The tests: upper bound on a6 and G(-1); G(x) must be squarefree; check on the small norm of prime ideals via the ZpX_primedec command; check on the number of real roots of the polynomial, in order to get the correct signature; irreducibility of the test polynomial; the squarefree part of the discriminant of the polynomial must be low */

							if(abs((a8^2)/a9-2*a7) >a7bound || abs(sm1)>P1B || !issquarefree(G)  || ZpX_primedec(G,2)[1,1]<3 || ZpX_primedec(G,3)[1,1]<2 || ZpX_primedec(G,5)[1,1]<2 || polsturm(G)!=r1 ||!polisirreducible(G) || abs(coredisc(poldisc(G)))> d, 					
		
							S8=S8-16*bonus;						/* If the test fails, skip to the next polynomial */
							a8=a8+2*bonus;
							s1=s1+2*bonus;
							sm1=sm1-2*bonus;
							G=G+2*bonus*x;
							next;
							);
						
							e++;					/* If the test succeeds, increase the number of survived polynomials, put G(x) in the list L and skip to the next 													polynomial */
							listput(L,G);	
							a8=a8+2*bonus;
							S8=S8-16*bonus;
							s1=s1+2*bonus;
							sm1=sm1-2*bonus;
							G=G+2*bonus*x;
						);	
												/* Here and in the following lines, whenever the inner cicle for S_{n+1} is over, decrease S_n, increase a_N  and 													restart the inner cycle. Go on until Sn < Ln*/
						S7=S7-7;					
						a7=a7+1;
					);		
					
					S6=S6-6;
					a6=a6+1;
				);
				
				S5=S5-5;
				a5=a5+1;
			);
			S4=S4-4;
			a4=a4+1;
		);

		
	t2=getabstime();								/* Take the time of computation for the given values of S2 and S3 */			
	print((t2-t1)/1000." seconds for S3 = "S3);


print("Number of polynomials created = " c);
print("Number of polynomials after the first test = " e);

t3=getabstime();					/* It begins to take the time of computation for the second test */								

LD=List();				/* List where we gather the polynomials of L which have small number field discriminant */					

for(i=1,#L,				/* Here we check if the polynomials of L have small number field discriminant */
	if(abs(nfdisc(L[i]))<d, 
		listput(LD,L[i]); 
	);
);

if(#LD!=0, 				/* If the second list LD contains polynomials, we reduce their number by computing the isomorphism classes of the number fields generated by them. For every number field 						found, we keep just one generating polynomial */
	LF=List();
	listput(LF,LD[1]);
	for(i=2,#LD, 
		counter=0; 
		for(j=1,#LF, 
			if(nfisisom(LD[i],LF[j])!=0,
				counter++; 
				break;
			); 
		); 
		if(counter==0, listput(LF,LD[i]);
		);
	);
print("Number of isomorphism classes found = " #LF);
);
t4=getabstime();
print((t4-t3)/60000." minutes employeed for the second test");
return(LF);
}


