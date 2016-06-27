start=time();
[samples,fs] = wavread("03.wav");
samples_length = length(samples);
 
p=4;
lp=4; #długosć pamięci dla estymatora
x=zeros(1,p); #wpółczynniki estymowanej funkcji
w=zeros(1,p); #regresory
yrp=zeros(lp); #pamięć dla RLS
yr=0; #wyjście r
yp=zeros(lp);
et=samples*0;
etx=samples*0;
tt=samples*0;
outsamples=samples*0;
etm=sqrt(samples(1)^2)
outsamples=[zeros(1,lp)]
tt(p-2)=start;
tt(p-1)=time()-start;
lambda=0.7
 
P=[];
for i=1:p
  Pp=[];
  for j=1:p
    if (i==j)
      Pp=[Pp 1];
    else
      Pp=[Pp 0];
    endif
  endfor
  P=[P;Pp];
endfor
999999*P

i=p
while i<samples_length #=p:samples_length
 
  y=samples(i);
 
  #estymacja:
  pom=outsamples(length(outsamples)-p+1:length(outsamples));
  x=pom;
  out=x*w';
  eta=y-out;
  et(i)=sqrt(eta^2);
  etx(i)=etm;
  if (i<100)
    outsamples(i)=y;
   start=time();
[samples,fs] = wavread("03.wav");
samples_length = length(samples);
 
p=4;
lp=4; #długosć pamięci dla estymatora
x=zeros(1,p); #wpółczynniki estymowanej funkcji
w=zeros(1,p); #regresory
yrp=zeros(lp); #pamięć dla RLS
yr=0; #wyjście r
yp=zeros(lp);
et=samples*0;
etx=samples*0;
tt=samples*0;
outsamples=samples*0;
etm=abs(samples(1))
outsamples=[zeros(1,lp)]
tt(p-2)=start;
tt(p-1)=time()-start;
lambda=0.7
 
P=[];
for i=1:p
  Pp=[];
  for j=1:p
    if (i==j)
      Pp=[Pp 1];
    else
      Pp=[Pp 0];
    endif
  endfor
  P=[P;Pp];
endfor
999999*P
 
for i=p:samples_length
 
  y=samples(i);
 
  #estymacja:
  pom=outsamples(length(outsamples)-p+1:length(outsamples));
  x=pom;
  out=x*w';
  eta=y-out;
  et(i)=abs(eta);
  etx(i)=etm;
  if (i<100)
    outsamples(i)=y;
   
    etm=(etm+abs(eta))/2
    #Obliczenia wzmocnienia RLS
      xp=x';
      Px=P*xp;
     
      k=(Px)/(1+x*Px);
      w=w+k'*eta;
      P=P-(Px*x*P)/(1+x*Px);
  else
    if abs(eta)>0.06
      outsamples(i)= samples(i-1);
    else
   
      etm=(etm+abs(eta))/2;
     
      outsamples(i)=y;
      #Obliczenia wzmocnienia RLS
      xp=x';
      Px=P*xp;
     
      k=(Px)/(lambda+x*Px);
      w=w+k'*eta;
      P=(P-(Px*x*P)/(1+x*Px))/lambda;
    endif
  endif
  tt(i)=time()-tt(i-1);
 endfor
 plot(et);
 title("eta")
 print -dpdf eta.pdf
 figure();
 plot(samples);
 title("oryginal");
 print -dpdf org.pdf
 figure()
 plot(outsamples);
 title("Wynik");
 print -dpdf out.pdf
 figure();
 plot(etx);
 print -dpdf etasr.pdf
 title("blad sredni");
 figure()
 plot(tt)
 title("czas wykonania")
 #wavwrite(outsamples,fs,'out.wav')
 en=time();
 en-start
    etm=(etm+sqrt(eta^2))/2;
    #Obliczenia wzmocnienia RLS
      xp=x';
      Px=P*xp;
     
      k=(Px)/(1+x*Px);
      w=w+k'*eta;
      P=P-(Px*x*P)/(1+x*Px);
  else
    if sqrt(eta^2)>3*etm
		% if i+5>samples_length
			% outsamples(i-1);
		% else
			% outsamples(i)= (samples(i-1)+samples(i+5))/6;
		% endif
		etq=sqrt(eta^2)
		j=1
		while (etq>3*etm && j<5)
			pomq=outsamples(length(outsamples)-p+1:length(outsamples));
			x=pomq;
			out=x*w';
			etq=sqrt((y-out)^2);
			%outsamples(i+j)=outsamples(i-1);
			j=j+1;
		end
		disp(length(outsamples))
		iint=(outsamples(i-1)+outsamples(i+j-1))/2;
		for il=1:j-1
			outsamples(i+il)=iint;
		endfor
		i=i+j-1;
    else
   
      etm=(etm+sqrt(eta^2))/2;
     
      outsamples(i)=y;
      #Obliczenia wzmocnienia RLS
      xp=x';
      Px=P*xp;
     
      k=(Px)/(lambda+x*Px);
      w=w+k'*eta;
      P=(P-(Px*x*P)/(1+x*Px))/lambda;
    endif
  endif
  #tt(i)=time()-tt(i-1);
  i=i+1;
 end
 plot(et);
 title("eta")
 print -dpdf eta.pdf
 figure();
 plot(samples);
 title("oryginal");
 print -dpdf org.pdf
 figure()
 plot(outsamples);
 title("Wynik");
 print -dpdf out.pdf
 figure();
 plot(etx);
 print -dpdf etasr.pdf
 title("blad sredni");
 figure()
 plot(tt)
 title("czas wykonania")
 wavwrite(outsamples,fs,'out.wav')
 who
 en=time();
 en-start