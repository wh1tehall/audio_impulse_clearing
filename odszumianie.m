###BEST OUTCOME UP TO THIS POINT - BASE TO REVERT###
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
etxs=samples*0;
tt=samples*0;
mt=samples*0;
outsamples=samples*0;
etm=sqrt(samples(1)^2)
ets=sqrt(samples(1)^2);
mt=0.03;
outsamples=[zeros(1,lp)]
tt(p-2)=start;
tt(p-1)=time()-start;
lambda=0.95
alf=0.99;
 
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
 
for i=p:300
 
  y=samples(i);
 
  #estymacja:
  pom=outsamples(length(outsamples)-p+1:length(outsamples));
  x=pom;
  out=x*w';
  eta=y-out;
  et(i)=eta;
  etx(i)=etm;
  etxs(i)=mt;
  
  outsamples(i)=y;
   
    mt=alf*mt+(1-alf)*abs(eta);
    etm=(etm+sqrt(eta^2))/2
	ets=(ets+((mt-eta)^2))/2;
    #Obliczenia wzmocnienia RLS
      xp=x';
      Px=P*xp;
     
      k=(Px)/(1+x*Px);
      w=w+k'*eta;
      P=P-(Px*x*P)/(1+x*Px);
endfor
for i=300:samples_length
      y=samples(i);
 
	  #estymacja:
	  pom=outsamples(length(outsamples)-p+1:length(outsamples));
	  x=pom;
	  out=x*w';
	  eta=y-out;
	  et(i)=eta;
	  etx(i)=etm;
	  etxs(i)=mt;
	  
    if eta^2>(3*mt)^2 && i+3<samples_length #eta^2>(0.02)^2
      outsamples(i)= 0.75*(outsamples(i-1)+0.25*samples(i+3));
    else
   
	  mt=alf*mt+(1-alf)*abs(eta);
      etm=(etm+sqrt(eta^2))/2;
      ets=(ets+((mt-eta)^2))/2;
      outsamples(i)=y;
      #Obliczenia wzmocnienia RLS
      xp=x';
      Px=P*xp;
     
      k=(Px)/(lambda+x*Px);
      w=w+k'*eta;
      P=(P-(Px*x*P)/(1+x*Px))/lambda;
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
 print -dpdf tt.pdf
 figure()
 plot(etxs)
 title("sigma")
 print -dpdf ets.pdf
 
 wavwrite(outsamples,fs,'out.wav')
 en=time();
 (en-start)/60