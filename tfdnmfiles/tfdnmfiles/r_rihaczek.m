function [tfr,t,f] = r_rihaczek(x,t,N,trace);

if (nargin == 0),
 error('At least 1 parameter required');
end;

[xrow,xcol] = size(x);
if (nargin == 1),
 t=1:xrow; N=xrow; trace=0;
elseif (nargin == 2),
 N=xrow; trace=0;
elseif (nargin == 3),
 trace = 0;
end;

if (N<0),
 error('N must be greater than zero');
end;
[trow,tcol] = size(t);
if (xcol==0) | (xcol>2),
 error('X must have one or two columns');
elseif (trow~=1),
 error('T must only have one row'); 
elseif (2^nextpow2(N)~=N & nargin==4),
 fprintf('For a faster computation, N should be a power of two\n');
end; 

tfr= zeros (N,tcol) ;  
if trace, disp('Rihaczek distribution'); end;
for icol=1:tcol,
 ti= t(icol); tau=-min([N-ti,xrow-ti]):(ti-1);
 indiceser= rem(N+tau,N)+1; indiceser=round(indiceser);
 if trace, disprog(icol,tcol,10); end;
 tfr(indiceser,icol)=x(ti,1)*conj(x(ti-tau,xcol));
end; 
if trace, fprintf('\n'); end;
tfr= fft(tfr); 

if (nargout==0),
 fprintf('if your goal is a graphical output, you should rather use tfrmh.\n');
 tfrqview(real(tfr),x,t,'tfrmh');
elseif (nargout==3),
 if rem(N,2)==0, 
  f=[0:N/2-1 -N/2:-1]'/N;
 else
  f=[0:(N-1)/2 -(N-1)/2:-1]'/N;  
 end;
end;

