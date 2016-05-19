function read_tds, fn
  if fn ne '' then $
  begin
   px=0
   py=0
   h=0.0
   k=0.0
   l=0.0
   Int=0.0
   i=0
   get_lun, lu
   openr, lu, fn
   data=fltarr(6)
   while not eof(lu) do $
   begin
     readf, lu, px,py, h, k, l, Int
     if i eq 0 then data=[px,py, h, k, l, Int] else $
     data=[[data],[px,py, h, k, l, Int]]
     i=i+1
   end
   close, lu
   free_lun, lu
  endif
  return, data
end


pro TDS
  ;C=[384.9,	246.7, 279.2, 82.1, 93.8, 96.5, 100.1, 99.2, 104.3] ; Olivine 8GPa
  ;C=[384.9,	246.7, 279.2, 82.1, 93.8, 96.5, 100.1, 99.2, 104.3] ; Olivine 8GPa
  ;---------------
  ;cubic
  C=[190.84, 190.84, 190.84, 80.0, 80.0, 80.0, 93.96, 93.96, 93.96] ; Silicon
  ;---------------

  ;C=[384.9,	279.2,246.7, 82.1, 93.8, 96.5, 100.1, 99.2, 104.3] ; Olivine 8GPa

  lp=   [  5.3584,   5.3584,   5.3584]
  ;lp=   [  9.9226, 5.8540,4.6833]

  fi=dialog_pickfile(/read, filter='*.txt', path='J:\Dropbox (UH Mineral Physics)\TDS\')
  data=read_tds(fi)
  NN=n_elements(data)/6
  Na=sqrt(NN)

  b1=fltarr(Na,Na)
  b2=fltarr(Na,Na)
  minX=min(data[0,*])
  minY=min(data[1,*])

  Q=round(data[2:4,NN/2])
  ;Q=[ -4.0,  0.0, 0.0]


  for i=0, NN-1 do $
  begin
     qp=data[2:4,i]
     p=qp-Q
     qp1=qp/lp
     p1=p/lp
     A1=(HydroDyn(C,p1))
     b1[data[0,i]-minX,data[1,i]-minY]=transpose(Qp1)#invert(A1)#(Qp1)
     b2[data[0,i]-minX,data[1,i]-minY]=data[5,i]
  endfor
  b3=(b2-median(b2))>0
  c1=congrid(b1, 400,400)<max(b1)/1.5
  c2=congrid(b3, 400,400)<max(b3)/10.

  window,1,  xsize=400, ysize=400, xpos=0, ypos=0
  ;tvscl, c1
  contour, c1, /isotropic;, LEVELS = [1.0e-4, 5.0e-4, 1.0e-3, 5.0e-3, 1.0e-2]
  window,2,  xsize=400, ysize=400, xpos=0, ypos=410
  contour, c2, /isotropic;, LEVELS = [1.0e-4, 5.0e-4, 1.0e-3, 5.0e-3, 1.0e-2]
  ;tvscl, (c2)

end

;-------------------------
function HydroDyn, C, p
; cubic case
C11=C[0]
C22=C[1]
C33=C[2]
C44=C[3]
C55=C[4]
C66=C[5]
C12=C[6]
C13=C[7]
C23=C[8]

lenp=1;(p[0]*p[0]+p[1]*p[1]+p[2]*p[2])

n1=p[0]/lenp
n2=p[1]/lenp
n3=p[2]/lenp

A=fltarr(3,3)

; modified for orthorhombic

A[0,0]=C11*n1^2 +C66*n2^2 + C55*n3^2
A[0,1]=(C12+C66)*n1*n2
A[0,2]=(C13+C55)*n3*n1

A[1,0]=A[0,1]
A[1,1]=C22*n2^2+C44*n3^2+C66*n1^2
A[1,2]=(C23+C44)*n2*n3

A[2,0]=A[0,2]
A[2,1]=A[1,2]
A[2,2]=C33*n3^2+C44*n2^2 + C55*n1^2

return, A

end
;-------------------------
function HydroDynG, C, p ; orthorhombic case
pp=fltarr(3,3)
for i=0, 2 do $
 for j=0, 2 do $
  pp[i,j]=p[i]*p[j]
  pv=[p[0,0],p[1,1],p[2,2],p[1,2],p[0,2],p[1,2]] ; Voigt notation p

C=[[C[0],C[6],C[7],0,    0,   0   ],$
   [C[6],C[1],C[8],0,    0,   0   ],$
   [C[7],C[8],C[2],0,    0,   0   ],$
   [0,   0,   0,   C[3], 0,   0   ],$
   [0,   0,   0,   0,    C[4],0   ],$
   [0,   0,   0,   0,    0,   C[5]]]


end