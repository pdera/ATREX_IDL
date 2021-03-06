function dac_profile, angle, r

; Corrections made on 2/14/2011

; Calculates the absorption correction for the cBN upstream seat
; with a conical opening and the upstream diamond (diamond contribution is negligible)
; The function describes I/Io changes as a function of omega angle.
; Input angle in degrees


   ;angle=angle-0.8
   an=angle/!radeg ; angle in radians


   ;==== DAC geometry parameters =====

 ;  diam_height   =   2.50000
 ;  CBN_height    =   5.800
 ;  CBN_bigdiam   =   4.000
 ;  CBN_smalldiam =   0.705

   ;diam_height   =   2.2000 ; przemek hpcat


   diam_height   =   1.8000 ; greg1 hpcat
   CBN_height    =   5.600
   CBN_bigdiam   =   4.000
   CBN_smalldiam =   0.705
   ;==================================

   diam_mu =  0.00462948
   CBN_mu =   0.1050460

   if n_params() gt 1 then $
   begin
    a=read_abs()

   ani=angle+a[5]
   an=abs(ani)/!radeg ; angle in radians

   diam_height   =   a[0]
   CBN_height    =   a[1]
   CBN_bigdiam   =   a[2]
   CBN_smalldiam =   a[3]
   CBN_mu        =   a[4]

   end



   ;CBN_mu =   0.1263460;greg3
   ;CBN_mu =   0.1313460;greg1 and 2
;   CBN_mu =   0.1063460 ; przemek hpcat
;   CBN_mu =   0.0863460 ; gsecars

   if an eq 0 then l_diam=diam_height else l_diam=abs(diam_height/cos(an)) ; path length through diamond

   cbn_cangle1=atan(0.5*CBN_smalldiam/diam_height)*!radeg
   cbn_cangle2=atan(0.5*CBN_bigdiam/(diam_height+CBN_height))*!radeg
   alpha=atan(0.5*(CBN_bigdiam-CBN_smalldiam)/CBN_height)*!radeg ; opening angle of the cbn cone
   beta=180.0-alpha
   gamma=180.0-abs(ani)-beta
   dh=(CBN_height+diam_height)-0.5*CBN_bigdiam/tan(alpha/!radeg)


   if abs(ani) le cbn_cangle1 then $ ; no cbn absorption
   begin
     l_CBN=0.0
   endif else $
   if abs(ani) lt cbn_cangle2 then $  ; inside of the cone
   begin
     l4=dh*sin(beta/!radeg)/sin(gamma/!radeg) ; total path length
     l_CBN=l4-l_diam
   endif else $
   begin
     l_CBN=CBN_height/cos(abs(an/2.)) 		    	                 ; full length of cbn
     print, ani, '-- full length'
   endelse

   att_diam=exp(-diam_mu*l_diam)
   att_cbn=exp(-CBN_mu*l_CBN)
   att=att_diam*att_cbn
   return, att
end

pro calculate_dac_profile


  xdata=fltarr(1000)
  ydata=fltarr(1000)
  om0=-30.0
  step=60.0/1000.0
  for i=0, 1000-1 do $
  begin
    om=om0+i*step
    ydata[i]=dac_profile(om)
    xdata[i]=om
  end
  m=min(ydata)
  plot, xdata, ydata;, /ynozero

  aby=fltarr(100)
  abx=fltarr(100)
  i=0
  a=0.0
  b=0.0

  ;openr, 2, 'c:\ASCII_test.txt'
  ;while not eof(2) do $
  ;begin
  ;  readf, 2, a, b
  ;  abx[i]=a
  ;  aby[i]=b
  ;  i=i+1
  ;endwhile
  ;close, 2
 ;  plot, abx[0:i-1]+89.5, aby[0:i-1]/max(aby[0:i-1]), color='00FF00'XL, psym=2
end


function read_profile
 ;re=dialog_pickfile(/read, filter='*.txt')
 re='C:\Documents and Settings\Przemek Dera\My Documents\Dropbox\manuscripts\calcite\Feb_2012\scans\ASCII_PD_calcitea.txt'
 openr, 2, re
 a=1.0
 b=1.0
 tab=[a,b]
 count=0
 tabs=tab
 while not(eof(2)) do $
 begin
    readf, 2, a, b
    tab[0]=a
    tab[1]=b
    if count eq 0 then tabs=tab else tabs=[[tabs],[tab]]
    count=count+1
 endwhile
 close, 2
 return, tabs
end

pro test_abs
 a=read_profile()
 a[0,*]=a[0,*]+90.
 window
 device, decomposed=1
 plot, a[0,*],a[1,*]
 s=n_elements(a)/2
 b=fltarr(s)
 for i=0, s-1 do $
 begin
   b[i]=dac_profile(a[0,i])*max(a[1,*])
 endfor
 oplot, a[0,*], b[*], color="FF00FF"XL


end

function generate_abs, omin, omax
 if omax lt omin then $
 begin
  a=omax
  omax=omin
  omin=a
 endif
 om=omin
 orange=omax-omin
 ostep=orange/181.0
 b=fltarr(2,181)
 for i=0, 180 do $
 begin
   b[0,i]=omin+i*ostep
   b[1,i]=dac_profile(b[0,i],1)
 endfor
 return, b
end