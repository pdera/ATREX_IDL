;// Routines for interaction with the main GUI which are used outside of
;// the main GUI event handler file

pro WID_MAR345_ext_commons

;--- PD change 07/29/2010

 @WID_GSE_ADA_COMMON

;--- PD change end

end

function use_compression
COMMON WID_MAR345_elements
  s=widget_info(WID_BUTTON_9b, /button_set) ;
 return, s
end

;-------------------------
function is_sequence_number_proper, fn
ON_IOERROR, erri
  a=strpos(fn, '_', /REVERSE_SEARCH)
  b=strpos(fn, '.', /REVERSE_SEARCH)
  lgt=strlen(fn)
  len=b-a
  numb=strmid(fn, lgt-a-2, len-1, /reverse_offset)
  ft='(I0'+strcompress(string(len-1),/remove_all)+')'
  numbi=string(long(numb), format=ft)
  if numb eq numbi then return, long(numb) else return, -1
  erri: return, -1
end
;-------------------------

function read_kappa_and_ttheta
COMMON WID_MAR345_elements
 widget_control, wid_text_36a, get_value=ka
 widget_control, wid_text_36b, get_value=tth
 return, [float(ka),float(tth)]
end

function read_compression
COMMON WID_MAR345_elements
 widget_control, wid_text_0a, get_value=ka
 return, float(ka)
end

pro set_kappa_and_ttheta, kt
COMMON WID_MAR345_elements
 widget_control, wid_text_36a, set_value=string(kt[0], format='(F10.5)')
 widget_control, wid_text_36b, set_value=string(kt[1], format='(F10.5)')
end

function read_spec_dat,fn
dn=fltarr(2)
if fn ne '' then $
begin
  a=0.0
  del=0.0
  nu=0.0
  aaa=''
  get_lun, lu
  openr, lu, fn
  i=0
  readf, lu, aaa
  while not eof(lu) do $
  begin
   readf, lu, a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,a,del,a,a,a,nu
   if i eq 0 then dn=[nu, del] else dn=[[dn],[nu, del]]
   i=i+1
  endwhile
  close, lu
  free_lun, lu
  print, del, nu
  endif
  return, dn
end


pro save_gonio, gonio, fn
if fn ne '' then $
begin
  get_lun, lu
  openw, lu, fn
  printf, lu, string(gonio[1], format='(F9.3)')+$
  string(gonio[3], format='(F9.3)')+$
  string(gonio[4], format='(F9.3)')+$
  string(gonio[5], format='(F9.3)')+$
  string(gonio[0], format='(F9.3)')+$
  string(gonio[2], format='(F9.3)')
  close, lu
  free_lun, lu
end
end

function read_gonio_from_file, fn
gonio=fltarr(6)
if fn ne '' then $
begin
  a=0.0
  get_lun, lu
  openr, lu, fn

  readf, lu, a1, a2, a3, a4, a5, a6
  gonio[1]=float(a1)
  gonio[3]=float(a2)
  gonio[4]=float(a3)
  gonio[5]=float(a4)
  gonio[0]=float(a5)
  gonio[2]=float(a6)
  close, lu
  free_lun, lu
end
return, gonio
end

function assign_cal
COMMON WID_MAR345_elements
 re0=widget_info(wid_button_20aa4, /button_set)
 return, re0
end

function read_gonio
COMMON WID_MAR345_elements
kp=read_kappa_and_ttheta()
widget_control, WID_TEXT_0z, get_value=a1
widget_control, WID_TEXT_1z, get_value=a2
widget_control, WID_TEXT_2z, get_value=a3
widget_control, WID_TEXT_3z, get_value=a4
return, [kp[0],kp[1],float(a1),float(a2),float(a3),float(a4)]
end

;---------------------

pro print_gonio, gonio
COMMON WID_MAR345_elements
set_kappa_and_ttheta, [gonio[0],gonio[1]]
widget_control, WID_TEXT_0z, set_value=string(gonio[2], format='(F9.3)')
widget_control, WID_TEXT_1z, set_value=string(gonio[3], format='(F9.3)')
widget_control, WID_TEXT_2z, set_value=string(gonio[4], format='(F9.3)')
widget_control, WID_TEXT_3z, set_value=string(gonio[5], format='(F9.3)')
end

;---------------------

pro update_hkl, hkl
COMMON WID_MAR345_elements
 widget_control, WID_LABEL_hkl_1, set_value=string(hkl[0], format='(F6.3)')
 widget_control, WID_LABEL_hkl_2, set_value=string(hkl[1], format='(F6.3)')
 widget_control, WID_LABEL_hkl_3, set_value=string(hkl[2], format='(F6.3)')
end

;---------------------

function assign_ub
COMMON WID_MAR345_elements
 re0=widget_info(wid_button_20aa5, /button_set)
 return, re0
end

function peak_filtering_settings
COMMON WID_MAR345_elements
  widget_control, aWID_TEXT_0, get_value=a0
  widget_control, aWID_TEXT_1, get_value=a1
  widget_control, aWID_TEXT_2, get_value=a2
  widget_control, aWID_TEXT_3, get_value=a3
  widget_control, aWID_TEXT_4, get_value=a4
  re0=widget_info(aWID_BUTTON_0, /button_set)
  re1=widget_info(aWID_BUTTON_1, /button_set)
  re2=widget_info(aWID_BUTTON_2, /button_set)
  re3=widget_info(aWID_BUTTON_3, /button_set)
  re4=widget_info(aWID_BUTTON_4, /button_set)
  return, [float(re0)*a0, float(re1)*a1, float(re2)*a2, float(re3)*a3, float(re4)*a4]

end

;---------------------------
function open_Cij, fname
  Cij=fltarr(6,6)
  get_lun, lu
  a1=0.0
  a2=0.0
  a3=0.0
  a4=0.0
  a5=0.0
  a6=0.0
  openr, lu, fname
  for i=0,5 do $
  begin
    readf, lu, a1, a2, a3, a4, a5, a6
    Cij[i,0:5]=[a1, a2, a3, a4, a5, a6]
  end
  Close, lu
  return, Cij
end
;---------------------------
function open_lp, fname
  lp=fltarr(6)
  get_lun, lu
  a1=0.0
  a2=0.0
  a3=0.0
  a4=0.0
  a5=0.0
  a6=0.0
  openr, lu, fname
  readf, lu, a1, a2, a3, a4, a5, a6
  Close, lu
  return, [a1, a2, a3, a4, a5, a6]
end
;---------------------------

pro print_Cij, cij, list
  text=strarr(6)
  for i=0,5 do $
    for j=0, 5 do text[i]=text[i]+string(Cij[i,j], format='(F10.4)')
  widget_control, list, set_value=text
end

;---------------------------

pro print_reference, lp, list
  text=strarr(4)
  for i=0,2 do $
    text[i]=text[i]+string(lp[i], format='(F10.4)')+string(lp[3+i], format='(F8.2)')
  V=v_from_lp(lp)
  text[3]=string(V, format='(F16.2)')
  widget_control, list, set_value=text
end

;---------------------------
pro print_strain, sig, list
  text=strarr(6)
  for i=0,2 do $
    for j=0, 2 do text[i]=text[i]+string(sig[i,j], format='(F10.6)')
  widget_control, list, set_value=text
end

;********************************** NEW CODE ****************
function read_om_rotation_dir
COMMON WID_MAR345_elements
   re=widget_info(WID_BUTTON_OMEGA_ROTATION, /button_set)
   ;---- PD change
   ;---- 6/23/2010
   ;---- changed default setting to GSECARS
   if re eq 1 then return, -1 else return, 1

   ;---- PD change end
end

;********************************** NEW CODE ****************

;-----------------------------

function gonio_zero_offset
COMMON WID_MAR345_elements
 widget_control, WID_TEXT_30aa, get_value=ka
 return, float(ka)
end

;-----------------------------

function read_image_settings, fname
COMMON WID_MAR345_elements
 fil=fname+'.txt'
 re=file_info(fil)
 if re.exists then $
 begin
   l1=''
   l2=''
   l3=''
   l4=''
   nums=fltarr(4)
   get_lun, un
   openr, un, fil
   readf, un, l1
   readf, un, l2
   readf, un, l3
   readf, un, l4
   close, un
   free_lun, un
   nums[0]=float(strmid(l1, 12,strlen(l1)-12))
   nums[1]=float(strmid(l2, 12,strlen(l2)-12))
   nums[2]=float(strmid(l3, 12,strlen(l3)-12))
   nums[3]=float(strmid(l4, 12,strlen(l4)-12))
   print, nums

   if read_om_rotation_dir() eq 1 then $
   begin
    WIDGET_CONTROL, WID_TEXT_30, SET_VALUE=string((nums[0]-gonio_zero_offset()), format='(F12.2)')
    WIDGET_CONTROL, WID_TEXT_31, SET_VALUE=string(nums[1], format='(F12.2)')
   end else $
   begin
    WIDGET_CONTROL, WID_TEXT_30, SET_VALUE=string( read_om_rotation_dir()*(nums[0]-gonio_zero_offset())-nums[1], format='(F12.2)')
    WIDGET_CONTROL, WID_TEXT_31, SET_VALUE=string(nums[1], format='(F12.2)')
   endelse
 end
end

;--------------------------------
function read_cal_type
COMMON WID_MAR345_elements
   return, widget_info(WID_BUTTON_refine_calt_p, /button_set)
end

;--------------------------------
function read_polar
COMMON WID_MAR345_elements
 widget_control, wid_text_14, get_value=uv
 return, float(uv)
end

;--------------------------------

pro update_progress, prog
COMMON WID_MAR345_elements
 ;widget_control, wid_text_8, set_value=string(prog*100.0, format='(I4)')+'%'
end



function read_overlap_limites
COMMON WID_MAR345_elements
 widget_control, wid_text_30r, get_value=ka
 widget_control, wid_text_31r, get_value=tth
 return, [float(ka),float(tth)]
end

;--------------------------------

function read_box_size, text1, text2
COMMON WID_MAR345_elements
 text1=wid_text_37
 text2=wid_text_37
 widget_control, text1, get_value=px
 widget_control, text2, get_value=py
 px=long(px)
 py=long(py)
 px=px[0]
 py=py[0]
 return, [px,py]
end


;--------------------------------

function hkl_labels_on
COMMON WID_MAR345_elements
 return, widget_info(wid_button_19, /button_set)
end

;--------------------------------

function whole_p_series
COMMON WID_MAR345_elements
 return, widget_info(WID_BUTTON_wholeseries, /button_set)
end



;-------------------------------

function omega_or_energy
COMMON WID_MAR345_elements
;re=widget_info(wid_button_7, /button_set)
re=widget_info(WID_DROPLIST_aab, /droplist_select)
return, re ; 1=om, 0=en
end

;-------------------------------

function update_on_browse
COMMON WID_MAR345_elements
re=widget_info(wid_button_23c, /button_set)
return, re ;
end

;-------------------------------

function fit_peaks_after_PS
COMMON WID_MAR345_elements
re=widget_info(wid_button_9a, /button_set)
return, re
end

;-------------------------------

function do_symmetric_correlation
 COMMON WID_MAR345_elements
 re=widget_info(wid_button_23b, /button_set)
 return, re
end

;-------------------------------

function I_corrections
COMMON WID_MAR345_elements
re1=widget_info(wid_button_51a, /button_set) ; dac abs
re2=widget_info(wid_button_51, /button_set) ; Lorenz
re3=widget_info(wid_button_52, /button_set) ; Polariz
return, [re1,re2, re3] ; 1=apply, 0=don't apply
end

;--------------------------------
function vert_polariz
COMMON WID_MAR345_elements
 s=widget_info(WID_BUTTON_52ab, /button_set) ;
 return, s
end

;===============================================================
;===============================================================
;===============================================================

pro peak_intensity_evaluation, pni, A, PCERROR, XY, axis, lbcgr

@COMMON_DATAS
@COMMON_DATAS2
;COMMON WID_MAR345_elements
;COMMON image_type_and_arrays
;common selections, selecting_status, select, excl, unselect, addpeak, mask, maskarr, zoom, unmask


       pt=opt->get_object()
       pt.peaks[pni].DetXY[0]= XY[0]+a[5]
       pt.peaks[pni].DetXY[1]= XY[1]+a[6]

       bs=[long(pt.peaks[pni].intssd[0]),long(pt.peaks[pni].intssd[0])]
       XX=Xinvec([bs[0]*2+1,bs[1]*2+1])
       yy=yinvec([bs[0]*2+1,bs[1]*2+1])
       p2=long(one2two(profile_function(XX, YY, A)))


       pic=oimage->get_zoomin(XY, bs, maskarr)

       ;---- if part of the box area is outside active circle, replace 0 with local background ------
       cc=which_pixels_in_mtx_outside_circle(bs[0], xy, arr1/2-2)
       c=where(cc eq 1)
       if c[0] ne -1 then pic[c]=replicate(lbcgr[xy[0],xy[1]], n_elements(c))

       ;---------------------------------------------------------------------------------------------


       attt=evaluate_fit_quality(pic, p2)
       if attt[0] gt 0 then $
       begin
          opt->select_peak, pni
       endif else $
       begin

       gonio=pt.peaks[pni].gonio
       kt=read_kappa_and_ttheta()
       gonio[1]=kt[1]
       xyz=oadetector->calculate_XYZ_from_pixels_mono(pt.peaks[pni].DetXY, gonio, wv)

     ;------------
       sd=oadetector->calculate_sd_from_pixels(pt.peaks[pni].DetXY, gonio)
       nu=ang_between_vecs([0.0,sd[1],sd[2]],[0.,-1.,0.])
       nu1=oadetector->get_nu_from_pix(pt.peaks[pni].DetXY)
       ;print, nu, nu1
       tth=get_tth_from_xyz(pt.peaks[pni].xyz)
       II=(A[1]*A[2]*A[3])*2.*!pi
       si=2.*!pi*(A[2]*A[3]*PCERROR[1]+A[1]*A[2]*PCERROR[3]+A[1]*A[3]*PCERROR[2])
       ty=I_corrections()

      ;-- Lorenz correction
      if ty[1] eq 1 then $
      begin
       ii=II/lorenz(xyz, gonio, wv)
       si=sI/lorenz(xyz, gonio, wv)
      endif

     ;-- polarization correction
      if ty[2] eq 1 then $
      begin
       P1=0.9
       mu=20.41
       h=0.005
       v=vert_polariz()*90.0+nu
       ii=II/(1.0-read_polar()*cos(2.0*v/!radeg));*cos(gonio[3]/!radeg)
       si=sI/(1.0-read_polar()*cos(2.0*v/!radeg));*cos(gonio[3]/!radeg)
      endif
      ;-- CBN absorption correction
      if ty[0] eq 1 then $
      begin
       om0= 0.0
       ii=II/dac_profile(pt.peaks[pni].gonio[axis]+om0,1)
       si=sI/dac_profile(pt.peaks[pni].gonio[axis]+om0,1)
      end

    ;-- omega polynomial correction
       aom=omega_correction(pt.peaks[pni].gonio[axis])
       ii=II*aom
       si=sI*aom

       pt.peaks[pni].energies[2]=A[2] ; width 2th
       pt.peaks[pni].energies[3]=A[3] ; width chi
       pt.peaks[pni].energies[4]=A[4] ; tilt
       pt.peaks[pni].energies[5]=A[0] ; background

       pt.peaks[pni].intAD[0]=ii
       pt.peaks[pni].intAD[1]=si
       opt->set_object, pt
       end

end
;----------------------------------



function cell_now_solution_n, n
@common_datas
fl=cell_now_dir+'cell_now.exe'
re=file_info(fl)
if not(re.exists) then begin
	Wid_cellnow_path_dlg
endif


fl=cell_now_dir+'cell_now.exe'
re=file_info(fl)
if re.exists eq 1 then $
begin


; write to last_directory.txt
save_last_directories
;dir='C:\Users\prze
dir=cell_now_dir
fil=dir+'cell_now_commands.txt
get_lun, ln
openw, ln, fil
printf, ln, 'xxx.p4p'
printf, ln, 'xxx._cn'
printf, ln, 'Enter'
printf, ln, '10'
printf, ln, '2 20'
printf, ln, '0.25'
a=string(n)
b=STRCOMPRESS(a, /REMOVE_ALL)
printf, ln,  b
printf, ln, '0.25'
printf, ln, 'A'
printf, ln, '1.p4p'
printf, ln, 'Q'
close, ln
free_lun, ln

fil = dir+'run_cellnow.cmd'
get_lun,ln
openw, ln, fil
printf, ln, '@ECHO OFF'
printf, ln, 'SETLOCAL EnableExtensions EnableDelayedExpansion'
printf, ln, 'REM script to automate indexing with cell_now'
printf, ln, 'REM 2015-02-14 pd'
printf, ln,  strmid(cell_now_dir,0,2)
printf, ln,  'cd "'+cell_now_dir+'"'
printf, ln, 'SET CELL="cell_now.exe"'
printf, ln, 'SET COMM="'+cell_now_dir+'cell_now_commands.txt"'
printf, ln, 'SET PROJ=xxx.p4p'
printf, ln, 'REM announce target'
printf, ln, 'REM announce target'
printf, ln, 'ECHO Processing %PROJ% ...'
printf, ln, 'ECHO %CELL%'
printf, ln, '%CELL% < %COMM%'
printf, ln, 'REM done'
printf, ln, 'ENDLOCAL'
printf, ln, 'ECHO ON'
close,ln
free_lun, ln
return, 0
endif else return, -1
end

function find_opx_cell, cells
  sz=size(cells)
  nm=sz[2]
  opx=-1
  cpx=-1
  for i=nm-1, 0, -1 do $
  begin
    if (abs(cells[3,i]-90.) lt 1.0) and (abs(cells[4,i]-90.) lt 1.0) and (abs(cells[5,i]-90.) lt 1.0) then opx=i
    if (abs(cells[3,i]-90.) lt 1.0) and (cells[4,i] gt 91.0) and (cells[4,i] lt 100.0) and (abs(cells[5,i]-90.) lt 1.0) then cpx=i
  endfor
  return, [opx,cpx]
end


;--------------------------------
function exclude_corners
COMMON WID_MAR345_elements
 re=WIDGET_INFO(WID_BUTTON_20aa3, /BUTTON_SET)
 return, re
end
;--------------------------------
function current_vs_series
COMMON WID_MAR345_elements
 re=WIDGET_INFO(WID_BUTTON_412b, /BUTTON_SET)
 return, re
end

;--------------------------------
pro read_cells_from_cellnow, lp, v, l, fom
	@COMMON_DATAS
   lp=fltarr(6)
   lp0=fltarr(6)
   V=0.0
   L=''
   V0=0.0
   L0=''
   fom= 0.0
   fom0 = 0.0
   ;dir='C:\Users\przemyslaw\Dropbox (UH Mineral Physics)\software\RSV_mSXD 2.5\'
   dir = cell_now_dir
   fil=dir+'xxx._cn'
   get_lun, ln
   openr, ln, fil
   re=''
   while strmid(re,0, 4) ne ' FOM' and not eof(ln) do $
   begin
    readf, ln, re
    ;print, re
   endwhile
    readf, ln, re
    ;print, re
   ;--- reading cell parameters sets
   readf, ln, re

   for i=0, 5 do $
   begin
     lp[i]=strmid(re, 18+i*8,8)
   endfor
     FOM=strmid(re, 11,5)
     V=strmid(re, 65,10)
     L=strmid(re, 75,3)
   num=1
   while re ne '' and not eof(ln) do $
   begin
    lp=[[lp],[lp0]]
    V=[V,V0]
    fom=[fom,fom0]
    L=[L,L0]
    readf, ln, re
    for i=0, 5 do lp[i, num]=strmid(re, 18+i*8,8)
    V[num]=strmid(re, 65,10)
    L[num]=strmid(re, 75,3)
    fom[num]=strmid(re,11,5)
    num=num+1
   endwhile
   close, ln
   free_lun, ln
    ;print, lp
   ;print, size(lp)
   lp=lp[*,0:num-2]
   l=l[0:num-2]
   v=v[0:num-2]
   fom=fom[0:num-1]
   print, lp
   print, v
   print, l
   print, fom

end

;---------------------------

function PS, om, axis, prog=prog
@COMMON_DATAS
@COMMON_DATAS2
COMMON WID_MAR345_elements
common status, PE_open, SIM_open
common calib_ref, zeroref

    opt->initialize
    read_ps_settings, a1,a2,a3,a4, a5

    ;a1, additive for gradient
    ;a2, max count
    ;a3, smooth window
    ;a4, min count
    ;a5, local background window

    capture_calibration, oadetector, wv

    widget_control, wid_text_6, get_value=ni
    widget_control, wid_text_7, get_value=i0
    widget_control, wid_text_4, get_value=om0
    widget_control, wid_text_5, get_value=omD

    i0=fix(i0)
    ni=fix(ni)
    om0=float(om0)
    omD=float(omD)
 ;   widget_control, wid_text_8, set_value=string(0, format='(I2)')+'/'+string(ni[0], format='(I2)')

    bx=read_box_size(wid_text_37,wid_text_37)

    im=oimage->get_object()
    lbcgr=oimage->calculate_local_background(0, a5)
    im1=smooth2(im.img-lbcgr, a3)

    if use_compression() eq 1 then $
    begin
     compression_factor=fix(read_compression())
     compression_factor=compression_factor[0]
     im2=congrid(im1, arr1/compression_factor, arr1/compression_factor)
     locmax, im2, a1, MASK=m,  IX=ix, IY=iy, VALUES=v, /NOEDGE
     print, 'number of local maxima found: ', n_elements(v)
     ix=ix*compression_factor
     iy=iy*compression_factor
    endif else $
    begin
     locmax, im1, a1, MASK=m,  IX=ix, IY=iy, VALUES=v, /NOEDGE
     print, 'number of local maxima found: ', n_elements(v)
    endelse


if n_elements(v) gt 0 then $
begin
    f=where(v gt a4)
    if f[0] ne -1 then $
    begin
     f1a=where(v[f] lt a2)
     if f1a[0] ne -1 then $
     begin
      gonio=fltarr(6)
      gonio[axis]=om

      if f1a[0] ne -1 then $
      for k=0, n_elements(f1a)-1 do $
      begin
        if oadetector->inside_detector_circle([IX[f[f1a[k]]],IY[f[f1a[k]]]], 12) then $
        begin
         ref_peak.detXY=[IX[f[f1a[k]]],IY[f[f1a[k]]]]
         ref_peak.intssd[0:1]=[bx[0],bx[0]]
         ref_peak.gonio=gonio
         opt->appendpeak, ref_peak
        endif
      end

      pt=opt->get_object()

      if  fit_peaks_after_PS() then fit_all_peaks, prog=prog

      update_peakno, opt->peakno()
  ;   opt->find_multiple_peak_copies
      opt->calculate_all_xyz_from_pix, oadetector, wv
      endif
     endif
    endif
    return, opt->peakno()
end
;--------------------------------------------------

function PS_Series
@COMMON_DATAS
@COMMON_DATAS2
COMMON WID_MAR345_elements
;COMMON image_type_and_arrays
;COMMON cont
;common selections
;common mcoordinates
;COMMON EX
;COMMON prof
;COMMON CLASS_peaktable_reference
;common status, PE_open, SIM_open
;common calib_ref, zeroref

    lp=fltarr(13)
    lpo=fltarr(13)
    lpc=fltarr(13)
    widget_control, wid_text_6, get_value=ni
    widget_control, wid_text_7, get_value=i0
    widget_control, wid_text_4, get_value=om0
    widget_control, wid_text_5, get_value=omD
    i0=fix(i0)
    ni=fix(ni)
    om0=float(om0)
    omD=float(omD)
   ; widget_control, wid_text_8, set_value=string(0, format='(I2)')+'/'+string(ni[0], format='(I2)')

    cgProgressBar = Obj_New("CGPROGRESSBAR", /Cancel, text='Progress in series')
    cgProgressBar -> Start


    for i=0, ni[0]-1 do $
    begin
     res.seq=i0[0]+i
     fn=generate_fname(res)
     res=analyse_fname(fn, dir, res.extno)
     widget_control, wid_text_9, set_value=res.name0
     oimage->load_image, fn, oadetector
     plot_image, oimage
     a=PS(om0+i*omD, 3, prog=0)
     plot_peaks, draw0, opt, arr1, arr2
     fn=dir+res.name0+'.pks'
     opt->write_object_to_file, fn, 1
     IF cgProgressBar -> CheckCancel() THEN $
     BEGIN
     	ok = Dialog_Message('Operation canceled')
     	cgProgressBar -> Destroy
     	RETURN, [1,1]
     ENDIF
     cgProgressBar -> Update, (float(i)/float(ni[0]))*100.0
   endfor
    cgProgressBar -> Destroy

    opt->find_multiple_peak_copies
    ;opt->recalculate_all_xyz, oadetector, wv
    update_peakno, opt->peakno()

    opt->calculate_all_xyz_from_pix, oadetector, wv
    plot_image, oimage
    plot_peaks, draw0, opt, arr1, arr2
    print_peak_list, opt, wid_list_3, oadetector, wv
    update_peakno, opt->peakno()



goto, ror

    opt->write_object_to_file, fn+'.pks'


    ;---------- determine UB matrix with

    cell_now_solution_n, 1
	dirs=cell_now_dir
    ;dirs='C:\Users\przemyslaw\Dropbox (UH Mineral Physics)\software\RSV_mSXD 2.5\'
    a=opt->save_p4p(dirs+'xxx.p4p')
    text='MYARG="'+dirs+'run_cellnow.cmd"'
    SETENV, text
    spawn, '%MYARG%'  , /LOG_OUTPUT


;    opx_cpx=find_opx_cell(read_cells_from_cellnow())
    if opx_cpx[0] gt 0 then $
    begin
      cell_now_solution_n, opx_cpx[0]+1
      spawn, '%MYARG%'  , /LOG_OUTPUT
    endif else $
    if opx_cpx[1] gt 0 then $
    begin
      cell_now_solution_n, opx_cpx[1]+1
      spawn, '%MYARG%'  , /LOG_OUTPUT
    endif
    pas:
    ;dirs='C:\Users\przemyslaw\Dropbox (UH Mineral Physics)\software\RSV_mSXD 2.5\'
	dirs=cell_now_dir
    ub=ReadUBfrom_p4p(dirs+'1.p4p')
    lp= lp_from_ub(UB)
    if lp[4] lt 89. or  lp[4] gt 91. then symm=12 else symm=2
    tol=0.15
    opt->select_indexable, ub, tol
    plot_peaks, draw0, opt, arr1, arr2
    opt->delete_selected
    lp=Refine_B_against_d(ub, opt, symm)
    opt->delete_selected
    lp=Refine_B_against_d(ub, opt, symm)
    opt->delete_selected
    if symm eq 2 then $
    lpo=Refine_B_against_d(ub, opt, symm) else $
    lpc=Refine_B_against_d(ub, opt, symm)
    plot_image, oimage
    plot_peaks, draw0, opt, arr1, arr2
    if symm eq 2 then $
    begin
      opt->write_object_to_file, fn+'_io.pks'
      if opx_cpx[1] gt -1 then $
      begin
        cell_now_solution_n, opx_cpx[1]+1
        spawn, '%MYARG%'  , /LOG_OUTPUT
        symm=12
        opt->select_indexable, ub, tol
    	plot_peaks, draw0, opt, arr1, arr2
    	opt->delete_selected
    	lp=Refine_B_against_d(ub, opt, symm)
    	opt->delete_selected
    	lp=Refine_B_against_d(ub, opt, symm)
    	opt->delete_selected
    	lpc=Refine_B_against_d(ub, opt, symm)
        opt->write_object_to_file, fn+'_ic.pks'
      endif
    endif else $
    begin
      opt->write_object_to_file, fn+'_ic.pks'
    endelse
    return, [[lpo],[lpc]]
    ror: return, [0,0]
end
;---------------------------

function read_cal_type
; 0 - fit2d style calibration
; 1 - Dera style calibration
COMMON WID_MAR345_elements
 rex=WIDGET_INFO(WID_BUTTON_refine_calt_s,/button_set)
return, rex
end
;---------------------------


function Refine_B_against_d, ub, optable1, sym

       ;sym=2 for orthorhombic

       ds=optable1->BUILD_d_list()
       hkls=optable1->BUILD_hkls()
       lp=lp_from_ub(ub)
       lp1=automatic_lp_refinement3(lp, ds, hkls, sym)
       lp=lp_from_ub(ub)
       ;change ub matrix
       b1=b_from_lp(lp1)
       b0=b_from_lp(lp)
       u0=ub ## invert(b0)
       ub=u0 ## b1
       lp1[0:5]=lp_from_ub(ub)
       print, 'refined cell parameters', lp1
       ddd=[-0.005, 0.005] ;$$$$$$$$$$$$$$$$$$$$$
       lp=lp_from_ub(ub)
       pt=optable1->get_object()
       N=pt.peakno
       ds=fltarr(N)
       dsc=fltarr(N)
       a=fltarr(N)
       ss=optable1->calculate_Ddd(lp)
       sel1=where(ss le ddd[0])
       sel2=where(ss ge ddd[1])
       if sel1[0] ne -1 then optable1->select_peaks, sel1
       if sel2[0] ne -1 then optable1->select_peaks, sel2
       print, n_elements(sel1)+n_elements(sel2), 'peaks selected'
       return, lp1

end

;---------------------------

function Refine_B_and_twist_against_d, ub, optable1, sym

       ;sym=2 for orthorhombic

       ds=optable1->BUILD_d_list()
       hkls=optable1->BUILD_hkls()

       lp=lp_from_ub(ub)
       lp1=automatic_lp_refinement_twist(lp, ds, hkls, sym)
       lp=lp_from_ub(ub)
       ;change ub matrix
       b1=b_from_lp(lp1)
       b0=b_from_lp(lp)
       u0=ub ## invert(b0)
       ub=u0 ## b1
       lp1[0:5]=lp_from_ub(ub)
       print, 'refined cell parameters', lp1
       ddd=[-0.005, 0.005] ;$$$$$$$$$$$$$$$$$$$$$
       lp=lp_from_ub(ub)
       pt=optable1->get_object()
       N=pt.peakno
       ds=fltarr(N)
       dsc=fltarr(N)
       a=fltarr(N)
       ss=optable1->calculate_Ddd(lp)
       sel1=where(ss le ddd[0])
       sel2=where(ss ge ddd[1])
       if sel1[0] ne -1 then optable1->select_peaks, sel1
       if sel2[0] ne -1 then optable1->select_peaks, sel2
       print, n_elements(sel1)+n_elements(sel2), 'peaks selected'
       return, lp1

end
;---------------------------
function ReadUBfrom_p4p, res
ub=fltarr(3,3)
	if res ne '' then $
	begin
 		FREE_LUN,2
 		OPENR, 2, res
 		str='       '
 		str1=''
 		str2=''
 		str3=''
 		; search for beginning of reflection block
 		while str ne 'CELLSD' and not eof(2) do $
   			readf, 2, str, format='(A6)'
 		if not eof(2) then $
 		begin
   			readf, 2, str1
			readf, 2, str2
			readf, 2, str3
     		CLOSE, 2
    		FREE_LUN,2
    		ub[0,0]=float(strmid(str1, 7,16))
    		ub[1,0]=float(strmid(str1, 23,16))
    		ub[2,0]=float(strmid(str1, 39,16))

    		ub[0,1]=float(strmid(str2, 7,16))
    		ub[1,1]=float(strmid(str2, 23,16))
    		ub[2,1]=float(strmid(str2, 39,16))

    		ub[0,2]=float(strmid(str3, 7,16))
    		ub[1,2]=float(strmid(str3, 23,16))
    		ub[2,2]=float(strmid(str3, 39,16))
 		endif
 		endif
 		return, ub
end
;--------------

function read_inversions
COMMON WID_MAR345_elements
 rex=WIDGET_INFO(WID_BUTTON_20aa,/button_set)
 rey=WIDGET_INFO(WID_BUTTON_20aa1,/button_set)
 tra=WIDGET_INFO(WID_BUTTON_20aa2,/button_set)
 res=0
 if rex eq 1 and rey eq 1 and tra eq 0 then res = 2 else $
 if rex eq 1 and rey eq 0 and tra eq 0 then res = 5 else $
 if rex eq 0 and rey eq 1 and tra eq 0 then res = 7 else $
 if rex eq 0 and rey eq 0 and tra eq 1 then res = 4 else $
 if rex eq 1 and rey eq 0 and tra eq 1 then res = 3 else $
 if rex eq 0 and rey eq 1 and tra eq 1 then res = 1 else $
 if rex eq 1 and rey eq 1 and tra eq 1 then res = 6
 return, res
end



;---------------------------



function get_dac_opening
COMMON WID_MAR345_elements
  WIDGET_CONTROL, WID_TEXT_opening,gET_VALUE=aa
  return, float(aa)
end


;------------------------------

function read_fcf_select_th

COMMON WID_MAR345_elements

 WIDGET_CONTROL, WID_text_fcf_1, GET_VALUE=a1
 WIDGET_CONTROL, WID_text_fcf_2, GET_VALUE=a2

 return, [float(a1),float(a2)]

end

function read_scale_select_th

COMMON WID_MAR345_elements

 WIDGET_CONTROL, WID_text_sel1, GET_VALUE=a1
 WIDGET_CONTROL, WID_text_sel2, GET_VALUE=a2

 return, [float(a1),float(a2)]

end
;---------------------------
function scale_spline

COMMON WID_MAR345_elements

 re=widget_info(WID_button_spline, /button_set)
 s=widget_info(WID_droplist_spline, /DROPLIST_SELECT) ;
 return, re*(s+5)
end

;----------------------

function fcf_spline

COMMON WID_MAR345_elements

 re=widget_info(WID_BUTTON_fcf_spline, /button_set)
 s=widget_info(WID_droplist_fcf_spline, /DROPLIST_SELECT) ;
 return, re*(s+5)
end

;---------------------------


function scale_poly
; reads a degree of polynomial to be used for scaling from the droplist in the scale tab
; available values are 2, 3, 4
COMMON WID_MAR345_elements
  s=widget_info(WID_droplist_poly, /DROPLIST_SELECT) ;
 return, s+2
end


pro print_scale_poly, coeff
@common_datas
 n=n_elements(coeff)
 poly_text=strarr(n)
 COMMON WID_MAR345_elements
 for i=0, n-1 do poly_text[i]=string(coeff[i], format='(F20.16)')
 widget_control, WID_list_poly, set_value=poly_text
end

function Lab6, dst, wh
@COMMON_DATAS

N=48
la=wv
;dst=200
d=fltarr(N)
tth=fltarr(N)
rad=fltarr(25)

d[0]=		4.1565	;	0	0	1
d[1]=		2.9391	;	0	1	1
d[2]=		2.3997	;	1	1	1
d[3]=		2.0782	;	0	0	2
d[4]=		1.8588	;	0	1	2
d[5]=		1.6969	;	1	1	2
d[6]=		1.4695	;	0	2	2
d[7]=		1.3855	;	0	0	3
d[8]=		1.3855	;	1	2	2
d[9]=		1.3144	;	0	1	3
d[10]=		1.2532	;	1	1	3
d[11]=		1.1999	;	2	2	2
d[12]=		1.1528	;	0	2	3
d[13]=		1.1109	;	1	2	3
d[14]=		1.0391	;	0	0	4
d[15]=		1.0081	;	0	1	4
d[16]=		1.0081	;	2	2	3
d[17]=		0.9797	;	1	1	4
d[18]=		0.9797	;	0	3	3
d[19]=		0.9536	;	1	3	3
d[20]=		0.9294	;	0	2	4
d[21]=		0.9070	;	1	2	4
d[22]=		0.8862	;	2	3	3
d[23]=		0.8484	;	2	2	4
d[24]=		0.8313	;	0	0	5
d[25]=		0.8313	;	0	3	4
d[26]=		0.8152	;	0	1	5
d[27]=		0.8152	;	1	3	4
d[28]=		0.7999	;	1	1	5
d[29]=		0.7999	;	3	3	3
d[30]=		0.7718	;	0	2	5
d[31]=		0.7718	;	2	3	4
d[32]=		0.7589	;	1	2	5
d[33]=		0.7348	;	0	4	4
d[34]=		0.7235	;	2	2	5
d[35]=		0.7235	;	1	4	4
d[36]=		0.7128	;	0	3	5
d[37]=		0.7128	;	3	3	4
d[38]=		0.7026	;	1	3	5
d[39]=		0.6927	;	0	0	6
d[40]=		0.6927	;	2	4	4
d[41]=		0.6833	;	0	1	6
d[42]=		0.6743	;	1	1	6
d[43]=		0.6743	;	2	3	5
d[44]=		0.6572	;	0	2	6
d[45]=		0.6491	;	1	2	6
d[46]=		0.6491	;	0	4	5
d[47]=		0.6491	;	3	4	4
;
for i=0, 47 do $
tth[i]=2.0*asin(la/(2.0*d[i]))*!radeg
rad=dst*tan(tth/!radeg)
if wh eq 0 then return, rad else return, d
end

function Neon, dst, wh
@COMMON_DATAS

N=11
la=wv
;dst=200
d=fltarr(11)
tth=fltarr(11)
rad=fltarr(11)



d[0] =	2.10560
d[1] =  1.82350
d[2] =  1.28941
d[3] = 	1.09961
d[4] =	1.05280
d[5] =	0.91175
d[6] =	0.83668
d[7] =	0.81549
d[8] =	0.74444
d[9] =	0.70187
d[10] =	0.64470

for i=0, 10 do $
tth[i]=2.0*asin(la/(2.0*d[i]))*!radeg
rad=dst*tan(tth/!radeg)
if wh eq 0 then return, rad else return, d
end

function CeO2, dst, wh
@COMMON_DATAS

N=25
la=wv
;dst=200
d=fltarr(25)
tth=fltarr(25)
rad=fltarr(25)

d[0] =3.12442
d[1] =2.70583
d[2] =1.91331
d[3] =1.63167
d[4] =1.56221
d[5] =1.35291
d[6] =1.24152
d[7] =1.21008
d[8] =1.10465
d[9] =1.04147
d[10]=0.95665
d[11]=0.91474
d[12]=0.90194
d[13]=0.85566
d[14]=0.82527
d[15]=0.81584
d[16]=0.78110
d[17]=0.75778
d[18]=0.75046
d[19]=0.72316
d[20]=0.70454
d[21]=0.67646
d[22]=0.66114
d[23]=0.65626
d[24]=0.63777
for i=0, 23 do $
tth[i]=2.0*asin(la/(2.0*d[i]))*!radeg
rad=dst*tan(tth/!radeg)
if wh eq 0 then return, rad else return, d
end

function CO2, dst, wh
@COMMON_DATAS

N=21
la=wv
;dst=200
d=fltarr(N)
tth=fltarr(N)
rad=fltarr(N)

d[0]=2.87261
d[1]=2.48775
d[2]=2.22511
d[3]=2.03124
d[4]=1.7591
d[5]=1.6585
d[6]=1.4363
d[7]=1.37996
d[8]=1.32976
d[9]=1.24388
d[10]=1.20674
d[11]=1.14146
d[12]=1.08574
d[13]=1.06078
d[14]=1.01562
d[15]=0.92393
d[16]=0.9084
d[17]=0.87955
d[18]=0.84101
d[19]=0.82925
d[20]=0.80713

for i=0, N-1 do $
tth[i]=2.0*asin(la/(2.0*d[i]))*!radeg
rad=dst*tan(tth/!radeg)
if wh eq 0 then return, rad else return, d
end

function which_calibrant, dst, no; choice of calibrant for powder based refinement of detector geometry
COMMON WID_MAR345_elements
  s=widget_info(WID_Droplist_36aa, /DROPLIST_SELECT) ;
  case s of
  0:a=ceo2(dst,no) ; CeO2
  1:a=LaB6(dst,no) ; LaB6
  2:a=Neon(dst,no) ; LaB6
  3:a=CO2(dst,no) ; LaB6
  endcase
 return, a
end

;WID_Droplist_36aa
function which_calibration
COMMON WID_MAR345_elements
  s=widget_info(WID_BUTTON_refine_cal_p, /button_set) ;
 return, s
end


function use_corners
COMMON WID_MAR345_elements
  s=widget_info(WID_BUTTON_33b, /button_set) ;
 return, s
end



function refine_twist
COMMON WID_MAR345_elements
  s=widget_info(WID_BUTTON_refine_twist, /button_set) ;
 return, s
end

function cal_IovS
COMMON WID_MAR345_elements
  WIDGET_CONTROL, WID_TEXT_36ac, Get_value=re
 return, float(re)
end



;--------------------------------
function dynamic_mask_on
COMMON WID_MAR345_elements
 s=widget_info(WID_BUTTON_dynamic_mask, /button_set) ;
 return, s
end

function prefitting
COMMON WID_MAR345_elements
 s=widget_info(WID_BUTTON_23a, /button_set) ;
 return, s
end




;-------------------------------------
function profile_function_name
COMMON WID_MAR345_elements
 s=widget_info(WID_button_413a, /button_set) ; pseudo-Voigt
 s1=widget_info(WID_button_412z, /button_set) ; Gauss1D
 if s eq 1 then return, 'Voigt2dwt' else $
 begin
 	if s1 eq 0 then return, 'Gaussian2dwt' else return, 'Gauss1D'
 endelse
end






;--------------------------------

function read_fcf_var
COMMON WID_MAR345_elements

 re=widget_info(WID_droplist_var1, /droplist_select)
 return, re+1
end

;--------------------------------
function get_laue_class


COMMON WID_MAR345_elements

 re=widget_info(WID_Droplist_4aa, /droplist_select)
 return, re
end

function which_PT
COMMON WID_MAR345_elements

 re=widget_info(WID_button_512, /button_set)
 return, re ; 1 for PT1 0 for PT2
end

function omega_correction, om
COMMON WID_MAR345_elements
;  re=widget_info(WID_BUTTON_51aomg,/button_set)
;  if re eq 0 then ret=[0,0,0] else $
;  begin
;   widget_control, WID_TEXT_14_da1om, get_value=a1
;   widget_control, WID_TEXT_14_da2om, get_value=a2
;   widget_control, WID_TEXT_14_da3om, get_value=a3
;   ret=[float(a1),float(a2),float(a3)]
;  endelse
;  b=ret[0]*om*om+ret[1]*om+ret[2]
  return, 1.0
end
;--------------------------------

pro print_R_int, text
COMMON WID_MAR345_elements

widget_control, WID_TEXT_37aa, set_value=string(text, format='(F7.4)')
widget_control, WID_text_rinta, set_value=string(text, format='(F7.4)')

end




function read_abs
COMMON WID_MAR345_elements
        widget_control, WID_TEXT_14_da1, get_value=va1
        widget_control, WID_TEXT_14_da2, get_value=va2
        widget_control, WID_TEXT_14_da3, get_value=va3
        widget_control, WID_TEXT_14_da4, get_value=va4
        widget_control, WID_TEXT_14_da5, get_value=va5
        widget_control, WID_TEXT_14_da6, get_value=va6
  return, [float(va1),float(va2),float(va3),float(va4),float(va5),float(va6)]
end
;--------------------------------
pro plot_image, oimage
COMMON WID_MAR345_elements
COMMON image_type_and_arrays, imt, arr1, arr2,cenx, ceny, rad, rad1
common selections, selecting_status, select, excl, unselect, addpeak, mask, maskarr, zoom, unmask

   im=oimage->get_object()
   data1=congrid(im.img[0:arr1-1,0:arr2-1]*maskarr, 600,600)
   widget_control, wid_text_24, get_value=imax
   widget_control, wid_text_25, get_value=imin
   imax=long(imax)
   imin=long(imin)
   imax=imax[0]
   imin=imin[0]
   data1=data1<imax
   data1=data1>imin

   wset, draw0
   device, decomposed=0
   tvscl, data1, true=0

        ;thisPostion = [0.1, 0.1, 0.9, 0.9]
        ;cgIMAGE, data1, POSITION=thisPosition, /KEEP_ASPECT_RATIO
        ;cgZImage, data1
  end
;--------------------------------
function get_md
COMMON WID_MAR345_elements
 widget_control, wid_text_39, get_value=mov
 return, reform(float(mov))
end
;--------------------------------
function get_alpha
COMMON WID_MAR345_elements
 widget_control, wid_text_36c, get_value=al
 return, reform(float(al))
end
;--------------------------------


function rotate_image
COMMON WID_MAR345_elements
 re=widget_info(WID_BUTTON_20aa,/button_set)
 return, re
end

function det_move
COMMON WID_MAR345_elements
 re=widget_info(wid_button_37,/button_set)
 return, re
end

;--------------------------------
pro plot_zoom, oimage
COMMON WID_MAR345_elements
COMMON image_type_and_arrays, imt, arr1, arr2,cenx, ceny, rad, rad1
common selections, selecting_status, select, excl, unselect, addpeak, mask, maskarr, zoom
   if cenx-rad le 0 then rad=cenx-1
   if ceny-rad1 le 0 then rad1=ceny-1
   if cenx+rad ge arr1 then rad=arr1-cenx-1
   if ceny+rad1 ge arr2 then rad1=arr2-ceny-1

   im=oimage->get_object()
   data1=congrid(im.img[cenx-rad:cenx+rad,ceny-rad1:ceny+rad1]*maskarr, 428,428)
   widget_control, wid_text_24, get_value=imax
   widget_control, wid_text_25, get_value=imin
   imax=long(imax)
   imin=long(imin)
   imax=imax[0]
   imin=imin[0]
   data1=data1<imax
   data1=data1>imin

   wset, draw3
   tvscl, data1, true=0
  end
;---------------------------
function read_box_change

; 0 is one peak
; 1 is all peaks

COMMON WID_MAR345_elements
 re=widget_info(wid_button_412, /button_set)
 if re eq 1 then return, 0 else return, 1
end

;--------------------------------

function read_color_scale
COMMON WID_MAR345_elements

   widget_control, wid_text_24, get_value=imax
   widget_control, wid_text_25, get_value=imin
   imax=long(imax)
   imin=long(imin)
   imax=imax[0]
   imin=imin[0]

   return, [imin, imax]
end

pro update_peakno, pn
COMMON WID_MAR345_elements
 widget_control, wid_text_10, set_value=string(pn, format='(I4)')
end
;-------------------------------
function omega_from_scan, filno
COMMON WID_MAR345_elements
   widget_control, wid_text_6, get_value=ni
   widget_control, wid_text_7, get_value=i0
   widget_control, wid_text_4, get_value=om0
   widget_control, wid_text_5, get_value=omD
   i0=fix(i0)
   ni=fix(ni)
   om0=float(om0)
   omD=float(omD)
   ene0=float(om0)
   eneD=float(omD)
   pos=filno-i0
   om=om0+pos*omD+omD/2.
   return, om
end

;-------------------------------

function im_seq_0
COMMON WID_MAR345_elements
   widget_control, wid_text_6, get_value=ni
   widget_control, wid_text_7, get_value=i0
   a=long(i0)
   return, a[0]
end

;-------------------------------

function im_seq_n
COMMON WID_MAR345_elements
   widget_control, wid_text_6, get_value=ni
   widget_control, wid_text_7, get_value=i0
   a=long(ni)
   return, a[0]
end


;********************************** NEW CODE ****************
pro WID_MAR345_Cleanup, WID_MAR345
@COMMON_DATAS
common closing, Wid_Image_simulation
common status, PE_open, SIM_open

if SIM_open eq 1 then wid_sim->destroy

end
;********************************** NEW CODE ****************
;-------------------------------
pro WID_MAR345_ext
end
;-------------------------------