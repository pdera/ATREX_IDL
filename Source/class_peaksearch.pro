
;--------------------------------------------------

function peaksearch_CLASS::get_object

COMMON CLASS_peaksearch_reference, pse

   pse.threshold = self.threshold
   pse.pbox      = self.pbox
   pse.bbox      = self.bbox
   pse.mindist   = self.mindist
   pse.peaktable = self.peaktable

   return, pse

end

;--------------------------------------------------

pro peaksearch_CLASS::set_object, pse

   self.threshold = pse.threshold
   self.pbox      = pse.pbox
   self.bbox      = pse.bbox
   self.mindist   = pse.mindist
   self.peaktable = pse.peaktable

end
;--------------------------------
function use_exclusions, butt
 re=widget_info(butt, /button_set)
 return, re
end

;---------------------------------------------------------------
FUNCTION EXCLUS, I, J,oadetector
COMMON EX, EXCLUSIONS
SZ=N_ELEMENTS(EXCLUSIONS)/2
 tth=oadetector->calculate_tth_from_pixels([I,J], [0,0,0,0,0,0])
 EX=0
 FOR K=0, SZ-1 DO IF ABS(TTH-EXCLUSIONS[0,K]) LT EXCLUSIONS[1,K] THEN RETURN, 1
 RETURN, 0
END


;---------------------------------------------------------

function estimate_local_background, nbinx, nbiny, img, thr, perc

s=size(img)
n0x=s[1]
n0y=s[2]

out_img=lonarr(n0x,n0y)
print, n0x, n0y
sx=fix(n0x/nbinx)
sy=fix(n0y/nbiny)

for i=0, nbinx-1 do $
begin
  for j=0, nbiny-1 do $
  begin
    x1=i*sx
    if i lt nbinx-1 then x2=(i+1)*sx else x2=n0x-1
    y1=j*sy
    if j lt nbiny-1 then y2=(j+1)*sy else y2=n0y-1
    box=img[x1:x2,y1:y2]
    ;h=histogram(box, locations=xs, min=thr, max=50000, nbins=1000)
    ;m=max(h,kk)
    ;out_img[x1:x2,y1:y2]=xs[kk]
    w=where(box gt 0)
    if w[0] ne -1 then out_img[x1:x2,y1:y2]=median(box[w])*perc else out_img[x1:x2,y1:y2]=thr*perc
  endfor ; j
endfor ; i
return, out_img
end


;---------------------------------------------------------

pro paint_rectangle, x0,x1,y0,y1, npixx, npixy
COMMON WID_MAR345_elements
 wset, draw0
 device, decomposed=1
 plots, float(x0)/float(npixx),float(y0)/float(npixy), /normal, color='FF0000'XL
 for i=y0, y1-1 do $
 begin
  plots, float(x1)/float(npixx),float(i)/float(npixy), /normal, color='FF0000'XL, /continue, linestyle=0
  plots, float(x0)/float(npixx),float(i+1)/float(npixy), /normal, color='FF0000'XL
 endfor
end

;---------------------------------------------------------

pro paint_peaks, tab, npixx, npixy
COMMON WID_MAR345_elements
wset, draw0
 device, decomposed=1
 plots, 0.,0.,/normal
 s=n_elements(tab)
 for i=0L, s-1 do $
 begin
   y=fix(tab[i]/float(npixy))
   x=fix(tab[i]-y*float(npixy))
   plots, float(x)/float(npixx),float(y)/float(npixy), /normal, color='FF00FF'XL, psym=3
  ; print, float(x)/float(npixx),float(y)/float(npixy)
 end
end




;========================================================================================
;========================================================================================



pro peaksearch_CLASS::execute, oimage, oad, excl, butt,  SHOW_PROGRESS_BAR=show_progress_bar
 COMMON CLASS_peaktable_reference, ref_peaktable, ref_peak
 pt=obj_new('CLASS_peaktable')
 pt->set_object, self.peaktable
 pt->initialize
 ptc=pt->get_object()
 ;----------- control parameters
 thr=self.threshold                       ; 100:       raw counts threshold for locating peaks
 max_peak_size=self.mindist               ; 10:        max allowed peak size with pixels above local background + Imin
 num_of_segments = [self.pbox,self.pbox]  ; [50.,50.]: number of segments in X and Y for local labckground estimation
 perc=self.bbox                           ; 1.0:       percent of median for background
 ;------------------------------
 t0=systime(/seconds)
 im=oimage->get_object()
 topX=long(im.sts.adet.nopixx/oimage.sts.binning)-1
 topY=long(im.sts.adet.nopixy/oimage.sts.binning)-1
 img1=congrid(im.img, 1000,1000)
 ff=where(img1 gt 10)
 ;h=histogram(img1, locations=xs, min=20, max=50000, nbins=1000)
 ;m=max(h,kk)
 ;print, 'histogram max: ', xs[kk]
;print, 'median: ', median(img1[ff])
; window
; plot, xs, h, xrange=[0,1000]

if ff[0] ne -1 then $
begin
 bg=estimate_local_background(num_of_segments[0], num_of_segments[1], img1,  median(img1[ff]), 1.0)
 w=where(img1-bg gt 100.)
 paint_peaks, w, 1000, 1000
 w0=n_elements(w)
 if show_progress_bar then $
 begin
 cgProgressBar = Obj_New("CGPROGRESSBAR", /Cancel)
 cgProgressBar -> Start
 endif
 oo:
 w=where(img1-bg gt thr)
 if w[0] ne -1 then $
 begin
    ;s=sort(img1[w])
    XY=ARRAY_INDICES([1000,1000], w[0], /DIMENSIONS)
    aa=grow_peak(img1, bg, xy[0],xy[1], thr, 1000, 1000)
    if max([aa[1]-aa[0], aa[3]-aa[2]]) lt 10 then $
    begin
     ref_peak.DetXY=[aa[4]*topX/1000.,aa[5]*topY/1000.]
     ref_peak.IntAD[0]=img1[aa[4],aa[5]]
     ref_peak.gonio=im.sts.gonio
     ref_peak.Selected[1]=1
     ref_peak.Selected[0]=0
     pt->Appendpeak, ref_peak
    endif
    ;plots, float(aa[4])/float(1000),float(aa[5])/float(1000), /normal, color='00FF00'XL, psym=1
    img1[aa[0]:aa[1],aa[2]:aa[3]]=0

    if show_progress_bar then $
    begin
     if cgProgressBar -> CheckCancel() THEN BEGIN
     ok = Dialog_Message('The user cancelled operation.')
     cgProgressBar -> Destroy
    RETURN
    endif
    ENDIF
    if show_progress_bar then cgProgressBar -> Update, 100.0-(n_elements(w)/float(w0))*100.0

    goto,oo
    endif
 end
 ptc=pt->get_object()
 self.peaktable=ptc
 obj_destroy, pt
 if show_progress_bar then $
 cgProgressBar -> Destroy
 ;print, 'Initial number of points:',w0
 print, 'Computation time: ',systime(/seconds)-t0
end



;---------------------------------------------------------------
pro CLASS_peaksearch

COMMON CLASS_peaktable_reference, ref_peaktable, ref_peak
COMMON CLASS_peaksearch_reference, pse

pse={peaksearch_CLASS, $
   threshold : 0.0,$
   pbox      : 0,$
   bbox      : 0,$
   mindist   : 0L,$
   peaktable : ref_peaktable}

end
