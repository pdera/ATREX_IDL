

; Last updated 05/02/2016 by PD
;	created list of methods
;   removed non class subroutines
;   removed all GUI-dependent calls
;   resolved all pre-definition problems

;   subroutines called with extra parameters:
;function adetector_class::tilt_mtx, cal
;function adetector_class::calculate_sd_from_pixels, pix, gonio, caltype
;function adetector_class::calculate_pixels_from_sd, sd, gonio, caltype
;function adetector_class::calculate_tth_from_pixels, pix, gonio, caltype
;function adetector_class::pixel_visible, XY, DACopening, XYZ_DAC, caltype
;pro adetector_class::visible_mask, mask, gonio, DACopening,wv, caltype
;function adetector_class::generate_one_peak, ubc, wv, pred, hkl, kt, caltype
;pro adetector_class::generate_all_peaks, ub, optable, wv, pred, exti, DAC_open, box, kt, rbc, caltype

;==============================================================
;---------------- Required libraries:
;   Crystallography
;   Vector_Math

;---------------- List of class methods:
; function adetector_class::tilt_mtx, caltype
; function adetector_class::calculate_sd_from_pixels, pix, gonio, caltype
; function adetector_class::calculate_pixels_from_sd, sd, gonio, caltype
; function adetector_class::calculate_tth_from_pixels, pix, gonio, caltype
; function adetector_class::pixel_visible, XY, DACopening, XYZ_DAC, caltype
; pro adetector_class::visible_mask, mask, gonio, DACopening,wv, caltype
; function adetector_class::generate_one_peak, ubc, wv, pred, hkl, kt, caltype
; pro adetector_class::generate_all_peaks, ub, optable, wv, pred, exti, DAC_open, box, kt, rbc, caltype
; pro adetector_class::generate_all_peaks, ub, optable, wv, pred, exti, DAC_open, box
; pro adetector_class::generate_peaks_laue, ub, optable, pred, en, exti, DAC_open
; function adetector_class::polarization, pix, p
; function adetector_class::get_nu_from_xyz, pix
; function adetector_class::calculate_pixels_from_xyz, xyz, gonio, caltype
; function adetector_class::calculate_XYZ_from_pixels_mono, pix, gonio, wv, caltype
; function adetector_class::generate_ds_ring, tth, kt, caltype
; function adetector_class::calculate_psi_angles, gonio, pix
; function adetector_class::get_nu_from_pix, pix
; function adetector_class::create_chi_bin_array,nbins, npo, chimax
; function adetector_class::create_tth_bin_array,nbins, npo, tthmax, caltype
; pro adetector_class::change_twist, newv
; pro adetector_class::set_values
; pro adetector_class::set_object, ad
; function adetector_class::get_object
; pro adetector_class::write_object_to_disk, filename
; pro adetector_class::read_object_from_disk, filename

;==============================================================
;-----------------------------------
;+
;NAME:
;   tilt_mtx
;
; PURPOSE:
;   This class method calculates a tilt matrix for the detector, which rotates vector 1 0 0 to
;   coincide with the detector normal. Two different definitions for the second tilt angle are available.
;   selectrion between these two definitions is dome using value of property self.self.ttheta0.
;   self.self.ttheta0=0 corresponds to fit2d
;
; CATEGORY:
;   Detector geometry
;
; CALLING SEQUENCE:
;   Result=object->tilt_mtx()
;
; INPUTS:
;   No formal parameters are accepted.
;   Second tilt angle definition is controlled by property self.self.ttheta0
;   self.self.ttheta0 = 0 is consistent witg fit2d
;
; OUTPUTS:
;   This function returns 3x3 dblarr rotation matrix
;
;
; COMMON BLOCKS:
;   Rota - used for rotation matrix algebra. Defined in vector_math
;
; MODIFICATION HISTORY:
;   May 19, 2016 Automated Documentation capability added
;-

function adetector_class::tilt_mtx
; using GUI interface to read cal type right now

     common Rota, Mtx
     caltype=self.ttheta0
     GenerateR, 3, self.tiltch
     omt=Mtx
     case caltype of
     0: begin
           GenerateR, 1, self.tiltom ; fit2d calibration, rotation defined about X axis
           cht=Mtx
     	   icht=invert(cht)
     	   tiltmtx=icht##omt##cht
        end
     1: begin
     	   GenerateR, 2, self.tiltom ; Dera calibration, rotation defined about Y axis
           cht=Mtx
     	   tiltmtx=omt##cht
     	end
	 else:
	 endcase
     tiltmtx1=dblarr(3,3)
     for i=0, 2 do for j=0, 2 do tiltmtx1[i,j]=double(tiltmtx[i,j])
     return, tiltmtx1

end

;--------------------------------------------------------------
;+
;NAME:
;   inside_detector_circle
;
; PURPOSE:
;	This function checks is a pixel with given coordinates is inside a detector center circle. It is assumed that the detector is square.
;
; CATEGORY:
;   Detector geometry
;
; CALLING SEQUENCE:
;   Result = object->inside_detector_circle(Pix, Marg)
;
; INPUTS:
;   Pix:    XY pixel coordinates
;   Marg:   Marigin size in pixels
;
; OUTPUTS:
;   This function returns 1 if pixel (Pix[0], Pix[1]) is inside center circle, and 0 is it is outside.
;
; COMMON BLOCKS:
;
; MODIFICATION HISTORY:
;   May 19, 2016 Automated Documentation capability added
;-

function adetector_class::inside_detector_circle, pix, marg
  XY=pix-[self.nopixx/2, self.nopixx/2]
  radius=min([self.nopixx/2, self.nopixy/2])
  R=sqrt(XY[0]^2+XY[1]^2)
  if R lt radius-marg then return, 1 else return, 0
end

;--------------------------------------------------------------
;+
;NAME:
;   calculate_sd_from_pixels
;
; PURPOSE:
;   This function calculates coordinates of the diffracted beam unit vector from pixel coorginates at given goniometer position.
;
; CATEGORY:
;   Detector geometry
;
; CALLING SEQUENCE:
;   Result = object->calculate_sd_from_pixels(Pix, Gonio)
;
; INPUTS:
;   Pix:    XY pixel coordinates
;   Gonio:  Goniometer position (6 angles, two of which desctibe detector position)
;      Gonio[0] = del
;      Gonio[1] = nu
;
; OUTPUTS:
;   This function returns coordinates of the diffracted beam unit vector
;
; COMMON BLOCKS:
;   Rota - used for rotation matrix algebra. Defined in vector_math
;
; MODIFICATION HISTORY:
;   May 19, 2016 Automated Documentation capability added
;-

function adetector_class::calculate_sd_from_pixels, pix, gonio

; Calculates the coordinates of the diffracted beam
; unit vectors from pixel coordinates

     caltype=self.ttheta0

     common Rota, Mtx

     xpix=pix[0]
     ypix=pix[1]

     vec1=[0.0,0.0,0.0]
     vec2=vec1
     vec3=vec1

    ; apply detector twist
     va=fltarr(3)
     va[1]=xpix-self.nopixx/2.0
     va[2]=ypix-self.nopixy/2.0

     GenerateR, 1,-self.angle
     chiA=Mtx
     vb=chiA##va

     vb[1]=vb[1]+self.nopixx/2.0
     vb[2]=vb[2]+self.nopixy/2.0

     ;--- calclate relative coordinares

     xrel=-(vb[1]-self.beamX)*self.psizeX
     yrel=(vb[2]-self.beamY)*self.psizeY

     ;--- create sd at 2theta 0
     vec1a=[0.0,0.0,0.0]
     vec1a[0]=0.0
     vec1a[1]=xrel
     vec1a[2]=yrel

     ;GenerateR, 1,-self.angle
     ;chiA=Mtx
     ;Vec1=chiA##vec1A

     Vec1=vec1A

    ;---- apply detector tilts

     mtx=self->tilt_mtx() ; Tilts are implemented here!!

     vec2=mtx##vec1

     ;---------------------

     vec2a=[0.0,0.0,0.0]

     vec2a[0]=self.dist
     vec2a[1]=0
     vec2a[2]=0
     vec3=vec2+vec2a

     ;--- include detector roatations ---
     GenerateR, 3,gonio[1]
     tth=Mtx
     vec3=tth ## (vec3)

     GenerateR, 2,gonio[0]
     nu=Mtx
     vec3=nu ## (vec3)

     ;--- normalize ----

     sd=transpose(vec3/vlength(vec3))
     return, sd

end

;----------------------------------------------------------------------
;+
;NAME:
;   calculate_pixels_from_sd
;
; PURPOSE:
;   This function calculates coordinates of the diffracted beam unit vector from pixel coorginates at given goniometer position.
;
; CATEGORY:
;   Detector geometry
;
; CALLING SEQUENCE:
;   Result = object->calculate_pixels_from_sd(Sd, Gonio)
;
; INPUTS:
;   Pix:    XY pixel coordinates
;   Gonio:  Goniometer position (6 angles, two of which desctibe detector position)
;      Gonio[0] = del
;      Gonio[1] = nu
;
; OUTPUTS:
;   This function returns coordinates of the diffracted beam unit vector
;
; COMMON BLOCKS:
;   Rota - used for rotation matrix algebra. Defined in vector_math
;
; MODIFICATION HISTORY:
;   May 19, 2016 Automated Documentation capability added
;-
function adetector_class::calculate_pixels_from_sd, sd, gonio

; calculates pixel coordinates from diffracted beam unit vector coordinates
; modified for GSECARS 11/30/2005
     caltype=self.ttheta0

     common Rota, Mtx

     Pi=acos(-1.0)
     vec1=[1.0,0.0,0.0]
     vec2=vec1
     vec3=vec1
     sd2=vec1
     pix=[0.0,0.0]


     sd2=sd/vlength(sd)

;--- include detector roatations ---

     GenerateR, 2,-gonio[0]
     nu=Mtx
     sd2=nu ## sd2

     GenerateR, 3,-gonio[1]
     tth=Mtx
     sd1=tth ## sd2

     ; calculate Bragg angle
     tth=ang_between_vecs(sd1,vec3)

     ;calculate detector normal vector
    ;---- apply detector tilts
     det=[-1.0,0.0,0.0]

     mtx=self->tilt_mtx() ; Tilts are applied here
     det1=mtx##det

     ;---------------------
     ; calculate point of intersection
     dv=line_plane_intersection(sd1, det1, [0.0,0.0,0.0], [self.dist,0.0,0.0])

     ;dvi=icht##iomt##cht##(dv+det*self.dist)
     dvi=invert(mtx)##(dv+det*self.dist)


     dv=dvi

     pix[0]=-dv[1]/self.psizeX+self.beamX
     pix[1]= dv[2]/self.psizeY+self.beamY

     ;---------------------- apply twist

     va=fltarr(3)
     va[1]=pix[0]-self.nopixx/2.0
     va[2]=pix[1]-self.nopixy/2.0


     GenerateR, 1,self.angle
     chiA=Mtx

;     dv=chiA##dvi
     vb=chiA##va
     pix=[vb[1]+self.nopixx/2.0,vb[2]+self.nopixy/2.0]

     ;pix[0]=-dv[1]/self.psizeX+self.beamX
     ;pix[1]= dv[2]/self.psizeY+self.beamY

     return, pix

end

;--------------------------------------------------------------------

function adetector_class::calculate_tth_from_pixels, pix, gonio

    caltype=self.ttheta0

    sd=self->calculate_sd_from_pixels(pix, gonio)
    s0=[1,0,0]
    return, ang_between_vecs(sd,s0)
end
;--------------------------------------------------------------------


function adetector_class::pixel_visible, XY, DACopening, XYZ_DAC
; currently this is done for zero goniometer position
; detector angles are determined in calculate_sd_from_pixels
; calculate_sd_from_pixels should be defined first
; XYZ_DAC is vector along dac axis
 caltype=self.ttheta0
 g=fltarr(6)
 Diffxyz=self->calculate_sd_from_pixels(XY, g)
 ang=ang_between_vecs(XYZ_DAC, Diffxyz)
 if abs(ang) le DACopening then return, 1 else return, 0

end

;----------------------------------------

pro adetector_class::visible_mask, mask, gonio, DACopening,wv
; generates a mask with 1 for visible and 0 for invisible, taking into account DAC opening

 common Rota, Mtx

 caltype=self.ttheta0

 GenerateR, 3, gonio[3]
 XYZ_DAC=Mtx # [1.,0.,0.]
 for i=0, self.nopixx-1 do $
  for j=0, self.nopixy-1 do $
  begin
     if mask[i,j] eq 1 then $
      if self->pixel_visible([i,j], DACopening,  XYZ_DAC) eq 0 then mask[i,j]=0
  endfor
end

;------------------------------------

function adetector_class::generate_one_peak, ubc, wv, pred, hkl, kt

; generates all peaks from crystal represented by orientation matrix UB
; ubc is the UB matrix
; pred is structure that carries prediction settings
; pred={  om_start : 0.0, $
;         om_range : 0.0, $
;         chi       : 0.0, $
;         d         : 0.0, $
;         h1        : 0, $
;         h2        : 0, $
;         k1        : 0, $
;         k2        : 0, $
;         l1        : 0, $
;         l2        : 0}
caltype=self.ttheta0

COMMON CLASS_peaktable_reference, ref_peaktable, ref_peak
 ;kt=read_kappa_and_ttheta()
 gonio=fltarr(6)
 gonio[1]=kt[1]
 gonio[4]=pred.chi
        ;hkl=[h,k,l]
        xyz=UBc ## hkl
        en=A_to_kev(wv)
        om=get_omega(en, xyz)
        ;om=get_omega_nonorthog(A_to_kev(wv), xyz, self.alpha)
        if om[0] ge pred.om_start and om[0] le pred.om_start+pred.om_range then $
        begin
          gonio[3]=om[0]
          pix=self->calculate_pixels_from_xyz(xyz, gonio)
          r=(pix-[self.beamx,self.beamy])
          r=sqrt(r[0]^2+r[1]^2)
          psi2=self->calculate_psi_angles(gonio, pix)
          if pix[0] gt 0 and pix[0] lt self.nopixx-1 then $
          begin
             ref_peak.detxy=pix
             ref_peak.gonio=gonio
             ref_peak.xyz=xyz
             ref_peak.hkl=hkl
          endif
        endif else if om[1] ge pred.om_start and om[1] le pred.om_start+pred.om_range then $
        begin
          gonio[3]=om[1]
          pix=self->calculate_pixels_from_xyz(xyz, gonio)
          r=(pix-[self.beamx,self.beamy])
          r=sqrt(r[0]^2+r[1]^2)
          psi2=self->calculate_psi_angles(gonio, pix)
          if pix[0] gt 0 and pix[0] lt self.nopixx-1  then $
          begin
             ref_peak.detxy=pix
             ref_peak.gonio=gonio
             ref_peak.xyz=xyz
             ref_peak.hkl=hkl
          endif
        endif
        return, ref_peak
end

;--------------------------------------------------------------------

pro adetector_class::generate_all_peaks, ub, optable, wv, pred, exti, DAC_open, box, kt, rbc
 ; generates all peaks from crystal represented by orientation matrix UB
 ; it is not using generate_one_peak

   caltype=self.ttheta0

COMMON CLASS_peaktable_reference, ref_peaktable, ref_peak
 ;kt=read_kappa_and_ttheta()
 gonio=fltarr(6)
 gonio[1]=kt[1] ; 2theta/Del
 gonio[0]=kt[0] ; nu
; gonio[4]=pred.chi ; this should be kappa/chi

 ;a=read_box_change()
 a=rbc
 ref_peak.intssd[0:1]=[a,a]

 for h=pred.h1, pred.h2 do $
 begin
   for k=pred.k1, pred.k2 do $
   begin
     for l=pred.l1, pred.l2 do $
     begin
        hkl=[h,k,l]
        if syst_extinction(exti, hkl) eq 1 then $
        begin
        xyz=UB ## hkl
        om=get_omega(A_to_kev(wv), xyz)

        if om[0]+360. ge pred.om_start and om[0]+360. le pred.om_start+pred.om_range then om[0]=om[0]+360.0
        if om[0]-360. ge pred.om_start and om[0]-360. le pred.om_start+pred.om_range then om[0]=om[0]-360.0

        if om[1]+360. ge pred.om_start and om[1]+360. le pred.om_start+pred.om_range then om[1]=om[1]+360.0
        if om[1]-360. ge pred.om_start and om[1]-360. le pred.om_start+pred.om_range then om[1]=om[1]-360.0

        ;om=get_omega_nonorthog(A_to_kev(wv), xyz, self.alpha)
        if om[0] ge pred.om_start and om[0] le pred.om_start+pred.om_range then $
        begin
          gonio[3]=om[0]
          pix=self->calculate_pixels_from_xyz(xyz, gonio)
          r=(pix-[self.beamx,self.beamy])
          r=sqrt(r[0]^2+r[1]^2)
          psi2=self->calculate_psi_angles(gonio, pix)
          if pix[0] gt 0 and pix[0] lt self.nopixx-1 and pix[1] gt 0 and pix[1] lt self.nopixy-1 and  (abs(psi2[1]) lt DAC_open) then $
          begin
             ref_peak.detxy=pix
             ref_peak.gonio=gonio
             ref_peak.xyz=xyz
             ref_peak.hkl=hkl
             ref_peak.intssd[0:1]=[box[0],box[0]]
             optable->appendpeak,ref_peak
          endif
        endif else if om[1] ge pred.om_start and om[1] le pred.om_start+pred.om_range then $
        begin
          gonio[3]=om[1]
          pix=self->calculate_pixels_from_xyz(xyz, gonio)
          r=(pix-[self.beamx,self.beamy])
          r=sqrt(r[0]^2+r[1]^2)
          psi2=self->calculate_psi_angles(gonio, pix)
          if pix[0] gt 0 and pix[0] lt self.nopixx-1 and (abs(psi2[1]) lt DAC_open) then $
          begin
             ref_peak.detxy=pix
             ref_peak.gonio=gonio
             ref_peak.xyz=xyz
             ref_peak.hkl=hkl
             ref_peak.intssd[0:1]=[box[0],box[0]]
             optable->appendpeak,ref_peak
          endif
        endif
        endif
     endfor
   endfor
 endfor


end

;--------------------------------------------------------------------

pro adetector_class::generate_peaks_laue, ub, optable, pred, en, exti, DAC_open
 ; generates all peaks from crystal represented by orientation matrix UB

COMMON CLASS_peaktable_reference, ref_peaktable, ref_peak

 gonio=fltarr(6)

 for h=pred.h1, pred.h2 do $
 begin
   for k=pred.k1, pred.k2 do $
   begin
     for l=pred.l1, pred.l2 do $
     begin
        if not (h eq 0 and k eq 0 and l eq 0) then $
        begin
        hkl=[h,k,l]
        if syst_extinction(exti, hkl) eq 1 then $
        begin
        xyz=UB ## hkl
        if vlength(xyz) ne 0 then $
        begin
        Ene=en_from_xyz(xyz)
        ;om=get_omega(A_to_kev(wv), xyz)
        if ene ge en[0] and  ene le en[1] then $
        begin
          pix=self->calculate_pixels_from_xyz(xyz, gonio)
          r=(pix-[self.beamx,self.beamy])
          r=sqrt(r[0]^2+r[1]^2)

          psi2=self->calculate_psi_angles(gonio, pix)
          if pix[0] gt 0 and pix[0] lt self.nopixx-1 and $
             pix[1] gt 0 and pix[1] lt self.nopixy-1 and $
              (abs(psi2[1]) lt DAC_open) then $
          begin
             ;print, hkl, ene, pix
             ref_peak.detxy=pix
             ref_peak.gonio=gonio
             ref_peak.xyz=xyz
             ref_peak.energies[0]=ene
             ref_peak.hkl=hkl
             optable->appendpeak,ref_peak
          endif
          endif
        endif
        endif
        endif
     endfor
   endfor
 endfor


end

;--------------------------------------------------------------

function adetector_class::polarization, pix, p

     xrel=-(pix[0]-self.beamX)*self.psizeX
     yrel= (pix[1]-self.beamY)*self.psizeY
     Ih=1.0
     Iv=(1.-p)/(1.+p)
     pol=(ih*abs(xrel)+iv*abs(yrel))/sqrt(xrel^2+yrel^2)
     return, pol
end


;--------------------------------------------------------------------

function adetector_class::get_nu_from_xyz, pix
; calculates the detector out of horiz. plane angle for a reflection from the Cartesian
; coordinates of reciprocal vector

  sd=[0.,pix[0],pix[1]]-[0.,self.beamx,self.beamx]
  nu=ang_between_vecs(sd,[0.,-1.,0.])
  if sd[2] lt 0.0 then return, nu else return, -nu
end

;--------------------------------------------------------------------

function adetector_class::calculate_pixels_from_xyz, xyz, gonio

; Calculates pixel coordinates from rec. vector coordinates, allows to apply further
; goniometer rotation gonio
; it is assumed that the goniometer rotation is relative to the original position of xyz
; this is used for recalc detxy
; angle gonio[0] is used to describe alpha - deviation of the omega axis from the vertical direction
; for now this is for an omega axis
     caltype=self.ttheta0

     sd=[0.0,0.0,0.0]
     pix=[0.0,0.0]
     xyz1=xyz/vlength(xyz)
     common Rota, Mtx

     GenerateR, 2, gonio[2]
     mu=Mtx

     GenerateR, 3, gonio[3]
     om=Mtx

     GenerateR, 1, gonio[4]
     ch=Mtx

     GenerateR, 3, gonio[5]
     ph=Mtx


     imu=invert(mu)
     iom=invert(om)
     ich=invert(ch)
     iph=invert(ph)

     hl=imu ## iom ## ich ## iph ## xyz1 ; already unit length

     sd=calculate_sd_from_hl(hl)

     pix=self->calculate_pixels_from_sd(sd, gonio)

     return, pix

 end

;--------------------------------------------------------------------

function adetector_class::calculate_XYZ_from_pixels_mono, pix, gonio, wv
; Calculates recip. versor coordinates from pixel coordinates
; (including non zero goinio position)
; Back rotation is applied, so the final vector coordinates are at zero
; goniometer position
; this is used for xyz from pixels

     common Rota, Mtx

     vec1=[0.0,0.0,0.0]
     vec2=vec1
     sd=vec1
     hl=vec1
     gonio0=fltarr(6)
     gonio0[0:1]=gonio[0:1]
     ;----------------------------
     sd=self->calculate_sd_from_pixels(pix, gonio0)
     ; sd has unit length
     ; it uses detector angles from gonio

     tth=self->calculate_tth_from_pixels(pix, gonio0)
     d=d_from_tth_and_en(tth, A_to_kev(wv))
     hl=(calculate_hl_from_sd(sd, wv))

     ; hl already has a proper length

     GenerateR, 2,-gonio[2]
     mu=Mtx
     GenerateR, 3,-gonio[3]
     om=Mtx
     GenerateR, 1,-gonio[4]
     ch=Mtx
     GenerateR, 3,-gonio[5]
     ph=Mtx

     vec2=ph # ch # om # mu # hl

     return, vec2

 end

;--------------------------------------------------------------------

function adetector_class::generate_ds_ring, tth, kt
; is using GUI right now

  caltype=self.ttheta0

  common Rota, Mtx

  vec=[1.,0.,0.]
  GenerateR, 3, tth
  sd=Mtx ## vec
  ;an=read_kappa_and_ttheta()
  an=kt
  ;gonio=[an[0],an[1],0.,0.,an[0],0.] ;
  gonio=[an[0],an[1],0.,0.,0.,0.] ; using nu and del
  pix=self->calculate_pixels_from_sd(sd, gonio)
  path=fltarr(360, 2)
  path[0,0:1]=pix
  for i=1, 359 do $
  begin
   GenerateR, 1, i*1.0
   vec2=Mtx ## sd
   pix=self->calculate_pixels_from_sd(vec2, gonio)
   path[i,0:1]=pix
  end
  return, path

end

;--------------------------------------------------------------------

function adetector_class::calculate_psi_angles, gonio, pix
; psi1 is the angle between incident beam and DAC axis,
; psi2 is the angle between diffracted beam and DAC axis

     caltype=self.ttheta0

     common Rota, Mtx
g=fltarr(6)
g[0:1]=gonio[0:1]
nn=n_elements(pix)/2
y=fltarr(2,nn)
for i=0,nn-1 do $
begin
 g[0:1]=gonio[0:1,i]
 sd=self->calculate_sd_from_pixels(pix[0:1,i], g)

 if gonio[3,i] lt -90 then $
 begin
    psi1=180.0+gonio[3,i]
    GenerateR, 3, -psi1
 endif else $
 if gonio[3,i] lt 90 then $
 begin
    psi1=gonio[3,i]
    GenerateR, 3, -psi1
 endif else $
 begin
    psi1=gonio[3,i]-180.0
    GenerateR, 3, -psi1
  endelse
 e1=[1.,0.,0.]
 e1=Mtx ## e1
 psi2=min([ang_between_vecs(sd,e1),ang_between_vecs(sd,-1.0*e1)])
 ;if psi2 gt 90. then psi2=psi2-90.
 y[0:1,i]=[psi1,psi2]
endfor
return, y
end

;------------------------------------------------------------------

function adetector_class::get_nu_from_pix, pix

;--- New class method replacing the old get_nu_from_pix
; calculates the detector out of horiz. plane angle for a reflection from the Cartesian
; coordinates of reciprocal vector

     caltype=self.ttheta0

  gonio=fltarr(6)
  nn=n_elements(pix)/2
  y=fltarr(nn)
  for i=0, nn-1 do $
  begin
   sd=self->calculate_sd_from_pixels(pix[0:1,i], gonio)
   sd[0]=0
   nu=ang_between_vecs(sd,[0.,-1.,0.])
   if sd[2] lt 0.0 then y[i]=nu else y[i]=-nu
 endfor
 return, y
end

;--------------------------------------------------------------------

function adetector_class::create_tth_bin_array,nbins, npo, tthmax
; generates 2theta bin array for integrating a powder pattern

     caltype=self.ttheta0

  binarr=fltarr(self.nopixx,self.nopixy)
  gonio=fltarr(6)
  spcx=fltarr(nbins)
  npo=lonarr(nbins)

  t1=self->calculate_tth_from_pixels([0,0], gonio)
  t2=self->calculate_tth_from_pixels([self.nopixx-1,0], gonio)
  t3=self->calculate_tth_from_pixels([0,self.nopixy-1], gonio)
  t4=self->calculate_tth_from_pixels([self.nopixx-1,self.nopixy-1], gonio)
  tthmax=max([t1,t2,t3,t4])
  dtth=tthmax/(nbins-2)
  for i=0, nbins-1 do spcx[i]=i*dtth

  for i=0, self.nopixx-1 do begin $
  print, i
  for j=0, self.nopixy-1 do $
  begin
       tth=self->calculate_tth_from_pixels([i,j], gonio)
       bin=long((tth)/dtth)+1
       binarr[i,j]=bin
       npo[bin]=npo[bin]+1
  endfor
  endfor
  return, binarr
end

;--------------------------------------------

function adetector_class::create_chi_bin_array,nbins, npo, chimax
; generates chi bin array for integrating a powder pattern
; get nu_from_pix should be replaced with class method

  binarr=fltarr(self.nopixx,self.nopixy)
  gonio=fltarr(6)
  spcx=fltarr(nbins)
  npo=lonarr(nbins)
  pix0=[self.beamX,self.beamY]
  chimax=180.
  chimin=-180.
  dtth=(chimax-chimin)/(nbins-2)
  for i=0, nbins-1 do spcx[i]=i*dtth

  for i=0, self.nopixx-1 do begin $
  print, i
  for j=0, self.nopixy-1 do $
  begin
       tth=get_nu_from_pix([i,j],pix0)
       bin=long((tth-chimin)/dtth)+1
       binarr[i,j]=bin
       npo[bin]=npo[bin]+1
  endfor
  endfor
  return, binarr
end

;----------------------------------------------------------

pro adetector_class::change_twist, newv
 self.angle  = newv
end

;--------------------------------------------------------------------

pro adetector_class::change_caltype, caltype
 self.ttheta0  = caltype
end

;--------------------------------------------------------------------

pro adetector_class::set_values

COMMON ADPV, dist,beamx,beamy,psizex,psizey,nopixx,nopixY,angle,omega0,ttheta0,tiltom,tiltch

  self.dist   = dist
  self.beamx  = beamx
  self.beamy  = beamy
  self.omega0  = omega0
  self.ttheta0 = ttheta0
  self.psizex = psizex
  self.psizey = psizey
  self.nopixx = nopixx
  self.nopixy = nopixy
  self.angle  = angle
  self.alpha  = alpha
  self.tiltom = tiltom
  self.tiltch = tiltch
  self.nu = nu
  self.del = del

end

;--------------------------------------------------------------------

pro adetector_class::set_object, ad

  self.dist   = ad.dist
  self.beamx  = ad.beamx
  self.beamy  = ad.beamy
  self.omega0 = ad.omega0
  self.ttheta0 = ad.ttheta0
  self.psizex = ad.psizex
  self.psizey = ad.psizey
  self.nopixx = ad.nopixx
  self.nopixy = ad.nopixy
  self.angle  = ad.angle
  self.alpha  = ad.alpha
  self.tiltom = ad.tiltom
  self.tiltch = ad.tiltch
  self.nu = ad.nu
  self.del = ad.del

end

;--------------------------------------------------------------------

function adetector_class::get_object

COMMON class_adetector_reference, ad

  ad.dist   = self.dist
  ad.beamx  = self.beamx
  ad.beamy  = self.beamy
  ad.omega0 = self.omega0
  ad.ttheta0 = self.ttheta0
  ad.psizex = self.psizex
  ad.psizey = self.psizey
  ad.nopixx = self.nopixx
  ad.nopixy = self.nopixy
  ad.angle  = self.angle
  ad.alpha  = self.alpha
  ad.tiltom = self.tiltom
  ad.tiltch = self.tiltch
  ad.nu = self.nu
  ad.del = self.del

  return, ad

end

;----------------------------------

pro adetector_class::write_object_to_disk, filename

COMMON class_adetector_reference, ad

     print, 'Writing detector parameters to:', filename

     ad=self->get_object()

     FREE_LUN, 2
     OPENW, 2, filename
     writeu, 2, ad
     CLOSE, 2
     FREE_LUN, 2

end

;----------------------------------

pro adetector_class::read_object_from_disk, filename

COMMON class_adetector_reference, ad

    print, 'Reading detector parameters from:', filename

     FREE_LUN, 2
     OPENR, 2, filename
     readu, 2, ad
     CLOSE, 2
     FREE_LUN, 2
     self->set_object, ad

end

;----------------------------------

pro class_adetector

COMMON class_adetector_reference, ad

ad={adetector_class, $
  dist   : 0.0,  $
  beamx  : 0.0,  $
  beamy  : 0.0,  $
  psizex : 0.0,     $
  psizey : 0.0,     $
  nopixx : 0,       $
  nopixY : 0,       $
  angle  : 0.0,     $
  alpha  : 0.0,     $
  omega0 : 0.0,     $
  ttheta0: 0.0,     $ ; this is calibration type
  tiltom : 0.0,     $
  tiltch : 0.0,	    $
  nu     : 0.0,     $
  del    : 0.0}

end