pro test_AD
COMMON CLASS_peaktable_reference, ref_peaktable, ref_peak
COMMON CLASS_Area_detector_parameters_reference, ref_adp
 crystallography
 CLASS_Area_detector_parameters
 Class_adetector
 oa=obj_new('adetector_class')
 print, '-- object instance created'

  ref=oa->get_object()

  ref.dist=200
  ref.beamx=50
  ref.beamy=50
  ref.psizex=0.1
  ref.psizey=0.1
  ref.nopixx=100
  ref.nopixy=100
  ref.tiltom=0
  ref.tiltch=0
  ref.nu=0
  ref.del=0

  oa->set_object, ref

  maskarr=lonarr(ref.nopixx,ref.nopixy)
 ;--------------------------------------
caltype=1
pix=[12,12]
xyz=[0.1,0.2,0.3]
gonio=fltarr(6)
sd=[1,0,0]
XY=pix
DACopening=40.
XYZ_DAC=[1,0,0]
ubc=fltarr(3,3)
ubc[0,0]=.1
ubc[1,1]=.1
ubc[2,2]=.1
wv=0.4
hkl=[1,2,3]
kt=[0,0]
p=0.9
tth=20.0

nbins=1000
npo=100
chimax=200
tthmax=20.0

pred={  om_start : 0.0, $
        om_range : 0.0, $
        chi       : 0.0, $
        d         : 0.0, $
        h1        : 0, $
        h2        : 0, $
        k1        : 0, $
        k2        : 0, $
        l1        : 0, $
        l2        : 0}

ref_peak = {CLASS_peak, $
                Stat     : 0,           $   ; reflection status
                HKL      : INTARR(3),   $   ; Miller indices
                XYZ      : FLTARR (3),  $   ; coordinates of the reciprocal sp. vector
                selected : INTARR(2),   $   ; 0-selcted, 1-visible
                DetXY    : FLTARR(2),   $   ; area detector pixel coordinates
                Gonio    : FLTARR(6),   $   ; goniometer settings for the Area detector
                GonioSS  : FLTARR(6),   $   ; setting angles for solid state detector
                nen      : 0,           $   ; number of different energy components
                energies : FLTARR(10),  $   ; energies
                IntAD    : FLTARR(2),   $   ; Intensity from area detector with e.s.d
                position : FLTARR(3),   $   ; Intensity from area detector with e.s.d
                IntSSD   : FLTARR(2),   $   ; Will now be used to store individualized peak fitting box size
                Adp      : ref_adp}         ; Area detector parameters

a=oa->tilt_mtx(caltype)
a=oa->calculate_sd_from_pixels(pix, gonio, caltype)
a=oa->calculate_pixels_from_sd(sd, gonio, caltype)
a=oa->calculate_tth_from_pixels(pix, gonio, caltype)
a=oa->pixel_visible(XY, DACopening, XYZ_DAC, caltype)
a=oa->generate_one_peak(ubc, wv, pred, hkl, kt, caltype)
a=oa->polarization(pix, p)
a=oa->calculate_pixels_from_xyz(xyz, gonio, caltype)
a=oa->calculate_XYZ_from_pixels_mono(pix, gonio, wv, caltype)
a=oa->generate_ds_ring(tth, kt, caltype)
a=oa->calculate_psi_angles(gonio, pix, caltype)
a=oa->get_nu_from_pix(pix, caltype)
a=oa->create_chi_bin_array(nbins, npo, chimax)
a=oa->create_tth_bin_array(nbins, npo, tthmax, caltype)
oa->visible_mask, maskarr, gonio, DACopening,wv, caltype

;----- the two routines below require peaktable object

;oa->generate_all_peaks, ub, optable, wv, pred, exti, DAC_open, box, kt, rbc, caltype
;oa->generate_peaks_laue, ub, optable, pred, en, exti, DAC_open

 ;--------------------------------------
 obj_destroy, oa
 print, '-- object instance destroyed'
end