pro WID_Peak_editor_cleanup, caller
common status, PE_open, SIM_open
   PE_open=0
end


;----------------------------------------------------------------
;----------------------------------------------------------------
;----------------------------------------------------------------


FUNCTION WID_Peak_editor::init


  self.widgets.WID_BASE_0 = Widget_Base(UNAME='WID_BASE_0'  $
      ,XOFFSET=5 ,YOFFSET=5 ,SCR_XSIZE=714 ,SCR_YSIZE=556  $
      ,TITLE='Peak editor' ,SPACE=3 ,XPAD=3 ,YPAD=3, /tlb_kill_request_events)


  self.widgets.WID_LABEL_0 = Widget_Label(self.widgets.WID_BASE_0, UNAME='WID_LABEL_0'  $
      ,XOFFSET=15 ,YOFFSET=16 ,/ALIGN_LEFT ,VALUE='Peak #')


  self.widgets.WID_TEXT_0 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_0' ,XOFFSET=60  $
      ,YOFFSET=10 ,SCR_XSIZE=54 ,SCR_YSIZE=21 ,XSIZE=20 ,YSIZE=1)


  self.widgets.WID_LABEL_1 = Widget_Label(self.widgets.WID_BASE_0, UNAME='WID_LABEL_1'  $
      ,XOFFSET=15 ,YOFFSET=52 ,/ALIGN_LEFT ,VALUE='hkl')


  self.widgets.WID_TEXT_2 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_2' ,XOFFSET=60  $
      ,YOFFSET=45 ,SCR_XSIZE=54 ,SCR_YSIZE=21 ,/EDITABLE ,XSIZE=20  $
      ,YSIZE=1)


  self.widgets.WID_TEXT_3 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_3'  $
      ,XOFFSET=119 ,YOFFSET=45 ,SCR_XSIZE=54 ,SCR_YSIZE=21 ,/EDITABLE  $
      ,XSIZE=20 ,YSIZE=1)


  self.widgets.WID_TEXT_4 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_4'  $
      ,XOFFSET=177 ,YOFFSET=45 ,SCR_XSIZE=54 ,SCR_YSIZE=21 ,/EDITABLE  $
      ,XSIZE=20 ,YSIZE=1)


  self.widgets.WID_TEXT_5 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_5'  $
      ,XOFFSET=119 ,YOFFSET=81 ,SCR_XSIZE=54 ,SCR_YSIZE=21 ,XSIZE=20  $
      ,YSIZE=1)


  self.widgets.WID_TEXT_6 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_6' ,XOFFSET=60  $
      ,YOFFSET=81 ,SCR_XSIZE=54 ,SCR_YSIZE=21 ,XSIZE=20 ,YSIZE=1)


  self.widgets.WID_LABEL_2 = Widget_Label(self.widgets.WID_BASE_0, UNAME='WID_LABEL_2'  $
      ,XOFFSET=15 ,YOFFSET=88 ,/ALIGN_LEFT ,VALUE='detXY')


  self.widgets.WID_TEXT_7 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_7'  $
      ,XOFFSET=177 ,YOFFSET=113 ,SCR_XSIZE=54 ,SCR_YSIZE=21 ,XSIZE=20  $
      ,YSIZE=1)


  self.widgets.WID_TEXT_8 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_8'  $
      ,XOFFSET=119 ,YOFFSET=113 ,SCR_XSIZE=54 ,SCR_YSIZE=21 ,XSIZE=20  $
      ,YSIZE=1)


  self.widgets.WID_TEXT_9 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_9' ,XOFFSET=60  $
      ,YOFFSET=113 ,SCR_XSIZE=54 ,SCR_YSIZE=21 ,XSIZE=20 ,YSIZE=1)


  self.widgets.WID_LABEL_3 = Widget_Label(self.widgets.WID_BASE_0, UNAME='WID_LABEL_3'  $
      ,XOFFSET=15 ,YOFFSET=120 ,/ALIGN_LEFT ,VALUE='xyz')


  self.widgets.WID_LABEL_4 = Widget_Label(self.widgets.WID_BASE_0, UNAME='WID_LABEL_4'  $
      ,XOFFSET=14 ,YOFFSET=151 ,/ALIGN_LEFT ,VALUE='d-spc')


  self.widgets.WID_TEXT_10 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_10'  $
      ,XOFFSET=59 ,YOFFSET=144 ,SCR_XSIZE=54 ,SCR_YSIZE=21 ,XSIZE=20  $
      ,YSIZE=1)


  self.widgets.WID_TEXT_11 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_11'  $
      ,XOFFSET=58 ,YOFFSET=171 ,SCR_XSIZE=54 ,SCR_YSIZE=21 ,/EDITABLE  $
      ,XSIZE=20 ,YSIZE=1)


  self.widgets.WID_LABEL_5 = Widget_Label(self.widgets.WID_BASE_0, UNAME='WID_LABEL_5'  $
      ,XOFFSET=13 ,YOFFSET=178 ,/ALIGN_LEFT ,VALUE='Energy')


  self.widgets.WID_TEXT_12 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_12'  $
      ,XOFFSET=175 ,YOFFSET=198 ,SCR_XSIZE=54 ,SCR_YSIZE=21  $
      ,/EDITABLE ,XSIZE=20 ,YSIZE=1)


  self.widgets.WID_TEXT_13 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_13'  $
      ,XOFFSET=117 ,YOFFSET=198 ,SCR_XSIZE=54 ,SCR_YSIZE=21  $
      ,/EDITABLE ,XSIZE=20 ,YSIZE=1)


  self.widgets.WID_TEXT_14 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_14'  $
      ,XOFFSET=58 ,YOFFSET=198 ,SCR_XSIZE=54 ,SCR_YSIZE=21 ,/EDITABLE  $
      ,XSIZE=20 ,YSIZE=1)


  self.widgets.WID_LABEL_6 = Widget_Label(self.widgets.WID_BASE_0, UNAME='WID_LABEL_6'  $
      ,XOFFSET=13 ,YOFFSET=205 ,/ALIGN_LEFT ,VALUE='Gonio')


  self.widgets.WID_TEXT_15 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_15'  $
      ,XOFFSET=58 ,YOFFSET=222 ,SCR_XSIZE=54 ,SCR_YSIZE=21 ,/EDITABLE  $
      ,XSIZE=20 ,YSIZE=1)


  self.widgets.WID_TEXT_16 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_16'  $
      ,XOFFSET=117 ,YOFFSET=222 ,SCR_XSIZE=54 ,SCR_YSIZE=21  $
      ,/EDITABLE ,XSIZE=20 ,YSIZE=1)


  self.widgets.WID_TEXT_17 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_17'  $
      ,XOFFSET=174 ,YOFFSET=222 ,SCR_XSIZE=54 ,SCR_YSIZE=21  $
      ,/EDITABLE ,XSIZE=20 ,YSIZE=1)


  self.widgets.WID_BUTTON_0 = Widget_Button(self.widgets.WID_BASE_0, UNAME='WID_BUTTON_0'  $
      ,XOFFSET=535 ,YOFFSET=208 ,SCR_XSIZE=75 ,SCR_YSIZE=43  $
      ,/ALIGN_CENTER ,VALUE='Update')


  self.widgets.WID_BUTTON_1 = Widget_Button(self.widgets.WID_BASE_0, UNAME='WID_BUTTON_1'  $
      ,XOFFSET=616 ,YOFFSET=208 ,SCR_XSIZE=75 ,SCR_YSIZE=43  $
      ,/ALIGN_CENTER ,VALUE='Close')


  self.widgets.WID_DRAW_0 = Widget_Draw(self.widgets.WID_BASE_0, UNAME='WID_DRAW_0' ,XOFFSET=9  $
      ,YOFFSET=271 ,SCR_XSIZE=688 ,SCR_YSIZE=245)


  self.widgets.WID_DRAW_1 = Widget_Draw(self.widgets.WID_BASE_0, UNAME='WID_DRAW_1'  $
      ,XOFFSET=511 ,YOFFSET=32 ,SCR_XSIZE=170 ,SCR_YSIZE=170)


  self.widgets.WID_BUTTON_2 = Widget_Button(self.widgets.WID_BASE_0, UNAME='WID_BUTTON_2'  $
      ,XOFFSET=237 ,YOFFSET=212 ,SCR_XSIZE=21 ,SCR_YSIZE=19  $
      ,/ALIGN_CENTER ,VALUE='-X')


  self.widgets.WID_BUTTON_3 = Widget_Button(self.widgets.WID_BASE_0, UNAME='WID_BUTTON_3'  $
      ,XOFFSET=272 ,YOFFSET=187 ,SCR_XSIZE=21 ,SCR_YSIZE=19  $
      ,/ALIGN_CENTER ,VALUE='+Y')


  self.widgets.WID_BUTTON_4 = Widget_Button(self.widgets.WID_BASE_0, UNAME='WID_BUTTON_4'  $
      ,XOFFSET=305 ,YOFFSET=210 ,SCR_XSIZE=21 ,SCR_YSIZE=20  $
      ,/ALIGN_CENTER ,VALUE='+X')


  self.widgets.WID_BUTTON_5 = Widget_Button(self.widgets.WID_BASE_0, UNAME='WID_BUTTON_5'  $
      ,XOFFSET=271 ,YOFFSET=236 ,SCR_XSIZE=21 ,SCR_YSIZE=18  $
      ,/ALIGN_CENTER ,VALUE='-Y')


  self.widgets.WID_DRAW_2 = Widget_Draw(self.widgets.WID_BASE_0, UNAME='WID_DRAW_2'  $
      ,XOFFSET=325 ,YOFFSET=33 ,SCR_XSIZE=170 ,SCR_YSIZE=170)


  self.widgets.WID_TEXT_18 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_18' ,FRAME=1  $
      ,XOFFSET=264 ,YOFFSET=211 ,SCR_XSIZE=36 ,SCR_YSIZE=21 ,XSIZE=20  $
      ,YSIZE=1)


  self.widgets.WID_BUTTON_6 = Widget_Button(self.widgets.WID_BASE_0, UNAME='WID_BUTTON_6'  $
      ,XOFFSET=340 ,YOFFSET=208 ,SCR_XSIZE=76 ,SCR_YSIZE=22  $
      ,/ALIGN_CENTER ,VALUE='Fit')


  self.widgets.WID_BUTTON_7 = Widget_Button(self.widgets.WID_BASE_0, UNAME='WID_BUTTON_7'  $
      ,XOFFSET=256 ,YOFFSET=10 ,SCR_XSIZE=21 ,SCR_YSIZE=19  $
      ,/ALIGN_CENTER ,VALUE='<')


  self.widgets.WID_BUTTON_8 = Widget_Button(self.widgets.WID_BASE_0, UNAME='WID_BUTTON_8'  $
      ,XOFFSET=284 ,YOFFSET=9 ,SCR_XSIZE=21 ,SCR_YSIZE=20  $
      ,/ALIGN_CENTER ,VALUE='>')


  self.widgets.WID_BUTTON_9 = Widget_Button(self.widgets.WID_BASE_0, UNAME='WID_BUTTON_9'  $
      ,XOFFSET=256 ,YOFFSET=34 ,SCR_XSIZE=48 ,SCR_YSIZE=24  $
      ,/ALIGN_CENTER ,VALUE='Prof.')


  self.widgets.WID_TEXT_19 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_19'  $
      ,XOFFSET=325 ,YOFFSET=7 ,SCR_XSIZE=272 ,SCR_YSIZE=21 ,XSIZE=20  $
      ,YSIZE=1)


  self.widgets.WID_TEXT_20 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_20'  $
      ,XOFFSET=606 ,YOFFSET=6 ,SCR_XSIZE=82 ,SCR_YSIZE=22 ,XSIZE=20  $
      ,YSIZE=1)


  self.widgets.WID_BUTTON_10 = Widget_Button(self.widgets.WID_BASE_0, UNAME='WID_BUTTON_10'  $
      ,XOFFSET=255 ,YOFFSET=63 ,SCR_XSIZE=49 ,SCR_YSIZE=24  $
      ,/ALIGN_CENTER ,VALUE='Print')


  self.widgets.WID_TEXT_21 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_21'  $
      ,XOFFSET=117 ,YOFFSET=246 ,SCR_XSIZE=54 ,SCR_YSIZE=21  $
      ,/EDITABLE ,XSIZE=20 ,YSIZE=1)


  self.widgets.WID_TEXT_22 = Widget_Text(self.widgets.WID_BASE_0, UNAME='WID_TEXT_22'  $
      ,XOFFSET=58 ,YOFFSET=246 ,SCR_XSIZE=54 ,SCR_YSIZE=21 ,/EDITABLE  $
      ,XSIZE=20 ,YSIZE=1)


  self.widgets.WID_LABEL_7 = Widget_Label(self.widgets.WID_BASE_0, UNAME='WID_LABEL_7'  $
      ,XOFFSET=13 ,YOFFSET=250 ,/ALIGN_LEFT ,VALUE='Intensity')

  Widget_Control, self.widgets.wid_base_0, set_uvalue=self

  Widget_Control, /REALIZE, self.widgets.WID_BASE_0

  widget_control, self.widgets.wid_text_18, editable=1
  widget_control, self.widgets.wid_text_18, set_value='2'
  widget_control, self.widgets.wid_draw_0, get_value=d0
  widget_control, self.widgets.wid_draw_1, get_value=d1
  widget_control, self.widgets.wid_draw_2, get_value=d2
  self.widgets.wid_draw_0=d0
  self.widgets.wid_draw_1=d1
  self.widgets.wid_draw_2=d2

  XManager, 'WID_Peak_editor', self.widgets.WID_BASE_0, /NO_BLOCK

  RETURN, 1

end

;--------------------------

pro WID_Peak_editor_event, ev
 widget_control, ev.top, get_uvalue=wid
 wid->event, ev
end

;--------------------------
pro WID_Peak_editor::hide
widget_control, self.widgets.wid_base_0, map=0
end

;--------------------------

pro WID_Peak_editor::show
widget_control, self.widgets.wid_base_0, map=1
end

;--------------------------

pro WID_Peak_editor::plot_peak, pn, pt, oimage

  common selections, selecting_status, select, excl, unselect, addpeak, mask, maskarr, zoom, unmask
  COMMON image_type_and_arrays, imt, arr1, arr2,cenx, ceny, rad, rad1
   ;bs=read_box_size()
   bs=[long(pt.peaks[pn].intssd[0]),long(pt.peaks[pn].intssd[0])]
   XY=pt.peaks[pn].detxy
   if not(XY[0]-bs[0] le 0 or $
           arr1 - XY[0] le bs[0] or $
           XY[1]-bs[1] le 0 or $
           arr2 - XY[1] le bs[1]) then $
   begin
    pic=oimage->get_zoomin(XY, bs, maskarr)
    pic1=congrid(pic, 170, 170)
    iscl=read_color_scale()
    pic1=pic1<iscl[1]
    pic1=pic1>iscl[0]
    wset, self.widgets.wid_draw_2
    tvscl, pic1
   endif else re=dialog_message('Peak too close to the edge of the image')
  end

;--------------------------

function WID_Peak_editor::read_peakshift_step
   widget_control, self.widgets.wid_text_18, get_value=pss
   return, float(pss)

end

;--------------------------

pro WID_Peak_editor::update_peak
  common local_pt,pt, pn
  COMMON CLASS_peaktable_reference, ref_peaktable, ref_peak

  pt=self.pt
  pn=self.pn

  if pn ge 0 and pn le pt.peakno-1 then $
  begin
  ref_peak=pt.peaks[pn]
  al=''
  widget_control, self.widgets.wid_text_2, get_value=al
  ref_peak.hkl[0]=long(al)
  widget_control, self.widgets.wid_text_3, get_value=al
  ref_peak.hkl[1]=long(al)
  widget_control, self.widgets.wid_text_4, get_value=al
  ref_peak.hkl[2]=long(al)
  widget_control, self.widgets.wid_text_11, get_value=al ; energy
  ref_peak.energies[0]=float(al)
  widget_control, self.widgets.wid_text_12, get_value=al
  ref_peak.gonio[2]=float(al)
  widget_control, self.widgets.wid_text_13, get_value=al
  ref_peak.gonio[1]=float(al)
  widget_control, self.widgets.wid_text_14, get_value=al
  ref_peak.gonio[0]=float(al)
  widget_control, self.widgets.wid_text_15, get_value=al
  ref_peak.gonio[3]=float(al)
  widget_control, self.widgets.wid_text_16, get_value=al
  ref_peak.gonio[4]=float(al)
  widget_control, self.widgets.wid_text_17, get_value=al
  ref_peak.gonio[5]=float(al)
  self.pt.peaks[pn]=ref_peak
  endif
end


;--------------------------
pro WID_Peak_editor::pe_print_filename, fm
  widget_control, self.widgets.wid_text_19, set_value=fm
end
;--------------------------
pro WID_Peak_editor::pe_print_omega, om
  widget_control, self.widgets.wid_text_20, set_value=string(om, format='(F6.1)')
end
;--------------------------
pro WID_Peak_editor::display_peak, pn, pt
  widget_control, self.widgets.wid_text_0, set_value=string(pn, format='(I4)')
  widget_control, self.widgets.wid_text_2, set_value=string(pt.peaks[pn].hkl[0], format='(I3)')
  widget_control, self.widgets.wid_text_3, set_value=string(pt.peaks[pn].hkl[1], format='(I3)')
  widget_control, self.widgets.wid_text_4, set_value=string(pt.peaks[pn].hkl[2], format='(I3)')
  widget_control, self.widgets.wid_text_5, set_value=string(pt.peaks[pn].detXY[1], format='(F7.2)')
  widget_control, self.widgets.wid_text_6, set_value=string(pt.peaks[pn].detXY[0], format='(F7.2)')
  widget_control, self.widgets.wid_text_7, set_value=string(pt.peaks[pn].xyz[2], format='(F7.4)')
  widget_control, self.widgets.wid_text_8, set_value=string(pt.peaks[pn].xyz[1], format='(F7.4)')
  widget_control, self.widgets.wid_text_9, set_value=string(pt.peaks[pn].xyz[0], format='(F7.4)')
  widget_control, self.widgets.wid_text_10, set_value=string(1./vlength(pt.peaks[pn].xyz), format='(F7.4)')
  widget_control, self.widgets.wid_text_11, set_value=string(pt.peaks[pn].energies[0], format='(F6.2)')
  widget_control, self.widgets.wid_text_12, set_value=string(pt.peaks[pn].gonio[2], format='(F10.2)')
  widget_control, self.widgets.wid_text_13, set_value=string(pt.peaks[pn].gonio[1], format='(F10.2)')
  widget_control, self.widgets.wid_text_14, set_value=string(pt.peaks[pn].gonio[0], format='(F10.2)')
  widget_control, self.widgets.wid_text_15, set_value=string(pt.peaks[pn].gonio[3], format='(F10.2)')
  widget_control, self.widgets.wid_text_16, set_value=string(pt.peaks[pn].gonio[4], format='(F10.2)')
  widget_control, self.widgets.wid_text_17, set_value=string(pt.peaks[pn].gonio[5], format='(F10.2)')
  widget_control, self.widgets.wid_text_22, set_value=string(pt.peaks[pn].intAD[0], format='(F10.2)')
  widget_control, self.widgets.wid_text_21, set_value=string(pt.peaks[pn].intAD[1], format='(F10.2)')
end


;---------------------------------------------------------------------------
;---------------------------------------------------------------------------
;---------------------------------------------------------------------------

pro WID_Peak_editor::event, ev

COMMON image_type_and_arrays, imt, arr1, arr2,cenx, ceny, rad, rad1
common selections, selecting_status, select, excl, unselect, addpeak, mask, maskarr, zoom, unmask
@COMMON_DATAS

@WID_profile_film_viewer_commons

common status, PE_open, SIM_open

widget_control, ev.id, get_uvalue=uv

pn=self.pn

widget_control, ev.top, get_uvalue=wid

if (tag_names(ev, /structure_name) eq 'WIDGET_KILL_REQUEST') then $
begin
  print, 'kill call captured'
  widget_control, self.widgets.wid_base_0, map=0
  PE_open=0
  return
endif

case ev.id of

 self.widgets.wid_button_0:$                        ; 'Update':$
  begin
    self->update_peak
    opt->set_object, self.pt
  end

;----

  self.widgets.wid_button_2:$                        ;'-X':$
  begin
   self.pt.peaks[pn].detxy[0]=self.pt.peaks[pn].detxy[0]-self->read_peakshift_step()
   self->plot_peak, self.pn, self.pt, oimage
   self->display_peak, self.pn, self.pt
  end

;----

  self.widgets.wid_button_3:$                        ;'+Y':$
  begin
   self.pt.peaks[pn].detxy[1]=self.pt.peaks[pn].detxy[1]+self->read_peakshift_step()
   self->plot_peak, self.pn, self.pt, oimage
   self->display_peak, self.pn, self.pt
   end

  self.widgets.wid_button_4:$                        ;'+X':$
  begin
   self.pt.peaks[pn].detxy[0]=self.pt.peaks[pn].detxy[0]+self->read_peakshift_step()
   self->plot_peak, self.pn, self.pt, oimage
   self->display_peak, self.pn, self.pt
   end

  self.widgets.wid_button_5:$                        ;'-Y':$
  begin
   self.pt.peaks[pn].detxy[1]=self.pt.peaks[pn].detxy[1]-self->read_peakshift_step()
   self->plot_peak, self.pn, self.pt, oimage
   self->display_peak, self.pn, self.pt
   end

  self.widgets.wid_button_6:$                        ;'fit':$
  begin
  pt=opt->get_object()
  bs=[long(pt.peaks[pn].intssd[0]),long(pt.peaks[pn].intssd[0])]
   ;bs=read_box_size()
   self.pt.peaks[pn].detxy[0]=self.pt.peaks[pn].detxy[0]+self->read_peakshift_step()
   XY=self.pt.peaks[pn].detxy
   if not(XY[0]-bs[0] le 0 or $
           arr1 - XY[0] le bs[0] or $
           XY[1]-bs[1] le 0 or $
           arr2 - XY[1] le bs[1]) then $
   begin
    pic=oimage->get_zoomin(XY, bs, maskarr)
    pic1=congrid(pic, 170, 170)
    wset, self.widgets.wid_draw_2
    tvscl, pic1
    self->display_peak, self.pn, self.pt
    Gauss = GAUSS2DFIT( pic, A, /TILT)
    pic2=congrid(Gauss, 170, 170)
    wset, self.widgets.wid_draw_1
    tvscl, pic2
    self.pt.peaks[pn].DetXY[0]= (A[4]-bs[0])+XY[0]
    self.pt.peaks[pn].DetXY[1]= (A[5]-bs[1])+XY[1]
  endif else re=dialog_message('Peak too close to the edge of the image')
  end


  self.widgets.wid_button_7:$                        ;'<':$
  begin
       if res.name0 ne '' and res.seq-1 ge f0 then $
     begin
       res.seq=res.seq-1
       fn=generate_fname(res)
       res=analyse_fname(fn, dir, 3)
       self->pe_print_filename, res.name0
       self->pe_print_omega, omega_from_scan(res.seq)
       oimage->load_image, fn, oadetector

       bs=[long(pt.peaks[pn].intssd[0]),long(pt.peaks[pn].intssd[0])]
       ;bs=read_box_size()
       self.pt.peaks[pn].detxy[1]=self.pt.peaks[pn].detxy[1]-2
       XY=self.pt.peaks[pn].detxy
       if not(XY[0]-bs[0] le 0 or $
           arr1 - XY[0] le bs[0] or $
           XY[1]-bs[1] le 0 or $
           arr2 - XY[1] le bs[1]) then $
      begin
       pic=oimage->get_zoomin(XY, bs, maskarr)
       pic1=congrid(pic, 170, 170)
       iscl=read_color_scale()
       pic1=pic1<iscl[1]
       pic1=pic1>iscl[0]
       wset, self.widgets.wid_draw_2
       tvscl, pic1
       self->display_peak, self.pn, self.pt
      end
      end
  end


  self.widgets.wid_button_8:$                        ;'>':$
begin
     if res.name0 ne '' and res.seq+1 le f1 then $
     begin
       res.seq=res.seq+1
       fn=generate_fname(res)
       res=analyse_fname(fn, dir, 3)
       self->pe_print_filename, res.name0
       self->pe_print_omega, omega_from_scan(res.seq)
       oimage->load_image, fn, oadetector
       bs=[long(pt.peaks[pn].intssd[0]),long(pt.peaks[pn].intssd[0])]
       ;bs=read_box_size()
       self.pt.peaks[pn].detxy[1]=self.pt.peaks[pn].detxy[1]-2
       XY=self.pt.peaks[pn].detxy
       if not(XY[0]-bs[0] le 0 or $
           arr1 - XY[0] le bs[0] or $
           XY[1]-bs[1] le 0 or $
           arr2 - XY[1] le bs[1]) then $
      begin
       pic=oimage->get_zoomin(XY, bs, maskarr)
       pic1=congrid(pic, 170, 170)
       iscl=read_color_scale()
       pic1=pic1<iscl[1]
       pic1=pic1>iscl[0]
       wset, self.widgets.wid_draw_2
       tvscl, pic1
       self->display_peak, self.pn, self.pt
      end
      end
 end

 self.widgets.wid_button_9:$                        ;'prof':$
  begin
       pt=opt->get_object()
       prof=fltarr(im_seq_n())
       b=fltarr(170)
       if res.name0 ne '' then  $
       begin
       for i=0, im_seq_n()-1 do $
       begin
        res.seq=i+im_seq_0()
        print,res.seq
        fn=generate_fname(res)
        res=analyse_fname(fn, dir, 3)
        oimage->load_image, fn, oadetector
        bs=[long(pt.peaks[pn].intssd[0]),long(pt.peaks[pn].intssd[0])]
        ;bs=read_box_size()
        XY=self.pt.peaks[pn].detxy
        if not(XY[0]-bs[0] le 0 or $
            arr1 - XY[0] le bs[0] or $
            XY[1]-bs[1] le 0 or $
            arr2 - XY[1] le bs[1]) then $
        begin
         pic=oimage->get_zoomin(XY, bs, maskarr)
         pic1=congrid(pic, 170, 170)
         prof[i]=total(pic)/n_elements(pic)
         if i eq 0 then film=pic1 else $
         begin
           film=[[film],[pic1]]
         end
        endif
        update_progress, float(i+1)/float(im_seq_n())
      endfor

      iscl=read_color_scale()
      film=film<iscl[1]
      film=film>iscl[0]

      for i=0, 169 do b[i]=max(film)
      for i=0, im_seq_n()-1 do $
      begin
        film[0:169,i*170]=b
      endfor
      wset, self.widgets.wid_draw_0
      plot, prof
      re=dialog_message('Omega profile computation complete')

      wid_film->show
      wid_film->resize_film, im_seq_n()
      wid_film->show
      wid_film->display, film, im_seq_n()
      endif
  end
  self.widgets.wid_button_10:$                        ;'print':$
  begin
    print, 'Print'
    print, prof

    re=dialog_pickfile(path=out_dir, default_extension=',txt',/write)
    if re ne '' then $
    begin
    get_lun, lu
    openw,lu,re
    for i=0, n_elements(prof)-1 do  printf,lu, prof[i]
    close, lu

    endif

  end

  self.widgets.wid_button_1:$                        ;'Close':$
  BEGIN
   widget_control, ev.top, map=0
   PE_open=0
  END
  else:
  endcase

end
;--------------------------



pro WID_Peak_editor::fillup , pn, pt

 COMMON image_type_and_arrays, imt, arr1, arr2,cenx, ceny, rad, rad1
  common datas


  self->display_peak, pn, pt
  self->plot_peak, pn, pt, oimage

  self->pe_print_filename, res.name0
  self->pe_print_omega, omega_from_scan(res.seq)
  self.pt=pt
  self.pn=pn

end
;--------------------------

pro Peak_editor, pl, opti
common datas
common local_pt,pt, pn
pn=pl
pt=opt->get_object()
  if pn ge 0 and pn le opt->peakno()-1 then $
  begin
   WID_Peak_editor
  end
end

;--------------------------

pro ed
end

;-------------------------------------------------------
;-------------------------------------------------------
;-------------------------------------------------------


pro WID_Peak_editor__define

COMMON CLASS_peaktable_reference, ref_peaktable, ref_peak

widgets={WID_Peak_editor_widgets, $
  WID_BASE_0 :0L, $
  WID_LABEL_0:0L, $
  WID_TEXT_0 :0L, $
  WID_LABEL_1 :0L, $
  WID_TEXT_2 :0L, $
  WID_TEXT_3 :0L, $
  WID_TEXT_4 :0L, $
  WID_TEXT_5 :0L, $
  WID_TEXT_6 :0L, $
  WID_LABEL_2 :0L, $
  WID_TEXT_7 :0L, $
  WID_TEXT_8 :0L, $
  WID_TEXT_9 :0L, $
  WID_LABEL_3 :0L, $
  WID_LABEL_4 :0L, $
  WID_TEXT_10 :0L, $
  WID_TEXT_11 :0L, $
  WID_LABEL_5 :0L, $
  WID_TEXT_12 :0L, $
  WID_TEXT_13 :0L, $
  WID_TEXT_14 :0L, $
  WID_LABEL_6 :0L, $
  WID_TEXT_15 :0L, $
  WID_TEXT_16 :0L, $
  WID_TEXT_17 :0L, $
  WID_BUTTON_0 :0L, $
  WID_BUTTON_1 :0L, $
  WID_DRAW_0 :0L, $
  WID_DRAW_1 :0L, $
  WID_BUTTON_2 :0L, $
  WID_BUTTON_3 :0L, $
  WID_BUTTON_4 :0L, $
  WID_BUTTON_5 :0L, $
  WID_DRAW_2 :0L, $
  WID_TEXT_18:0L, $
  WID_BUTTON_6 :0L, $
  WID_BUTTON_7:0L, $
  WID_BUTTON_8 :0L, $
  WID_BUTTON_9 :0L, $
  WID_TEXT_19 :0L, $
  WID_TEXT_20 :0L, $
  WID_BUTTON_10 :0L, $
  WID_TEXT_21 :0L, $
  WID_TEXT_22 :0L, $
  WID_LABEL_7 :0L}


  WID_Peak_editor={WID_Peak_editor, widgets: widgets, pt: ref_peaktable, pn:0L}

end

