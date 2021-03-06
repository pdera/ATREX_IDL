pro WID_profile_film_viewer_event, ev
   widget_control, ev.top, get_uvalue=wid
   wid->event, ev

end

;-------------------------------------------------------------

function WID_profile_film_viewer::resize_film, n
 if n gt 0 then $
 widget_control,self.widgets.WID_DRAW_0,XSIZE=fix(n)*170
end

;-------------------------------------------------------------

function WID_profile_film_viewer::init

 self.widgets.WID_profile_film_viewer = Widget_Base(UNAME='WID_profile_film_viewer' $
      ,XOFFSET=5 ,YOFFSET=5  $
      ,SCR_XSIZE=873 ,SCR_YSIZE=300 ,TITLE='Profile film viewer'  $
      ,SPACE=3 ,XPAD=3 ,YPAD=3)

 self.widgets.WID_DRAW_0 = Widget_Draw(self.widgets.WID_profile_film_viewer,  $
      UNAME='WID_DRAW_0' ,XOFFSET=5 ,YOFFSET=6 ,SCR_XSIZE=850, SCR_YSIZE=200  $
      ,/SCROLL ,XSIZE=1700 ,YSIZE=170)

  self.widgets.WID_BUTTON_1 = Widget_Button(self.widgets.WID_profile_film_viewer,  $
      UNAME='WID_BUTTON_1' ,XOFFSET=10 ,YOFFSET=220 ,SCR_XSIZE=60  $
      ,SCR_YSIZE=29 ,/ALIGN_CENTER ,VALUE='Save')


  self.widgets.WID_BUTTON_2 = Widget_Button(self.widgets.WID_profile_film_viewer,  $
      UNAME='WID_BUTTON_2' ,XOFFSET=72 ,YOFFSET=220 ,SCR_XSIZE=60  $
      ,SCR_YSIZE=29 ,/ALIGN_CENTER ,VALUE='Close')

  Widget_Control, /REALIZE, self.widgets.WID_profile_film_viewer

  widget_control, self.widgets.wid_draw_0, get_value=st
  self.widgets.wid_draw_0=st
  widget_Control, self.widgets.WID_profile_film_viewer, set_uvalue=self

  XManager, 'WID_profile_film_viewer', self.widgets.WID_profile_film_viewer, /NO_BLOCK

end

;------------------------------------------------------------

pro WID_profile_film_viewer::event, ev

if (tag_names(ev, /structure_name) eq 'WIDGET_KILL_REQUEST') then $
begin
  print, 'kill call captured'
  widget_control, self.widgets.WID_profile_film_viewer, map=0
  return
endif

case ev.id of

  self.widgets.WID_BUTTON_1:$                ;'save':$
    begin
      wset, self.widgets.WID_DRAW_0
      img=tvrd(0,0,xsize,170)
      fn=dialog_pickfile(/write, filter='*.jpg', default_extension='.jpg', path=out_dir)
      if fn ne '' then write_jpeg, fn, img, quality=100
    end
  self.widgets.WID_BUTTON_2:$                ;'close':$
    begin
         widget_control, self.widgets.WID_profile_film_viewer, map=0
    end
  else:
  endcase

end

;------------------------------------------------------------
;-------------------------------------------------------------
;-------------------------------------------------------------

pro WID_profile_film_viewer__define

  widgets={WID_profile_film_viewer_widgets, $
  WID_profile_film_viewer :0L,$
  WID_DRAW_0 :0L,$
  WID_BUTTON_1 :0L,$
  WID_BUTTON_2:0L}

  WID_profile_film_viewer={WID_profile_film_viewer, widgets:widgets}

end
