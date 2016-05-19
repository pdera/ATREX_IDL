; 
; IDL Widget Interface Procedures. This Code is automatically 
;     generated and should not be modified.

; 
; Generated on:	04/26/2016 15:51.25
; 
pro WID_Cij_MainWindow_event, Event

  wTarget = (widget_info(Event.id,/NAME) eq 'TREE' ?  $
      widget_info(Event.id, /tree_root) : event.id)


  wWidget =  Event.top

  case wTarget of

    Widget_Info(wWidget, FIND_BY_UNAME='WID_Cij_MainWindow'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_KILL_REQUEST' )then $
        OnKill, Event
    end
    Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_Openfile'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        Cij_Openfile, Event
    end
    Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_Close'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        Cij_CloseWindow, Event
    end
    Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_CalculateStress'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        Cij_CalculateStress, Event
    end
    Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_CalculateStrain'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        Cij_CalculateStrain, Event
    end
    Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_ExportResults'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        Cij_ExportResults, Event
    end
    Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_ref_cubic'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        OnReferenceCubic, Event
    end
    Widget_Info(wWidget, FIND_BY_UNAME='WID_BUTTON_ref_custom'): begin
      if( Tag_Names(Event, /STRUCTURE_NAME) eq 'WIDGET_BUTTON' )then $
        OnReferenceCustom, Event
    end
    else:
  endcase

end

pro WID_Cij_MainWindow, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_

  Resolve_Routine, 'Wid_Cij_eventcb',/COMPILE_FULL_FILE  ; Load event callback routines
  
  WID_Cij_MainWindow = Widget_Base( GROUP_LEADER=wGroup,  $
      UNAME='WID_Cij_MainWindow' ,XOFFSET=5 ,YOFFSET=5 ,SCR_XSIZE=926  $
      ,SCR_YSIZE=496 ,NOTIFY_REALIZE='OnRealize'  $
      ,KILL_NOTIFY='OnDestroy' ,/TLB_KILL_REQUEST_EVENTS ,TITLE='Cij'  $
      ,SPACE=3 ,XPAD=3 ,YPAD=3)

  
  WID_LIST_Cij = Widget_List(WID_Cij_MainWindow, UNAME='WID_LIST_Cij'  $
      ,XOFFSET=29 ,YOFFSET=51 ,SCR_XSIZE=342 ,SCR_YSIZE=187 ,XSIZE=11  $
      ,YSIZE=2)

  
  WID_BUTTON_Openfile = Widget_Button(WID_Cij_MainWindow,  $
      UNAME='WID_BUTTON_Openfile' ,XOFFSET=542 ,YOFFSET=257  $
      ,SCR_XSIZE=110 ,SCR_YSIZE=30 ,/ALIGN_CENTER ,VALUE='Open cij'+ $
      ' file')

  
  WID_BUTTON_Close = Widget_Button(WID_Cij_MainWindow,  $
      UNAME='WID_BUTTON_Close' ,XOFFSET=778 ,YOFFSET=257  $
      ,SCR_XSIZE=110 ,SCR_YSIZE=30 ,/ALIGN_CENTER ,VALUE='Close')

  
  WID_LIST_strain = Widget_List(WID_Cij_MainWindow,  $
      UNAME='WID_LIST_strain' ,XOFFSET=381 ,YOFFSET=51 ,SCR_XSIZE=248  $
      ,SCR_YSIZE=187 ,XSIZE=11 ,YSIZE=2)

  
  WID_LIST_Stress = Widget_List(WID_Cij_MainWindow,  $
      UNAME='WID_LIST_Stress' ,XOFFSET=638 ,YOFFSET=50 ,SCR_XSIZE=248  $
      ,SCR_YSIZE=187 ,XSIZE=11 ,YSIZE=2)

  
  WID_LABEL_0 = Widget_Label(WID_Cij_MainWindow, UNAME='WID_LABEL_0'  $
      ,XOFFSET=29 ,YOFFSET=31 ,/ALIGN_LEFT ,VALUE='Cij')

  
  WID_LABEL_1 = Widget_Label(WID_Cij_MainWindow, UNAME='WID_LABEL_1'  $
      ,XOFFSET=385 ,YOFFSET=30 ,/ALIGN_LEFT ,VALUE='Strain deviator')

  
  WID_LABEL_stress = Widget_Label(WID_Cij_MainWindow,  $
      UNAME='WID_LABEL_stress' ,XOFFSET=642 ,YOFFSET=26 ,/ALIGN_LEFT  $
      ,VALUE='Stress deviator')

  
  WID_BUTTON_CalculateStress = Widget_Button(WID_Cij_MainWindow,  $
      UNAME='WID_BUTTON_CalculateStress' ,XOFFSET=658 ,YOFFSET=293  $
      ,SCR_XSIZE=110 ,SCR_YSIZE=30 ,/ALIGN_CENTER ,VALUE='Calculate'+ $
      ' stress')

  
  WID_BUTTON_CalculateStrain = Widget_Button(WID_Cij_MainWindow,  $
      UNAME='WID_BUTTON_CalculateStrain' ,XOFFSET=658 ,YOFFSET=257  $
      ,SCR_XSIZE=110 ,SCR_YSIZE=30 ,/ALIGN_CENTER ,VALUE='Calculate'+ $
      ' strain')

  
  WID_BUTTON_ExportResults = Widget_Button(WID_Cij_MainWindow,  $
      UNAME='WID_BUTTON_ExportResults' ,XOFFSET=542 ,YOFFSET=293  $
      ,SCR_XSIZE=110 ,SCR_YSIZE=30 ,/ALIGN_CENTER ,VALUE='Export'+ $
      ' results')

  
  WID_LABEL_2 = Widget_Label(WID_Cij_MainWindow, UNAME='WID_LABEL_2'  $
      ,XOFFSET=29 ,YOFFSET=246 ,/ALIGN_LEFT ,VALUE='Reference state')

  
  WID_LIST_reference = Widget_List(WID_Cij_MainWindow,  $
      UNAME='WID_LIST_reference' ,XOFFSET=29 ,YOFFSET=265  $
      ,SCR_XSIZE=195 ,SCR_YSIZE=149 ,XSIZE=11 ,YSIZE=2)

  
  WID_BASE_0 = Widget_Base(WID_Cij_MainWindow, UNAME='WID_BASE_0'  $
      ,XOFFSET=232 ,YOFFSET=266 ,TITLE='IDL' ,COLUMN=1 ,/EXCLUSIVE)

  
  WID_BUTTON_ref_cubic = Widget_Button(WID_BASE_0,  $
      UNAME='WID_BUTTON_ref_cubic' ,/ALIGN_LEFT ,VALUE='Cubic'+ $
      ' iso-volume')

  
  WID_BUTTON_ref_custom = Widget_Button(WID_BASE_0,  $
      UNAME='WID_BUTTON_ref_custom' ,/ALIGN_LEFT ,VALUE='Custom')

  Widget_Control, /REALIZE, WID_Cij_MainWindow

  XManager, 'WID_Cij_MainWindow', WID_Cij_MainWindow, /NO_BLOCK  ,CLEANUP='OnDestroy'  

end
; 
; Empty stub procedure used for autoloading.
; 
pro Wid_Cij, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_
  WID_Cij_MainWindow, GROUP_LEADER=wGroup, _EXTRA=_VWBExtra_
end
