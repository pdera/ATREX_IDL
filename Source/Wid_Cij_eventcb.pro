;
; IDL Event Callback Procedures
; Wid_Cij_eventcb
;
; Generated on:	01/13/2016 20:23.46
;
;
; Empty stub procedure used for autoloading.
;
pro Wid_Cij_eventcb
@COMMON_DATAS

end
;-----------------------------------------------------------------
; Activate Button Callback Procedure.
; Argument:
;   Event structure:
;
;   {WIDGET_BUTTON, ID:0L, TOP:0L, HANDLER:0L, SELECT:0}
;
;   ID is the widget ID of the component generating the event. TOP is
;       the widget ID of the top level widget containing ID. HANDLER
;       contains the widget ID of the widget associated with the
;       handler routine.

;   SELECT is set to 1 if the button was set, and 0 if released.
;       Normal buttons do not generate events when released, so
;       SELECT will always be 1. However, toggle buttons (created by
;       parenting a button to an exclusive or non-exclusive base)
;       return separate events for the set and release actions.

;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro Cij_Openfile, Event
@COMMON_DATAS

 Print, 'Open File'
  fname=dialog_pickfile(FILTER='*.cij', /READ, path=out_dir)
  if fname ne '' then $
     begin
        Cij=open_Cij(fname)
        print, cij
        list=widget_info(event.top, find_by_uname='WID_LIST_Cij')
        print_Cij, Cij,list
     end
end
;-----------------------------------------------------------------
; Activate Button Callback Procedure.
; Argument:
;   Event structure:
;
;   {WIDGET_BUTTON, ID:0L, TOP:0L, HANDLER:0L, SELECT:0}
;
;   ID is the widget ID of the component generating the event. TOP is
;       the widget ID of the top level widget containing ID. HANDLER
;       contains the widget ID of the widget associated with the
;       handler routine.

;   SELECT is set to 1 if the button was set, and 0 if released.
;       Normal buttons do not generate events when released, so
;       SELECT will always be 1. However, toggle buttons (created by
;       parenting a button to an exclusive or non-exclusive base)
;       return separate events for the set and release actions.

;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro Cij_CloseWindow, Event
  Print, 'Close Window'
  	widget_control, event.top,/destroy
end
;-----------------------------------------------------------------
; Activate Button Callback Procedure.
; Argument:
;   Event structure:
;
;   {WIDGET_BUTTON, ID:0L, TOP:0L, HANDLER:0L, SELECT:0}
;
;   ID is the widget ID of the component generating the event. TOP is
;       the widget ID of the top level widget containing ID. HANDLER
;       contains the widget ID of the widget associated with the
;       handler routine.

;   SELECT is set to 1 if the button was set, and 0 if released.
;       Normal buttons do not generate events when released, so
;       SELECT will always be 1. However, toggle buttons (created by
;       parenting a button to an exclusive or non-exclusive base)
;       return separate events for the set and release actions.

;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro Cij_CalculateStress, Event
  @COMMON_DATAS
  if total(Cij) gt 0.0 then $
  begin
  sig=calculate_strain(ub, lp_ref)
  list=widget_info(event.top, find_by_uname='WID_LIST_strain')
  print_strain, sig, list

  sig=calculate_stress(ub,cij, lp_ref)
  list=widget_info(event.top, find_by_uname='WID_LIST_Stress')
  print_strain, sig, list
  endif else  re=dialog_message('Read Cij first')


end
;-----------------------------------------------------------------
; Activate Button Callback Procedure.
; Argument:
;   Event structure:
;
;   {WIDGET_BUTTON, ID:0L, TOP:0L, HANDLER:0L, SELECT:0}
;
;   ID is the widget ID of the component generating the event. TOP is
;       the widget ID of the top level widget containing ID. HANDLER
;       contains the widget ID of the widget associated with the
;       handler routine.

;   SELECT is set to 1 if the button was set, and 0 if released.
;       Normal buttons do not generate events when released, so
;       SELECT will always be 1. However, toggle buttons (created by
;       parenting a button to an exclusive or non-exclusive base)
;       return separate events for the set and release actions.

;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro Cij_CalculateStrain, Event
@COMMON_DATAS
  if total(cij) gt 0.0 then $
  begin
  sig=calculate_strain(ub, lp_ref)
  list=widget_info(event.top, find_by_uname='WID_LIST_strain')
  print_strain, sig, list
  list_ref=widget_info(event.top, find_by_uname='WID_LIST_reference')
  print_reference, lp_ref, list_ref
  end else re=dialog_message('Read Cij first')

end
;-----------------------------------------------------------------
; Activate Button Callback Procedure.
; Argument:
;   Event structure:
;
;   {WIDGET_BUTTON, ID:0L, TOP:0L, HANDLER:0L, SELECT:0}
;
;   ID is the widget ID of the component generating the event. TOP is
;       the widget ID of the top level widget containing ID. HANDLER
;       contains the widget ID of the widget associated with the
;       handler routine.

;   SELECT is set to 1 if the button was set, and 0 if released.
;       Normal buttons do not generate events when released, so
;       SELECT will always be 1. However, toggle buttons (created by
;       parenting a button to an exclusive or non-exclusive base)
;       return separate events for the set and release actions.

;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro Cij_ExportResults, Event

end
;-----------------------------------------------------------------
; Notify Realize Callback Procedure.
; Argument:
;   wWidget - ID number of specific widget.
;
;
;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro OnRealize, wWidget
@COMMON_DATAS
@WID_GSE_ADA_COMMON
 	print, 'On realize'
 	    Cij_window=wWidget
        cubic=widget_info(wWidget, find_by_uname='WID_BUTTON_ref_cubic')
        widget_control,cubic, set_button=1
        list=widget_info(wWidget, find_by_uname='WID_LIST_Cij')
        print, list
        print, ub
        print_Cij, Cij,list

        V0=V_from_UB(ub)
    	a0=V0^(1.0/3)
    	lp_ref=[a0,a0,a0,90.,90.,90.]
    	list_ref=widget_info(wWidget, find_by_uname='WID_LIST_reference')
    	print_reference, lp_ref, list_ref

end
;-----------------------------------------------------------------
; Activate Button Callback Procedure.
; Argument:
;   Event structure:
;
;   {WIDGET_BUTTON, ID:0L, TOP:0L, HANDLER:0L, SELECT:0}
;
;   ID is the widget ID of the component generating the event. TOP is
;       the widget ID of the top level widget containing ID. HANDLER
;       contains the widget ID of the widget associated with the
;       handler routine.

;   SELECT is set to 1 if the button was set, and 0 if released.
;       Normal buttons do not generate events when released, so
;       SELECT will always be 1. However, toggle buttons (created by
;       parenting a button to an exclusive or non-exclusive base)
;       return separate events for the set and release actions.

;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro OnReferenceCustom, Event
@COMMON_DATAS
  custom=widget_info(event.top, find_by_uname='WID_BUTTON_ref_custom')
  re=widget_info(custom, /button_set)
  if re eq 1 then $
  begin
    Print, 'OnReferenceCustom'
    fname=dialog_pickfile(FILTER='*.lp', /READ, path=out_dir)
    if fname ne '' then $
    begin
       list_ref=widget_info(event.top, find_by_uname='WID_LIST_reference')
       print_reference, lp_ref, list_ref
       lp_ref=open_lp(fname)
       V0=V_from_lp(lp_ref)
       V1=V_from_ub(ub)
       if abs(V0-V1) gt 0.01 then re=dialog_message('Reference unit cell has volume different than from ub matrix')
    endif
  end
end
;-----------------------------------------------------------------
; Activate Button Callback Procedure.
; Argument:
;   Event structure:
;
;   {WIDGET_BUTTON, ID:0L, TOP:0L, HANDLER:0L, SELECT:0}
;
;   ID is the widget ID of the component generating the event. TOP is
;       the widget ID of the top level widget containing ID. HANDLER
;       contains the widget ID of the widget associated with the
;       handler routine.

;   SELECT is set to 1 if the button was set, and 0 if released.
;       Normal buttons do not generate events when released, so
;       SELECT will always be 1. However, toggle buttons (created by
;       parenting a button to an exclusive or non-exclusive base)
;       return separate events for the set and release actions.

;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro OnReferenceCubic, Event
@COMMON_DATAS
  custom=widget_info(event.top, find_by_uname='WID_BUTTON_ref_custom')
  re=widget_info(custom, /button_set)
  if re eq 0 then $
  begin
    Print, 'OnReferenceCubic'
    V0=V_from_UB(ub)
    a0=V0^(1.0/3)
    lp_ref=[a0,a0,a0,90.,90.,90.]
    list_ref=widget_info(event.top, find_by_uname='WID_LIST_reference')
    print_reference, lp_ref, list_ref
  end
end
;-----------------------------------------------------------------
; Kill Notify Callback Procedure.
; Argument:
;   wWidget - ID number of specific widget.
;
;
;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro OnDestroy, wWidget
@WID_GSE_ADA_COMMON
 Cij_window=-1
 print, 'OnDestroy'
end
;-----------------------------------------------------------------
; TLB_KILL_REQUEST_EVENTS Callback Procedure.
; Argument:
;   Event structure:
;
;   {WIDGET_KILL_REQUEST, ID:0L, TOP:0L, HANDLER:0L }
;
;   ID is the widget ID of the component generating the event. TOP is
;       the widget ID of the top level widget containing ID. HANDLER
;       contains the widget ID of the widget associated with the
;       handler routine.

;
;   Retrieve the IDs of other widgets in the widget hierarchy using
;       id=widget_info(Event.top, FIND_BY_UNAME=name)

;-----------------------------------------------------------------
pro OnKill, Event
@WID_GSE_ADA_COMMON
 Cij_window=-1
  print, 'OnKill'
end
