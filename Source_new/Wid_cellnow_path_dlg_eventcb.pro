;
; IDL Event Callback Procedures
; Wid_cellnow_path_dlg_eventcb
;
; Generated on:	05/19/2015 16:43.42
;
;
; Empty stub procedure used for autoloading.
;
pro Wid_cellnow_path_dlg_eventcb
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
pro Browse_CellNowDir, Event
	file = DIALOG_PICKFILE (title='Select Path to Cell_now.exe', filter='cell_now.exe', get_path=path)
	id= widget_info (event.top,FIND_BY_UNAME='WID_TEXT_CellnowPath')
	widget_control, id, set_value=path

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
pro CloseUp, Event
	@common_datas
	id=widget_info( event.top, find_by_uname='WID_TEXT_CellnowPath')
	widget_control, id, get_value=path
	re=file_info(path+'cell_now.exe')
	if re.exists then  $
	begin
	   cell_now_dir = path
       widget_control, event.top,/destroy
	endif else mi=dialog_message('cell_now.exe is not present in the selected directory')
end

pro Cancel_action, Event
        widget_control, event.top,/destroy
end