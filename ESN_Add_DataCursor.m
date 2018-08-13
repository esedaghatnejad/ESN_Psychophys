function ESN_Add_DataCursor(h_plot_, index_)

x_data_ = get(h_plot_,'XData')';
y_data_ = get(h_plot_,'YData');

if ((index_>length(x_data_))||(index_<0)); index_ = mod(index_, length(x_data_)); end;

pos = [x_data_(index_) y_data_(index_) 0];

hDataCursorMgr = datacursormode(ancestor(h_plot_,'figure'));
hDatatip = createDatatip(hDataCursorMgr, h_plot_);

hDatatip.Cursor.DataIndex = index_;
hDatatip.Cursor.Position = pos;

updateDataCursors(hDataCursorMgr);

end