import tkinter
import tkinter.filedialog
import inspect
import scipy.io

def _get_all_vars():
    for f in inspect.stack():
        if '__name__' in f[0].f_locals.keys():
            if f[0].f_locals['__name__'] == '__main__':
                return f[0].f_locals
    return None 

def export_variables():   
    v = _get_all_vars()
    assert (v is not None)
    n = v.keys()
    
    top = tkinter.Tk()  
    def done_fcn():
        #print("done")
        selected = listbox.curselection()
        selected_names = []
        if len(selected) > 0:
            mdict = {}
            for si in selected:
                curr_name = listbox.get(si)
                selected_names.append(curr_name)
                mdict[curr_name] = v[curr_name]
        
        #print(mdict)
        outputfile = tkinter.filedialog.asksaveasfilename(filetypes=[("MATLAB file", "*.mat")],defaultextension='.mat')
        if outputfile != "":
            print(outputfile)
            scipy.io.savemat(outputfile,mdict)
            print("Written", selected_names)
        top.destroy()
        
    top.geometry("200x250")  
    lbl = tkinter.Label(top,text = "Select variables to export")  
    listbox = tkinter.Listbox(top,selectmode=tkinter.MULTIPLE)  
    for vi in n:
        #print(vi)
        if vi.startswith("_"):
            continue
        if (vi in ['scipy','np','plt','simos','quit','In','Out','exit','get_ipython','tkinter','qu']):
            continue
        # Exclude modules, to get type, let's just use an imported module
        if type(v[vi]) == type(tkinter):
            continue
        listbox.insert(tkinter.END,vi)  
   
    btn = tkinter.Button(top, text = "done", command = done_fcn)  
    top.attributes('-topmost', True)
    top.update()
    lbl.pack()  
    listbox.pack()  
    btn.pack()  
    top.mainloop() 
