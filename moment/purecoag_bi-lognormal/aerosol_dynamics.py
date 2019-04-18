#file: aerosol_dynamics.py
import sys
import os
import random
import math

try:
    import Tkinter as tk
    import ttk
    import tkMessageBox as tkmb
    tk.messagebox = tkmb
except:
    import tkinter as tk
    import tkinter.messagebox
    import tkinter.ttk as ttk

import pygubu

CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))

class AerosolDynamics:

    def __init__(self, master):

        self.master = master

        self.about_dialog = None
        
        self.builder = builder = pygubu.Builder()
        builder.add_from_file(os.path.join(CURRENT_DIR, 'aerosol_dynamics.ui'))
        builder.add_resource_path(os.path.join(CURRENT_DIR, 'img'))
        
        self.mainwindow = b.get_object('mainwindow', self.master)
        self.mainmenu = b.get_object('mainmenu', self.mainwindow)
        self.btn_menu = b.get_object('btn_menu')
        
        self.mainwindow.protocol("WM_DELETE_WINDOW", self.quit)
        
        builder.connect_callbacks(self)

    def on_mainmenu_action(self, option_id=None):
        if option_id == 'mhelp_about':
            self.show_about_dialog()
        if option_id == 'mfile_quit':
            self.mainwindow.quit()
    
    def show_about_dialog(self):
        if self.about_dialog is None:
            dialog = self.builder.get_object('dlg_about', self.mainwindow)
            self.about_dialog = dialog
            
            def dialog_btnclose_clicked():
                dialog.close()
            
            btnclose = self.builder.get_object('btn_about_close')
            btnclose['command'] = dialog_btnclose_clicked
            
            dialog.run()
        else:
            self.about_dialog.show()
    
    def quit(self, event=None):
        self.mainwindow.quit()

    def run(self):
        self.mainwindow.mainloop()

if __name__ == '__main__':
    root = tk.Tk()
    app = AerosolDynamics(root)
    app.run()

