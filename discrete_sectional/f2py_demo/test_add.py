import add_fortran_mod
from tkinter import *
from tkinter import ttk

class ModelApp:

    def __init__(self, master):

        self.x = ttk.Label(master, text = "X")
        self.x.grid(row = 0, column = 0, columnspan = 1)
        self.x_entry = ttk.Entry(master, width = 10)
        self.x_entry.grid(row = 0, column = 1, columnspan = 1)

        self.y = ttk.Label(master, text = "Y")
        self.y.grid(row = 0, column = 2, columnspan = 1)
        self.y_entry = ttk.Entry(master, width = 10)
        self.y_entry.grid(row = 0, column = 3, columnspan = 1)
        
        self.add = ttk.Button(master, text = "Add", command = self.addRun)
        self.add.grid(row = 1, column = 0, columnspan = 2)

        self.sum = ttk.Label(master)
        self.sum.grid(row = 1, column = 2, columnspan = 2)

    def addRun(self):

        print(add_fortran_mod.mod.__doc__)
        add_fortran_mod.mod.x = self.x_entry.get()
        add_fortran_mod.mod.y = self.y_entry.get()
        add_fortran_mod.mod.addtwonumbers()
        s = add_fortran_mod.mod.z
        self.sum.config(text = 'The sum is ' + str(s))
            
def main():            
    
    root = Tk()
    model_app = ModelApp(root)
    root.mainloop()
    
if __name__ == "__main__":
    main()
