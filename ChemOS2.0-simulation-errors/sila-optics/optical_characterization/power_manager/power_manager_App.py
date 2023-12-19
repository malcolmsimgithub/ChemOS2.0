try:
    from power_manager import PowerManager
except:
    from .power_manager import PowerManager

import tkinter as tk
from tkinter import messagebox, ttk

column_width = 200

class PowerManagerApp(tk.Frame):

    def __init__(self, master = None):

        self.PUSB = PowerManager()
        self.dev_num = max([val['PowerUSB'] for val in self.PUSB.dev_list.values()]) + 1
        # self.dev_num = self.PUSB.PUSB.getdevicenumber() # number of PowerUSB

        super().__init__(master, height = 200, width =column_width * self.dev_num)
        self.pack()

        self.master.title("PowerManager")

        self.color_on = '#00ff00'
        self.color_off = '#ff0000'

        self.dev_state = {}

        for key in self.PUSB.dev_list.keys():
            self.dev_state[key] = tk.IntVar()
            self.dev_state[key].set(self.PUSB.read_state(key))

        self.create_widgets()


    def create_widgets(self):

        x_offset, y_offset = 20, 0

        self.status_label = {}
        self.button_on = {}
        self.button_off = {}

        for key, val in self.PUSB.dev_list.items():

            if val['port']  == 1:      
                text = "PowerUSB_%s" %val['PowerUSB']
                label = tk.Label(self,  text = text, anchor = 'w').place(x = x_offset + val['PowerUSB'] * column_width, y = y_offset + 20)
            
            self.status_label[key] = tk.Label(self, borderwidth = 3, relief = 'sunken', bg=self.color_off, width = 2)
            self.status_label[key].place(x = x_offset + val['PowerUSB'] * column_width, y = y_offset + val['port']* 40 + 20)
            if self.dev_state[key].get() == 1:
                self.status_label[key]['bg'] = self.color_on

            label = tk.Label(self,  text = key, anchor = 'w').place(x = x_offset + val['PowerUSB'] * column_width + 30, y = y_offset + val['port']* 40 + 20)

            self.button_on[key] = tk.Button(self, text = 'on')
            self.button_on[key].place(x = x_offset + val['PowerUSB'] * column_width + 90, y = y_offset + val['port']* 40 + 20, height = 25, width = 25)
            self.button_on[key].bind("<ButtonPress>", self.button_on_click)

            self.button_off[key] = tk.Button(self, text = 'off')
            self.button_off[key].place(x = x_offset + val['PowerUSB'] * column_width + 120, y = y_offset + val['port']* 40 + 20, height = 25, width = 25)
            self.button_off[key].bind("<ButtonPress>", self.button_off_click)


    def button_on_click(self, event):
        if event.widget.cget("state") == 'normal':
            device_num = self.get_device(event)
            dev = list(self.PUSB.dev_list.keys())[device_num]
            self.PUSB.set_state(dev, 'on')
            self.status_label[dev]['bg'] = self.color_on
            self.dev_state[dev].set(1)


    def button_off_click(self, event):
        if event.widget.cget("state") == 'normal':
            device_num = self.get_device(event)
            dev = list(self.PUSB.dev_list.keys())[device_num]
            self.PUSB.set_state(dev, 'off')
            self.status_label[dev]['bg'] = self.color_off
            self.dev_state[dev].set(0)

    def get_device(self, event):
        num = event.widget.winfo_name()
        num = ''.join(filter(lambda i: i.isdigit(), num)) 
        if num.isdigit() == False:
            num = 1
        return (int(num)-1)//2


root = tk.Tk()
main = PowerManagerApp(master=root)
main.mainloop()