from tkinter import *
import Prime_Edit as PE
import Base_Edit as BE

window = Tk()
window.title("Single Base Editing Helper Tool")

def main():
    if clicked.get() == "Prime Editing Tool":
        PE.main(FASTAEntry.get(), mutationEntry.get(), filenameEntry.get(), pamEntry.get(), PBSlengthEntry.get())

    elif clicked.get() == "Base Editing Tool":
        pass
        # BE.main()


canvas = Canvas(window, height = 200, width = 600)
canvas.pack()

frame = Frame(window,relief = 'groove')
frame.place(relx = 0.1, rely = 0.1, relwidth = 0.8, relheight = 0.8)

welcome = Label(frame, text = "Welcome to the Single Base Editing Helper Tool", fg = "Black")
welcome.pack(side = "top")

options = ["--Select tool to use--", "Prime Editing Tool", "Base Editing Tool"]

clicked = StringVar()
clicked.set(options[0])
drop = OptionMenu(window, clicked, *options)
drop.pack(side = "top")

FASTAEntry = Entry(frame, width = 50)
FASTAEntry.pack(side = "top")
FASTAEntry.insert(0, "Please enter the DNA sequence")

mutationEntry = Entry(frame, width = 50)
mutationEntry.pack(side = "top")
mutationEntry.insert(0, "Please enter the desired mutation")

pamEntry = Entry(frame, width = 50)
pamEntry.pack(side = "top")
pamEntry.insert(0, "Please enter the desired PAM sequence if applicable")

PBSlengthEntry = Entry(frame, width = 50)
PBSlengthEntry.pack(side = "top")
PBSlengthEntry.insert(0, "Please enter PBS Length (from 7-17) if applicable")

filenameEntry = Entry(frame, width = 50)
filenameEntry.pack(side = "top")
filenameEntry.insert(0, "Please enter the desired .txt filename")

enterButton = Button(window, text = "Start", padx = 10, pady = 5, fg = "Black", bg = "gray", command = main)
enterButton.pack()

window.mainloop()
