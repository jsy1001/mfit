#!/usr/bin/env python
#
# $Id: fitgui.py,v 1.1 2003/06/13 17:25:43 jsy1001 Exp $

"""Graphical user interface for clfit.

Usage: fitgui &

"""


import string, sys, os, tempfile, popen2, fcntl, time
from Tkinter import *
from ScrolledText import ScrolledText
import tkFileDialog

_revision = string.split("$Revision: 1.1 $")[1]


class GUI:
    """Graphical user interface to clfit.

    Data attributes that may sensibly be used externally:

    exe (string) -- path to clfit executable

    device (string) -- PGPLOT device, passed to clfit

    fileName (Tkinter.StringVar) -- path to OI-FITS/Mapping Data file
    (use set() & get() methods)

    ChangeFileButton (Tkinter.Button) -- brings up data file dialog
    box (may wish to disable)
    
    initialdir (string) -- initial directory for file dialog boxes

    preFitCallback -- function to call just before running clfit. If this
    returns false, clfit won't be run

    postFitCallback -- function to call after running clfit (if successful)
    
    """

    def __init__(self, parent, dismissCommand=None):
        """Constructor.

        Arguments are:

        parent -- parent window

        dismissCommand -- function to call when Dismiss button clicked

        """
        # Initialise data attributes
        self.parent = parent
        self.exe = 'clfit'
        self.device = '/xserv'
        self.fileName = StringVar()
        self.fileName.set('(unset)')
        self.initialdir = os.getcwd()
        self.preFitCallback = None
        self.postFitCallback = None
        self.calErr = StringVar()
        self.calErr.set('0.0')
        self.cwl = StringVar()
        self.bw = StringVar()
        self.nofit = IntVar()
        self.nofit.set(0)
        self.plots = ['No plot', 'vis2', 't3amp', 't3phi']
        self.selPlot = StringVar()
        self.selPlot.set(self.plots[1])
        self.plotFrom = StringVar()
        self.plotTo = StringVar()

        # Initialise GUI elements
        fileFrame = Frame(parent)
        fileFrame.pack(side=TOP)
        Label(fileFrame, text='Data file:').pack(side=LEFT)
        Label(fileFrame, textvariable=self.fileName).pack(side=LEFT)
        self.ChangeFileButton = Button(fileFrame, text='Change',
                                       command=self._ChangeFileName)
        self.ChangeFileButton.pack(side=LEFT)
        calErrFrame = Frame(parent)
        calErrFrame.pack(side=TOP, fill=X, pady=4)
        Label(calErrFrame, text='Calibration Error (extra frac. error in system vis.)').pack(side=LEFT, anchor=W)
        Entry(calErrFrame, textvariable=self.calErr, width=5).pack(
            side=LEFT, anchor=W, padx=4)
        wbFrame = Frame(parent)
        wbFrame.pack(side=TOP, fill=X, pady=4)
        Label(wbFrame, text='Waveband:').pack(side=LEFT, anchor=W)
        Entry(wbFrame, textvariable=self.cwl, width=5).pack(side=LEFT,
                                                            anchor=W, padx=4)
        Entry(wbFrame, textvariable=self.bw, width=5).pack(side=LEFT,
                                                           anchor=W, padx=4)
        Label(parent, text='Model:').pack(side=TOP, anchor=W)
        self.ModelText = ScrolledText(parent, height=19, width=40,
                                      font=('Helvetica', 10))
        self.ModelText.pack(side=TOP, expand=1, fill=BOTH)
        midFrame1 = Frame(parent)
        midFrame1.pack(side=TOP, fill=X, pady=4)
        Label(midFrame1, text='Plot:').pack(side=LEFT, anchor=NW)
        plotFrame = Frame(midFrame1)
        plotFrame.pack(side=LEFT)
        for p in self.plots:
            Radiobutton(plotFrame, text=p, variable=self.selPlot,
                        value=p).pack(side=TOP, anchor=W)
        Label(midFrame1, text='From:').pack(side=LEFT)
        Entry(midFrame1, textvariable=self.plotFrom, width=5).pack(side=LEFT)
        Label(midFrame1, text='To:').pack(side=LEFT)
        Entry(midFrame1, textvariable=self.plotTo, width=5).pack(side=LEFT)
        Button(midFrame1, text='Go', command=self.Go).pack(side=RIGHT,
                                                           anchor=NE, padx=4)
        Button(midFrame1, text='Save model',
               command=self.SaveModel).pack(side=RIGHT, anchor=NE, padx=4)
        Button(midFrame1, text='Load model',
               command=self.LoadModel).pack(side=RIGHT, anchor=NE, padx=4)
        Checkbutton(parent, text="Don't fit (report goodness-of-fit only)",
                    variable=self.nofit).pack(side=TOP, anchor=W, pady=4)
        midFrame2 = Frame(parent)
        midFrame2.pack(side=TOP, fill=X)
        Label(midFrame2, text='Results:').pack(side=LEFT, anchor=SW)
        if dismissCommand is None: dismissCommand = parent.quit
        Button(midFrame2, text='Dismiss',
               command=dismissCommand).pack(side=RIGHT, padx=4, pady=4)
        Button(midFrame2, text='Clear results',
               command=self.ClearResults).pack(side=RIGHT, padx=4, pady=4)
        self.Results = ScrolledText(parent, height=31, width=90,
                                    font=('Courier', 10), state=DISABLED)
        self.Results.tag_config('result', foreground='#1e90ff') # dodger blue
        self.Results.tag_config('commentary', foreground='#ff8c00') # dark orange
        self.Results.tag_config('error', foreground='#8b0000') # dark red
        self.Results.pack(side=TOP, expand=1, fill=BOTH)

    def LoadModel(self):
        """Get filename and read model from file."""
        fileName = tkFileDialog.askopenfilename(parent=self.parent,
            initialdir=self.initialdir,
            filetypes=[('mfit model files','*.model'),
                       ('All files','*')])
        if fileName != '': self.ReadModel(fileName)

    def ReadModel(self, fileName):
        """Read model from file."""
        try:
            fil = file(fileName, 'r')
            text = fil.read()
            fil.close()
        except IOError, (errNo, errStr):
            self.ShowResult('Error reading %s: %s\n' % (fileName, errStr),
                            'error')
        else:
            self.SetModel(text)

    def SetModel(self, text):
        """Set model text."""
        self.ModelText.delete(1.0, END)
        self.ModelText.insert(END, text)

    def ClearResults(self):
        """Clear results window."""
        self.Results.configure(state=NORMAL)
        self.Results.delete(1.0, END)
        self.Results.configure(state=DISABLED)

    def SaveModel(self):
        """Get filename and write model to file."""
        fileName = tkFileDialog.asksaveasfilename(parent=self.parent,
            initialdir=self.initialdir,
            filetypes=[('mfit model files','*.model'),
                       ('All files','*')])
        if fileName != '': self.WriteModel(fileName)

    def WriteModel(self, fileName):
        """Write model text to file."""
        fil = file(fileName, 'w')
        fil.write(self.ModelText.get(1.0, END))
        fil.close()

    def _ChangeFileName(self):
        newName = tkFileDialog.askopenfilename(
            parent=self.parent,
            initialdir=self.initialdir,
            title='Choose data file for fit',
            filetypes=[('(OI-)FITS files','*.fits'), ('COAST Mapping Data files', '*.mapdat'),
                       ('All files','*')])
        if newName != '': self.fileName.set(newName)

    def ShowResult(self, text, tag):
        """Display text in 'Results'."""
        self.Results.configure(state=NORMAL)
        self.Results.insert(END, text, tag)
        self.Results.yview(END)
        self.Results.configure(state=DISABLED)
        self.parent.update_idletasks()

    def Go(self):
        """Fit the current model."""
        # Execute pre-callback
        if callable(self.preFitCallback):
            if not self.preFitCallback(): return
        # Write model text to tempfile
        self._tempName = tempfile.mktemp('.model')
        self.WriteModel(self._tempName)
        # Run clfit so we can grab its output
        optText = ''
        try:
            c = float(self.calErr.get())
        except ValueError:
            pass
        else:
            optText += ' --calerr %.3f' % c
        try:
            cwl = float(self.cwl.get())
            bw = float(self.bw.get())
        except ValueError:
            pass
        else:
            optText += ' --waveband %.2f %.2f' % (cwl, bw)
        p = self.selPlot.get()
        if p != self.plots[0]: # not 'No plot'
            try:
                xmin = float(self.plotFrom.get())
                xmax = float(self.plotTo.get())
            except ValueError:
                optText += ' --plot %s' % p
            else:
                optText += ' --zoomplot %s %.3f %.3f' % (p, xmin, xmax)
        if self.nofit.get(): optText += ' --nofit'
        command = '%s%s --device %s %s %s' % (
            self.exe, optText, self.device, self.fileName.get(),
            self._tempName)
        self.ShowResult('Running %s:\n' % command, tag='commentary')
        self._popen4 = popen2.Popen4(command, 10)
        fd = self._popen4.fromchild.fileno() # stdout & stderr of child process
        oldflags = fcntl.fcntl(fd, fcntl.F_GETFL)
        fcntl.fcntl(fd, fcntl.F_SETFL, oldflags|os.O_NONBLOCK)
        tkinter.createfilehandler(fd, READABLE, self._DisplayOutput)

    def _DisplayOutput(self, fd, mask):
        """Read and display output from child process."""
        self.ShowResult(self._popen4.fromchild.read(), 'result')
        status = self._popen4.poll()
        if status != -1:
            # child process completed
            time.sleep(0.5)
            self.ShowResult(self._popen4.fromchild.read(), 'result')
            if os.WEXITSTATUS(status) != 0:
                self.ShowResult('Subprocess exited with code %d\n' \
                                % os.WEXITSTATUS(status), 'error')
            else:
                self.ShowResult('Subprocess exited normally\n', 'commentary')
            tkinter.deletefilehandler(fd)
            os.remove(self._tempName)
            # Execute post-callback
            if callable(self.postFitCallback):
                self.postFitCallback()
    

def _main():
    """Main routine."""

    # Find out if Python Megawidgets are installed
    try:
        import Pmw
    except ImportError:
        havePmw = 0
    else:
        havePmw = 1

    if len(sys.argv) == 1:
        # No command-line arguments - run graphical user interface
        root = Tk()
        if havePmw:
            Pmw.initialise(root, fontScheme='pmw1')
        root.title('fitgui %s' % _revision)
        main = GUI(root)
        #main.fileName.set('/net/oberon/home/jsy1001/reductions/alp_ori_02/alp_ori_0203-04_corr2.mapdat')
        #main.ReadModel('/net/oberon/home/jsy1001/reductions/alp_ori_02/alpOri.model')
        root.mainloop()
    else:
        # Too many arguments
        print "%s %s\n" % (sys.argv[0], _revision)
        print __doc__
        sys.exit(2)

if __name__ == '__main__':
    _main()
 
# Local Variables:
# mode: python
# End:
