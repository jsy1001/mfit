#!/usr/bin/env python
#
# $Id: fitgui.py,v 1.11 2006/08/18 16:22:19 jsy1001 Exp $

"""Graphical user interface for clfit.

Usage: fitgui &

"""


import string, sys, os, tempfile, popen2, fcntl, time
from Tkinter import *
from ScrolledText import ScrolledText
import tkFileDialog

_revision = string.split("$Revision: 1.11 $")[1]


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
        self.wlmin = StringVar()
        self.wlmax = StringVar()
        self.target_id = StringVar()
        self.nofit = IntVar()
        self.nofit.set(0)
        self.plots = ['No plot', 'uv', 'vis2', 't3amp', 't3phi',
                      'vis2wl', 't3ampwl', 't3phiwl',
                      'vis2mjd', 't3ampmjd', 't3phimjd', 'post', 'mpost',
                      'post2d', 'mpost2d']
        self.selPlot = StringVar()
        self.selPlot.set(self.plots[1])
        self.plotXIndex = StringVar()
        self.plotXIndex.set('1')
        self.plotYIndex = StringVar()
        self.plotYIndex.set('2')
        self.plotXFrom = StringVar()
        self.plotXTo = StringVar()
        self.plotYFrom = StringVar()
        self.plotYTo = StringVar()
        self.margErr = IntVar()
        self.margErr.set(0)
        self.margErrVar = StringVar()
        self.margErrVar.set('1')

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
        Label(wbFrame, text='or Wavelength range:').pack(side=LEFT, anchor=W)
        Entry(wbFrame, textvariable=self.wlmin,
              width=5).pack(side=LEFT, anchor=W, padx=4)
        Entry(wbFrame, textvariable=self.wlmax,
              width=5).pack(side=LEFT, anchor=W, padx=4)
        targetFrame = Frame(parent)
        targetFrame.pack(side=TOP, fill=X, pady=4)
        Label(targetFrame, text='TARGET_ID (blank to use 1st in OI_TARGET table):').pack(side=LEFT, anchor=W)
        Entry(targetFrame, textvariable=self.target_id,
              width=5).pack(side=LEFT, anchor=W, padx=4)
        Label(parent, text='Model:').pack(side=TOP, anchor=W)
        self.ModelText = ScrolledText(parent, height=19, width=40,
                                      font=('Helvetica', 10))
        self.ModelText.pack(side=TOP, expand=1, fill=BOTH)
        midFrame1 = Frame(parent)
        midFrame1.pack(side=TOP, fill=X, pady=4)
        Label(midFrame1, text='Plot:').pack(side=LEFT, anchor=NW)
        plotFrame = Frame(midFrame1)
        plotFrame.pack(side=LEFT)
        for i in range(len(self.plots)):
            p = self.plots[i]
            Radiobutton(plotFrame, text=p, variable=self.selPlot,
                        value=p).grid(row=(i+1)/2, column=(i+1)%2, sticky=W)
        Entry(plotFrame, textvariable=self.plotXIndex,
              width=3).grid(row=len(self.plots)/2-1, column=2)
        Entry(plotFrame, textvariable=self.plotYIndex,
              width=3).grid(row=len(self.plots)/2, column=2)
        rangeFrame = Frame(midFrame1)
        rangeFrame.pack(side=LEFT)
        Label(rangeFrame, text='X From:').grid(row=0, column=0, sticky=E)
        Entry(rangeFrame, textvariable=self.plotXFrom, width=5).grid(row=0, column=1)
        Label(rangeFrame, text='To:').grid(row=0, column=2)
        Entry(rangeFrame, textvariable=self.plotXTo, width=5).grid(row=0, column=3)
        Label(rangeFrame, text='[(m)post2d only] Y From:').grid(row=1, column=0, sticky=E)
        Entry(rangeFrame, textvariable=self.plotYFrom, width=5).grid(row=1, column=1)
        Label(rangeFrame, text='To:').grid(row=1, column=2)
        Entry(rangeFrame, textvariable=self.plotYTo, width=5).grid(row=1, column=3)
        Button(midFrame1, text='Go', command=self.Go).pack(side=RIGHT,
                                                           anchor=NE, padx=4)
        Button(midFrame1, text='Save model',
               command=self.SaveModel).pack(side=RIGHT, anchor=NE, padx=4)
        Button(midFrame1, text='Load model',
               command=self.LoadModel).pack(side=RIGHT, anchor=NE, padx=4)
        midFrame2 = Frame(parent)
        midFrame2.pack(side=TOP, fill=X, pady=4)
        Checkbutton(midFrame2, text="Don't fit (report goodness-of-fit only)",
                    variable=self.nofit).pack(side=LEFT, anchor=W, padx=8)
        Entry(midFrame2, textvariable=self.margErrVar,
              width=5).pack(side=LEFT, anchor=W)
        Checkbutton(midFrame2, text="Error bar by marginalising",
                    variable=self.margErr).pack(side=LEFT, anchor=W)
        midFrame3 = Frame(parent)
        midFrame3.pack(side=TOP, fill=X)
        Label(midFrame3, text='Results:').pack(side=LEFT, anchor=SW)
        if dismissCommand is None: dismissCommand = parent.quit
        Button(midFrame3, text='Dismiss',
               command=dismissCommand).pack(side=RIGHT, padx=4, pady=4)
        Button(midFrame3, text='Clear results',
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
            filetypes=[('(OI-)FITS files','*fits'), ('COAST Mapping Data files', '*.mapdat'),
                       ('wbCalib / nbCalib files','*calib'), ('All files','*')])
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
            try:
                wlmin = float(self.wlmin.get())
                wlmax = float(self.wlmax.get())
                optText += ' --waverange %.2f %.2f' % (wlmin, wlmax)
            except ValueError:
                pass # don't specify waveband(s)
        else:
            optText += ' --waveband %.2f %.2f' % (cwl, bw)
        try:
            target_id = int(self.target_id.get())
            optText += ' --target_id %d' % target_id
        except ValueError:
            pass
        p = self.selPlot.get()
        if p != self.plots[0]: # not 'No plot'
            if p == 'post' or p == 'mpost':
                try:
                    index = int(self.plotXIndex.get())
                except ValueError:
                    index = 1
            if p == 'post2d' or p == 'mpost2d':
                try:
                    indx = (int(self.plotXIndex.get()),
                            int(self.plotYIndex.get()))
                except ValueError:
                    indx = (1, 2)
            try:
                xmin = float(self.plotXFrom.get())
                xmax = float(self.plotXTo.get())
                if p[-2:] == '2d':
                    ymin = float(self.plotYFrom.get())
                    ymax = float(self.plotYTo.get())
            except ValueError:
                if p == 'post' or p == 'mpost':
                    optText += ' --plot %s %d' % (p, index)
                elif p == 'post2d' or p == 'mpost2d':
                    optText += ' --plot %s %d %d' % (p, indx[0], indx[1])
                else:
                    optText += ' --plot %s' % p
            else:
                if p == 'post' or p == 'mpost':
                    optText += ' --zoomplot %s %d %.4f %.4f' % (p, index, xmin, xmax)
                elif p == 'post2d' or p == 'mpost2d':
                    optText += ' --zoomplot %s %d %d %.4f %.4f  %.4f %.4f' % \
                               (p, indx[0], indx[1], xmin, xmax, ymin, ymax)
                else:
                    optText += ' --zoomplot %s %.4f %.4f' % (p, xmin, xmax)
        if self.nofit.get(): optText += ' --nofit'
        if self.margErr.get():
            optText += ' --margerr %s' % (self.margErrVar.get())
        command = '%s%s --device %s %s %s' % (
            self.exe, optText, self.device, self.fileName.get(),
            self._tempName)
        self.ShowResult('Running %s:\n' % command, tag='commentary')
        self._popen4 = popen2.Popen4(command, 10)
        fd = self._popen4.fromchild.fileno() # stdout & stderr of child process
        oldflags = fcntl.fcntl(fd, fcntl.F_GETFL)
        fcntl.fcntl(fd, fcntl.F_SETFL, oldflags|os.O_NONBLOCK)
        self.parent.tk.createfilehandler(fd, READABLE, self._DisplayOutput)

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
            self.parent.tk.deletefilehandler(fd)
            os.remove(self._tempName)
            # Execute post-callback
            if callable(self.postFitCallback):
                self.postFitCallback()
    

def _main(altExe=None):
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
        if altExe is not None: main.exe = altExe
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
