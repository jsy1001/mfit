#!/usr/bin/env python

# Copyright (C) 2003-2018, 2022 John Young, Matthew Worsley
#
# This file is part of mfit.
#
# Mfit is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/ .

"""Graphical user interface for clfit.

Usage: fitgui &

"""

import sys
import os
import tempfile
import time
from subprocess import Popen, PIPE, STDOUT
import tkinter as tk
import tkinter.filedialog as tk_filedialog
from tkinter.scrolledtext import ScrolledText

_revision = "$Revision: 1.14 $".split()[1]


class GUI:
    """Graphical user interface to clfit.

    Data attributes that may sensibly be used externally:

    exe (string) -- path to clfit executable

    device (string) -- PGPLOT device, passed to clfit

    fileName (tkinter.StringVar) -- path to OI-FITS/Mapping Data file
    (use set() & get() methods)

    ChangeFileButton (tkinter.Button) -- brings up data file dialog
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
        self.exe = "clfit"
        self.device = "/xserv"
        self.fileName = tk.StringVar()
        self.fileName.set("(unset)")
        self.initialdir = os.getcwd()
        self.preFitCallback = None
        self.postFitCallback = None
        self.calErr = tk.StringVar()
        self.calErr.set("0.0")
        self.cwl = tk.StringVar()
        self.bw = tk.StringVar()
        self.wlmin = tk.StringVar()
        self.wlmax = tk.StringVar()
        self.target_id = tk.StringVar()
        self.nofit = tk.IntVar()
        self.nofit.set(0)
        self.plots = [
            "No plot",
            "uv",
            "vis2",
            "t3amp",
            "t3phi",
            "vis2-wl",
            "t3amp-wl",
            "t3phi-wl",
            "vis2-mjd",
            "t3amp-mjd",
            "t3phi-mjd",
            "vis2-st",
            "t3amp-st",
            "t3phi-st",
            "vis2-pa",
            "post",
            "mpost",
            "post2d",
            "mpost2d",
        ]
        self.selPlot = tk.StringVar()
        self.selPlot.set(self.plots[1])
        self.plotXIndex = tk.StringVar()
        self.plotXIndex.set("1")
        self.plotYIndex = tk.StringVar()
        self.plotYIndex.set("2")
        self.plotXFrom = tk.StringVar()
        self.plotXTo = tk.StringVar()
        self.plotYFrom = tk.StringVar()
        self.plotYTo = tk.StringVar()
        self.margErr = tk.IntVar()
        self.margErr.set(0)
        self.margErrVar = tk.StringVar()
        self.margErrVar.set("1")

        # Initialise GUI elements
        fileFrame = tk.Frame(parent)
        fileFrame.pack(side=tk.TOP)
        tk.Label(fileFrame, text="Data file:").pack(side=tk.LEFT)
        tk.Label(fileFrame, textvariable=self.fileName).pack(side=tk.LEFT)
        self.ChangeFileButton = tk.Button(
            fileFrame, text="Change", command=self._ChangeFileName
        )
        self.ChangeFileButton.pack(side=tk.LEFT)
        calErrFrame = tk.Frame(parent)
        calErrFrame.pack(side=tk.TOP, fill=tk.X, pady=4)
        tk.Label(
            calErrFrame, text="Calibration Error (extra frac. error in system vis.)"
        ).pack(side=tk.LEFT, anchor=tk.W)
        tk.Entry(calErrFrame, textvariable=self.calErr, width=5).pack(
            side=tk.LEFT, anchor=tk.W, padx=4
        )
        wbFrame = tk.Frame(parent)
        wbFrame.pack(side=tk.TOP, fill=tk.X, pady=4)
        tk.Label(wbFrame, text="Waveband:").pack(side=tk.LEFT, anchor=tk.W)
        tk.Entry(wbFrame, textvariable=self.cwl, width=5).pack(
            side=tk.LEFT, anchor=tk.W, padx=4
        )
        tk.Entry(wbFrame, textvariable=self.bw, width=5).pack(
            side=tk.LEFT, anchor=tk.W, padx=4
        )
        tk.Label(wbFrame, text="or Wavelength range:").pack(side=tk.LEFT, anchor=tk.W)
        tk.Entry(wbFrame, textvariable=self.wlmin, width=5).pack(
            side=tk.LEFT, anchor=tk.W, padx=4
        )
        tk.Entry(wbFrame, textvariable=self.wlmax, width=5).pack(
            side=tk.LEFT, anchor=tk.W, padx=4
        )
        targetFrame = tk.Frame(parent)
        targetFrame.pack(side=tk.TOP, fill=tk.X, pady=4)
        tk.Label(
            targetFrame, text="TARGET_ID (blank to use 1st in OI_TARGET table):"
        ).pack(side=tk.LEFT, anchor=tk.W)
        tk.Entry(targetFrame, textvariable=self.target_id, width=5).pack(
            side=tk.LEFT, anchor=tk.W, padx=4
        )
        tk.Label(parent, text="Model:").pack(side=tk.TOP, anchor=tk.W)
        self.ModelText = ScrolledText(
            parent, height=19, width=40, font=("Helvetica", 10)
        )
        self.ModelText.pack(side=tk.TOP, expand=1, fill=tk.BOTH)
        midFrame1 = tk.Frame(parent)
        midFrame1.pack(side=tk.TOP, fill=tk.X, pady=4)
        tk.Label(midFrame1, text="Plot:").pack(side=tk.LEFT, anchor=tk.NW)
        plotFrame = tk.Frame(midFrame1)
        plotFrame.pack(side=tk.LEFT)
        ncol = 3
        for i in range(len(self.plots)):
            p = self.plots[i]
            tk.Radiobutton(plotFrame, text=p, variable=self.selPlot, value=p).grid(
                row=int((i + 1) / ncol), column=(i + 1) % ncol, sticky=tk.W
            )
        tk.Entry(plotFrame, textvariable=self.plotXIndex, width=3).grid(
            row=int(len(self.plots) / ncol) - 1, column=ncol
        )
        tk.Entry(plotFrame, textvariable=self.plotYIndex, width=3).grid(
            row=int(len(self.plots) / ncol), column=ncol
        )
        rangeFrame = tk.Frame(midFrame1)
        rangeFrame.pack(side=tk.LEFT)
        tk.Label(rangeFrame, text="X From:").grid(row=0, column=0, sticky=tk.E)
        tk.Entry(rangeFrame, textvariable=self.plotXFrom, width=5).grid(row=0, column=1)
        tk.Label(rangeFrame, text="To:").grid(row=0, column=2)
        tk.Entry(rangeFrame, textvariable=self.plotXTo, width=5).grid(row=0, column=3)
        tk.Label(rangeFrame, text="Y From:").grid(row=1, column=0, sticky=tk.E)
        tk.Entry(rangeFrame, textvariable=self.plotYFrom, width=5).grid(row=1, column=1)
        tk.Label(rangeFrame, text="To:").grid(row=1, column=2)
        tk.Entry(rangeFrame, textvariable=self.plotYTo, width=5).grid(row=1, column=3)
        tk.Label(rangeFrame, text="[Y for (m)post2d only]").grid(row=2, columnspan=4)
        tk.Button(midFrame1, text="Go", command=self.Go).pack(
            side=tk.RIGHT, anchor=tk.NE, padx=4
        )
        tk.Button(midFrame1, text="Save model", command=self.SaveModel).pack(
            side=tk.RIGHT, anchor=tk.NE, padx=4
        )
        tk.Button(midFrame1, text="Load model", command=self.LoadModel).pack(
            side=tk.RIGHT, anchor=tk.NE, padx=4
        )
        midFrame2 = tk.Frame(parent)
        midFrame2.pack(side=tk.TOP, fill=tk.X, pady=4)
        tk.Checkbutton(
            midFrame2,
            text="Don't fit (report goodness-of-fit only)",
            variable=self.nofit,
        ).pack(side=tk.LEFT, anchor=tk.W, padx=8)
        tk.Entry(midFrame2, textvariable=self.margErrVar, width=5).pack(
            side=tk.LEFT, anchor=tk.W
        )
        tk.Checkbutton(
            midFrame2, text="Error bar by marginalising", variable=self.margErr
        ).pack(side=tk.LEFT, anchor=tk.W)
        midFrame3 = tk.Frame(parent)
        midFrame3.pack(side=tk.TOP, fill=tk.X)
        tk.Label(midFrame3, text="Results:").pack(side=tk.LEFT, anchor=tk.SW)
        if dismissCommand is None:
            dismissCommand = parent.quit
        tk.Button(midFrame3, text="Dismiss", command=dismissCommand).pack(
            side=tk.RIGHT, padx=4, pady=4
        )
        tk.Button(midFrame3, text="Clear results", command=self.ClearResults).pack(
            side=tk.RIGHT, padx=4, pady=4
        )
        self.Results = ScrolledText(
            parent, height=31, width=90, font=("Courier", 10), state=tk.DISABLED
        )
        self.Results.tag_config("result", foreground="#1e90ff")  # dodger blue
        self.Results.tag_config("commentary", foreground="#ff8c00")  # dark orange
        self.Results.tag_config("error", foreground="#8b0000")  # dark red
        self.Results.pack(side=tk.TOP, expand=1, fill=tk.BOTH)

    def LoadModel(self):
        """Get filename and read model from file."""
        fileName = tk_filedialog.askopenfilename(
            parent=self.parent,
            initialdir=self.initialdir,
            filetypes=[("mfit model files", "*.model"), ("All files", "*")],
        )
        if fileName != "":
            self.ReadModel(fileName)

    def ReadModel(self, fileName):
        """Read model from file."""
        try:
            fil = open(fileName, "r")
            text = fil.read()
            fil.close()
        except IOError as e:
            (errNo, errStr) = e.args
            self.ShowResult("Error reading %s: %s\n" % (fileName, errStr), "error")
        else:
            self.SetModel(text)

    def SetModel(self, text):
        """Set model text."""
        self.ModelText.delete(1.0, tk.END)
        self.ModelText.insert(tk.END, text)

    def ClearResults(self):
        """Clear results window."""
        self.Results.configure(state=tk.NORMAL)
        self.Results.delete(1.0, tk.END)
        self.Results.configure(state=tk.DISABLED)

    def SaveModel(self):
        """Get filename and write model to file."""
        fileName = tk_filedialog.asksaveasfilename(
            parent=self.parent,
            initialdir=self.initialdir,
            filetypes=[("mfit model files", "*.model"), ("All files", "*")],
        )
        if fileName != "":
            self.WriteModel(fileName)

    def WriteModel(self, fileName):
        """Write model text to file."""
        fil = open(fileName, "w")
        fil.write(self.ModelText.get(1.0, tk.END))
        fil.close()

    def _ChangeFileName(self):
        newName = tk.filedialog.askopenfilename(
            parent=self.parent,
            initialdir=self.initialdir,
            title="Choose data file for fit",
            filetypes=[
                ("(OI-)FITS files", "*fits"),
                ("COAST Mapping Data files", "*.mapdat"),
                ("wbCalib / nbCalib files", "*calib"),
                ("All files", "*"),
            ],
        )
        if newName != "":
            self.fileName.set(newName)

    def ShowResult(self, text, tag):
        """Display text in 'Results'."""
        self.Results.configure(state=tk.NORMAL)
        self.Results.insert(tk.END, text, tag)
        self.Results.yview(tk.END)
        self.Results.configure(state=tk.DISABLED)
        self.parent.update_idletasks()

    def Go(self):
        """Fit the current model."""
        # Execute pre-callback
        if callable(self.preFitCallback):
            if not self.preFitCallback():
                return
        # Write model text to tempfile
        self._tempName = tempfile.mktemp(".model")
        self.WriteModel(self._tempName)
        # Run clfit so we can grab its output
        args = ["nice", self.exe, "--device", self.device]
        try:
            c = float(self.calErr.get())
        except ValueError:
            pass
        else:
            args += f"--calerr {c:.3f}".split()
        try:
            cwl = float(self.cwl.get())
            bw = float(self.bw.get())
        except ValueError:
            try:
                wlmin = float(self.wlmin.get())
                wlmax = float(self.wlmax.get())
                args += f"--waverange {wlmin:.2f} {wlmax:.2f}".split()
            except ValueError:
                pass  # don't specify waveband(s)
        else:
            args += f"--waveband {cwl:.2f} {bw:.2f}".split()
        try:
            target_id = int(self.target_id.get())
            args += f"--target_id {target_id}".split()
        except ValueError:
            pass
        p = self.selPlot.get()
        if p != self.plots[0]:  # not 'No plot'
            if p == "post" or p == "mpost":
                try:
                    index = int(self.plotXIndex.get())
                except ValueError:
                    index = 1
            if p == "post2d" or p == "mpost2d":
                try:
                    indx = (int(self.plotXIndex.get()), int(self.plotYIndex.get()))
                except ValueError:
                    indx = (1, 2)
            try:
                xmin = float(self.plotXFrom.get())
                xmax = float(self.plotXTo.get())
                if p[-2:] == "2d":
                    ymin = float(self.plotYFrom.get())
                    ymax = float(self.plotYTo.get())
            except ValueError:
                if p == "post" or p == "mpost":
                    args += f"--plot {p} {index}".split()
                elif p == "post2d" or p == "mpost2d":
                    args += f"--plot {p} {indx[0]} {indx[1]}".split()
                else:
                    args += f"--plot {p}".split()
            else:
                if p == "post" or p == "mpost":
                    args += f"--zoomplot {p} {index} {xmin} {xmax}".split()
                elif p == "post2d" or p == "mpost2d":
                    args += (
                        f"--zoomplot {p} {indx[0]} {indx[1]} {xmin:.4f} {xmax:.4f}"
                        f" {ymin:.4f} {ymax:.4f}".split()
                    )
                else:
                    args += f"--zoomplot {p} {xmin:.4f} {xmax:.4f}".split()
        if self.nofit.get():
            args += ["--nofit"]
        if self.margErr.get():
            args += f"--margerr {self.margErrVar.get()}".split()
        args += [self.fileName.get(), self._tempName]
        self.ShowResult("Running %s:\n" % " ".join(args), tag="commentary")

        # https://gitpress.io/u/1282/tkinter-read-async-subprocess-output
        self._proc = Popen(args, bufsize=10, stdout=PIPE, stderr=STDOUT, close_fds=True)
        self.parent.createfilehandler(
            self._proc.stdout, tk.READABLE, self._HandleChildOutput
        )

    def _HandleChildOutput(self, fd, mask):
        """Handle output from child process."""
        self._ShowOutput()
        status = self._proc.poll()
        if status != -1:
            # child process completed
            time.sleep(0.5)  # wait for last output
            self._ShowOutput()
            if os.WEXITSTATUS(status) != 0:
                self.ShowResult(
                    "Subprocess exited with code %d\n" % os.WEXITSTATUS(status), "error"
                )
            else:
                self.ShowResult("Subprocess exited normally\n", "commentary")
            self.parent.deletefilehandler(self._proc.stdout)
            os.remove(self._tempName)
            # Execute post-callback
            if callable(self.postFitCallback):
                self.postFitCallback()

    def _ShowOutput(self):
        """Read child process output and pass to ShowResult()."""
        try:
            result = self._proc.stdout.read()
        except IOError as e:
            (errNo, errMsg) = e.args
            self.ShowResult("I/O Error %d: %s\n" % (errNo, errMsg), "commentary")
        else:
            self.ShowResult(result, "result")


def _main(altExe=None):
    """Main routine."""
    if len(sys.argv) == 1:
        # No command-line arguments - run graphical user interface
        root = tk.Tk()
        root.title("fitgui %s" % _revision)
        main = GUI(root)
        main.fileName.set("test.oifits")
        main.ReadModel("test.model")
        if altExe is not None:
            main.exe = altExe
        root.mainloop()
    else:
        # Too many arguments
        print("%s %s\n" % (sys.argv[0], _revision))
        print(__doc__)
        sys.exit(2)


if __name__ == "__main__":
    _main()
