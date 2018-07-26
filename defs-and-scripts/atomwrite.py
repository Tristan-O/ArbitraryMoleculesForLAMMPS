#!/usr/bin/env python

#===== atomwrite.py =====
#Routines for writing configuration coordinates
#to popular molecular file formats.

import os
import numpy as np
import gzip


def MinImage(Pos, L):
    """Returns a new Pos array with minimum-imaged positions."""
    return Pos - L * np.round_(Pos / L)


def PdbStr(Pos, L = None, AtomNames = ["C"], ModelNum = 1):
    """Gets a Pdb string.
Input:
    Pos: (N,3) numpy array of atomic positions
    L: scalar or vector of box lengths (or None to skip min-imaging)
    AtomNames: list of atom names; duplicated as necessary to length N
"""
    N = len(Pos)
    #make a new AtomNames the same length as the position array
    #by repeating it over and over
    AtomNames = AtomNames * (N / len(AtomNames) + 1)
    AtomNames = AtomNames[:N]        
    #minimum image the positions
    if not L is None:
        Pos = MinImage(Pos, L)
    #make the pdb header
    s = "MODEL     %-4i\n" % ModelNum
    #add the coordinates to the pdb string
    for (i, (x,y,z)) in enumerate(Pos):
        an = AtomNames[i].strip()
        s += "HETATM%5i %4s %3s  %4i    %8.3f%8.3f%8.3f                     %2s  \n" % (
             i+1, an.ljust(4)[:4], "SYS", i+1, x, y, z, an.rjust(2)[:2])
    #make the footer
    s += "TER\nENDMDL\n"
    return s


def WritePdb(Filename, Pos, L = None, AtomNames = ["C"], First = False):
    """Writes a position array to a pdb file.
Input:
    Filename: string filename of file to write
    Pos: (N,3) numpy array of atomic positions
    L: scalar or vector of box lengths (or None to skip min-imaging)
    AtomNames: list of atom names; duplicated as necessary to length N
    First: True will overwrite or start a new Pdb file; False will append
"""
    #check to see if the file exists;
    #if so, and we are appending, get the model number
    ModelNum = 1
    if not First and os.path.isfile(Filename):
        for line in file(Filename, "r"):
            if line.startswith("MODEL"):
                ModelNum = int(line[10:].strip())
        ModelNum += 1
    #get the data
    s = PdbStr(Pos, L, AtomNames, ModelNum)
    #write to the file
    if First:
        file(Filename, "w").write(s)
    else:
        file(Filename, "a").write(s)


class pdbfile:
    def __init__(self, FileName, L = None, AtomNames = ["C"],
                 Compressed = True):
        """Creates a new class for writing to a pdb file.
        Input:
        Filename: string filename of file to write
        L: scalar or vector of box lengths (or None to skip min-imaging)
        AtomNames: list of atom names; duplicated as necessary to length N
        Compressed: True will write to a gzipped file (default True)
        """
        self.L = L
        self.AtomNames = AtomNames
        self.ModelNum = 1
        if Compressed:
            self.FileName = FileName + ".gz"
            self.fobj = gzip.GzipFile(self.FileName, "w")
        else:
            self.FileName = FileName
            self.fobj = file(self.FileName, "w")

    def write(self, Pos):
        """Writes positions to a Pdb file.
        Input:
        Pos: (N,3) numpy array of atomic positions
        """
        s = PdbStr(Pos, self.L, self.AtomNames, self.ModelNum)
        self.fobj.write(s)
        self.ModelNum += 1

    def close(self):
        """Closes Pdb file object."""
        self.fobj.close()

        
            

    
        