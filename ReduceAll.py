__author__ = 'kgulliks'

import Reduce_Components
from collections import defaultdict
import warnings
import os
import sys


class Files:
    def __init__(self, logfilename):
        """
        This class reads in the logfile, and stores each file
        in an appropriate dictionary with the key as
        the type of file (off or on). If the file is a
        standard or object file, the dictionary will contain
        a further set of keys for the name of the object


        :param logfilename:
        :defines: self.Flats
                  self.Arcs
                  self.Standards
                  self.Targets
        :return: None
        """

        self.Standards = defaultdict(lambda: defaultdict(list))
        self.Targets = defaultdict(lambda: defaultdict(list))
        self.Flats = defaultdict(list)
        self.Arcs = defaultdict(list)
        self.WorkingDirectory = os.path.dirname(logfilename)

        infile = open(logfilename)
        lines = infile.readlines()[2:]
        infile.close()
        for line in lines:
            segments = line.split(",")
            filename = segments[0]
            objname = segments[4]
            objtype = segments[5]
            frametype = segments[6]

            # Assign this file to the appropriate dictionary
            if objtype.lower() == "std":
                self.Standards[objname][frametype].append("%s/%s" %(self.WorkingDirectory, filename))
            elif objtype.lower() == "tar":
                self.Targets[objname][frametype].append("%s/%s" %(self.WorkingDirectory, filename))
            elif objtype.lower() == "flat":
                self.Flats[frametype].append("%s/%s" %(self.WorkingDirectory, filename))
            elif objtype.lower() == "arc":
                self.Arcs[frametype].append("%s/%s" %(self.WorkingDirectory, filename))
            else:
                warnings.warn("Unknown object type found: %s" % objtype)

        # Check to make sure there are the right number of on/off frames
        # The arcs often don't have off frames and we just use off flats
        self.CheckOnOff(self.Flats)
        if len(self.Arcs["OFF"]) == 0:
            self.Arcs["OFF"] == []
        self.CheckOnOff(self.Arcs, insert_flats=True)

        # Do the same for A/B frames for the objects
        if len(self.Targets) > 0:
            self.CheckAB(self.Targets)
        if len(self.Standards) > 0:
            self.CheckAB(self.Standards)


    def CheckOnOff(self, dictionary, insert_flats=False):
        """
        This method checks that the given dictionary has the right
        number of on/off frames. If insert_flats=True, it will
        insert off flats if there are not enough of them.

        If there is a mismatch, it will warn the user and return False

        :param dictionary to check:
        :param insert_flats (optional, default False):
        :return: bool
        """
        if len(dictionary["OFF"]) == len(dictionary["ON"]):
            return True
        elif len(dictionary["ON"]) > len(dictionary["OFF"]) and insert_flats:
            i = 0
            while len(dictionary["ON"]) > len(dictionary["OFF"]):
                dictionary["OFF"].append(self.Flats["OFF"][i])
                i += 1
            return True
        else:
            warnings.warn("Different number of off and on frames!")
            return False

    def CheckAB(self, dictionary):
        """
        This method checks that the given dictionary has the same
        number of A and B frames.

        If there is a mismatch, it will warn the user and return False

        :param dictionary:
        :return: bool
        """
        starnames = dictionary.keys()
        retval = True
        for starname in starnames:
            if len(dictionary[starname]["A"]) == len(dictionary[starname]["B"]):
                return True
            else:
                warnings.warn("Star %s has an incorrect number of on/off frames!" % starname)



if __name__ == "__main__":
    # Read in the logfile (you must give that)
    files = Files(sys.argv[1])

    #Figure out what band we are working in
    band = files.WorkingDirectory[-1]

    # Reduce the flats
    print "Reducing Flats..."
    Reduce_Components.FlatReduce(files, band)


