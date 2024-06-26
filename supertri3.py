#!/usr/bin/env python


import copy
import getopt
import os
import re
import sys
import warnings
from functools import reduce
from collections import OrderedDict

__version__ = "$Revision: 59 $"
__licence__ = """
supertri3.py

Copyright (C) 2007 Anne Ropiquet, Blaise Li and Alexandre Hassanin

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
"""


MAX_TAX_NAME = 0

def oui(question):
    return (input("\n" + question + " (y/n) ") in
            ["O", "o", "Oui", "oui", "OUI", "Y", "y", "Yes", "yes", "YES", "J", "j", "Ja", "ja", "JA", "D", "d", "Da",
             "da", "DA"])


def attend():
    pause = input("\npress enter to continue...")
    return 0


def askfile(question, list_of_files=None):
    while True:
        file = input(question)
        if list_of_files is None:
            return file
        elif file not in list_of_files:
            print(f"You must chose among: {' '.join(list_of_files)}.")
        else:
            return file


def mac2x(string):
    # This fct is not used anymore i think
    """mac2x(string) returns a string corresponding to <string> but with the mac-style ends of line translated to
    unix-style ends of line. """
    new_string = str.replace("\r", "\n", string)
    return new_string


###########
# Objects #
###########

class NodeIndex(float):
    """
    Represents a node support index

    Attributes
    ----------
    initval: float
        Initial value of the node support
    """


    format = "%.2f"

    def __init__(self, initval):
        float.__init__(initval)

    def __str__(self):
        return self.__class__.format % float(self)

    def reformat(self, newformat):
        """This function changes globally the format for all indexes of the same class as self."""
        assert isinstance(newformat, str), "newformat must be a string"
        self.__class__.format = newformat

    def __add__(self, other):
        return self.__class__(float(self) + float(other))

    def __sub__(self, other):
        return self.__class__(float(self) - float(other))

    def __mul__(self, other):
        return self.__class__(float(self) * float(other))

    def __div__(self, other):
        return self.__class__(float(self) / float(other))


class Bootstrap(NodeIndex):
    """
    This object represents a node bootstrap percentage.

    Attributes
    ----------
    initval: float
        Initial Boostrap value

    """
    format = "%.0f"

    def __init__(self, initval):
        NodeIndex.__init__(self, initval)

    def __str__(self):
        # perc = 100*self
        return f"{(int(self * 100)):}"  # Bootstrap values are displayed as percentages.


class PosteriorP(NodeIndex):
    """
    This object represents a node posterior probability.
    """
    format = "%.2f"

    def __init__(self, initval):
        NodeIndex.__init__(self, initval)


class Bipartition(object):
    """
    This object represents a bipartition.

    """

    def __init__(self, inclade, outclade, brlen=None):
        self.datasets = set([])  # initiate empty set out of a list (an iterable))
        self.weights = {}
        self.weight = NodeIndex(0)
        self.support = NodeIndex(0)
        self.meansup = NodeIndex(0)
        self.occurences = 0
        self.repro = 0
        assert isinstance(inclade, int) and isinstance(outclade, int)
        assert inclade & outclade == 0, "the two parts of a bipartition cannot have shared taxa"
        self.id = frozenset([inclade, outclade])
        self.inclade = inclade
        self.outclade = outclade
        self.brlen = brlen  # optional, indicates the length of the branch separating the two parts of the bipart
        self.domain = inclade + outclade
        self.quartets = None

    def __str__(self):
        return f"Bipart ({self.inclade},{self.outclade})"

    def compatible(self, other):
        """self.compatible(other) returns True if self and other can be found in the same tree."""
        assert isinstance(other, (Bipartition, Clade, list)), \
            "Must be compared with a clade or a bipartition or a list of such objects."
        if isinstance(other, list):
            return self.allcompat(other)
        assert self.domain == other.domain, "Cannot state about compatibility if the validity domains are not the same."
        if self.inclade == 0 or other.inclade == 0:
            return True
        a = self.inclade & other.inclade
        b = self.outclade & other.inclade
        c = self.inclade & other.outclade
        d = self.outclade & other.outclade
        if isinstance(self, Clade) and isinstance(other, Clade):
            return a == 0 or b == 0 or c == 0
        else:
            return a == 0 or b == 0 or c == 0 or d == 0

    def allcompat(self, list_of_others):
        """self.allcompat(list_of_others) returns True if self is compatible with all the elements of the list
        <listofothers>. """
        assert isinstance(list_of_others, list), "allcompat only works with lists."
        result = True
        i = 0
        while result:  # No need to continue once one incompatibility has been found.
            if i == len(list_of_others):  # All members of the list have been tested.
                return result
            result = result and self.compatible(list_of_others[i])
            i += 1
        return result

    def update(self, w, newdata=None):
        """self.update(w, newdata) increases the weight of the bipartition and records the new data that produced the
        bipatition. """
        assert isinstance(w, (int, float)) and w >= 0, "The weight is supposed to be a positive number."
        # A NodeIndex object is also a float.
        assert isinstance(newdata,
                          list) or newdata is None, "If there is new data implied, it should be given as a list of " \
                                                    "datasets names. "
        if newdata is not None:
            for m in newdata:
                self.weights[m] = w + self.weights.get(m, 0)
        self.weight = w + self.weight
        self.support = w + self.support
        if newdata is not None:
            for m in newdata:
                if m in self.datasets:
                    raise Exception(f"dataset {m} already taken into account")
                else:
                    self.datasets.add(m)

    def clade(self, root):
        """self.clade(root) returns a clade corresponding to the rooting of self by the taxon <root>, given as a
        power of 2. """
        if root & self.outclade == root:
            clad = Clade(self.inclade, self.outclade, root)
        elif root & self.inclade == root:
            clad = Clade(self.outclade, self.inclade, root)
        else:
            raise Exception("Cannot root with a taxon not placed inside neiter outside the bipartition.")
        for attr in self.__dict__:
            if attr not in ['id', 'inclade', 'outclade', 'domain']:
                clad.__dict__[attr] = copy.deepcopy(self.__dict__[attr])
        return clad


class Clade(Bipartition):
    """
       "This object represents a clade, that is an oriented bipartition.
    It also has a dictionnary of subclades, to optionally specify its internal topology.

       Attributes
       ----------
       inclade :
           smth
       outclade: smth
           smth
       brlen: float
           Represents branch length of the bipartition
       root: float
           Root of the tree

       """

    def __init__(self, inclade, outclade, root=None, brlen=None):
        Bipartition.__init__(self, inclade, outclade, brlen)
        self.id = (self.inclade, self.outclade)
        assert (root is None) or (
                root & outclade == root), "If explicit, the root should be one of the taxa outside the clade."
        self.root = root  # optional attribute that indicates with respect to which taxon the clade is rooted
        self.subclades = {}
        self.tgf = None  # description of the clade in tgf format
        self.triplets = None  # We don't necessarily need to know the triplets implied by the clade. Use
        # self.findtriplets for that.

    def bipart(self):
        """self.bipart() returns a bipartition corresponding to the unrooting of self."""
        bip = Bipartition(self.inclade, self.outclade, self.brlen)
        for attr in self.__dict__:
            if attr not in ['id', 'inclade', 'outclade', 'root', 'domain', 'subclades', 'tgf', 'triplets']:
                bip.__dict__[attr] = copy.deepcopy(self.__dict__[attr])
        return bip

    # Beware of the implication of this method on the use of the "in" operator.
    def __contains__(self, other):
        """self.__contains__(other) returns True if Clade <other> can be a subclade of self, False otherwise.
                Empty clades are considered included in any clade that shares the same validity domain. Inclusion implies compatibility."""

        assert isinstance(other, Clade), "To state about inclusion, the object must be a Clade type object."
        if other.domain == self.domain:
            inter = other.inclade & self.inclade
            # Verification
            if other.inclade == inter:
                assert self.compatible(other), "other is in self, so it should be compatible."
            return other.inclade == inter
        else:
            return False

    # All these functions are used instead of the "cmp" function in the python 2 version
    # Python 3 does not support __cmp__ anymore, every comparaison operation has to be hard coded
    def __lt__(self, other):
        if self.id == other.id:
            return 0
        if self.domain == other.domain:
            return self.inclade < other.inclade

    def __le__(self, other):
        if self.id == other.id:
            return 0
        if self.domain == other.domain:
            return self.inclade <= other.inclade

    def __eq__(self, other):
        if self.id == other.id:
            return True
        if self.domain == other.domain:
            return False
        else:
            raise ValueError("""To be compared, two clades must have the same validity domain,
                                         that is, define a partition of the same set of taxa.""")

    def __ne__(self, other):
        if self.id == other.id:
            return 0
        if self.domain == other.domain:
            return self.inclade != other.inclade

    def __gt__(self, other):
        if self.id == other.id:
            return 0
        if self.domain == other.domain:
            return self.inclade > other.inclade

    def __ge__(self, other):
        if self.id == other.id:
            return 0
        if self.domain == other.domain:
            return self.inclade >= other.inclade

    def __str__(self):
        return f"Clade {self.inclade},{self.outclade})"

    def addsubclades(self, subclades):
        """self.addsubclades(subclades) recursively includes the clades in the list <subclades> into self.
                This method should not be called with a partial list of subclades, since all taxa of the clade not
                present in the subclades are considered being placed in an unresolved basal position (a "rake")."""

        assert isinstance(subclades, list), "This method only accepts a list of subclades as argument."
        rejected = []
        to_include = []
        to_keep = []
        for subclade in subclades:
            assert isinstance(subclade, Clade), "The list should only contain Clade type objects."
            assert subclade.domain == self.domain, "It doesn't make sense to combine clades with different validity domains."
            if subclade != self:
                if subclade in self:
                    to_include.append(subclade)
                elif self.compatible(subclade):
                    to_keep.append(subclade)
                else:
                    rejected.append(subclade)
        to_include.sort()

        while to_include:
            candidate = to_include.pop()
            possible = True
            for subclade in self.subclades:
                if not self.subclades[subclade].compatible(candidate):
                    rejected.append(candidate)
                    possible = False
                    break
            if possible:
                self.subclades[candidate.id] = candidate
                to_include = self.subclades[candidate.id].addsubclades(to_include)
            else:
                print(f"A clade that was not compatible with {self}'s subclades was found.")

        invalue = 0
        for subclade in self.subclades:
            invalue += self.subclades[subclade].inclade
        raketaxa = self.inclade - invalue
        for taxa in decomposition(raketaxa):  # decomposition is defined later in main program
            assert taxa != self.root, "The root ca,,ot be included as a subclade"
            self.subclades[(taxa, self.domain - taxa)] = Clade(taxa, self.domain - taxa, self.root)

        # to_keep.sort()
        return sorted(to_keep)

    def set_support_type(self, type):
        """self.set_support_type(type) recursively sets the support type to the class <self>."""
        self.support.__class__ = type
        for subclade in self.subclades.values():
            subclade.set_support_type(type)

    def find_clade(self, c):
        """self.find_clade(c) tries to find if clade c is a subclade of self or of one of its subclades, and recursively.
         It returns a copy of the clade found, if there is one, None otherwise."""
        # The cause of the error is here, but why ?
        assert isinstance(c, Clade), "Can only tell inclusion for Clade objects."
        try:
            assert c.domain == self.domain, "Can only tell inclusion for clades with the same validity domain."
        except AssertionError:
            warnings.warn("Default for clades with different validity domains is to consider that one is not included "
                          "in the other. Don't panic, this is just a little warning; the execution is probably going "
                          "on...")
            return None
        if c == self:
            return copy.deepcopy(self)

        detected = None
        for s in self.subclades.values():
            if detected:
                return detected  # no need to search further
            detected = s.find_clade(c)

        return detected

    def is_supported_by(self, tree):
        """self.is_supported_by(tree) returns the support value of the clade of <tree> "equivalent" to the clade self
        if there is one, 0.0 if there is not, and None if this is not relevant. If there are missing taxa for <tree>,
        the support of a clade is actually the support of its restriction to the validity domain of <tree>. It is not
        relevant to assess the support when the minimal set of taxa necessary to "recognize" the clade self is
        missing from <tree>. """

        for subclade in self.subclades.values():
            if subclade.inclade & tree.domain == 0:
                return None  # The tree is not relevant because it doesn't have the essential taxa;
                # at least one of its subclades has no possible equivalent in the restricted domain.
        sup_inclade = self.inclade & tree.domain
        sup_outclade = self.outclade & tree.domain
        supporting = tree.find_clade(Clade(sup_inclade, sup_outclade))
        if supporting:
            return supporting.support
        return tree.support.__class__(0)

    def compute_indices(self, source_trees, bipartitions):
        """self.compute_indices(sourcetrees) computes the mean support self.meansup and the reproducibility self.repro
        using the trees in the list <sourcetrees> and the bipartitions in <bipartitions>. It proceeds recursively for
        all subclades of self. """

        self.meansup = Bootstrap(self.meansup)
        assert not (self.support or self.occurences), "index computation must start from zero."
        # all bipartitions considered successively
        for bipart in bipartitions.values():
            relevant = True
            for subclade in self.subclades.values():
                if subclade.inclade & bipart.domain == 0:
                    relevant = False
                    break
            if relevant:
                if bipart.inclade == self.inclade & bipart.domain:
                    self.support = (self.support + bipart.support)
                elif bipart.outclade == self.inclade & bipart.domain:
                    self.support += (self.support + bipart.support)
        relevant = 0

        for tree in source_trees:
            supp = self.is_supported_by(tree)  # this happens bc we never have supp
            if supp is not None:
                relevant += 1
                if supp:
                    self.occurences += 1
        self.repro = self.occurences / float(relevant)
        self.meansup = Bootstrap(self.support / relevant)

        for subclade in self.subclades:
            self.subclades[subclade].compute_indices(source_trees, bipartitions)

    def parenthesize(self, converter, index=""):
        """self.parenthesize(converter, index) returns the parenthesized description of self.
                <converter> is a dictionnary that translates from powers of 2 to taxon names.
                If <index> is specified, the attribute whose name is <index> will be written as node labels."""
        assert isinstance(index, str), "If defined, index should be provided as a string."
        if self.inclade in converter:
            if self.brlen:
                return_str = f"{converter[self.inclade]}:{self.brlen}"
                return return_str
                # return converter[self.inclade] + ":%s" % self.brlen
            return converter[self.inclade]

        else:
            arguments = []
            for subclade in self.subclades.values():
                arguments.append((subclade, converter, index))
            if index:
                if self.brlen:
                    fstr = f"({''.join([Clade.parenthesize(a[0], a[1], a[2]) for a in arguments])})"
                    return "(" + ",".join([Clade.parenthesize(a[0], a[1], a[2]) for a in arguments]) + \
                           f'){self.__dict__.get(index, "")}:{self.brlen}'
                return "(" + ",".join([Clade.parenthesize(a[0], a[1], a[2]) for a in arguments]) + \
                       f'){(self.__dict__.get(index, ""))}'
                # ")%s" % (self.__dict__.get(index, ""))
            if self.brlen:
                return "(" + ",".join([Clade.parenthesize(a[0], a[1], a[2]) for a in arguments]) + f"):{self.brlen}"
            return "(" + ",".join([Clade.parenthesize(a[0], a[1], a[2]) for a in arguments]) + ")"

    def gettgf(self, converter, indent):
        """self.gettgf(converter, indent) sets self.tgf and returns it.
                <converter> is the dictionnary that translates the powers of 2 into taxa names.
                <indent> is just to improve readability of the tgf format."""
        if self.tgf is not None:
            return self.tgf
        else:
            if self.inclade in converter:  # It's a terminal taxon.
                if self.brlen is not None:
                    # TODO include { to fstr (with a backslash probably
                    # self.tgf = f"{indent* ' '}\\len{self.brlen} \\r{converter[self.brlen]},\n"
                    self.tgf = indent * " " + "\\len{%f} \\r{%s},\n" % (self.brlen, converter[self.inclade])
                else:
                    # self.tgf = f"{indent* ' '}\\r{converter[self.inclade]},\n"
                    self.tgf = indent * " " + "\\r{%s},\n" % converter[self.inclade]
            else:  # The clade has subclades (at least its terminal taxa).
                if self.brlen is not None:
                    # self.tgf = f"{indent* ' '}\\len{self.brlen} \\u1{self.meansup} \\u2{self.repro}"
                    self.tgf = indent * " " + "\\len{%f} \\u1{%s} \\u2{%s}(\n" % (self.brlen, self.meansup, self.repro)
                else:
                    # self.tgf = f"{indent* ' '}\\u1{self.meansup} \\"
                    self.tgf = indent * " " + "\\u1{%s} \\u2{%s}(\n" % (self.meansup, self.repro)
                indent = indent + 1  # The indentation is increased before writing the subclades.
                for s in self.subclades.values():
                    self.tgf = self.tgf + s.gettgf(converter, indent)
                self.tgf = self.tgf[:-2] + "\n"  # The comma following the last subclade is removed.
                indent = indent - 1  # The indentation is decreased after writing the subclades.
                self.tgf = self.tgf + indent * " " + "),\n"  # The clade is closed.
            return self.tgf

    def next_tree(self, converter, index="", to_root=True):
        """self.next_tree(converter, index) returns the nexus representation of self.
                <converter> is a dictionnary that translates from powers of 2 to taxon names.
                If <index> is specified, the attribute whose name is <index> will be written as node labels.
                If <toroot> is True, the root of the tree must be explicitly added to the nexus parenthesised form."""

        if index == "meansup":
            self.meansup = Bootstrap(self.meansup)
            tmpfmt = self.meansup.format
            self.__dict__[index].reformat(tmpfmt)  # somehow this returns a float and we should have a str?
        else:
            tmpfmt = "%2f"

        if not to_root:
            newick = f"tree {index} = ({self.parenthesize(converter, index)});"
            if index == "meansup":
                self.__dict__[index].reformat(tmpfmt)
            return newick
        newick = "tree %s = (%s," % (index, converter[self.root]) + self.parenthesize(converter, index) + ");"
        # TODO fstr
        if index == "meansup":
            self.__dict__[index].reformat(tmpfmt)
        return newick

    def tgftree(self, converter):

        # TODO : do smth abt that massive string :)
        tgf = "\\begindef\n\\paper{a0}\n\\style{r}{italic}{12}\n\\style{u1}{plain}{11}\n\\style{u2}{plain}{11}\n\\separator{ / }\n\\enddef\n\\label{root}(\n"
        indent = 1
        for subclade in self.subclades.values():
            tgf = tgf + subclade.gettgf(converter, indent)
        indent -= 1
        tgf = tgf[:-2] + "\n)\n"
        return tgf


class BipartSet(dict):
    """
      This object represents a set of bipartitions.
    It is a dictionnary of bipartitions of the type Bipart.
    bipartitions should be added using the method add_bipart()

       """

    def __init__(self):
        dict.__init__(self)
        self.tax2val = {}

    def __str__(self):
        name = "{"
        for b in self:
            name = name + self[b].__str__() + ","  # TODO format this string better
        return name[:-1] + "}"

    def set_converter(self, dic):
        assert isinstance(dic,
                          dict), "The converter is supposed to be a dictionnary converting taxon names into powers of 2."
        self.tax2val = dic

    def add_bipart(self, bipart):
        """self.add_bipart(bip) adds bipartition <bip> to self if it is new or updates self according to the new
        informations carried by <bip>. """
        assert isinstance(bipart, Bipartition), "Only Bipartition objects can be added to a BipartSet object."

        if bipart.id in self:
            self[bipart.id].update(bipart.weight, list(bipart.datasets))
        else:
            self[bipart.id] = bipart

    def matrix_line(self, taxon, converter=None):
        """self.matrix_line(taxon,converter) returns a nexus matrix representation line for taxon <taxon> with "1",
        "0", or "?" for the columns representing the bipartitions of self. "1" is when the taxon is present in the
        bipartition, 0 when it's outside and "?" when the information is missing. converter is the dictionnary used
        to convert taxa to values as powers of 2. """
        assert isinstance(taxon, str), "Taxon names should be strings."

        def valid_name(tax):
            """valid_name(taxon) returns the conversion of the name <taxon> to a string valid in PAUP. Spaces are
            replaced with underscores ("_") and if other non-alphanumeric characters are met, the string is placed
            between single quotes. """
            special = False # True once special characters have been met
            valid = ""
            for c in tax:
                if c.isalnum():
                    valid += c
                elif c.isspace() or c == "_":
                    valid += "_"
                else:
                    valid += c
                    special = True
            if special:
                valid = "'" + valid + "'"  # TODO : fstr
            return valid

        line = valid_name(taxon) + (3 + maxtaxname - len(valid_name(taxon))) * " "  # maxtaxname = gobal var in main fct
        if converter is None:
            converter = self.tax2val

        assert isinstance(converter, dict), "The converter is supposed to be a dictionnary converting taxon names " \
                                            "into powers of 2. "
        val = converter[taxon]
        cols = list(self.keys())  # turned into list because there is no .sort for dict keys anymore
        #sorted(cols)
        for col in cols:

            # This should be correct, yet somehow we have different values apparently ?
            for marker in self[col].weights:
                if self[col].inclade & val:
                    line += "1"
                elif self[col].outclade & val:
                    line += "0"
                else:
                    line += "?"
        return line

    def weight_set(self, set_name):
        """self.weight_set() returns a nexus weight set corresponding to the weights of the bipartitions of self.
                <setname is the name to give to the set>"""
        assert isinstance(set_name, str), "The name of a weight set should be a string."
        cols = list(self.keys())
        cols.sort()

        def weights(column):
            """weights(column) returns the list of the pairs (weight, marker) corresponding the bipartition whose
            column is <column>. """
            weights = []
            for marker in self[column].weights:
                weights.append((self[column].weights[marker], marker))
            return weights

        wtset = f"wtset {set_name} = "
        indent = " " * len(wtset)
        wtset += "\n"
        i = 1

        for col in cols:
            for w in weights(col):
                wtset = wtset + indent + "%s: %s [%s],\n" % (w[0], i, w[1])
                i = i + 1
        wtset = wtset[:len(indent)] + wtset[2 * len(indent) + 1:-2] + ";\n"
        return wtset

    def nex_matrix_rep(self, taxons, converter=None, weight=True, set_name=None):
        """self.nexmatrixrep(weight, setname) returns the matrix representation of the bipartitions stored in self.
                <taxons> is the list of the taxon names that are supposed to be in the bipartitions.
                <converter> is the dictionnary converting from taxon names to powers of 2.
                If <weight> is true, an assumptions bolck is added, defining the weights of the columns of the matrix.
                <setname> is the name to give to the weights set."""

        assert isinstance(taxons, list), "Taxons should be a list of taxon names."
        if converter is None:
            converter = self.tax2val
        assert isinstance(converter, dict), \
            "The converter is supposed to be a dictionary converting taxon names into powers of 2."
        if set_name is None:
            set_name = "weights"
        assert isinstance(set_name, str), "The name of a weight set should be a string"
        nchar = 0
        for bipart in self:
            for marker in self[bipart].weights:
                nchar += 1
        matrix = """#Nexus
[produced by supertri.py %s]
begin data;
    dimensions ntax=%s nchar=%s;
    format datatype=standard missing=?;
    matrix
""" % (__version__, len(taxons), nchar)  # TODO: make an fstring out of this ?
        for taxon in taxons:
            matrix += (self.matrix_line(taxon, converter) + "\n")
        matrix += ";\nend;\n"
        if weight:
            matrix = matrix + "\nbegin assumptions;\n" + self.weight_set(set_name) + "end;"  # TODO fstring
        return matrix


# FUNCTION DEFINITION STARTS HERE ------------------------------------------------------------------------------------


def decomposition(integer):
    """decomposition(integer) returns the list of the powers of 2 composing <integer>."""
    powers = []
    sum = 0
    puiss = 1
    while sum != integer:
        if puiss & integer:
            powers.append(puiss & integer)
            sum = sum + powers[-1]
        puiss *= 2
    return powers


def read_bootlog(bootlog, ntaxs):
    """read_bootlog(ntaxs) reads the PAUP bootstrap log file <bootlog> that contains the description of bipartitions
    in the form of "." and "*", with values associated to these bipartitions (for example, bootstrap supports). This
    function returns a list of pairs representing the bipartitions read in <bootlog>. and their associated values. """

    print(f"reading {bootlog} ...")
    file = open(bootlog, 'r')  # TODO with open ?
    lines = file.read().split("\n")
    file.close()
    selection = []
    for line in lines:
        goodline = re.match("^[\.\*]+.*$", re.sub("\\[.*\\]", "", line))
        if goodline is not None:
            selection.append(goodline.group())
    bipparts = []
    part_index = 0
    i = len(selection) - 1
    ref_len = len(selection[i].split()[0])
    tot_len = ref_len
    print(tot_len)
    while len(selection[i].split()) >= 2:
        line_len = len(selection[i].split()[0])
        assert line_len == ref_len, "The lines describing bipartitions are supposed to be of equal lengths if they " \
                                    "are in the same block. "
        bipparts.append([selection[i].split()[0], selection[i].split()[-1]])
        i -= 1
        if i < 0:
            break  # avoids Index Error
    n_biparts = len(bipparts)
    if i >= 0:
        bip_index = 0
        assert len(selection[
                       i].split()) == 1, "Only the last block of the bipartitions description may have a support " \
                                         "value following stars and points after a space. "
        ref_len = len(selection[i].split()[0])
        tot_len = tot_len + ref_len
    while i >= 0:  # There are parts of line to add.
        line_len = len(selection[i].split()[0])
        assert len(selection[
                       i].split()) == 1, "Only the last block of the bipartitions description may have a support " \
                                         "value following stars and points after a space. "
        assert line_len == ref_len, "The lines describing bipartitions are supposed to be of equal lengths if they are in the same block."
        bipparts[bip_index][0:0] = [
            selection[i].split()[0]]  # The new part of line is inserted as first item of the list of parts of the line.
        i = i - 1  # previous line in selection
        if i < 0:
            break  # to avoid an IndexError
        if bip_index == n_biparts - 1:
            if i >= 0:  # There will be other parts of line to add.
                bip_index = 0
                assert len(selection[i].split()) == 1, "Only the last block of the bipartitions description may have " \
                                                       "a support value following stars and points after a space. "
                ref_len = len(selection[i].split()[0])
                tot_len += ref_len
        else:
            bip_index += 1

    bipartitions = []
    for bip_line in bipparts:
        bipartitions.append(("".join(bip_line[:-1]), float(re.sub("\%", "", bip_line[-1])) / 100))
    bipartitions.reverse()
    return bipartitions


def read_biparts(bipart_file):
    """read_biparts(partsfile) reads the file <partsfile> that contains the description of bipartitions in the form of
    "." and "*", with values associated to these bipartitions (for example, posterior probabilities). This function
    returns a list of pairs representing the bipartitions read in <partsfile>. and their associated values. """

    bipartitions = []
    fich = open(bipart_file, 'r')
    lines = fich.read().split("\n")
    fich.close()
    for line in lines:
        bipline = re.match("^.*[\.\*]+.*$", re.sub("\\[.*\\]", "", line))
        if bipline is not None:
            fields = bipline.group().split()
            bipartitions.append((re.findall("[\.\*]+", bipline.group())[0], float(fields[2])))
    return bipartitions


def val_stars(l, missing_vals=[]):
    """val_stars(line,missing) reads a line <l> composed of "." and "*" representing a bipartition. This function
    returns the value of the set of taxa represented by a "*" and the value of the set of taxa represented by a ".",
    assuming that taxa whose values are in the list <missing_vals> are not represented in <l>. """
    puiss = 1
    val = 0
    compl = 0
    for c in l:
        if c == "?":
            raise Exception(
                "Question marks are not authorized yet; missing taxa must be missing for all bipartitions produced by "
                "the dataset and be declared in a .abs file.")
        assert c in ["*", "."], "A taxon must be in or out of the bipartition, or be declared as missing from the " \
                                "dataset in a .abs file. "
        while puiss in missing_vals:
            puiss = 2 * puiss  # skipping missing taxa
        if c == "*":
            val = val + puiss
        else:
            compl = compl + puiss
        puiss = puiss * 2
    return val, compl


def readnextrees(treefile, converter):
    """readnextrees(file, taxons) reads the trees in nexus format that are written in the file <treefile>. These trees must use only taxa that are keys of the dictionnary <converter>.
        <converter> is a dictionnary that translates from taxon names to powers of 2.
        This function returns a list of Clade objects that correspond to the trees, without a specified root."""
    with open(treefile, 'r', encoding='latin1') as fich:
        lines = fich.read().split("\n")
    trees = []
    for line in lines:
        treeline = re.match(".*tree.*\\(.*\\);", re.sub("\\[.*\\]", "", line))  # Comments are removed from the line.
        if treeline is not None:
            parenth = re.sub("\\);[^()]*$", ")", re.sub("^[^()]*\\(", "(", treeline.group()))
            domain = 0
            for t in converter:
                if t in parenth:
                    domain = domain + converter[t]
            trees.append(read_parenth(parenth, converter, domain))
    return trees


def regtaxname(taxon):
    """regtaxname(taxon) returns a compiled regular expression that should match <taxon> as a taxon name in a newick tree,
    that is <taxon> between separators like ":", "(", "," or ")". """
    reg = re.compile("\\W%s\\W" % taxon)  # The taxon name should be between 2 non-word characters.
    return reg


def read_brlen(parenth, i):
    """read_brlen(parenth, i) returns the value of the branch length written in the parenthesized clade description
    <parenth> just after position <i> in this string. If there is no branch length, it return None. This function
    also retruns the position in <parenth> where the branch length ends. """
    if i + 1 < len(parenth):
        if parenth[i + 1] in [",", ")"]:
            return None, i
        if parenth[i + 1] == ":":
            m = re.search(":\d*\.*\d*", parenth[i + 1:])
            assert m is not None, "In a newick tree, a colon should be followed by a branch length value."
            return float(m.group()[1:]), i + m.end()
        raise NotImplementedError("This program can not read trees with node labels or similar ornementation.")
    return None, i


def read_terminal(parenth, i, converter):
    """read_terminal(parenth, i, converter) reads a terminal taxon and its possible branch length in the parenthesized
    clade <parenth> starting at positin <i> in this string. This function returns the value of the terminal taxon as
    obtained by the use of the dictionnary <converter>, the value of the branch length if there is one, or None,
    and the position in <parenth> corresponding to the last character of the taxon or of its branch length. """

    r = re.compile("\w+:*\d*\.*\d*")  # pattern matching a taxon name possibly followed by a branch length
    # r = re.compile("\w+\(*\d*\)*:*\d*\.*\d*") # pattern matching a taxon name followed by a branch length
    m = r.match(parenth, i)  # use match rather than search: search tilts too soon
    if m is not None:
        pair = m.group().split(":")
        if len(pair) == 2:
            return converter[pair[0]], float(pair[1]), m.end()
        else:
            return converter[pair[0]], None, m.end()
    else:
        return None, None, i


def read_parenth(parenth, converter, domain, root=None, brlen=None):
    """read_parenth(parenth, converter, domain, root, brlen) builds the Clade object that corresponds to the newick format tree <parenth>. This tree must only use taxa that are keys in the dictionnary <converter>.
        <converter> is a dictionnary that translates from taxon names to powers of 2.
        <domain> is the sum of the values of the taxa on which the clade is defined; its validity domain.
        If <root> is specified, it will be used as explicit root for the clade.
        If <brlen> is provided, it will be recorded in the Clade object."""
    assert parenth[0] == "(" and parenth[-1] == ")", "A newick clade must be delimited by parentheses."
    if root is not None:
        assert root & domain, "The root must be in the validity domain."
    subclades = []
    inclade = 0
    found = False
    for t in converter:
        r = regtaxname(t)
        if r.search(parenth) is not None:
            inclade = inclade + converter[t]
    if found:
        print(f"The value of {parenth} is now {inclade}")
    clade = Clade(inclade, domain - inclade, root, brlen)
    parlevel = 1
    i = 0
    while parlevel != 0:
        i = i + 1
        if parlevel == 1 and parenth[i] not in ["(", ")", ","]:  # A first-level terminal taxon is begining.
            termval, brlen, i = read_terminal(parenth, i, converter)
            if termval is not None:
                subclades.append(Clade(termval, domain - termval, root, brlen))
        if parenth[i] == "(":
            parlevel += 1
            if parlevel == 2:  # A first-level sub clade has just been opened.
                ibegin = i
        if parenth[i] == ")":
            parlevel -= 1
            if parlevel == 1:  # A first-level subclade has just been closed.
                iclose = i
                subparenth = parenth[ibegin:iclose + 1]
                brlen, i = read_brlen(parenth, i)  # override brlen with value for new subclade
                subclades.append(read_parenth(subparenth, converter, domain, root, brlen))
    clade.addsubclades(subclades)
    return clade


def get_root(interactive, taxons, converter, missingvalues=[]):
    """get_root(interactive, missingvalue) returns the name of the root taxon, chosen by the user if <interactive> is True, or Arbitrarily chosen otherwise.
        <taxons> is the list of taxon names.
        <converter> is a dictionnary translating from taxon names to powers of 2.
        <missingvalues> is a list of the values of the taxa that cannot be chosen as root."""
    assert isinstance(taxons, list), "taxons is supposed to be a list of taxon names."
    assert isinstance(converter, dict), "converter should be a dictionnary converting from taxon names to powers of 2."
    root = False
    if interactive:
        while not root:
            for t in taxons:
                if converter[t] not in missingvalues:
                    sys.stdout.write("%s " % t)
            r = input("\nPlease choose a root for this dataset\n")
            if r in taxons and converter[r] not in missingvalues:
                root = r
    i = 0
    while not root:
        try:
            r = taxons[i]
            if converter[r] not in missingvalues:
                root = r
                sys.stdout.write(f"root arbitrarily set to {r}\n")
            i = i + 1
        except IndexError as e:
            print("All taxa seem to be missing; no valid root can be chosen.", e)
    return root


def sort_clades(clades, attr):
    """sort_clades(clades, attr) returns the sorted list corresponding to the list <clades>, sorted according to the
    attribute <attr> of the clades in the list. The clades with the highest value of attr are placed on top of the
    list. """
    # implementation taken from module func2 of package p4
    def pairing(clade, attr=attr):
        """pairing(clade, attr) returns the pair clade.attr, attr."""
        return getattr(clade, attr), clade

    pairs = [pairing(clade) for clade in clades]
    pairs.sort()

    def strip(pair):
        """strip(pair) returns the second element of a pair."""
        return pair[1]

    cmplist = [strip(pair) for pair in pairs]
    return cmplist


def consensusclades(clades, type):
    """consensusclades(clades, type) returns a selection of the clades in list <clades>.
        The selection is made to obtain a consensus of the type specified by <type>."""
    sorted_clades = sort_clades(clades, "weight")
    sorted_clades.reverse()  # The clades with the highest weight come first.
    selection = []
    i = 0
    while sorted_clades[i].weight == 1:
        selection.append(sorted_clades[i])
        i += 1
        if i == len(sorted_clades):
            break
    if type == "strict":
        return selection  # else try adding more clades to the selection
    if i < len(sorted_clades):
        while sorted_clades[i].weight > 0.5:
            selection.append(sorted_clades[i])
            i += 1
            if i == len(clades):
                break
    if type == "majrule":
        return selection
    if i < len(sorted_clades):
        while sorted_clades[i]:
            if sorted_clades[i].allcompat(
                    selection):  # If clades with the same weight are not compatible, this greedy consensus
                # arbitrarily selects the clade that appears first; the result is thus dependent on the way the sort
                # method works.  # TODO: is this maybe related to the sorting problem? :')
                selection.append(sorted_clades[i])
            i += 1
            if i == len(sorted_clades):
                break
    if type == "greedy":
        return selection
    else:
        # TODO: fstring
        raise NotImplementedError("%s consensus not implemented.\n"
                                  "Possible types are strict, majrule and greedy.") % type


################
# Main program #
################

def main():
    sys.stdout.write(
        "\nSuperTRI\n%s\n\nA python script to help assembling supertree matrices with clade weights\n\nDistributed "
        "under the GNU General Public Licence\n\n" % __version__[1:-1])
    # TODO: alternatives for psyco? not used in python3

    try:
        opts, args = getopt.getopt(sys.argv[1:], "io:",
                                   ["help", "licence", "datasets=", "intree=", "root=", "suffix=", "taxlist="])
    except ImportError:
        print(__doc__)
        sys.exit(1)
    greedy = True
    interactive = False
    intreefiles = []
    markfich = None
    outfich = None
    r = None
    suffix = None
    taxfich = None
    for o, a in opts:  # TODO : use argparse ?
        if o == "-i":
            interactive = True
        if o == "--help":
            print(__doc__)
        if o == "--licence":
            print(__licence__)
            sys.exit()
        if o == "--datasets":
            markfich = a
        if o == "-o":
            outfich = a
        if o == "--intree":
            intreefiles.append(a)
        if o == "--root":
            r = a
        if o == "--suffix":
            suffix = a
        if o == "--taxlist":
            taxfich = a
    files = os.listdir(".")
    toopen = []
    marks = set()
    if interactive:
        if not markfich:
            if oui("Do you want to enter the names of the markers used in the analyses ?"):
                nmarks = int(input("How many different markers are there ? "))
                for i in range(1, nmarks + 1):
                    marks.add(input(f"What name do you use for marker {i}"))

            else:
                sys.stdout.write("Then, you'll have to provide a file with the names.\n")
                markfich = askfile("""
                Please enter the name of a file containing the names of the markers you used in the analyses.
                The names must contain no whitespaces and must be separated by whitespaces.
                """, files)

    if markfich:  # TODO: with open
        fich = open(markfich, 'r')
        for m in fich.read().split():
            marks.add(m)
        fich.close()
        nmarks = len(marks)
        marks = sorted(marks)

    elif len(marks) == 0:
        print(__doc__)
        raise Exception("No dataset names provided. Either use option -i or provide a file with option --datasets.")

    sys.stdout.write(f"\nThere are {nmarks} markers:\n")
    sys.stdout.write(", ".join(marks) + "\n")
    if interactive:
        attend()
        # obtaining the taxon names
        if not taxfich:
            taxfich = askfile("""
        Please enter the name of a file containing the names of the taxa in the order they were in the matrix you used to do the analyses.
        The names must contain no whitespaces and must be separated by whitespaces.
        """, files)
        if not suffix:
            suffix = input(
                "Please enter the suffix of the files containing the results to be read [.parts/.log]\n")
            if suffix == "":
                suffix = ".parts"  # default value, for a use with MrBayes
    try:
        if suffix[0] != ".":
            suffix = "." + suffix
    except:
        print("No file extension specified. Defaulting to MrBayes .parts file.")
        suffix = ".parts"

    for f in files:
        if f.endswith(suffix):
            if interactive:
                if oui("Is file %s supposed to be read ?" % f):
                    toopen.append(f)
            else:  # The default is to take files matching exactly marker + suffix
                # else: # The default is to take all files matching
                for m in marks:
                    if f == m + suffix:
                        toopen.append(f)
                        break

    sys.stdout.write(f"\nTaxon names will be read in {taxfich}\n")
    if taxfich is not None:
        with open(taxfich, 'r') as filin:
            taxons = filin.read().split()
    else:
        print(__doc__)
        raise Exception("No taxon names provided. Either use option -i or provide a file with option --taxlist.")
    sys.stdout.write(f"The {len(taxons)} taxa that have been read are, in this order:\n")
    sys.stdout.write(", ".join(taxons) + "\n")

    # Building dictionnaries to convert taxon names to powers of 2 and vice-versa
    puiss = 1
    val2tax = {}  # to convert from powers of 2 to names
    tax2val = {}  # to convert from names to powers of 2

    global maxtaxname
    maxtaxname = 0  # used to align the supertree matrix

    for t in taxons:
        if len(t) > maxtaxname:
            maxtaxname = len(t)
        val2tax[puiss] = t
        tax2val[t] = puiss
        puiss = 2 * puiss
    tot_val = puiss - 1  # value of the set of all taxa
    # A taxon is uniquely defined by its value as a power of 2
    sys.stdout.write("\nThe files in which the bipartions and their weights will be read are:\n")
    sys.stdout.write(", ".join(toopen) + "\n")
    # dictionnary whose keys are datasets names and the associated values are the trees produced by the dataset
    datasetsdict = {}
    B = BipartSet()
    B.set_converter(tax2val)
    # Loop over the files (and hence, to the markers):
    for file in toopen:
        if suffix == ".parts":
            new_biparts = read_biparts(file)
        else:
            try:
                new_biparts = read_bootlog(file, len(taxons))
            except:
                print("Is %s really a PAUP bootstrap log file ?\nTrying to read it as MrBayes .parts file." % file)
                try:
                    new_biparts = read_biparts(file)
                except:
                    print("For the moment, only PAUP log files or MrBayes .parts files are allowed.")
                    sys.exit(1)
        abs_files = []
        data_markers = []
        for mark in marks:
            if (mark + suffix == file) or (mark in file and interactive):
                data_markers.append(mark)
                while len(data_markers) == 2:
                    print("Name conflict. Which marker corresponds to file %s ?\n" % file)
                    try:
                        num = int(input("1) %s\n2) %s\n" % (data_markers[0], data_markers[1]))) - 1
                    except ValueError:
                        num = 2
                    if num in [0, 1]:
                        data_markers.remove(data_markers[1 - num])
                    else:
                        print("Wrong number. Try again!")
                if mark + ".abs" in files:
                    abs_files.append(mark + ".abs")
                    sys.stdout.write(
                        "%s was found. Expecting to find the names of the taxa missing for dataset %s in it.\n" % (
                            m + ".abs", m))
                    if interactive:
                        attend()
        if len(data_markers) == 0:
            raise Exception("The file %s cannot be assigned to one of the declared datasets. Check the spelling.\n"
                            "The declared datasets were: %s" % (file, ", ".join(marks)))  # TODO fstring ????
        assert len(data_markers) == 1, "No partial combinations are allowed."
        assert len(abs_files) < 2, "There should not be more than one file specifying missing taxa for a particular " \
                                   "dataset. "
        clades = []
        maj_clades = []
        root = False
        if len(abs_files) == 1:
            fich = open(abs_files[0], 'r')
            items = fich.read().split()  # list of the values corresponding to the taxa read in absfile[0]
            missingvalues = [tax2val.get(item) for item in items]
            # total value of the missing taxa, obtained by summing via reduce() function
            missingtot = reduce(lambda x, y: x + y, missingvalues)
            fich.close()
            if r in taxons and tax2val[r] not in missingvalues:  # r is the root proposed by the user.
                root = r
            if not root:
                root = get_root(interactive, taxons, tax2val, missingvalues)
            majruletree = Clade(tot_val - (tax2val[root] + missingtot), tax2val[root], tax2val[root])
            # The tree must be defined on a domain not comprising the missing values.
            # The sum of the missing values is therefore substracted from the tentative outclade value.

            # Set the support type:
            if suffix == ".log":
                majruletree.weight = Bootstrap(0)
                majruletree.support = Bootstrap(0)
            elif suffix == ".parts":
                majruletree.weight = PosteriorP(0)
                majruletree.support = PosteriorP(0)
            else:
                majruletree.weight = NodeIndex(0)
                majruletree.support = NodeIndex(0)
            # Loop over bipartitions read in the file:
            for b in new_biparts:
                inclade, outclade = val_stars(b[0], missingvalues)
                bip = Bipartition(inclade, outclade)
                if suffix == ".log":
                    # The bootstrap proportion for bip is set as its support and as its weight.
                    bip.update(Bootstrap(b[1]), data_markers)
                elif suffix == ".parts":
                    # The posterior probability for bip is set as its support and as its weight.
                    bip.update(PosteriorP(b[1]), data_markers)
                else:
                    bip.update(NodeIndex(b[1]), data_markers)
                clades.append(bip.clade(tax2val[root]))
                # The clade that is made inherits from bip's bootstrap proportion.
                assert bip.weight == clades[-1].weight, "The clade has not the same bootstrap value as his mother " \
                                                        "bipartition. This should not happen. "
                # Majority-rule selection:
                if bip.weight > 0.5:
                    maj_clades.append(bip.clade(tax2val[root]))
                B.add_bipart(bip)  # If B already contains bip, the existing record is updated.
        else:
            if r in taxons:
                root = r
            if not root:
                root = get_root(interactive, taxons, tax2val)
            # After defining the root, a tree is initiated with the "total clade": the one comprising all the taxa:
            majruletree = Clade(tot_val - tax2val[root], tax2val[root], tax2val[root])
            # Set the support type:
            if suffix == ".log":
                majruletree.weight = Bootstrap(0)
                majruletree.support = Bootstrap(0)
            elif suffix == ".parts":
                majruletree.weight = PosteriorP(0)
                majruletree.support = PosteriorP(0)
            else:
                majruletree.weight = NodeIndex(0)
                majruletree.support = NodeIndex(0)
            # Loop over bipartitions read in the file:
            for b in new_biparts:
                inclade, outclade = val_stars(b[0])
                bip = Bipartition(inclade, outclade)
                if suffix == ".log":
                    bip.update(Bootstrap(b[1]), data_markers)
                elif suffix == ".parts":
                    bip.update(PosteriorP(b[1]), data_markers)
                else:
                    bip.update(NodeIndex(b[1]), data_markers)
                clades.append(
                    bip.clade(tax2val[root]))  # The clade that is made inherits from bip's bootstrap proportion.
                assert bip.weight == clades[-1].weight, "The clade has not the same support value as his mother " \
                                                        "bipartition. This should not happen. "
                # Majority-rule selection:
                if bip.weight > 0.5:
                    maj_clades.append(bip.clade(tax2val[root]))
                B.add_bipart(bip)
        clades.sort()
        maj_clades.sort()  # majority-rule consensus; contains only the clades appearing in more than half the profile
        if greedy:  # majclades is overriden, it was there for a test version
            maj_clades = consensusclades(clades, "greedy")
        majruletree.addsubclades(maj_clades)
        datasetsdict[data_markers[0]] = majruletree
    if suffix == ".parts":
        tmpfmt = "%.2f"
        PosteriorP(0).reformat("%s")
        matrix = B.nex_matrix_rep(taxons, set_name="Cumulated_posterior_probabilities")
        PosteriorP(0).reformat(tmpfmt)
    elif suffix == ".log":
        tmpfmt = "%.0f"
        Bootstrap(0).reformat("%s")
        matrix = B.nex_matrix_rep(taxons, set_name="Cumulated_bootstrap_proportions")
        Bootstrap(0).reformat(tmpfmt)
    else:
        tmpfmt = "%.2f"
        NodeIndex(0).reformat("%s")
        matrix = B.nex_matrix_rep(taxons)
        NodeIndex(0).reformat(tmpfmt)

    if interactive:
        confirm = False
        while not confirm:
            outfich = askfile("Please enter a name for the file to save the matrix representation. ")
            if outfich in files:
                confirm = oui("File %s already exists. Do you want to overwrite it ?")
            else:
                confirm = True
    try:
        fich = open(outfich, 'w')
        fich.write(matrix)
        fich.close()
        sys.stdout.write("\nFile %s written\n" % outfich)
    except:
        filename = "".join(marks) + (suffix == ".parts") * ".mrpp" + (suffix == ".log") * ".mrbp" + ".nex"
        sys.stdout.write("\nThe matrix will be written in file %s\n" % filename)
        fich = open(filename, 'w')
        fich.write(matrix)
        fich.close()

    ###########################################################################
    # Reading trees, calculating the indices, and reporting them on the trees #
    ###########################################################################

    for tree_file in intreefiles:
        out_file_name = tree_file + ".withindices"
        output = "#Nexus\n\n"
        new_trees = readnextrees(tree_file, tax2val)
        for tree in new_trees:
            if suffix == ".log":
                tree.set_support_type(Bootstrap(0).__class__)
            elif suffix == ".parts":
                tree.set_support_type(PosteriorP(0).__class__)
            else:
                tree.set_support_type(NodeIndex(0).__class__)
            tree.compute_indices(datasetsdict.values(), B)  # maybe pb is here
            output += "begin trees;\n"
            output += tree.next_tree(val2tax, "meansup", False) + "\n"
            output += "end;\n"
            output += "begin trees;\n"
            output += tree.next_tree(val2tax, "repro", False) + "\n"
            output += "end;\n"
        sys.stdout.write(output)
        with open(out_file_name, "w") as filout:
            filout.write(output)
        sys.stdout.write("\nThese trees have been written in %s\n" % out_file_name)
        tree_number = 0
        for tree in new_trees:
            fich = open(tree_file + ".%s.tgf" % tree_number, 'w')
            fich.write(tree.tgftree(val2tax))
            fich.close()
            sys.stdout.write(
                "A .tgf version with mean support and reproducibility has been written in " + tree_file + ".%s.tgf\n" % tree_number)
            tree_number += 1


if __name__ == "__main__":
    main()