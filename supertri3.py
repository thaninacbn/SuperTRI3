import copy
import getopt
import os
import re
import warnings


def oui(question):
    return (input("\n" + question + " (y/n) ") in
            ["O", "o", "Oui", "oui", "OUI", "Y", "y", "Yes", "yes", "YES", "J", "j", "Ja", "ja", "JA", "D", "d", "Da",
             "da", "DA"])


def attend():
    pause = input("\npress enter to continue...")
    return 0


def askfile(question, list_of_files=None):
    return 0


class NodeIndex(float):
    format = "%.2f"

    def __init__(self, initval):
        float.__init__(initval)

    def __str__(self):
        return self.__class__.format % float(self)

    def reformat(self, newformat):
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
    """This object represents a node bootstrap percentage."""
    format = "%.0f"

    def __init__(self, initval):
        NodeIndex.__init__(self, initval)

    def __str__(self):
        return self.__class__.format % (100 * float(self))  # Bootstrap values are displayed as percentages.


class PosteriorP(NodeIndex):
    """This object represents a node posterior probability."""
    format = "%.2f"

    def __init__(self, initval):
        NodeIndex.__init__(self, initval)


class Bipartition(object):
    """put useful docstrings here"""

    def __init__(self, inclade, outclade, brlen=None):
        self.datasets = set([])  # initiate empty set out of a list (an iterable))
        self.weights = {}
        self.weight = NodeIndex(0)
        self.support = NodeIndex  # this is useless no?
        self.meansup = NodeIndex(0)  # again isnt this useless?
        self.occurences = 0
        self.repro = 0
        assert isinstance(inclade, int) and isinstance(outclade, int)
        assert inclade and outclade == 0, "the two parts of a bipartition cannot have shared taxa"
        self.id = frozenset([inclade, outclade])
        self.inclade = inclade
        self.outclade = outclade
        self.brlen = brlen
        self.domain = inclade + outclade
        self.quartets = None

    def __str__(self):
        return f"Bipart ({self.inclade},{self.outclade})"

    def compatible(self, other):
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
        """self.allcompat(listofothers) returns True if self is compatible with all the elements of the list
        <listofothers>. """
        assert isinstance(list_of_others, list), "allcompat only works with lists."
        result = True
        i = 0
        # a voir si on peut pas changer ça en une boucle for ?
        while result:  # No need to continue once one incompatibility has been found.
            # while result and listofothers[i]:
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
                # self.weights[m] = self.weights.get(m, 0) + w
        ## The order of the terms to be added might matter; the expression above needs to be re-checked.
        self.weight = w + self.weight
        # adding in this order, self.weight should inherit the type of w
        # print self.weight.__class__
        self.support = w + self.support  # TODO why is this like this????
        # adding in this order, self.support should inherit the type of w
        # print self.support.__class__
        # if w > 0.5:
        #    self.occurrences = self.occurrences + 1
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
    """This object represents a clade, that is an oriented bipartition.
    It also has a dictionnary of subclades, to optionally specify its internal topology."""

    def __init__(self, inclade, outclade, root=None, brlen=None):
        Bipartition.__init__(self, inclade, outclade, brlen)
        self.id = (self.inclade,
                   self.outclade)  # id is a tuple instead of a frozenset; this allows the distinction between a
        # clade and its complementary.
        assert (root is None) or (
                root & outclade == root), "If explicit, the root should be one of the taxa outside the clade."
        self.root = root  # optional attribute that indicates with respect to which taxon the clade is rooted
        self.subclades = {}  # # Peut-etre faudrait-il definir une sous-classe de Clade, avec la topologie definie au
        # lieu d'ajouter cet attribut. Mais que mettre comme id a un tel clade ?
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

        ## Beware of the implication of this method on the use of the "in" operator.

    def __contains__(self, other):
        """self.__contains__(other) returns True if Clade <other> can be a subclade of self, False otherwise.
        Empty clades are considered included in any clade that shares the same validity domain. Inclusion implies compatibility."""
        assert isinstance(other, Clade), "To state about inclusion, the object must be a Clade type object."
        if other.domain == self.domain:
            inter = other.inclade & self.inclade
            ## Verification
            if other.inclade == inter:
                assert self.compatible(other), "other is in self, so it should be compatible."
            return other.inclade == inter
        else:
            return False

    def __cmp__(self, other):
        """self.__cmp__(other) is called when <, ==, or > are called to compare self and other. Based on this __cmp__
        implementation, the order relation between clades with the same validity domain that are compatible will
        place clades included in other clades before these other clades. """
        if self.id == other.id:
            return 0  # Clades are equal if they have the same outside and the same inside
        if self.domain == other.domain:
            return cmp(self.inclade, other.inclade)
        else:
            raise ValueError("To be compared, two clades must have the same validity domain, that is, concern the "
                             "same set of taxa.")

    def __str__(self):
        return f"Clade {self.inclade},{self.outclade})"

    def addsubclades(self, subclades):

        """self.addsubclades(subclades) recursively includes the clades in the list <subclades> into self.
        This method should not be called with a partial list of subclades, since all taxa of the clade not present in the subclades are considered being placed in an unresolved basal position (a "rake")."""

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

    