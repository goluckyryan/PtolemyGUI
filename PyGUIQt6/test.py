#!/usr/bin/python3

from IAEANuclearData import IsotopeClass

iso = IsotopeClass()

iso.PrintIso('16O')

pd = iso.GetExList('16O', 10)

print(pd)