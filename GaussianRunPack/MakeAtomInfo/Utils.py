from collections import OrderedDict

Multiplicity = OrderedDict([
  ("H" ,2),
  ("He",1),
  ("Li",2),
  ("Be",1),
  ("B" ,2),
  ("C" ,3),
  ("N" ,4),
  ("O" ,3),
  ("F" ,2),
  ("Ne",1),
  ("Na",2),
  ("S",3),
  ("Cl",2),
  ("Si" ,3),
  ("P" ,4),
  ("Br" ,2),
  ("I" ,2)
])

Basis = [
  "LANL2DZ",
  "STO-3G",
  "3-21G",
  "6-31G",
  "6-311G", 
  "3-21G*",
  "3-21+G*",
  "6-31G*", 
  "6-311G**",
  "6-31G**",
  "6-31+G*",
  "6-31+G**",
  "6-311+G*",
  "6-311+G**"
]

Functional = [
  "BLYP",
  "B3LYP", 
  "X3LYP",
  "LC-BLYP",
  "CAM-B3LYP"
]
