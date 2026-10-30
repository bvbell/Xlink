# Load the inspection set with the macrocycle bond drawn.
#   pymol view.pml
# Rosetta writes no connectivity record for the cyclisation, so without these
# bond commands every design renders as an open chain. The geometry is fine:
# the closure C-N measures 1.33 A at 116.2/121.7 deg in every amide design.

load apj/good/apj_i1__medium__amide__DGLN-GLY__c00.pdb, apj_i1__medium__amide__DGLN-GLY__c00
bond apj_i1__medium__amide__DGLN-GLY__c00 and chain P and resi 12 and name C, apj_i1__medium__amide__DGLN-GLY__c00 and chain P and resi 1 and name N
load apj/good/apj_i1__long__amide__DTYR-GLY__c00.pdb, apj_i1__long__amide__DTYR-GLY__c00
bond apj_i1__long__amide__DTYR-GLY__c00 and chain P and resi 14 and name C, apj_i1__long__amide__DTYR-GLY__c00 and chain P and resi 1 and name N
load apj/good/apj__short__amide__ASN-GLY__c01.pdb, apj__short__amide__ASN-GLY__c01
bond apj__short__amide__ASN-GLY__c01 and chain P and resi 10 and name C, apj__short__amide__ASN-GLY__c01 and chain P and resi 1 and name N
load apj/bad/apj_i1__medium__amide__DILE-GLY__c01.pdb, apj_i1__medium__amide__DILE-GLY__c01
bond apj_i1__medium__amide__DILE-GLY__c01 and chain P and resi 12 and name C, apj_i1__medium__amide__DILE-GLY__c01 and chain P and resi 1 and name N
load apj/bad/apj_i1__medium__amide__DALA-GLY__c01.pdb, apj_i1__medium__amide__DALA-GLY__c01
bond apj_i1__medium__amide__DALA-GLY__c01 and chain P and resi 12 and name C, apj_i1__medium__amide__DALA-GLY__c01 and chain P and resi 1 and name N
load apj/bad/apj_i1__medium__amide__DASP-GLY__c03.pdb, apj_i1__medium__amide__DASP-GLY__c03
bond apj_i1__medium__amide__DASP-GLY__c03 and chain P and resi 12 and name C, apj_i1__medium__amide__DASP-GLY__c03 and chain P and resi 1 and name N
load apj/good/apj__short__ssapp__CYD-DCYD__c02.pdb, apj__short__ssapp__CYD-DCYD__c02
# apj__short__ssapp__CYD-DCYD__c02 closes through S-S; PyMOL draws disulfides automatically
load apj/good/apj__short__ssapp__DCYD-CYD__c01.pdb, apj__short__ssapp__DCYD-CYD__c01
# apj__short__ssapp__DCYD-CYD__c01 closes through S-S; PyMOL draws disulfides automatically
load apj/good/apj__medium__ssapp__DCYD-CYD__c03.pdb, apj__medium__ssapp__DCYD-CYD__c03
# apj__medium__ssapp__DCYD-CYD__c03 closes through S-S; PyMOL draws disulfides automatically
load apj/bad/apj__medium__ssapp__DCYD-CYD__c04.pdb, apj__medium__ssapp__DCYD-CYD__c04
# apj__medium__ssapp__DCYD-CYD__c04 closes through S-S; PyMOL draws disulfides automatically
load apj/bad/apj__long__ssapp__CYD-CYD__c01.pdb, apj__long__ssapp__CYD-CYD__c01
# apj__long__ssapp__CYD-CYD__c01 closes through S-S; PyMOL draws disulfides automatically
load apj/bad/apj__medium__ssapp__DCYD-CYD__c02.pdb, apj__medium__ssapp__DCYD-CYD__c02
# apj__medium__ssapp__DCYD-CYD__c02 closes through S-S; PyMOL draws disulfides automatically
load uor/good/uor__long__amide__DALA-GLY__c01.pdb, uor__long__amide__DALA-GLY__c01
bond uor__long__amide__DALA-GLY__c01 and chain P and resi 16 and name C, uor__long__amide__DALA-GLY__c01 and chain P and resi 1 and name N
load uor/good/uor_i1__short__amide__DTYR-GLY__c02.pdb, uor_i1__short__amide__DTYR-GLY__c02
bond uor_i1__short__amide__DTYR-GLY__c02 and chain P and resi 10 and name C, uor_i1__short__amide__DTYR-GLY__c02 and chain P and resi 1 and name N
load uor/good/uor_i1__short__amide__DTYR-GLY__c01.pdb, uor_i1__short__amide__DTYR-GLY__c01
bond uor_i1__short__amide__DTYR-GLY__c01 and chain P and resi 10 and name C, uor_i1__short__amide__DTYR-GLY__c01 and chain P and resi 1 and name N
load uor/bad/uor_i1__long__amide__DALA-GLY__c00.pdb, uor_i1__long__amide__DALA-GLY__c00
bond uor_i1__long__amide__DALA-GLY__c00 and chain P and resi 16 and name C, uor_i1__long__amide__DALA-GLY__c00 and chain P and resi 1 and name N
load uor/bad/uor__medium__amide__DPRO-GLY__c00.pdb, uor__medium__amide__DPRO-GLY__c00
bond uor__medium__amide__DPRO-GLY__c00 and chain P and resi 12 and name C, uor__medium__amide__DPRO-GLY__c00 and chain P and resi 1 and name N
load uor/bad/uor__short__amide__DPRO-GLY__c02.pdb, uor__short__amide__DPRO-GLY__c02
bond uor__short__amide__DPRO-GLY__c02 and chain P and resi 10 and name C, uor__short__amide__DPRO-GLY__c02 and chain P and resi 1 and name N
load uor/good/uor__short__ssapp__CYD-CYD__c00.pdb, uor__short__ssapp__CYD-CYD__c00
# uor__short__ssapp__CYD-CYD__c00 closes through S-S; PyMOL draws disulfides automatically
load uor/good/uor__short__ssmut__CYS-CYS__c00.pdb, uor__short__ssmut__CYS-CYS__c00
# uor__short__ssmut__CYS-CYS__c00 closes through S-S; PyMOL draws disulfides automatically
load uor/good/uor__short__ssapp__CYD-CYD__c02.pdb, uor__short__ssapp__CYD-CYD__c02
# uor__short__ssapp__CYD-CYD__c02 closes through S-S; PyMOL draws disulfides automatically
load uor/bad/uor__short__ssapp__CYD-CYD__c03.pdb, uor__short__ssapp__CYD-CYD__c03
# uor__short__ssapp__CYD-CYD__c03 closes through S-S; PyMOL draws disulfides automatically
load uor/bad/uor__medium__ssapp__DCYD-DCYD__c00.pdb, uor__medium__ssapp__DCYD-DCYD__c00
# uor__medium__ssapp__DCYD-DCYD__c00 closes through S-S; PyMOL draws disulfides automatically
load uor/bad/uor__short__ssapp__DCYD-CYD__c01.pdb, uor__short__ssapp__DCYD-CYD__c01
# uor__short__ssapp__DCYD-CYD__c01 closes through S-S; PyMOL draws disulfides automatically

hide everything
show cartoon, chain B
show sticks, chain P
color grey80, chain B
color cyan, chain P
set cartoon_transparency, 0.6
zoom chain P
