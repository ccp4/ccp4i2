"""
     python/ui/selection_protocols.py: CCP4MG Molecular Graphics Program
     Copyright (C) 2001-2008 University of York, CCLRC
     Copyright (C) 2009 University of York

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""
generic_selection = {
        'all'            :  [None,'/*/*/*,*/*:*'],
        'peptide'	:  [None,'/*/*/(%amino_acid%)'] ,
        'amino_acid'	:  [None,'/*/*/(%amino_acid%)'] ,
	'nucleic'	: [None,'/*/*/(%nucleic_acid%)'],
	'saccharide'	: [None,'/*/*/(%saccharide%)'],
	'solvent'	: [None,'/*/*/(%solvent%)'],
	'solute'	: [None,'/*/*/(%solute%)'],
	'catrace'	:  [None,'/*/*/(%amino_acid%)/CA:*'] ,
        'trace'	        :  [None,'/*/*/(%nucleic_acid%)/P:*'] ,
	'main'		:  [None,'/*/*/(%amino_acid%)/%main_chain%:*'] ,
	'peptide_link'	:  [None,'/*/*/(%amino_acid%)/%peptide_link%:*'] ,
	'backbone'	:  [None,'/*/*/(%nucleic_acid%)/%backbone_atoms%:*'],
	'sugar_ring'	:  [None,'/*/*/(%nucleic_acid%)/%sugar_ring_atoms%:*'],
        'CA+side'       :  ['and',[None,'amino_acid'],['not',[None,'/*/*/*.*/%not_CA_side%:*']]] ,
	'side'		:  ['and',[None,'amino_acid'],['not', [None,'main']]],
	'base'		:  ['and',[None,'/*/*/(%nucleic_acid%)'],['not', [None,'/*/*/*,*/%not_base%:*']]],
        'polar_atoms'   :  [None,'[N,O,P]'],
        'metal'         :  [None, '[%metal%]'],
        'allnotsolv'     :   ['excl',[None,'all'],['or',[None,'solvent'],[None,'solute']]],
        'helix'         : [None, 'SEC == helix'],
        'alpha-helix'   : [None, 'SEC == helix'],
        'strand'        : [None, 'SEC == strand'],
        'beta-strand'        : [None, 'SEC == strand'],
        '3turn'         : [None,'SEC == 3turn'],
        '4turn'         : [None,'SEC == 4turn'],
        '5turn'         : [None,'SEC == 5turn'],
        'beta-bulge'    : [None,'SEC == bulge'],
        'bulge'         : [None,'SEC == bulge']
        }

metal = { 'ORDER' : ['AR','CA','CD','CO','CR','CU','FE','GD','HG','IR','K','MG','MN','NA','UR','ZN'] }

solvent = {
    'ORDER' : ['HOH','H2O','WAT','SOL','TIP','DOD','D2O', 'DUM'] }
    
dummy = {
    'ORDER' : ['DUM'] }

solute = {
    'ORDER' : ['SUL','SO4','NO3','MOH','EOH','GOL', 'EDO','PO4','CL','ATP','PHO'],
    'SUL' : 'Sulphate' }

main_chain = {
    'ORDER' : ['N','CA','C','O','H','HA' ] }

peptide_link = {
    'ORDER' : ['N','C','O','H' ] }

backbone_atoms = {
    'ORDER' : ['C3*','O3*','H3*','C4*','C5*','1H5*','O5*','2H5*','H5T','P','O1P','O2P',"C3'","O3'","H3'","C4'","C5'","1H5'","O5'","2H5'",'H7','P','OP1','OP2',"HO5'","HO3'" ] }

dummy_atoms = { 'ORDER' : ['DUM'] }

not_base = {
    'ORDER' : ['C3*','O3*','H3*','C5*','1H5*','O5*','2H5*','H5T','P','O1P','O2P',"C3'","O3'","H3'","C5'","1H5'","O5'","2H5'","H7",'P','OP1','OP2',"HO5'","HO3'" ] }

sugar_ring_atoms = {
    'ORDER' : ['C1*','H1*','C2*','2H2*','1H2*','C3*','O3*','H3*','C4*','H4*', \
               'C5*','1H5*','2H5*','O5*','H5T','O4*','O2*',
               "C1'","H1'","C2'","2H2'","1H2'","C3'","O3'","H3'","C4'","H4'", \
               "C5'","1H5'","2H5'","O5'","H5T'","O4'","O2'"] }

not_CA_side = {
    'ORDER' : ['N','C','O','H','HA' ] }

saccharide = {
    'ORDER' :  ['GLC','MAN','NAG','RAM','RIB','XYS'], \
    'GLC' : 'Glucose',
    'MAN' : 'Mannose',
    'NAG' : 'N-acetyl-Glucosamine',
    'RAM' : 'Rhamnose',
    'RIP' : 'Ribose',
    'XYS' : 'Xylose' 
    }

amino_acid = {
    'ORDER' :  [ 'GLY', 'ALA', 'VAL', 'PRO', 'SER', 'THR' , \
                'LEU', 'ILE', 'CYS', 'ASP', 'GLU', 'ASN', 'GLN', \
                'ARG', 'LYS',  'MET', 'MSE', 'HIS', 'PHE', 'TYR', 'TRP','UNK'], \
    'GLY' : 'Glycine',
    'ALA' : 'Alanine',
    'VAL' : 'Valine',
    'PRO' : 'Proline',
    'SER' : 'Serine',
    'THR' : 'Threonine',
    'LEU' : 'Leucine',
    'ILE' : 'Isoleucine',
    'CYS' : 'Cysteine',
    'ASP' : 'Aspartic acid',
    'GLU' : 'Glutamic acid',
    'ASN' : 'Asparagine',
    'GLN' : 'Glutamine',
    'ARG' : 'Arginine',
    'LYS' : 'Lysine',
    'MET' : 'Methionine',
    'MSE' : 'Se-Methionine',
    'HIS' : 'Histidine',
    'PHE' : 'Phenylalanine',
    'TYR' : 'Tyrosine',
    'TRP' : 'Tryptophan',
    'UNK' : 'Unknown'
    }


amino_acid_groups = { \
     'Small' : [ 'GLY', 'ALA', 'VAL', 'PRO', 'SER', 'THR' ], \
     'Large' : [ 'MET', 'MSE', 'HIS', 'PHE', 'TYR', 'TRP' ], 
     'Hydrophobic': [ 'ALA', 'VAL', 'PRO', 'LEU', 'ILE', 'MET', \
                      'PHE','TRP' ],
     'Hydrophilic' : ['ASP','GLU','TYR','ASN','GLN','THR','SER', \
                        'CYS','HIS','LYS','ARG' ], 
     'Acidic': [ 'ASP','GLU' ],
     'Basic' : [ 'HIS','LYS','ARG']
    }

nucleic_acid = {
     'ORDER' : ['DA','DC','DG','DT','Ad','Cd','Gd','Td','A', 'G', 'I' , 'U', 'T', 'C','ADE','CYT','GUA','INO', 'THY','URA', 'YG','PSU'], \
    'DA' : 'Adenine',
    'DC' : 'Cytosine',
    'DG' : 'Guanine',
    'DT' : 'Thymine',
    'Ad' : 'Adenine',
    'Cd' : 'Cytosine',
    'Gd' : 'Guanine',
    'Td' : 'Thymine',
    'I' : 'Inosine',
    'A' : 'Adenine',
    'C' : 'Cytosine',
    'G' : 'Guanine',
    'I' : 'Inosine',
    'T' : 'Thymine',
    'U' : 'Uracil',
    'ADE' : 'Adenine',
    'CYT' : 'Cytosine',
    'GUA' : 'Guanine',
    'INO' : 'Inosine',
    'THY' : 'Thymine',
    'URA' : 'Uracil',
    'YG': 'Wybutosine',
    'PSU' :  'Pseudouridine'
    }

nucleic_acid_groups = { \
     'Pyrimidines': [ 'Cd','Td','DC','DT','U', 'T', 'C', 'PSU', 'URA', 'THY', 'CYT' ] , \
    'Purines': [ 'Ad','Gd','DA','DG','A', 'G', 'I' , 'ADE', 'GUA' ]  \
    }


#element = {
#    'ORDER' : [ 'C', 'N','O','H','S','P','CA'],
#    'C' : 'carbon' ,
#    'N' : 'nitrogen' ,
#    'O' : 'oxygen' ,
#    'H' : 'hydrogen' ,
#    'S' : 'sulphur' ,
#    'P' : 'phosporus',
#    'CA' : 'calcium'
#    }
element = {
    'ORDER' : [ ' C', ' N',' O'],
    ' C' : ' C' ,
    ' N' : ' N' ,
    ' O' : ' O' ,
    }

element_groups = { }

nucleic_acid_code = { 'GLY' : 'G',
    'DA' : 'A',
    'DC' : 'C',
    'DG' : 'G',
    'DT' : 'T',
    'Ad' : 'A',
    'Cd' : 'C',
    'Gd' : 'G',
    'Td' : 'T',
    'I' : 'I',
    'A' : 'A',
    'C' : 'C',
    'G' : 'G',
    'I' : 'I',
    'T' : 'T',
    'U' : 'U',
    'ADE' : 'A',
    'CYT' : 'C',
    'GUA' : 'G',
    'INO' : 'I',
    'THY' : 'T',
    'URA' : 'U',
    'YG': 'W',
    'PSU' :  'Q'

}

amino_acid_code = { 'GLY' : 'G',
		    'PRO' : 'P',
		    'ALA' : 'A',
		    'VAL' : 'V',
		    'LEU' : 'L',
		    'ILE' : 'I',
		    'MET' : 'M',
		    'MSE' : 'M',
		    'CYS' : 'C',
		    'PHE' : 'F',
		    'TYR' : 'Y',
		    'TRP' : 'W',
		    'HIS' : 'H',
		    'LYS' : 'K',
            'GLN' : 'Q',
            'ASN' : 'N',
		    'ARG' : 'R',
		    'GLU' : 'E',
		    'ASP' : 'D',
		    'SER' : 'S',
		    'THR' : 'T',
            'YPO' : 'Y',
            'TPO' : 'T'  }
