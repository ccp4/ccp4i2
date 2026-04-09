#!/usr/bin/env python

###########################################
#
# Dictionary Accumulator
# 2018-2021
# 29/09/2021
#
# Authors: Rob Nicholls and Garib Murshudov
#
###########################################

import sys
import argparse
from gemmi import cif

PRINT_CONTENT_TO_COUT = False
CREATE_NEW_CIF = True
DISPLAY_FINAL_METADATA = False

DUPLICATE_MODIFICATION_FAIL = False
DUPLICATE_LINK_FAIL = False
REQUIRE_ID_EQUALS_TLC = False
ALLOW_MISSING_DATA_BLOCKS = True

def check_identical(block_id1,block_id2,cif1,cif2):
   data_block1 = cif1.find_block(block_id1)
   data_block2 = cif2.find_block(block_id2)
   if data_block1.get_mmcif_category_names() != data_block2.get_mmcif_category_names():
      return False
   for loop_name in data_block1.get_mmcif_category_names():
      category1 = data_block1.find_mmcif_category(loop_name)
      category2 = data_block2.find_mmcif_category(loop_name)
      if len(category1) != len(category2):
         return False
      for i in range(1,len(category1)):
         if tuple(category1[i]) != tuple(category2[i]):
            return False
   return True


def accumulate(input_cif_paths,fileout,opt={}):
 try:
   if not 'unknown' in opt:
      opt['unknown'] = False
    
    
   if PRINT_CONTENT_TO_COUT:
      for path in input_cif_paths:
         doc = cif.read_file(path)
         for block in doc:
            print("BLOCK: ",block.name)
            for loop in block.get_mmcif_category_names():
               print("LOOP: ",loop)
               category = block.find_mmcif_category(loop)
               print("COLS: ")
               for column_id in category.tags:
                  for row in block.find_loop(column_id):
                     print(column_id,row)
               print("ROWS:")
               for row in category:
                  print(row)
               print("")
            print("--- END OF BLOCK ---")
         print("--- END OF FILE ---")

   if CREATE_NEW_CIF:
      print("#######################################")
      print("###### Accumulating dictionaries ######")
      print("#######################################")
      
      # Read input CIF files
      data = []
      for path in input_cif_paths:
         print("Reading file: "+path)
         data.append({"name":path,"data":cif.read_file(path),"comp":[],"mod":[],"link":[]})

      chemcomp_list = []   # lists comps to be included, as well as reference to cif data blocks
      blocks_added = {"comp_list":0,"comp":0,"link_list":0,"link":0,"mod_list":0,"mod":0,"unknown":0}

      USE_COMPLIST = 0
      USE_MODLIST = 0
      USE_LINKLIST = 0

      for el in data:
         input_cif = el["data"]
         INVALID_CIF = 1
         if input_cif.find_block("comp_list"):
            if len(input_cif.find_block("comp_list").get_mmcif_category_names()) > 0:
               USE_COMPLIST = 1
               INVALID_CIF = 0
         if input_cif.find_block("mod_list"):
            if len(input_cif.find_block("mod_list").get_mmcif_category_names()) > 0:
               USE_MODLIST = 1
               INVALID_CIF = 0
         if input_cif.find_block("link_list"):
            if len(input_cif.find_block("link_list").get_mmcif_category_names()) > 0:
               USE_LINKLIST = 1
               INVALID_CIF = 0
         if INVALID_CIF:
            raise Exception("Input CIF does not contain any comp_list mod_list or link_list blocks: {}".format(el["name"]))

      new_cif = cif.Document()

      if USE_COMPLIST:
         added_comps = []
         complist_block = new_cif.add_new_block("comp_list")
         blocks_added["comp_list"] += 1
         required_fields = ["id","three_letter_code","name","group","number_atoms_all","number_atoms_nh"]
         optional_fields = ["desc_level"]
         complist_loop = complist_block.init_mmcif_loop("_chem_comp",required_fields+optional_fields)
         chemcomp_category = complist_block.find_mmcif_category("_chem_comp")
         # Process comp_list from all input CIFs, and get list of all code-CIF pairs that will be combined
         for el in data:
            input_cif = el["data"]
            block = input_cif.find_block("comp_list")
            if block:
               if len(block.get_mmcif_category_names()) == 0:
                  continue
               category = block.find_mmcif_category("_chem_comp")
               print("Preparing to add "+str(len(category))+" components from "+el["name"])

               for (id,three_letter_code) in block.find(['_chem_comp.id','_chem_comp.three_letter_code']):
                  if REQUIRE_ID_EQUALS_TLC:
                     if id != three_letter_code:
                        raise Exception("Comp ID not equal to three letter code: {} {}".format(id,three_letter_code))
                  if id in added_comps:
                     for item in chemcomp_list:
                        if item[0] == id:
                           item.append(input_cif)
                     print("Duplicate code encountered: {}. Only the first instance will be retained.".format(id))
                  else:
                     added_comps.append(id)
                     row = category.find_row(id)
                     # Get rows that match the output tags
                     field_dict = {}
                     for i in range(0,len(row)):
                        field_dict[category.tags.get(i)] = row.get(i)

                     new_row = []
                     for item in ["_chem_comp."+field for field in required_fields]:
                        if item in field_dict.keys():
                           new_row.append(field_dict[item])
                           field_dict.pop(item)
                        else:
                           raise Exception("Incompatible comp_list._chem_comp tags.\nRequired tags: "+" ".join(required_fields)+"\nOptional tags: "+" ".join(optional_fields)+"\nThe following tag is missing: "+item)
                     for item in ["_chem_comp."+field for field in optional_fields]:
                        if item in field_dict.keys():
                           new_row.append(field_dict[item])
                           field_dict.pop(item)
                        else:
                           new_row.append(".")
            
                     # in future, we may want to expand table rather than requiring identical tags (columns)
                     for item in field_dict:
                        print("Warning: unknown tag will be ignored: "+item)
            
                     complist_loop.add_row(new_row)
                     chemcomp_list.append([id,input_cif])
                  el["comp"].append(id)
               
      if USE_LINKLIST:
         linklist_block = new_cif.add_new_block("link_list")
         blocks_added["link_list"] += 1
         linklist_loop = linklist_block.init_mmcif_loop("_chem_link",["id","comp_id_1","mod_id_1","group_comp_1","comp_id_2","mod_id_2","group_comp_2","name"])
         chemlink_category = linklist_block.find_mmcif_category("_chem_link")
         for el in data:
            input_cif = el["data"]
            block = input_cif.find_block("link_list")
            if block:
               if len(block.get_mmcif_category_names()) == 0:
                  continue
               category = block.find_mmcif_category("_chem_link")
               print("Preparing to add "+str(len(category))+" links from "+el["name"])
               if not (set(category.tags) - set(chemlink_category.tags)):
                  for id in block.find_loop('_chem_link.id'):
                     row = category.find_row(id)
                     el["link"].append({"id":id,"row":row})
               else:
                  print(chemlink_category.tags)
                  print(category.tags)
                  raise Exception("Incompatible link_list._chem_link tags")

      if USE_MODLIST:
         modlist_block = new_cif.add_new_block("mod_list")
         blocks_added["mod_list"] += 1
         modlist_loop = modlist_block.init_mmcif_loop("_chem_mod",["id","name","comp_id","group_id"])
         chemmod_category = modlist_block.find_mmcif_category("_chem_mod")
         for el in data:
            input_cif = el["data"]
            block = input_cif.find_block("mod_list")
            if block:
               if len(block.get_mmcif_category_names()) == 0:
                  continue
               category = block.find_mmcif_category("_chem_mod")
               print("Preparing to add "+str(len(category))+" modifications from "+el["name"])
               if not (set(category.tags) - set(chemmod_category.tags)):
                  for id in block.find_loop('_chem_mod.id'):
                     IS_DUPLICATE = 0
                     for new_id in modlist_block.find_loop('_chem_mod.id'):
                        if id==new_id:
                           IS_DUPLICATE = 1
                           break
                     row = category.find_row(id)
                     if IS_DUPLICATE:
                        if DUPLICATE_MODIFICATION_FAIL:
                           raise Exception("Duplicate modification encountered: {}".format(id))
                     el["mod"].append({"id":id,"row":row})
               else:
                  print(chemmod_category.tags)
                  print(category.tags)
                  raise Exception("Incompatible mod_list._chem_mod tags")

# --- COMPONENTS ---

      if USE_COMPLIST:
         # Perform checks:
         # (1) Atom IDs are unique
         # (2) Duplications are compatible
         #     At present just checks for same atom_ids, but in future there will be an option for graph matching
         for code_cif in chemcomp_list:
            data_block_name = "comp_{}".format(code_cif[0])
            block1 = code_cif[1].find_block(data_block_name)
            if block1:
               atom_ids1 = block1.find_loop('_chem_comp_atom.atom_id')
               atom_set1 = set(atom_ids1)
               if len(atom_ids1) != len(atom_set1):
                  raise Exception("Non-unique atom_id found in comp_id: {}".format(code_cif[0]))
               for i in range(2,len(code_cif)):
                  block2 = code_cif[i].find_block(data_block_name)
                  atom_ids2 = block2.find_loop('_chem_comp_atom.atom_id')
                  atom_set2 = set(atom_ids2)
                  if atom_set1 != atom_set2:
                     raise Exception("Different atomic compositions in instances of duplicate comp_id: {}".format(code_cif[0]))
            else:
               if not ALLOW_MISSING_DATA_BLOCKS:
                  raise Exception("Missing data block: "+data_block_name)

         # Add all data_comp_ blocks
         for code_cif in chemcomp_list:
            data_block_name = "comp_{}".format(code_cif[0])
            data_block = code_cif[1].find_block(data_block_name)
            if data_block:
               print("Adding component to dictionary: {}".format(code_cif[0]))
               new_block = new_cif.add_new_block(data_block_name)
               blocks_added["comp"] += 1
               for loop in data_block.get_mmcif_category_names():
                  category = data_block.find_mmcif_category(loop)
                  new_loop = new_block.init_loop(loop,[tag[len(loop):] for tag in category.tags])
                  for row in category:
                     new_loop.add_row(row)
            else:
               if not ALLOW_MISSING_DATA_BLOCKS:
                  raise Exception("Missing data block: "+data_block_name)

# --- MODIFICATIONS ---

      if USE_MODLIST:

         def is_integer(x):
            try:
               int(x)
               return True
            except ValueError:
               return False

         def get_modified_code(code):
            arr = code.split('mod')
            if len(arr) > 1:
               if(is_integer(arr[-1])):
                  arr[-1] = str(int(arr[-1]) + 1)
                  return "mod".join(arr)
            arr = code.split('-')
            if len(arr) > 1:
               if(is_integer(arr[-1])):
                  arr[-1] = str(int(arr[-1]) + 1)
                  return "-".join(arr)
            return code+"-1"

         chemmod_list_renamed = []
         for el in data:
            mods = el["mod"]
            for i in range(0,len(mods)):
               NEW_ENTRY = True
               IS_UNIQUE = True
               mod_id = mods[i]["id"]
               for entry in chemmod_list_renamed:
                  if entry[2] == mod_id:
                     NEW_ENTRY = False
                     if check_identical("mod_{}".format(entry[0]),"mod_{}".format(mod_id),entry[1],el["data"]):
                        print("Skipping identical modification: {}".format(entry[2]))
                        IS_UNIQUE = False
                        break
               if NEW_ENTRY:
                  chemmod_list_renamed.append([mod_id,el["data"],mod_id,mods[i]["row"]])
               elif IS_UNIQUE:
                  new_name = get_modified_code(mod_id)
                  IS_VALID = False
                  while not IS_VALID:
                     IS_VALID = True
                     for entry in chemmod_list_renamed:
                        if new_name == entry[2]:
                           new_name = get_modified_code(entry[2])
                           IS_VALID = False
                           break
                  print("Renaming duplicated modification name {} to: {}".format(mod_id,new_name))
                  chemmod_list_renamed.append([mod_id,el["data"],new_name,mods[i]["row"]])
                  mods[i]["new_id"] = new_name

         # Add row to data_mod_list
         for code_cif in chemmod_list_renamed:
            code_cif[3][0] = code_cif[2]
            modlist_loop.add_row(code_cif[3])

         # Add all data_mod_ blocks
         for code_cif in chemmod_list_renamed:
            data_block_name = "mod_{}".format(code_cif[0])
            data_block_name_new = "mod_{}".format(code_cif[2])
            data_block = code_cif[1].find_block(data_block_name)
            if(data_block):
               print("Adding modification to dictionary: {}".format(code_cif[2]))
               new_block = new_cif.add_new_block(data_block_name_new)
               blocks_added["mod"] += 1
               for loop in data_block.get_mmcif_category_names():
                  category = data_block.find_mmcif_category(loop)
                  new_loop = new_block.init_loop(loop,[tag[len(loop):] for tag in category.tags])
                  if code_cif[0] != code_cif[2]:   # Modify renamed mod_id in table
                     idx = -1
                     for i in range(0,len(category.tags)):
                        if str(category.tags[i].split('.')[-1]) == 'mod_id':
                           idx = i
                           break
                     if idx >= 0:
                        for row in category:
                           row[idx] = code_cif[2]
                  for row in category:
                     new_loop.add_row(row)
            else:
               raise Exception("Input CIF does not contain {} block".format(data_block_name))

# --- LINKS ---

      if USE_LINKLIST:

         if USE_MODLIST:
            # Update links with new modification IDs
            for el in data:
               for link in el["link"]:
                  for mod in el["mod"]:
                     for i in [2,5]:
                        if link["row"][i] == mod["id"]:
                           if "new_id" in mod:
                              link["row"][i] = mod["new_id"]

         chemlink_list = []
         for el in data:
            for link in el["link"]:
               NEW_LINK = True
               for existing_link in chemlink_list:
                  if link["id"] == existing_link["id"]:
                     link_block_id = "link_{}".format(link["id"])
                     if check_identical(link_block_id,link_block_id,el["data"],existing_link["data"]["data"]):
                        if list(link["row"]) == list(existing_link["link"]["row"]):
                           print("Skipping identical link: {} from {}".format(link["id"],el["name"]))
                           NEW_LINK = False
                           break
                        print(" ".join(list(link["row"])))
                        print(" ".join(list(existing_link["link"]["row"])))
                        raise Exception("Identical links encountered with inconsistent renamed modifications.")
                     raise Exception("Non-identical links encountered with the same ID: {}".format(link["id"]))
               if NEW_LINK:
                  chemlink_list.append({"id":link["id"],"data":el,"link":link})

#for link in chemlink_list:
#            print("LINK: {} {} {}".format(link["id"],link["data"]["name"],link["data"]["data"]))
#            print("  {}".format(link["link"]["row"]))

         # Add rows to the link_list block and add all data_link_ blocks
         for link in chemlink_list:
            linklist_loop.add_row(link["link"]["row"])
            data_block_name = "link_{}".format(link["id"])
            data_block = link["data"]["data"].find_block(data_block_name)
            if(data_block):
               print("Adding link to dictionary: {}".format(link["id"]))
               new_block = new_cif.add_new_block(data_block_name)
               blocks_added["link"] += 1
               for loop in data_block.get_mmcif_category_names():
                  category = data_block.find_mmcif_category(loop)
                  new_loop = new_block.init_loop(loop,[tag[len(loop):] for tag in category.tags])
                  for row in category:
                     new_loop.add_row(row)
            else:
               raise Exception("Input CIF does not contain {} block".format(data_block_name))

# --- UNKNOWN ---

      unknown_blocks = []
      for el in data:
         input_cif = el["data"]
         for block in input_cif:
            if not new_cif.find_block(block.name):
               print("Unknown block '"+block.name+"' found in "+el["name"])
               unknown_blocks.append(block)

      if opt['unknown']:
         for block in unknown_blocks:
            if new_cif.find_block(block.name):
               raise Exception("Cowardly refusing to continue - unknown block with the same name '"+block.name+"' found in multiple input CIF files")
            print("Adding unknown block to dictionary: "+block.name)
            new_block = new_cif.add_new_block(block.name)
            blocks_added["unknown"] += 1
            for loop in block.get_mmcif_category_names():
               category = block.find_mmcif_category(loop)
               new_loop = new_block.init_loop(loop,[tag[len(loop):] for tag in category.tags])
               for row in category:
                  new_loop.add_row(row)
      else:
         print("Unknown blocks will not be copied to the output CIF dictionary.")

# --- OUTPUT ---

      if DISPLAY_FINAL_METADATA:
         for el in data:
            print("DATA:")
            for obj in el:
               print("  {} : {}".format(obj,el[obj]))

      print("CIF files combined")
      print("List blocks included:")
      if blocks_added["comp_list"]:
         print("  "+"comp_list")
      if blocks_added["mod_list"]:
         print("  "+"mod_list")
      if blocks_added["link_list"]:
         print("  "+"link_list")
      print("Data blocks included:")
      if blocks_added["comp"] > 0:
         print("  Components:  \t\t"+str(blocks_added["comp"]))
      if blocks_added["mod"] > 0:
         print("  Modifications:\t"+str(blocks_added["mod"]))
      if blocks_added["link"] > 0:
         print("  Links:       \t\t"+str(blocks_added["link"]))
      if blocks_added["unknown"] > 0:
         print("Unknown blocks included:  "+str(blocks_added["unknown"]))


      new_cif.write_file(fileout)
      print("Output written to: %s" % fileout)

 except Exception as e:
   print("Error: %s" % e)
   print("Cannot continue - program terminated.")
   sys.exit(1)

 print("Completed.")

if __name__ == '__main__':
   parser = argparse.ArgumentParser()
   parser.add_argument('CIFs', nargs='+', help="One or more CIF dictionaries to be merged")
   parser.add_argument('-o','--output', dest='fileout', required=False, help="Output accumulated dictionary file", default='output.cif')
   parser.add_argument('-u','--unknown', dest='unknown', action='store_true', help="Include unknown blocks - will fail if unknown blocks with the same name are found in multiple input CIF files")
   args = parser.parse_args()
   
   opt = {}
   if args.unknown:
      opt['unknown'] = True
   
   accumulate(args.CIFs,args.fileout,opt)

