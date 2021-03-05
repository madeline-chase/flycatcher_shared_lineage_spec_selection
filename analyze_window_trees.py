import ete3
from ete3 import Tree
from sys import argv


# function to convert a list of nwk trees into a ETE tree objects
def multiTrees(lines):
  treeList = []
  for line in lines:
    if len(line) > 1:
      tree = Tree(line.rstrip())
      treeList.append(tree)
  return treeList


input = argv[1]
scaff = argv[2]
#function to test whether a specific set of individuals is monophyletic
def areMono(tree,leaves):
  ancestor = tree.get_common_ancestor(leaves)
  if len(ancestor.get_leaf_names()) == len(leaves):
    return True
  else:
    return False


# get trees

treeFile = open(input, "rU")

treeFileLines = treeFile.readlines()
treeList = multiTrees(treeFileLines)


# Define taxa

coll = "93F24", "93F26", "93F30", "93F32", "93F34", "93F35", "93F42", "93F44", "93F45", "93F47", "93F54", "93F56", "93F59", "93F74", "93F75", "93F77", "93F82", "93F88", "93F89", "93F90", "93F92", "93F93", "93F94", "93M25", "93M27", "93M28", "93M29", "93M36", "93M38", "93M39", "93M40", "93M41", "93M46", "93M53", "93M55", "93M58", "93M71", "93M72", "93M73", "93M78", "93M79", "93M80", "93M81", "93M83", "93M84", "93M86", "93M91", "15F129", "15F130", "15F131", "15F135", "15F142", "15F143", "15F145", "15F149", "15F151", "15F17", "15F18", "15F21", "15F22", "15F23", "15F24", "15F25", "15F29", "15F447", "15F448", "15F450", "15F453", "REF","15F457", "15F459", "15F460", "15M153", "15M155", "15M158", "15M160", "15M161", "15M162", "15M163", "15M201", "15M202", "15M203", "15M204", "15M207", "15M468", "15M469", "15M475", "15M477", "15M49", "15M537", "15M568", "15M571", "15M573", "15M589", "15M684", "15M724"
pied = "SP_11_M", "SP_12_M", "SP_13_M", "SP_14_M", "SP_15_M", "SP_16_F", "SP_SV1_M", "SP_SV2_M", "SP_SV4_M", "SP_SV5_M", "SP_SV7_F"
redb = "Sample_31", "Sample_32", "Sample_4", "Sample_48", "Sample_49", "Sample_51", "Sample_52", "Sample_53", "Sample_54", "Sample_55", "Sample_56", "Sample_77", "Sample_82", "Sample_83", "Sample_92"
taig = "Sample_1", "Sample_10", "Sample_11", "Sample_12", "Sample_13", "Sample_16", "Sample_17", "Sample_18", "Sample_2", "Sample_20", "Sample_21", "Sample_22", "Sample_23", "Sample_24", "Sample_25", "Sample_26", "Sample_28", "Sample_29", "Sample_3", "Sample_30", "Sample_33", "Sample_34", "Sample_35", "Sample_36", "Sample_37", "Sample_38", "Sample_39", "Sample_40", "Sample_41", "Sample_42", "Sample_44", "Sample_45", "Sample_5", "Sample_57", "Sample_58", "Sample_59", "Sample_6", "Sample_60", "Sample_61", "Sample_63", "Sample_64", "Sample_66", "Sample_67", "Sample_68", "Sample_69", "Sample_7", "Sample_70", "Sample_71", "Sample_72", "Sample_73", "Sample_74", "Sample_75", "Sample_76", "Sample_78", "Sample_79", "Sample_8", "Sample_80", "Sample_81", "Sample_86", "Sample_87", "Sample_88", "Sample_89", "Sample_9", "Sample_90", "Sample_91"



for x in range(len(treeList)):
    t = treeList[x]
    if not "F_hyperythra_X" in t:
        print(scaff, x+1, 'Outgroup_missing')
        out = 'outgroup_missing'
    elif len(t.get_leaf_names()) < 188:
        print(scaff, x+1, "Sample_missing")
        out = sample_missing
    else:
        t.set_outgroup(t&"F_hyperythra_X")
        if areMono(t, coll) and areMono(t, pied) and areMono(t, taig) and areMono(t, redb):
            if areMono(t, coll+pied) and areMono(t, taig+redb):
                print(scaff,x+1, "Species_tree_topology")
                out = "species_tree_topology"
            elif areMono(t, coll+pied):
                if areMono(t, coll+pied+taig):
                    out = "coll_pied_sister_plus_taig"
                    print(scaff, x+1, out)
                elif areMono(t, coll+pied+redb):
                    out = "coll_pied_sister_plus_redb"
                    print(scaff, x+1, out)
                else:
                    out = "Diff_topology_cp_group_all_mono"
                    print(scaff, x+1, out)
            else:
                out = "all_mono_diff_top"
                print(scaff, x+1, out)
        elif areMono(t, taig) and areMono(t, redb) and areMono(t, coll+pied):
            if areMono(t,redb+taig):
                out = "Coll_pied_not_mono_RT_sisters"
                print(scaff, x+1, out)
            elif areMono(t, taig+coll+pied) and areMono(t, redb):
                out = "taig_sister_to_cp_mixed"
                print(scaff, x+1, out)
            elif areMono(t, redb+coll+pied) and areMono(t, taig):
                out = "redb_sister_to_cp_mixed"
                print(scaff, x+1, out)
            else:
                out = 'Diff_topology'
                print(scaff, x+1, out)

        elif areMono(t, taig+redb) and areMono(t, coll) and areMono(t, pied):
            if areMono(t, coll+pied):
                out = "sister_groups_RT_not_recip_CP_recip"
                print(scaff, x+1, out)
            else:
                out = 'sister_groups_RT_not_recip_CP_diff'
                print(scaff, x+1, out)
        elif areMono(t, taig+redb) and areMono(t, coll+pied):
            out = 'sister_pairs_nested_groups'
            print(scaff, x+1, out)
        elif areMono(t, coll+pied) and areMono(t,taig) and not areMono(t,redb):
            out = 'cp_pair_par_not_mono'
            print(scaff, x+1, out)
        elif areMono(t, coll+pied) and not areMono(t, taig) and areMono(t, redb):
            out = 'cp_pair_taig_not_mono'
            print(scaff, x+1, out)
        elif areMono(t, coll+pied) and not areMono(t, taig) and not areMono(t, redb):
            out = 'cp_pair_taig_rb_both_not_mono'
            print(scaff, x+1, out)
        elif areMono(t, redb+taig) and not areMono(t, coll) and areMono(t, pied):
            out = 'rt_mono_group_pied_mono_coll_not_mono'
            print(scaff, x+1, out)
        elif not areMono(t, coll) and not areMono(t, pied) and not areMono(t, taig) and not areMono(t, redb):
            out = 'nobody_mono'
            print(scaff, x+1, out)
        else:
            out = "Diff_topology"
            print(scaff, x+1, out)
