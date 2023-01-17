#Treemmer_v0.3

#Copyright 2019 Fabrizio Menardo


#   This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Dependencies:
#  ete3 
#  joblib
#  numpy
#  matplotlib	

# Tremmer_v0.3 is compatible with joblib 0.13.1 and 0.13.2




# If you use Treemmer for your research, please cite:
# Treemmer: a tool to reduce large phylogenetic datasets with minimal loss of diversity. Menardo et. al (2018),BMC Bioinformatics 19:164. https://doi.org/10.1186/s12859-018-2164-8

from collections import defaultdict
from joblib import Parallel, delayed
from ete3 import Tree
import sys
import random
import operator
import argparse
import csv



############################################################			define arg type float 0 < X > 1		###############################################################

def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]"%(x,))
    return (x)


##########################################			FIND LEAVES NEIGHBORS OF A LEAF (2 NODE OF DISTANCE MAX) and calc DISTANCE			#######################	

def find_N(t,leaf):
	dlist ={}
	parent= leaf.up
	dist_parent=leaf.dist
	flag=0

	if arguments.verbose==3:
		print ("leaf findN at iteration:	" + str(counter))
		print (leaf)
		print ("parent findN at iteration:	" + str(counter))
		print (parent)
		print (parent.get_children())

	sister_flag=0
	for n in range(0,len(parent.get_children())):				##this for loop start from parent and climb up max two nodes, if it finds leaves calculate the distances, 
		if parent.is_root():
			flag=1
		#	break			#this  would stuck the algorithm into a infinite loop when the  root is  polytomy
		if arguments.verbose==3:							
			print ("children	" + str(n))
			print (parent.children[n])
					
		if (parent.children[n].is_leaf()):						# search at one node of distance
			if (parent.children[n] != leaf):
				DIS = leaf.get_distance(parent.children[n])
				dlist.update({leaf.name + "," +parent.children[n].name : DIS})
				flag=flag+1
				if arguments.verbose==3:
					print (leaf.name + "," +parent.children[n].name + str(DIS) + "have one node of distance")
		else:	
			if flag == 0:									
				if arguments.verbose==3:					#going up, search at two nodes of distance
					print ("going up, brother is node")
				
				temp_dlist={}
				for nn in range(0,len(parent.children[n].get_children())):
					if (parent.children[n].children[nn].is_leaf()):
						DIS = leaf.get_distance(parent.children[n].children[nn])
						temp_dlist.update({leaf.name + "," +parent.children[n].children[nn].name : DIS})
						sister_flag=sister_flag +1
					
				
	if ((sister_flag==1) and (flag==0)):						#collect results at two nodes of distance only if there are no leaves that are closer
		dlist.update(temp_dlist)
		if arguments.verbose==3:
			print (str(temp_dlist) + "	are not sister taxa, but neighbours first is leaf, second is upper neighbor")


	
	if (flag == 0):			#### this means that the leaf has no neighbors at one node of dist
		parent=parent.up 		#### therefore I climb the tree down towards the root of one more step and look for leaves
		multi_flag=0
		if arguments.verbose==3:
			print ("going down")		
			print ("gran parent")		
			print (parent)
		temp_dlist={}
		for n in range(0,len(parent.get_children())):		#this for loop start from gran parent and climb up max one nodes, if it finds leaves calculate the distances, 
			if parent.is_root():
				break										
			if (parent.children[n].is_leaf()):
				DIS = leaf.get_distance(parent.children[n])
				multi_flag = multi_flag+1
				temp_dlist.update({leaf.name + "," +parent.children[n].name : DIS})
		if multi_flag==1:					# this is to deal with polytomies 
			dlist.update(temp_dlist)
			if arguments.verbose==3:
				print (leaf.name + "," +parent.children[n].name + str(DIS) + "	are not sister taxa, but neighbours first is leaf, second is neighbor of downstair (towards root)")
	
	return (dlist)
		
##########################################				Check if leaf is protected			#######################

def check_protected(leaf_to_prune):
	

	if arguments.verbose > 1:
		print ("checking " + str(leaf_to_prune))
	warning=0	
	for tag in  dict_meta[leaf_to_prune]:      # loop thru all the tags of the selected leaf
		tag_counter=0

		for k,v in dict_meta.items():		#loop thru all the leaves,tag 
			flag=0
			for value in v:					
				if (str(value) == str(tag)):		#count the leaves with the tag 
					flag=1
			if flag==1:				# if the leaf has several time the same tag (mistake in input) it counts only one
				tag_counter = tag_counter +1
		if arguments.verbose > 1:
			print (str(leaf_to_prune)+ " " + str(tag) + " "+ str(tag_counter))

		if (arguments.list_meta_count):				# if -lmc option
			if (int(dict_meta_count[tag] == [])):		# if the tag is not present in list_meta_count
				dict_meta_count[tag] = 0					
			if tag_counter <= int(dict_meta_count[tag]):
				warning=1
		else:
			if (arguments.meta_count):				# if -mc option

				if (int(tag_counter) <= int(arguments.meta_count)):
					warning=1
		if arguments.verbose > 1:
			print ("warning = "+ str(warning))		
	return (warning)			

##########################################		IDENTIFY  LEAF TO PRUNE RANDOM			#######################


def find_leaf_to_prune_random(leaves):					
	random.shuffle(leaves)
	flag=0
	for leaf in leaves:
		warning=0
		if (arguments.list_meta and dict_meta[leaf.name]):			#check if leaf is protected
			warning=check_protected(leaf.name)

		if warning == 0:
			flag=1
			break

	if flag == 0:
		leaf_to_pr="stop,"
	else:
		leaf_to_pr=leaf.name
		
				
	return(leaf_to_pr)

##########################################		IDENTIFY  LEAF TO PRUNE			#######################


def find_leaf_to_prune(dlist):					#parse the list with all neighbor pairs and distances, find the closest pair and select the leaf
	warning=1
	while warning != 0:
		if (len(dlist) == 0):
			leaf_to_prune = "stop,"
			break

		min_key=min(dlist, key=dlist.get)
		min_val=dlist[min_key]


		d_min={}
		d_min.update({min_key:min_val})
	

	
		pair_unsplit= str(random.choice(list(d_min)))
		
		pair=pair_unsplit.split(",")
		leaf1 = t.search_nodes(name=pair[0])[0]
		leaf2 = t.search_nodes(name=pair[1])[0]
	
		if  (leaf1.dist > leaf2.dist):
			if (arguments.leaves_pair == 1):
				leaf_to_prune = leaf2.name
				leaf_to_keep = leaf1.name
			if (arguments.leaves_pair == 0):
				leaf_to_prune = leaf1.name
				leaf_to_keep = leaf2.name
	
		if  (leaf1.dist < leaf2.dist):
			if (arguments.leaves_pair == 1):		
				leaf_to_prune = leaf1.name
				leaf_to_keep = leaf2.name
			if (arguments.leaves_pair == 0):		
				leaf_to_prune = leaf2.name
				leaf_to_keep = leaf1.name
				
	
		if  ((leaf1.dist == leaf2.dist) or (arguments.leaves_pair ==2)):
			leaf_to_prune = random.choice(list(pair))			#this select the leaf at random within the couple
			for leaf in pair:
				if leaf != leaf_to_prune:
					leaf_to_keep = leaf

		if arguments.verbose > 1:		
			print ("leaf_to_check  " + str(leaf_to_prune))

		if (arguments.list_meta and dict_meta[leaf_to_prune]):			#check if leaf is protected
			warning=check_protected(leaf_to_prune)
	
			if warning == 1:

				if arguments.verbose > 1:
					print (" checking the neighbour")
					print ("leaf_to_check " + str(leaf_to_keep))
				warning=check_protected(leaf_to_keep)			# if leaf is protected I check the sister
				if warning == 0 :
					leaf_to_prune=leaf_to_keep

			
		else: warning = 0
		if warning ==1:
			del dlist[pair_unsplit]						#if both leaves of the pair are protected => delete the pair from dlist and make another cycle
			

	return (leaf_to_prune)





##########################################				PRUNE LEAF FROM TREE			#######################


def prune_t(leaf_to_prune,tree):

	G = tree.search_nodes(name=leaf_to_prune)[0]		
	parent= G.up
	dist_parent=G.dist

	if (len(parent.get_children()) == 2):
	

		if  parent.children[0] != G:
			parent.children[0].dist = parent.children[0].dist + parent.dist		

		if parent.children[1] != G:
			parent.children[1].dist = parent.children[1].dist + parent.dist	
			
	G.detach()

	if (len(parent.get_children()) == 1):
		parent.delete()		# after pruning the remaining branch will be like this ---/---leaf_name. I delete useless node keeping the b length


	return (tree)


####################################################################	calculate Tree length ##########################################################3

def calculate_TL(t):
	tree_length=0
	for n in t.traverse():
		tree_length=tree_length+n.dist
	tot_TL = tree_length
	return(tot_TL)


##########################################		PRUNE LEAF FROM MATRIX		#######################
def prune_dist_matrix(dlist,leaf_to_prune):
	key_del=[]
	for k, v in dlist.items():

		(one,two)=k.split(",")
		if ((one == leaf_to_prune) or (two == leaf_to_prune)):
			key_del.append(k)

	for KK in key_del:
		del dlist[KK]
	return (dlist)

##########################################		parallel loop		#######################
def parallel_loop(t,leaves,i):
	n=i
	DLIST_temp={}
	while n < len(leaves):
		N_list=find_N(t,leaves[n])
		n=n+arguments.cpu    			#n of  threads
		if N_list:
			DLIST_temp.update(N_list)
	return (DLIST_temp)


##########################################		write output with stop option		#######################

def write_stop(t,output1,output2):
	F=open(output1,"w")
	F.write(t.write())
	F.close()
	leaves = t.get_leaves()
	list_names=[]
	for leaf in leaves:
		list_names.append(leaf.name) 	
	F=open(output2,"w")
	F.write("\n".join(list_names))
	F.close()

##########################################		read list leaf_name,tag		#######################

def read_list_meta (path_to_list_meta):

	dict_meta = defaultdict(list)
	with open(path_to_list_meta, 'r') as f:
		reader = csv.reader(f)
		list_meta = list(reader)
	list_meta=filter(None, list_meta)
   		 
	for taxa_name, tag in list_meta:	
		dict_meta[taxa_name].append(tag)


	return (dict_meta)

##########################################		read list tag,number 	#######################

def read_list_tags (path_to_list_tag):
	dict_tag= defaultdict(list)
	with open(path_to_list_tag, 'r') as f:
    		reader = csv.reader(f)
    		list_tag = list(reader)
	list_tag=filter(None, list_tag)

	for tag, count in list_tag:	
		dict_tag[tag] = count


	return (dict_tag)

##########################################		select clade 	#######################


def select_clade (tip1,tip2):
	temp=[]
	ancestor = t.get_common_ancestor(tip1,tip2)
	write_stop(ancestor,arguments.INFILE+"_subtree",arguments.INFILE+"_subtree_leaves_list")

	sys.exit(0)

##########################################		select all 	#######################

def select_all (t):
	temp=[]
	all_leaves = t.get_leaves()

	for all_leaf in all_leaves:
		temp.append(all_leaf.name)

	F=open(arguments.INFILE+"_leaf_names_all","w")
	F.write("\n".join(temp))
	F.close()

	sys.exit(0)

#################################################  make plot  ##########################################################
def make_plot ():
	import numpy as np
	import matplotlib.pyplot as plt
	from matplotlib.ticker import MaxNLocator
	ax = plt.figure().gca()
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))	
	plt.scatter(x, y, s= 2, c= 'black')
	plt.xlim(ori_length,0)
	plt.ylim(-0.02,1.02)
	plt.xlabel('Number of leaves')
	plt.ylabel('Relative tree length')
	#plt.savefig(arguments.INFILE+'_res_'+ str(arguments.resolution)+'_TLD.png')
	plt.savefig(arguments.INFILE+'_res_'+ str(arguments.resolution)+'_TLD.pdf')	




############################################################		arguments and command line menu			###########################################33

parser = argparse.ArgumentParser(usage="Treemmer_v0.3.py INFILE [options (-h to see all options)]")

parser.add_argument('INFILE',type=str,help='path to the newick tree')
parser.add_argument('-X','--stop_at_X_leaves', metavar='X', default='0', help='Output reduced tree with X leaves. If multiple values are given Treemmer will produce multiple reduced datsets in the same run' , type =int, nargs='*')
parser.add_argument('-RTL','--stop_at_RTL', metavar='0-1', default='0', help='Output reduced tree with the specified RTL. If multiple values are given Treemmer will produce multiple reduced datsets in the same run', type =restricted_float, nargs='*')
parser.add_argument('-r','--resolution', metavar='INT', default=1,help='number of leaves to prune at each iteration (default: 1)',type =int, nargs='?')
parser.add_argument('-p','--solve_polytomies',help='resolve polytomies at random (default: FALSE)',action='store_true',default =False)
parser.add_argument('-pr','--prune_random',help='prune random leaves (default: FALSE)',action='store_true',default =False)
parser.add_argument('-lp','--leaves_pair', metavar='0,1,2', default=2,help='After the pair of leaves with the smallest distance is dentified Treemmer prunes: 0: the longest leaf\n1: the shortest leaf\n2: random choice (default: 2)',type =int, nargs='?')
parser.add_argument('-np','--no_plot',help='do not load matplotlib and plot (default: FALSE)',action='store_true',default =False)
parser.add_argument('-fp','--fine_plot',help='when --resolution > 1, plot RTL vs n leaves every time a leaf is pruned  (default: FALSE => plot every X leaves (X = -r))',action='store_true',default =False)
parser.add_argument('-c','--cpu', metavar='INT', default=1,help='number of cpu to use (default: 1)',type =int, nargs='?')
parser.add_argument('-lm','--list_meta', metavar='path/to/file', default="",help='path to file with metainformation. Format for each line: "leaf_name,tag". Leaves can appear mutiple times with different tags, or not appear at all',type =str, nargs='?')
parser.add_argument('-mc','--meta_count', metavar='INT', default= '0' ,help='if the -lm option is active -mc defines the minimum number of leaves that will be kept for each category defined in the metainformation file (default = 0)',type =int, nargs='?')
parser.add_argument('-lmc','--list_meta_count', metavar='path/to/file', default="",help='path to file. Format for each line: "tag,number", this option is alternative to -mc and allows to specify the different minimum number of leaves that shuld be retained for different categories',type =str, nargs='?')
parser.add_argument('-v' ,'--verbose', metavar='0,1,2', default='1', help='0: silent (almost), 1: show progress, 2: print tree at each iteration, 3: only for testing (findN), 4: only for testing (prune_t) (default: 1)', type =int, nargs='?',choices=[0,1,2,3,4])
parser.add_argument('-sc1' ,'--select_clade_1', metavar='leaf_name', default='', help='use together with -sc2. Treemmer will identify the smallest monophyletic clade including two specified leaves and output a list of leaves belonging to this clade. This can be usefull to prepare the --list_meta input file in case you want to prune only leaves belonging (or not belonging) to a certain clade', type =str, nargs='?')
parser.add_argument('-sc2' ,'--select_clade_2', metavar='leaf_name', default='', help='use together with -sc1. Treemmer will identify the smallest monophyletic clade including two specified leaves and output a list of leaves belonging to this clade. This can be useful to prepare the --list_meta input file in case you want to prune only leaves belonging (or not belonging) to a certain clade', type =str, nargs='?')
parser.add_argument('-sa' ,'--select_all', default= False, help='output the list of leaf names in the input tree and exit', action='store_true')
parser.add_argument('-pa' ,'--plot_always', default= False, help='output the RTL plot with the smallest tree defined by the -X or -RTL option', action='store_true')
parser.add_argument('-pc' ,'--plot_complete', default= False, help='plot the complete RTL plot and file when the -X or -RTL options are specified ', action='store_true')
parser.add_argument('-sX','--switch_at_X', metavar='sX', default=1, help='Treemmer will start normally and switch to random subsampling when the tree has less than sX leaves. This option can be used with -sRTL, Treemmer will change behaviour as soon as one of the two criteria is met' , type =int, nargs='?')
parser.add_argument('-sRTL','--switch_at_RTL', metavar='0-1', default=0, help='Treemmer will start normally and switch to random subsampling when the tree is shorter than sRTL.  This option can be used with -sX, Treemmer will change behaviour as soon as one of the two criteria is met', type =restricted_float, nargs='?')
arguments = parser.parse_args()

if ((not (arguments.stop_at_RTL)) and (not(arguments.stop_at_X_leaves))):
	arguments.plot_complete = True

#	raise argparse.ArgumentTypeError("-X and -RTL are mutually exclusive options")

######   SOFTWARE STARTS

t = Tree(arguments.INFILE,format=1)

if arguments.solve_polytomies:
	t.resolve_polytomy()

if ((arguments.select_clade_1) and (arguments.select_clade_2)):						#select clade routine
	select_clade(arguments.select_clade_1,arguments.select_clade_2)

if (arguments.select_all):										#select all routine
	select_all(t)

if (arguments.list_meta):
	dict_meta = read_list_meta(arguments.list_meta)

if (arguments.list_meta_count):
	dict_meta_count = read_list_tags(arguments.list_meta_count)


if arguments.verbose > 0:
												# print progress on standard output
	print ("N of taxa in tree is : "+ str(len(t)))
	
	if arguments.solve_polytomies:
		print ("\nPolytomies will be solved at random")
	else:
		print ("\nPolytomies will be kept")
	
	if arguments.prune_random:
		print ("\nA random leaf is pruned at each iteration")

	if (arguments.leaves_pair == 0):
		print ("\nAfter the pair of leaves with the smallest distance is dentified Treemmer will prune the longest of the two leaves")

	if (arguments.leaves_pair == 1):
		print ("\nAfter the pair of leaves with the smallest distance is dentified Treemmer will prune the shortest of the two leaves")

	if (arguments.leaves_pair == 2):
		print ("\nAfter the pair of leaves with the smallest distance is dentified Treemmer will prune one of the two leaves picked at random")
	
	if (arguments.switch_at_X !=1):
		print ("\nWhen the tree is reduced to " + str(arguments.switch_at_X) + " leaves, Treemmer will switch behaviour and start pruning random leaves")

	if arguments.switch_at_RTL:
		print ("\nWhen the tree is reduced to " + str(arguments.switch_at_RTL) + " of the original tree length, Treemmer will switch behaviour and start pruning random leaves")

	if arguments.stop_at_X_leaves:
		print ("\nTreemmer will reduce the tree to " + str(arguments.stop_at_X_leaves) + " leaves")
	
	if arguments.stop_at_RTL:
		print ("\nTreemmer will reduce the tree to " + str(arguments.stop_at_RTL) + " of the original tree length")
	else: 
		print ("\nTreemmer will calculate the tree length decay")

	if arguments.list_meta:
		print ("\nsome leaves are protected by the -lm options and will not be pruned based on what specified with -mc or -lmc")  

	
	
	print ("\nTreemmer will prune " + str(arguments.resolution) + " leaves at each iteration")
	print ("\nTreemmer will use " + str(arguments.cpu) + " cpu(s)")

x=[]
y=[]	
counter =0
output=[]
stop=0
TOT_TL=calculate_TL(t)
ori_length = len(t)
output.append ('1	' + str(len(t))) 			#append first point to the output with RTL = 1 (before starting pruning)################################
x.append(ori_length)
y.append(1)
leaves = t.get_leaves()

sys.setrecursionlimit(50000)
leaf_names=[]
leaf_to_p=""

if (arguments.list_meta):				############# update the dictionary of taxa_names with tags, only taxa in the tree stays are kept in the dict
	for leaf in leaves:					    	
		leaf_names.append(leaf.name) 			    
	dict_meta_new={key: dict_meta[key] for key in leaf_names}   
	dict_meta=dict_meta_new

while (len(t) > 3):								#################### Main loop ################################
	counter = counter +1
	leaves = t.get_leaves()
	DLIST={}

	if arguments.verbose > 0:
		print ("\niter		" + str(counter))
		if arguments.verbose > 1:
			print ("\ncalculating distances\n")

	
	DLIST = Parallel(n_jobs=arguments.cpu)(delayed(parallel_loop)(t,leaves,i) for i in range(0,arguments.cpu))	#loop all leaves and find neighbours, report pairs and distances
	result = {}

	for d in DLIST:								#when running in parallel DLIST is a dict of dicts, this for loop merge them all in one
		result.update(d)

	DLIST=result

	if arguments.verbose > 1:
		print (DLIST)
		print ("\npruning\n")

	for r in range (1,arguments.resolution+1):    				#resolution loop (find leaf to prune, prune it, update matrix r times)


		if arguments.list_meta:						
			if ((len(DLIST)==0) and (len(t)>4)):
				leaf_to_p="stop,"

		if leaf_to_p != "stop,":
			if ((len(DLIST)<1) or (len(t) < 4)):
				break

		if arguments.prune_random:

			(leaf_to_p)= find_leaf_to_prune_random(leaves)		#find leaf to prune,  protections (from -lm option) are embedded in the function
		else:
			(leaf_to_p) = find_leaf_to_prune(DLIST)			#find leaf to prune,  protections (from -lm option) are embedded in the function


		if (leaf_to_p == "stop,"):
			if r == 1:				# if r > 1 some taxa might have been no considered because surrounded by already pruned leaves
				print ("WARNING: all remaining leaves are protected by the -lm option, outputting the results at current iteration")
			if r > 1:
				leaf_to_p = ""
				break				#if r > 1 make another cycle recalculating distances and maybe find some more leaves to prune 

		if (leaf_to_p != "stop,"):
	
			leaf_to_prune = t.search_nodes(name=leaf_to_p)[0]		
			t = prune_t(leaf_to_p,t)				#do the tree pruning
			leaves = t.get_leaves()

			if (arguments.list_meta):		####### update the dictionary of taxa_names with tags, only taxa in the tree are kept in the dict

				leaf_names=[]
				for leaf in leaves:					    	
					leaf_names.append(leaf.name) 			    
				dict_meta_new={key: dict_meta[key] for key in leaf_names}   
				dict_meta=dict_meta_new

			TL= calculate_TL(t)
			if not arguments.prune_random:	
				DLIST=prune_dist_matrix(DLIST,leaf_to_p)		#### purge the distance list of all pairs that have the pruned leaf 
			rel_TL=TL/TOT_TL
		
		#################################  		OUTPUT 		##########################################################		

		if ((arguments.fine_plot) and (leaf_to_p != "stop,")):							# plot point in rtld after every leaf independently of -r 
			output.append (str(rel_TL) + '	' + str(len(t)))
			length=len(t)	
			x.append(length)
			y.append(rel_TL)

		if (arguments.switch_at_X >= len(t)):
			arguments.prune_random = True
			
		if (arguments.switch_at_RTL >= rel_TL):
			arguments.prune_random = True
			
		if arguments.stop_at_X_leaves:										# if stop criterium is met (X) ==> output
			if ((max(arguments.stop_at_X_leaves) >= len(t)) or (leaf_to_p == "stop,")):
				output1=arguments.INFILE+"_trimmed_tree_X_" + str(max(arguments.stop_at_X_leaves))
				output2=arguments.INFILE+"_trimmed_list_X_" + str(max(arguments.stop_at_X_leaves))
				write_stop(t,output1,output2)
				arguments.stop_at_X_leaves.remove(max(arguments.stop_at_X_leaves))
				#stop=1
				#break
			
		if arguments.stop_at_RTL:										# if stop criterium is met (RTL) ==> output
			if ((max(arguments.stop_at_RTL) >= rel_TL ) or (leaf_to_p == "stop,")):
				output1=arguments.INFILE+"_trimmed_tree_RTL_" + str(max(arguments.stop_at_RTL))
				output2=arguments.INFILE+"_trimmed_list_RTL_" + str(max(arguments.stop_at_RTL)) 
				write_stop(t,output1,output2)
				arguments.stop_at_RTL.remove(max(arguments.stop_at_RTL))
				#stop=1
				#break

		if (leaf_to_p == "stop,"):
			break

		if arguments.verbose > 1:										# print progress to standard output

			print ("\n ITERATION RESOLUTION:	" + str(r))
			print ("leaf to prune:\n" + str(leaf_to_p) + "	" + str(leaf_to_prune.dist))
			print ("\n new tree")
			print (t)
			print ("\nRTL :	" + str(rel_TL) + " N_seq:	" +str(len(t)))
			print ("\nnew matrix\n")
			print (DLIST)


	if  ((not(arguments.stop_at_RTL)) or  (len(arguments.stop_at_RTL) ==0)):
		if ( (not(arguments.stop_at_X_leaves)) or (len(arguments.stop_at_X_leaves) ==0)):
			if not arguments.plot_complete:
				stop = 1

	if (stop ==1):
		if arguments.verbose > 0:
			print ("\nRTL :	" + str(rel_TL) + " N_seq:	" +str(len(t)))
		break	

	if (leaf_to_p == "stop,"):
		if arguments.verbose > 0:
			print ("\nRTL :	" + str(rel_TL) + " N_seq:	" +str(len(t)))
		break

	if not (arguments.fine_plot):											# normal plot (with -fp = FALSE)
		output.append (str(rel_TL) + '	' + str(len(t)))
		length=len(t)	
		x.append(length)
		y.append(rel_TL)
		
	if arguments.verbose > 0:
		print ("\nRTL :	" + str(rel_TL) + " N_seq:	" +str(len(t)))
		 	

if ((stop == 0) or (arguments.plot_always)):														# create file for plot of rltd
	F=open(arguments.INFILE+"_res_"+ str(arguments.resolution) + "_LD","w")
	F.write("\n".join(output))

	if not arguments.no_plot:
		make_plot()

if (arguments.verbose > 0):
	print ("\n If you use Treemmer, please cite:\n\n \"Treemmer: a tool to reduce large phylogenetic datasets with minimal loss of diversity\" Menardo et. al., BMC Bioinformatics (2018) 19:164\n\n")

