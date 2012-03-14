#!/usr/bin/python

##############################################################################################
import os
os.environ['HOME'] = "./" #required for write rights for matplotlib in www-user environment	 #	
##############################################################################################

import matplotlib as plt
plt.use("Agg")
import warnings
import cPickle
import sys
import argparse
import pylab
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import collections
import numpy
import h5py
import rpy2.rlike.container as rlc



#Global vars :
#-------------

path_to_use = './'

#-------------


#################################################### Functions ############################################

#We do not want this depreciated warning
def fxn():
	warnings.warn("deprecated", DeprecationWarning)
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	fxn()
	pass
#
def import_sample_data(_file,sep='\t',verbose=False, line_to_skip=0):
	"""docstring for import_sample_data
	import data from tab separated file containing sample description. 
	"""
	f=open(_file,'r')
	t=f.readlines()
	f.close()
	data=[]
	for sample in t[line_to_skip:]:
		data.append( [i.strip() for i in sample.strip().split(sep)] )
		if verbose:
			print '#%s#'% data[-1] 
	return data
	pass
#


def beeswarm_adv(data, names,out='beeswarm.pdf', celltypes_to_use=None, title='',plot_median=False, use_mean=False, use_log=False):
	"""docstring for beeswarm
	plots a beeswram plot.
	data is an array or values
	names is the corresponding names for each data values.
	
	"""
	import pylab
	if use_mean:
		meam_method = numpy.mean
	else:
		meam_method = numpy.median
		
	in_gene = title
	graphics=importr('graphics')
	beeswarm= importr('beeswarm')
	grdevices = importr('grDevices')
	classes2use_order=[]
	classes2use={}
	for index,classe in enumerate(names):
		if classe not in classes2use:
			classes2use_order.append(classe)
			classes2use[classe]=[]
			classes2use[classe].append(index)
		else:
			classes2use[classe].append(index)


	#reduce the set
	reduced_set={}
	reduced_order=[]
	reduced_names=[]
	reduced_data=numpy.array([])
	reduced_colors=[]
	try:
		for i in celltypes:
			reduced_set[i]=classes2use[i]
	except:
		print head
		print 'You are calling for cell types %s not in the dataset.'%(i)
		sys.exit(42)
	classes2use = reduced_set 

	for index,i in enumerate(classes2use_order):
		if reduced_set.has_key(i):
			reduced_order.append(i)
			reduced_colors.append(colors[index])
	classes2use_order = reduced_order
	
	for i,j in enumerate(names):
		if classes2use.has_key(j):
			reduced_names.append(j)
			reduced_data = numpy.append(reduced_data, data[i])
	
	#model:
	nn=robjects.StrVector(numpy.array(names))
	fcdf=robjects.FloatVector(numpy.array(data))
	cool=robjects.DataFrame({'Expression': fcdf , 'Class':nn})
	Classes= robjects.r.factor(cool.rx2(2), levels = robjects.r.unique(cool.rx2(2)), ordered = 1)
	robjects.globalenv['Classes']=Classes
	
	a=beeswarm.beeswarm( robjects.r("Expression ~ Classes") , data = cool, method = "center",pch=1, bg=robjects.r("rainbow(6)") ,col = robjects.r("rainbow(16)") , las =2 )
	opt={'main':in_gene}
	graphics.title( **opt)
	grdevices.dev_off()
	
	x=numpy.array(a[0])
	y=numpy.array(a[1])
	c=numpy.array(a[3])
	
	fig=pylab.figure()
	
	if use_log:
		pylab.yscale('log', basey=2)
		pylab.ylabel('Expression (log2)')
		
	if plot_median:
		for	 i in range(len(classes2use_order)):
			z= numpy.array(classes2use[classes2use_order[i]]) 
			m=meam_method(numpy.exp2(y[z]))
			pylab.plot([i+1-.3, i+1+.3],[m,m],'b',linewidth=1.5, alpha=.7)
#		pylab.plot([1,2],[0,0],'b',linewidth=1.5,alpha=.5)
	pylab.xticks(range(len(classes2use_order)+1), numpy.concatenate([[''],classes2use_order],axis=0), rotation=90,size=9)
	fig.autofmt_xdate(bottom=0.18)	
		
	
	nn=robjects.StrVector(numpy.array(reduced_names))
	fcdf=robjects.FloatVector(numpy.array(reduced_data))
	cool=robjects.DataFrame({'Expression': fcdf , 'Class':nn})
	Classes= robjects.r.factor(cool.rx2(2), levels = robjects.r.unique(cool.rx2(2)), ordered = 1)
	robjects.globalenv['Classes']=Classes
	
	a=beeswarm.beeswarm( robjects.r("Expression ~ Classes") , data = cool, method = "center",pch=1, bg=robjects.r("rainbow(6)") ,col = robjects.StrVector(reduced_colors) , las =2 )
	x=numpy.array(a[0])
	y=numpy.array(a[1])
	c=numpy.array(a[3]) #r-colors
	
	for i,j in enumerate(x): 
		pylab.plot(x[i],numpy.exp2(y[i]),'o',color=c[i][0:-2],alpha=.7) # R adds FF at the end of the color string, which is bad.#But I do not use R-colors,
		
	pylab.title(in_gene)

	if type(out)==list:
		for i in out:
			pylab.savefig(i)
	else:
		pylab.savefig(out)

	pylab.close()
	return reduced_set
	
	
	
def beeswarm(data, names,out='beeswarm.pdf',title='', plot_median=False, use_mean=False, use_log=False):
	"""docstring for beeswarm
	plots a beeswram plot.
	data is an array or values
	names is the corresponding names for each data values.
	
	"""
	import pylab
	if use_mean:
		meam_method = numpy.mean
	else:
		meam_method = numpy.median
		
	in_gene = title
	graphics=importr('graphics')
	beeswarm= importr('beeswarm')
	grdevices = importr('grDevices')
	classes2use_order=[]
	classes2use={}
	for index,classe in enumerate(names):
		if classe not in classes2use:
			classes2use_order.append(classe)
			classes2use[classe]=[]
			classes2use[classe].append(index)
		else:
			classes2use[classe].append(index)
	
	nn=robjects.StrVector(numpy.array(names))
	fcdf=robjects.FloatVector(numpy.array(data))
	cool=robjects.DataFrame({'Expression': fcdf , 'Class':nn})
	Classes= robjects.r.factor(cool.rx2(2), levels = robjects.r.unique(cool.rx2(2)), ordered = 1)
	robjects.globalenv['Classes']=Classes
	
	a=beeswarm.beeswarm( robjects.r("Expression ~ Classes") , data = cool, method = "center",pch=1, bg=robjects.r("rainbow(6)") ,col = robjects.StrVector(colors) , las =2 )
	opt={'main':in_gene}
	graphics.title( **opt)
	grdevices.dev_off()
	
	x=numpy.array(a[0])
	y=numpy.array(a[1])
	c=numpy.array(a[3]) #r-colors
	
	fig=pylab.figure()
	
	if use_log:
		pylab.yscale('log', basey=2)
		pylab.ylabel('Expression (log2)')
	
	for i,j in enumerate(x): 
		pylab.plot(x[i],numpy.exp2(y[i]),'o',color=c[i][0:-2],alpha=.7) # R adds FF at the end of the color string, which is bad.
		
	
	if plot_median:
		for	 i in range(len(classes2use_order)):
			z= numpy.array(classes2use[classes2use_order[i]]) 
			m=meam_method(numpy.exp2(y[z]))
			pylab.plot([i+1-.3, i+1+.3],[m,m],'b',linewidth=1.5, alpha=.7)
#		pylab.plot([1,2],[0,0],'b',linewidth=1.5,alpha=.5)
	pylab.xticks(range(len(classes2use_order)+1), numpy.concatenate([[''],classes2use_order],axis=0), rotation=90,size=9)
	fig.autofmt_xdate(bottom=0.18)	
		
	pylab.title(in_gene)
	
	if type(out)==list:
		for i in out:
			pylab.savefig(i)
	else:
		pylab.savefig(out)
		
	pylab.close()
	return classes2use



def pickle_load(in_file):
	"""docstring for pickle_load
	load a pickled variable into memory.
	input is the file name and path.
	"""
	x=None
	if os.path.exists(in_file):
		pkl_file = open(in_file, 'rb')
		x=cPickle.load(pkl_file)
		pkl_file.close()
	else:
		print 'File %s could not be found ! ' % (in_file)
	
	return x
	pass

def find_probes(gene,organism):
	"""docstring for find_probes
		loads the annotation DBI R package converted by the generate_annotation_dics.py script.
		loops through all possible annotations, and the fuction hopefully returns the correct probes ;)
	"""
	gene=gene.lower()
	
	if 0: #slow method.
		p=None
		if organism == 'human':
			annot_to_probe = pickle_load(path_to_use + '/beeswarn_data/h_annot_to_probe.pkl')
		else:
			annot_to_probe = pickle_load(path_to_use + '/beeswarn_data/m_annot_to_probe.pkl')
			
			
		for k in annot_to_probe:
			p = annot_to_probe[k].get(gene)
			if p != None:
				if not WEB_opt:
					print 'found in '+k
				break
		
	else: #fast method
		if organism == 'human':
			r=h5py.File(path_to_use + '/beeswarn_data/h_annot_to_probe.hdf5','r')
		else:
			r=h5py.File(path_to_use + '/beeswarn_data/m_annot_to_probe.hdf5','r')
		
		root= r['root']
		p=None
		for k in root:
			if gene in root[k]:
				p = root[k].get(gene)[...]
			if p != None:
				if not WEB_opt:
					print 'found in '+k
				break
		r.close()
	return p
	pass
	
	
def make_raw_outputFile(data, probe_names, sample_names, adv=False, use_log=False, celltypes_to_use=None, filename='rawData.txt'):
	'''docstring for make_raw_outputFile
		prints the data to a file. If the advancedmode is on, then only samples/celltypes given, is printed, by reducing the input data.
		input: data, probe_names and sample names must be provided with same reduced index syntax, so that only probes wanted in file is given, an in constitent index order.
	'''
	
	outfile=open(filename, 'w')
	
	if adv:
		temp_index=0
		reduced_samplenames=sample_names
		for i in sample_names:
			if not i in celltypes_to_use:
				data = numpy.delete(data, temp_index, 1)
				reduced_samplenames = numpy.delete(reduced_samplenames, temp_index, 0)
			else:
				temp_index += 1
		
		sample_names=reduced_samplenames

	
	outfile.write(',')
	outfile.write(', '.join(sample_names))
	outfile.write('\n')
	for index,i in enumerate(probe_names):
		outfile.write(i)
		outfile.write(', ')
		for y in data[index]:
			outfile.write('%s,'%(y))
		outfile.write('\n')


#################################################### End of functions ########################################
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Plots gene expression accross the hematopoietic system.')
	parser.add_argument('-i','--input', help='Gene 1 to be plotted.')
	parser.add_argument('-j','--input2', help='Gene 2 to be plotted.')#, default='NULL')
	parser.add_argument('-f','--foldchange', action='store_const', const=True, default=False, help='Plot foldchange of gene against normal counterpart.')
	parser.add_argument('data_seq',help='dummy argument - remove when common2 has been updated', default=argparse.SUPPRESS )
	parser.add_argument('-o','--organism', type=str,choices= ['human','mouse'] , default='human', help='use human or mouse dataset.'  )
	parser.add_argument('-l','--data', help='Data kind to be plotted.', choices= ['leukemia','normal','both'], default='both')
	parser.add_argument('-lg1','--log2gene1', action='store_const', const=True, default=False, help='Give output of gene 1 in log2 transformed axis')
	parser.add_argument('-lg2','--log2gene2', action='store_const', const=True, default=False, help='Give output of gene 2 in log2 transformed axis')

#advanced options

	parser.add_argument('-adv','--advanced', action='store_const', const=True, default=False, help='Used advanced option to clear all default celltypes in output. Cells must be selected manually') 

#cell types
	parser.add_argument("-HSC_BM","--HSC_BM", dest='celltypes', action="append_const", const="HSC_BM", default=argparse.SUPPRESS, help="Hematopoietic") 
	parser.add_argument("-Early","--Early", dest='celltypes', action="append_const", const="early HPC_BM", default=argparse.SUPPRESS, help="-Early") 
	parser.add_argument("-CMP","--CMP", dest='celltypes', action="append_const", const="CMP", default=argparse.SUPPRESS, help="Common") 
	parser.add_argument("-GMP","--GMP", dest='celltypes', action="append_const", const="GMP", default=argparse.SUPPRESS, help="Granulocyte") 
	parser.add_argument("-MEP","--MEP", dest='celltypes', action="append_const", const="MEP", default=argparse.SUPPRESS, help="Megakaryocyte-erythroid") 
	parser.add_argument("-PM_BM","--PM_BM", dest='celltypes', action="append_const", const="PM_BM", default=argparse.SUPPRESS, help="Promyelocyte") 
	parser.add_argument("-MY_BM","--MY_BM", dest='celltypes', action="append_const", const="MY_BM", default=argparse.SUPPRESS, help="Myelocyte") 
	parser.add_argument("-PMN_BM","--PMN_BM", dest='celltypes', action="append_const", const="PMN_BM", default=argparse.SUPPRESS, help="Polymorphonuclear") 
	parser.add_argument("-PMN_PB","--PMN_PB", dest='celltypes', action="append_const", const="PMN_PB", default=argparse.SUPPRESS, help="Polymorphonuclear") 
	parser.add_argument("-cd14+monocytes","--cd14+monocytes", dest='celltypes', action="append_const", const="cd14+ monocytes", default=argparse.SUPPRESS, help="cd14+") 
	parser.add_argument("-AMLI_ETO","--AMLI_ETO", dest='celltypes', action="append_const", const="AMLI_ETO", default=argparse.SUPPRESS, help="AML") 
	parser.add_argument("-APL","--APL", dest='celltypes', action="append_const", const="APL", default=argparse.SUPPRESS, help="AML") 
	parser.add_argument("-inv16","--inv16", dest='celltypes', action="append_const", const="AML with inv(16)/t(16;16)", default=argparse.SUPPRESS, help="AML") 
	parser.add_argument("-t11q23","--t11q23", dest='celltypes', action="append_const", const="AML with t(11q23)/MLL", default=argparse.SUPPRESS, help="AML") 
	parser.add_argument("-LT_HSC","--LT_HSC", dest='celltypes', action="append_const", const="LT-HSC", default=argparse.SUPPRESS, help="Long") 
	parser.add_argument("-ST_HSC","--ST_HSC", dest='celltypes', action="append_const", const="ST-HSC", default=argparse.SUPPRESS, help="Short") 
	parser.add_argument("-LMPP","--LMPP", dest='celltypes', action="append_const", const="LMPP", default=argparse.SUPPRESS, help="Lymphoid-primed") 
	parser.add_argument("-CLP","--CLP", dest='celltypes', action="append_const", const="CLP", default=argparse.SUPPRESS, help="Common") 
	parser.add_argument("-ETP","--ETP", dest='celltypes', action="append_const", const="ETP", default=argparse.SUPPRESS, help="Early") 
	parser.add_argument("-ProB","--ProB", dest='celltypes', action="append_const", const="ProB", default=argparse.SUPPRESS, help="Pro-B") 
	parser.add_argument("-PreB","--PreB", dest='celltypes', action="append_const", const="PreB", default=argparse.SUPPRESS, help="Pre-B") 
	parser.add_argument("-IgM+SP","--IgM+SP", dest='celltypes', action="append_const", const="IgM+SP", default=argparse.SUPPRESS, help="Immunoglobulin") 
	parser.add_argument("-CD4","--CD4", dest='celltypes', action="append_const", const="CD4", default=argparse.SUPPRESS, help="CD4") 
	parser.add_argument("-NKmature","--NKmature", dest='celltypes', action="append_const", const="NKmature", default=argparse.SUPPRESS, help="Mature") 
	parser.add_argument("-MkE","--MkE", dest='celltypes', action="append_const", const="MkE", default=argparse.SUPPRESS, help="Megakaryocyte") 
	parser.add_argument("-MkP","--MkP", dest='celltypes', action="append_const", const="MkP", default=argparse.SUPPRESS, help="Megakaryocyte") 
	parser.add_argument("-PreCFUE","--PreCFUE", dest='celltypes', action="append_const", const="PreCFUE", default=argparse.SUPPRESS, help="Pre-colony-forming") 
	parser.add_argument("-CFUE","--CFUE", dest='celltypes', action="append_const", const="CFUE", default=argparse.SUPPRESS, help="Colony-forming") 
	parser.add_argument("-ProE","--ProE", dest='celltypes', action="append_const", const="ProE", default=argparse.SUPPRESS, help="Erythroid")
	

	args = parser.parse_args()

	in_gene = args.input
#hack for second gene.
	in_gene2 = args.input2
	if in_gene2 == 'NULL':
		in_gene2 = None
	
	log2gene1 = args.log2gene1
	log2gene2 = args.log2gene2

	fold_change = args.foldchange
	organism = args.organism
	data_to_use = args.data
	
	adv_mode = args.advanced
	if adv_mode:
		try:
			celltypes = args.celltypes
		except:
			print "When using advanced mode each celltype must be given as argument"
			sys.exit(42)
	else:
		celltypes=[]
	
	
						# CONFIG FILES #
##################################################################################	
##################################################################################	
	#read in configuration file:
#	config = import_sample_data(path_to_use+'JARID2.ini')
	config = import_sample_data('/webdata/servers.binf.ku.dk/programs/shs/JARID2_w.ini') # WEB server.
	
	path_to_use = config[0][1]	#first line is the path
	WEB_opt = int(config[1][1]) #second line is the web option.
##################################################################################	
	colors=['#FAAFBEFF','#808000FF','#FFF8C6FF','#8E35EFFF','#A52A2AFF', '#FFFFFFFF','#FF0000FF', '#008000FF', '#C0C0C0FF', '#0000FFFF', '#FFFF00FF', '#FF00FFFF', '#8AFB17FF','#FFA500FF', '#000000FF', '#00FFFFFF']
	py_colors=['#FAAFBE','#808000','#FFF8C6','#8E35EF','#A52A2A', '#FFFFFF','#FF0000', '#008000', '#C0C0C0', '#0000FF', '#FFFF00', '#FF00FF', '#8AFB17','#FFA500', '#000000', '#00FFFF']

######################################################
	
	# load appropriate data file: 
	
	if organism == 'mouse':
		all_data = pickle_load( path_to_use + '/beeswarn_data/nl_mouse_data.pkl')
	else: # for human there are three choices : ['leukemia','normal','both'], if the fc flag is activated, then only leukemia are plotted.
		if not fold_change:
			if data_to_use == 'leukemia' :
				all_data = pickle_load(path_to_use + '/beeswarn_data/all_data_expr.pkl')
			elif data_to_use == 'normal' :
				all_data = pickle_load(path_to_use + '/beeswarn_data/nl_human_data.pkl')
			else:
				merged_data = pickle_load( path_to_use + '/beeswarn_data/nl_human_data.pkl')	   # NOTA : nl_human_data will be updated as soon as the Cancer vs normal paper is out.
																								   
				all_data = pickle_load(path_to_use + '/beeswarn_data/all_data_expr.pkl')
				
				all_data['colnames'] = numpy.concatenate([all_data['colnames'],merged_data['colnames']],axis=0)
				all_data['data'] = numpy.concatenate([all_data['data'],merged_data['data']],axis=1)
				
		else:
			all_data = pickle_load(path_to_use + '/beeswarn_data/all_data_fc.pkl')

###########################################################
	# HTML header for output
	head='''<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
		<html xmlns="http://www.w3.org/1999/xhtml" lang="en">
		<head>
		<title> Servers.binf.ku.dk | SHS </title>
		<meta name="language" content="en" />
		<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1" />
		<style type="text/css">
		body {
		font-family: 'ArialMT', 'Arial', 'sans-serif';
		font-size: 14px;
		color: #3E3535;
		background-color: #FFFFFF;
		line-height: 20px;
		max-width: 630px;
		margin: 20px auto 20px auto;
		padding: 0px;
		}

		h1.title {font-size: 36px;}
		h2.title {font-size: 36px;}
		pre		 {font-size: 12px;}
		boring	 {font-size: 12px}
		</style>	
		</head>
		<body>
		<br><p style="font-size:150%;font-weight:bold;text-align:center">Welcome to the <span style="color:silver;font-size:250%">HemaExplorer</span></p><p style="text-align:center"><i><FONT COLOR="#000000"> --- See the expression of your favourite gene in the hematopoietic system --- </FONT></i></p>
		<!-- ##### ##### HEADER END ##### ##### -->'''
#########################################################
	# get list of probes for the different genes.
	
	probes_names = find_probes(in_gene, organism)
	
	if probes_names == None:
		print head
		print '%s could not be found on array. Maybe you are using an ambiguous gene name. If the gene is in database you might find an alias to your gene name at <a href="http://www.genecards.org/">Genecards</a>. <br><br> To get more details about genes in our database see our <a href="http://servers.binf.ku.dk/shs/help.php">help page</a>. <br><br> <a href="http://servers.binf.ku.dk/shs/">Go back</a>'%(in_gene)
		sys.exit(2)
	
	probes_to_use=[]
	for g in probes_names:
		try:
			probes_to_use.append(all_data['rownames'].tolist().index(g))
		except Exception,  e:
			pass
	
	if in_gene2 != None:
		probes_names2 = find_probes(in_gene2, organism)
		
		if probes_names2 == None:
			print head
			print '%s could not be found on array. Maybe you are using an ambiguous gene name. If the gene is in database you might find an alias to your gene name at <a href="http://www.genecards.org/">Genecards</a>. <br><br> To get more details about genes in our database see our <a href="http://servers.binf.ku.dk/shs/help.php">help page</a>.<br><br> <a href="http://servers.binf.ku.dk/shs/">Go back</a>'%(in_gene2)
			sys.exit(2)
			
		probes_to_use2=[]
		for g in probes_names2:
			try:
				probes_to_use2.append(all_data['rownames'].tolist().index(g))
			except Exception,  e:
				pass
		
#########################################################
#make raw datafile to user - as a complement
	if in_gene2==None:
		make_raw_outputFile(all_data['data'][probes_to_use,:], probe_names=all_data['rownames'][probes_to_use], sample_names=all_data['colnames'], adv=adv_mode, use_log = log2gene1, celltypes_to_use=celltypes, filename=in_gene+'_rawData.txt')
	else:
		make_raw_outputFile(all_data['data'][probes_to_use,:], probe_names=all_data['rownames'][probes_to_use], sample_names=all_data['colnames'], adv=adv_mode, use_log = log2gene1, celltypes_to_use=celltypes, filename=in_gene+'_rawData.txt')
		make_raw_outputFile(all_data['data'][probes_to_use2,:], probe_names=all_data['rownames'][probes_to_use2], sample_names=all_data['colnames'], adv=adv_mode, use_log = log2gene2, celltypes_to_use=celltypes, filename=in_gene2+'_rawData.txt')
		

#########################################################	
	
	if in_gene2 == None : # just plot the data
#		print probes_to_use
		data=[]
		if len(probes_to_use) == 0:
			print head
			print '%s could not be found on array. Maybe you are using an ambiguous gene name. If the gene is in database you might find an alias to your gene name at <a href="http://www.genecards.org/">Genecards</a>. <br><br> To get more details about genes in our database see our <a href="http://servers.binf.ku.dk/shs/help.php">help page</a>.<br><br> <a href="http://servers.binf.ku.dk/shs/">Go back</a>'%(in_gene)
			sys.exit(999)

		a = all_data['data'][probes_to_use,:]
		maxes = numpy.argmax(numpy.abs(a), axis=0)
		for index,i in enumerate(maxes):
			d = all_data['data'][probes_to_use,index][i]
			data.append(d)
			
		if organism == 'mouse':
			genename = in_gene.capitalize()
		else: #it is human
			genename = in_gene.upper()
		
		if fold_change:                                    
			if adv_mode:
				classes2use = beeswarm_adv(numpy.array(data),all_data['colnames'], celltypes_to_use=celltypes, title = genename, out= [in_gene+'_fc.pdf', in_gene+'_fc.png'], plot_median=True, use_log=log2gene1) 
			else:
				classes2use = beeswarm(numpy.array(data),all_data['colnames'], title = genename, out= [in_gene+'_fc.pdf', in_gene+'_fc.png'], plot_median=True, use_log=log2gene1) 
		else:                                              
			if adv_mode:
				classes2use = beeswarm_adv(numpy.array(data),all_data['colnames'], celltypes_to_use=celltypes, title = genename, out= [in_gene+'.pdf', in_gene+'.png'], plot_median=True, use_log=log2gene1)
			else:
				classes2use = beeswarm(numpy.array(data),all_data['colnames'], title = genename, out= [in_gene+'.pdf', in_gene+'.png'], plot_median=True, use_log=log2gene1)

	
	else:	   # here we plot the correlation
		import matplotlib.pyplot as plt
		fc=[]
		fc2=[]
		all_names=all_data['colnames']
		if not fold_change:
			for i,s in enumerate(all_data['colnames']):
				maxi=0
				for j in all_data['data'][probes_to_use,i]:
					if j>0:
						if numpy.abs(numpy.exp2(j)) > maxi:
							maxi=numpy.abs(numpy.exp2(j))
							signe=numpy.sign(j)
					else:
						if numpy.abs(numpy.exp2(j)) > maxi:
							maxi=numpy.abs(numpy.exp2(-j))
							signe=numpy.sign(j)
				fc.append(maxi*signe)
				maxi=0
				for j in all_data['data'][probes_to_use2,i]:
					if j>0:
						if numpy.abs(numpy.exp2(j)) > maxi:
							maxi=numpy.abs(numpy.exp2(j))
							signe=numpy.sign(j)
					else:
						if numpy.abs(numpy.exp2(j)) > maxi:
							maxi=numpy.abs(numpy.exp2(-j))
							signe=numpy.sign(j)
				fc2.append(maxi*signe)

		else:
			for i,s in enumerate(all_data['colnames']):
				maxi=0
				for j in all_data['data'][probes_to_use,i]:
					if numpy.abs(j) > maxi:
						maxi=numpy.abs(j)
						signe=numpy.sign(j)
				fc.append(maxi*signe)

				maxi=0
				for j in all_data['data'][probes_to_use2,i]:
					if numpy.abs(j) > maxi:
						maxi=numpy.abs(j)
						signe=numpy.sign(j)
				fc2.append(maxi*signe)
#		print 'using %d probeset(s) for %s'%(len(probes_to_use),in_gene)
#		print 'using %d probeset(s) for %s'%(len(probes_to_use2),in_gene2)

		markers = ['s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ,'s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ,'s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ,'s' , 'o' , '^' , '>' , 'v' , '<' , 'd' , 'p' , 'h' , '8' ]
		
		classes2use_order=[]
		classes2use={}
		
		for index,classe in enumerate(all_names):
			if classe not in classes2use:
				classes2use_order.append(classe)
				classes2use[classe]=[]
				classes2use[classe].append(index)
			else:
				classes2use[classe].append(index)


		if adv_mode:
			reduced_set={}
			reduced_order=[]
			try:
				for i in celltypes:
					reduced_set[i]=classes2use[i]
			except:
				print head
				print 'You are calling for cell types %s not in the dataset.'%(i)
				sys.exit(42)
			classes2use = reduced_set 

			for i in classes2use_order:
				if reduced_set.has_key(i):
					reduced_order.append(i)
			classes2use_order=reduced_order

		#compute corr coef line.
		lin_fit = numpy.polyfit(numpy.array(fc),numpy.array(fc2),1) # this returns the coef of the polynomial fit
		corr_coef = numpy.corrcoef(numpy.array(fc),numpy.array(fc2))[0][1] # this is R
		line_x = numpy.linspace(numpy.array(fc).min(),numpy.array(fc).max()) # this is to have some points to actually draw the line. 

		plots=[]
		for i,classe in enumerate(classes2use):
			a=plt.plot(numpy.array(fc)[classes2use[classe]],numpy.array(fc2)[classes2use[classe]],'o',alpha=.5, color=py_colors[i], marker=markers[i], label=classe)
			plots.append(a)
		if not log2gene1 and not log2gene2:
			plots.append( plt.plot(line_x, line_x*lin_fit[0] + lin_fit[1] , '--b', label='$R^2$ = %.2f'%(corr_coef*corr_coef) )) #we append the plot of the line here
		kwarg={'size':6 }
		
		plt.legend(loc='upper right', prop=kwarg)
	
		if log2gene1:
			plt.xscale('log', basex=2)
			plt.xlabel('Expression (log2) %s'%in_gene)
			plt.xlim(xmax=plt.xlim()[1]*1.1)
		else:
			plt.xlabel('%s'%in_gene)
			plt.xlim(xmax=plt.xlim()[1]*1.05) #make room for the legend
		
		if log2gene2:
			plt.yscale('log', basey=2)
			plt.ylabel('Expression (log2) %s'%in_gene2)
			plt.ylim(ymax=plt.ylim()[1]*1.1) #make room for the legend
		else:
			
			plt.ylabel('%s'%in_gene2)
			plt.ylim(ymax=plt.ylim()[1]*1.05) #make room for the legend
		
		
		
		if organism == 'mouse':
			genename1=in_gene.capitalize()
			genename2=in_gene2.capitalize()
		else: #it is human
			genename1= in_gene.upper()
			genename2= in_gene2.upper()
		
		if not fold_change:
			plt.title('%s vs %s'%(genename1,genename2))
		else:
			plt.title('%s vs %s (versus normal counterpart)'%(genename1,genename2))


		if not fold_change:
			pylab.savefig('%s.pdf'%in_gene)
			pylab.savefig('%s.png'%in_gene)
		else:
			pylab.savefig('%s_fc.png'%in_gene)
			pylab.savefig('%s_fc.pdf'%in_gene)





##################################################################################	
#MAKE HTML text

	abbreviations = {'cd14+ monocytes':'CD14 positive Monocytes','HSC_BM':'Hematopoietic stem cells from bone marrow','early HPC_BM':'Hematopoietic progenitor cells from bone marrow','CMP':'Common myeloid progenitor cell','GMP':'Granulocyte monocyte progenitors','MEP':'Megakaryocyte-erythroid progenitor cell','PM_BM':'Promyelocyte from bone marrow','MY_BM':'Myelocyte from bone marrow','PMN_BM':'Polymorphonuclear cells from bone marrow','PMN_PB':'Polymorphonuclear cells from peripheral blood','AMLI_ETO':'AML with t(8;21)','APL':'AML with t(15;17)','AML with inv(16)/t(16;16)':'AML with inv(16)/t(16;16)','AML with t(11q23)/MLL':'AML with t(11q23)/MLL','LT-HSC':'Long term Hematopoietic stem cell','ST-HSC':'Short term Hematopoietic stem cell','LMPP':'Lymphoid-primed multipotential progenitors','CLP':'Common lymphoid progenitor cells','ETP':'Early T-cell progenitor','ProB':'Pro-B cell','PreB':'Pre-B cell','IgM+SP':'Immunoglobulin M positive side population cells','CD4':'CD4 cells','NKmature':'Mature natural killer cells','GMP':'Granulocyte monocyte progenitors','MkE':'Megakaryocyte erythroid precursors','MkP':'Megakaryocyte precursor','PreCFUE':'Pre-colony-forming unit erythroid cells','CFUE':'Colony-forming unit erythroid cells','ProE':'Erythroid progenitor cells'}


	singletxt='''	<b>Single gene lookup</b><br>
		Each dot in the plot corresponds the expression of '''+in_gene+''' in a microarray. Horizontal lines represent the median expression for each class of cells. '''
	
	singletxt_fc='''	<b>Single gene lookup</b><br>
	Each dot in the plot corresponds the foldchange of '''+in_gene+''' in a microarray. Horizontal lines represent the median expression for each class of cells. '''

	if in_gene2 is not None:
		cortxt='''<b>Correlation</b><br>
		Each dot in the plot represent a microarray experiment with different
		cell types, as specified in the legend of the plot. Expression for gene '''+in_gene+''' is given on the x-axis and expression for 
		'''+in_gene2+''' is given on the y-axis. '''

		cortxt_fc='''<b>Correlation</b><br>
		Each dot in the plot represent a microarray experiment with different
		cell types, as specified in the legend of the plot. The foldchange for '''+in_gene+''' is given on the x-axis and foldchange for '''+in_gene2+''' is given on the y-axis. '''

		
	cor_loglog=''' x- and y-axis are in log2 scale.'''
	cor_lg1='''x-axis is in log2 scale.'''
	cor_lg2='''y-axis is in log2 scale.'''
	single_log='''Expression is given on y-axis on a log2 scale.'''
	
	
	
	
	if WEB_opt: # some display
		if in_gene is not None:
			print head
			
			#print plot
			if not fold_change:
				print '<img src=\"%s.png\" align=center>' % in_gene
			else:
				print '<img src=\"%s_fc.png\" align=center>' % in_gene
				
			#print "what is on plot"
			if fold_change:
				if in_gene2 is None:
					print singletxt_fc
				else:
					print cortxt_fc
			else: 
				if in_gene2 is None:
					print singletxt
				else:
					print cortxt

			#print info on log2 axis - in case user haven't noted from the plot...
			if in_gene2 is not None:#correlations
				if log2gene2 and log2gene1: 
					print cor_loglog
					print 'The linear Pearson product-moment correlation coefficient is: R<sup>2</sup> = %.2f <br>'%(corr_coef*corr_coef)
				elif log2gene2:
					print cor_lg2
				elif log2gene1:
					print cor_lg2
				else: #no log = print add line
					print 'The stippled line represent a theoretical perfect correlation and the R<sup>2</sup> value for the correlation is given in the legend.<br>'
			elif in_gene2 is None: #singles
				if log2gene1:
					print single_log				

			#print abbreviation table
			print '<br><br><br>Abrieviations:<br> <table border="0">'
			for classe in classes2use:
				print '<tr> <td>', classe, '</td> <td>', abbreviations[classe], '</td> </tr>'
			print '</table>'

				
			if not fold_change:
				print '<br><a target=\"_blank\" href=\"%s.pdf\" title=\"\">Get plot in pdf format</a>' % in_gene		
			else:
				print '<a target=\"_blank\" href=\"%s_fc.pdf\" title=\"\">Get plot in pdf format</a>' % in_gene
				
		else:
			print head
			print "Please provide a gene!"
			sys.exit(2)
	
	else:
		if in_gene is not None:
			print 'gene is #%s#' % in_gene
			if in_gene2 is not None:
				print 'gene 2 is #%s#' % in_gene2
		else:
			print "Please provide a gene!"
