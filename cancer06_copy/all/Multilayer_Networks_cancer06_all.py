#!/usr/bin/python
# -*- coding: utf-8 -*-

import matplotlib as plt
from pymnet import *
import pymnet
from itertools import islice
import math

snp_name = ["rs2819348","rs11644424","rs11257188","rs1432679","rs3773651","rs3757318","rs10443215","rs2075570","rs580384",
"rs1493354","rs11654964","rs2032224","rs11196174","rs594206","rs13387042","rs3844412","rs11059635","rs2974935","rs2075555",
"rs6788895","rs2290203","rs7904519"]

gene_name = ["LMOD1", "TGFBR2", "BCL9", "TMEM132C", "COL1A1", "TCF7L2", "PFKFB3", "MTX1", "CADPS", "KLF4", 
"PRC1", "CDKN2B.AS1", "STXBP6", "TMEM132E", "CDH13", "SIAH2", "FPR3", "EBF1", "SLC52A3", "CBX8", "RBMS3", 
"CEP83", "CHST9", "RPUSD4", "APOBEC3A", "PAX9", "WASF2", "ROPN1L", "LINC.PINT", "GRIN2B", "MAML2", "PTPN4", 
"PDE4D", "MSRA", "GCKR", "GALNT2", "TAB2", "CASC8", "TNNT3", "TOX3", "ELL", "DENND1A", "DNAJC1", "TPD52", 
"PBX1", "MIR4757", "CD80", "ADGRE2", "ZNF462", "GPC5", "CDON"]#,'CCDC170', 'RNU6-830P', 'LOC105378740', 'LOC100505878', 'RPL10AP9', 'LOC105371740', 'SERPINB5', 'SERPINB12', 'LOC105373874', 'LOC101928278', 'LOC101928392', 'PRC1-AS1']

mirna_name = ["hsa.mir.139", "hsa.mir.21", "hsa.mir.183", "hsa.mir.96", "hsa.mir.190b", "hsa.mir.6507", "hsa.mir.204", 
"hsa.mir.107", "hsa.mir.1262", "hsa.mir.195", "hsa.mir.379", "hsa.mir.10b", "hsa.mir.484", "hsa.mir.337", 
"hsa.mir.129.2", "hsa.mir.6892", "hsa.mir.374a", "hsa.mir.584", "hsa.mir.548y", "hsa.mir.99a", "hsa.mir.1248", 
"hsa.mir.423", "hsa.mir.145", "hsa.mir.486.2", "hsa.mir.1301", "hsa.mir.486.1", "hsa.mir.744", "hsa.mir.142", 
"hsa.mir.592", "hsa.mir.129.1", "hsa.mir.130b", "hsa.mir.877", "hsa.mir.203a", "hsa.mir.1287", "hsa.mir.3614", 
"hsa.mir.342", "hsa.mir.10a", "hsa.mir.7705", "hsa.mir.511", "hsa.mir.378c", "hsa.mir.206", "hsa.mir.103a.1", 
"hsa.mir.1247", "hsa.mir.1185.2", "hsa.mir.3174", "hsa.mir.1295a", "hsa.mir.141", "hsa.mir.155", "hsa.mir.493", 
"hsa.mir.215"]

protein_name = ["Bax", "GSK3.alpha.beta", "E.Cadherin", "Rab11", "Caveolin.1", "Collagen_VI", 
"c.Myc", "PKC.alpha", "GAPDH", "P.Cadherin", "PDK1", "XBP1.G.C", "CDK1", "c.Kit", "Syk", 
"p27", "Shc_pY317", "mTOR", "AR", "Ku80", "Rictor", "MSH6", "14.3.3_zeta", "PI3K.p85", 
"Chk2_pT68", "Bap1.c.4", "Lck", "STAT3_pY705", "S6_pS240_S244", "ATM", "Src_pY416", 
"14.3.3_epsilon", "Bak", "N.Ras", "Myosin.IIa_pS1943", "LKB1", "YB.1", "p70S6K", 
"Claudin.7", "Annexin.1", "4E.BP1_pT37_T46", "MIG.6", "FOXO3a"]

net=pymnet.MultilayerNetwork(aspects=1, fullyInterconnected=False)
#---------------------------------------------------------------------------------------------------------------------
# get line
def isGeneRnaProteinlist(name):
	if name in snp_name: return name, 'SNP'
	if name in gene_name: return name, 'Gene'
	if name in mirna_name: return name, 'miRNA'
	if name in protein_name: return name, 'Protein'


MIC06 = open('MIC0.6.csv','r')
for line in islice(MIC06, 1, None): # ignore 1st line
	name1 = line.split(',')[0]
	name2 = line.split(',')[1]
	value = line.split(',')[2]
#1	if (name1 in gene_name and name2 in snp_name) or (name1 in gene_name and name2 in gene_name) or (name1 in mirna_name and name2 in mirna_name) or (name1 in protein_name and name2 in protein_name): # MIC0.6_2.csv
#		net[isGeneRnaProteinlist(name1)][isGeneRnaProteinlist(name2)] = float(value)
	net[isGeneRnaProteinlist(name1)][isGeneRnaProteinlist(name2)] = float(value)
#	if (name1 in gene_name and name2 in snp_name) or (name1 in gene_name and name2 in mirna_name) or (name1 in gene_name and name2 in protein_name):# or (name1 in mirna_name and name2 in protein_name):
#		net[isGeneRnaProteinlist(name1)][isGeneRnaProteinlist(name2)] = float(value)
MIC06.close()

#---------------------------------------------------------
MIC06 = open('MIC0.6.csv','r')
snp_unique, gene_unique, mirna_unique, protein_unique = [], [], [], []
for line in islice(MIC06, 1, None): # ignore 1st line
	name1 = line.split(',')[0]
	name2 = line.split(',')[1]
	if name1 in snp_name and name1 not in snp_unique: snp_unique.append(name1)
	if name2 in snp_name and name2 not in snp_unique: snp_unique.append(name2)
	if name1 in gene_name and name1 not in gene_unique: gene_unique.append(name1)
	if name2 in gene_name and name2 not in gene_unique: gene_unique.append(name2)
	if name1 in mirna_name and name1 not in mirna_unique: mirna_unique.append(name1)
	if name2 in mirna_name and name2 not in mirna_unique: mirna_unique.append(name2)
	if name1 in protein_name and name1 not in protein_unique: protein_unique.append(name1)
	if name2 in protein_name and name2 not in protein_unique: protein_unique.append(name2)

MIC06.close()


MIC06 = open('MIC0.6.csv','r')
sgmp = {}
for line in islice(MIC06, 1, None): # ignore 1st line
	key = line.split(',')[0]
	value = line.split(',')[1]
	sgmp.setdefault(key,[]).append(value)
	sgmp.setdefault(value,[]).append(key)

MIC06.close()

'''
MIC06 = open('MIC0.6.csv','r')
snp2gene = []
for line in islice(MIC06, 1, None): # ignore 1st line
	name1 = line.split(',')[0]
	name2 = line.split(',')[1]
	if name2 in snp_unique:
		snp2gene.append(name1)

MIC06.close()


for i in snp2gene:
	for j in sgmp[i]:
		if j in protein_name:
			for k in sgmp[i]:
				net[isGeneRnaProteinlist(i)][isGeneRnaProteinlist(k)]=1
			break

for i in mirna_unique:
	for j in sgmp[i]:
		if j in gene_name:
			for k in sgmp[i]:
				net[isGeneRnaProteinlist(i)][isGeneRnaProteinlist(k)]=1
			break


for i in protein_unique:
	for j in sgmp[i]:
		if j in mirna_name:
			for k in sgmp[i]:
				net[isGeneRnaProteinlist(i)][isGeneRnaProteinlist(k)]=1
			break
'''
#---------------------------------------------------------------------------------------------------------------------
# get nodeColorDict && nodeLabelColorDict
def nodeColor(name):
	if name in snp_name: return {(name, 'SNP'):'#00B2EE'}
	if name in gene_name: return {(name, 'Gene'):'g'}
	if name in mirna_name: return {(name, 'miRNA'):'#FF8C00'}
	if name in protein_name: return {(name, 'Protein'):'k'}

MIC06 = open('MIC0.6.csv','r')
nodeColorDict0, nodeLabelColorDict0 = {}, {}
for line in islice(MIC06, 1, None): # ignore 1st line
	name1 = line.split(',')[0]
	name2 = line.split(',')[1]
	value = line.split(',')[2]
	nodeColorDict0 = dict(nodeColorDict0.items() + nodeColor(name1).items() + nodeColor(name2).items())
	nodeLabelColorDict0 = dict(nodeLabelColorDict0.items() + nodeColor(name1).items() + nodeColor(name2).items())
MIC06.close()

#---------------------------------------------------------------------------------------------------------------------

# get edgeColorDict
def edgeColor(name1, name2):
	if name1 in gene_name and name2 in gene_name:
		return {((name1, 'Gene'), (name2, 'Gene')):'#FF3030'}
	if name1 in mirna_name and name2 in mirna_name:
		return {((name1, 'miRNA'), (name2, 'miRNA')):'#FF3030'}
	if name1 in protein_name and name2 in protein_name:
		return {((name1, 'Protein'), (name2, 'Protein')):'#FF3030'}
	if name1 in gene_name and name2 in snp_name:
		return {((name1, 'Gene'), (name2, 'SNP')):'#BF3EFF'}
	if name1 in gene_name and name2 in mirna_name:
		return {((name1, 'Gene'), (name2, 'miRNA')):'#F08080'}
	if name1 in gene_name and name2 in protein_name:
		return {((name1, 'Gene'), (name2, 'Protein')):'#008000'}
	if name1 in mirna_name and name2 in protein_name:
		return {((name1, 'miRNA'), (name2, 'Protein')):'#0099CC'}

MIC06 = open('MIC0.6.csv','r')
edgeColorDict0 = {}
for line in islice(MIC06, 1, None): # ignore 1st line
	name1 = line.split(',')[0]
	name2 = line.split(',')[1]
	value = line.split(',')[2]
	edgeColorDict0 = dict(edgeColorDict0.items() + edgeColor(name1, name2).items())
MIC06.close()

#---------------------------------------------------------------------------------------------------------------------
# get nodeSizeDict
def nodeSize(name, count):
	if name in snp_name: return {(name, 'SNP'): count}
	if name in gene_name: return {(name, 'Gene'): count}
	if name in mirna_name: return {(name, 'miRNA'): count}
	if name in protein_name: return {(name, 'Protein'): count}

MIC06 = open('MIC0.6.csv','r')
nodeCountDict = {}
for line in islice(MIC06, 1, None): # ignore 1st line
	name1 = line.split(',')[0]
	name2 = line.split(',')[1]
	if name1 not in nodeCountDict.keys():
		nodeCountDict[name1] = 1
	else:
		nodeCountDict[name1] += 1
	if name2 not in nodeCountDict.keys():
		nodeCountDict[name2] = 1
	else:
		nodeCountDict[name2] += 1

MIC06.close()
nodeSizeDict0 = {}
nodeLabelSizeDict0 = {}
for name in nodeCountDict.keys():
	nodeSizeDict0 = dict(nodeSizeDict0.items() + nodeSize(name, nodeCountDict[name]/40.0).items())
	nodeLabelSizeDict0 = dict(nodeLabelSizeDict0.items() + nodeSize(name, 5).items())

#---------------------------------------------------------------------------------------------------------------------
# get edge width
def edgeWidth(name1, name2, name3):
	if name1 in gene_name and name2 in gene_name:
		return {((name1, 'Gene'), (name2, 'Gene')):name3}
	if name1 in mirna_name and name2 in mirna_name:
		return {((name1, 'miRNA'), (name2, 'miRNA')):name3}
	if name1 in protein_name and name2 in protein_name:
		return {((name1, 'Protein'), (name2, 'Protein')):name3}
	if name1 in gene_name and name2 in snp_name:
		return {((name1, 'Gene'), (name2, 'SNP')):name3}
	if name1 in gene_name and name2 in mirna_name:
		return {((name1, 'Gene'), (name2, 'miRNA')):name3*1.5}
	if name1 in gene_name and name2 in protein_name:
		return {((name1, 'Gene'), (name2, 'Protein')):name3*1.5}
	if name1 in mirna_name and name2 in protein_name:
		return {((name1, 'miRNA'), (name2, 'Protein')):name3*1.5}

MIC06 = open('MIC0.6.csv','r')
edgeWidthDict0 = {}
for line in islice(MIC06, 1, None): # ignore 1st line
	name1 = line.split(',')[0]
	name2 = line.split(',')[1]
	value = float(line.split(',')[2])
	edgeWidthDict0 = dict(edgeWidthDict0.items() + edgeWidth(name1, name2, value).items())
MIC06.close()

#---------------------------------------------------------------------------------------------------------------------
# get coords, must be simultaneously
def Coords(gmp_unique):
	Count = len(gmp_unique)
	dis = 1.6/9 #1.6/int(math.sqrt(Count))
	dis2 = 1.6/8
	Dict = {}
	for i in range(Count):
		zheng = i//9 #(int(math.sqrt(Count))+1)
		yu = i%9 #(int(math.sqrt(Count))+1)
		zheng2 = i//8
		yu2 = i%8
		if set(gmp_unique).issubset(set(snp_name)): Dict[(gmp_unique[i], 'SNP')] = (-0.8+yu*dis, 0.8-zheng*dis)
		if set(gmp_unique).issubset(set(gene_name)): Dict[(gmp_unique[i], 'Gene')] = (-0.8+yu*dis, 0.8-zheng*dis)
		if set(gmp_unique).issubset(set(mirna_name)): Dict[(gmp_unique[i], 'miRNA')] = (-0.8+yu2*dis2, 0.8-zheng2*dis2)
		if set(gmp_unique).issubset(set(protein_name)): Dict[(gmp_unique[i], 'Protein')] = (-0.8+yu*dis, 0.8-zheng*dis)
	return Dict

nodeCoords0 = dict(Coords(snp_unique).items() + Coords(gene_unique).items() + Coords(mirna_unique).items() + Coords(protein_unique).items())

#---------------------------------------------------------------------------------------------------------------------
# main
fig=pymnet.draw(net, show=False, azim=-60, elev=25, layergap=0.75,
						layout='shell',
#						layershape='circle',
						autoscale=True,
						defaultNodeLabelAlpha=0.75,
						defaultLayerAlpha=0.25,

						# color
						nodeColorDict = nodeColorDict0,
						nodeLabelColorDict = nodeColorDict0,
						edgeColorDict = edgeColorDict0,

						# size
						nodeSizeDict = nodeSizeDict0,
						nodeSizeRule={'scalecoeff': 0.2, 'rule': 'scaled'},
						nodeLabelSizeDict = nodeLabelSizeDict0,

						# width
						edgeWidthDict = edgeWidthDict0,
						defaultEdgeAlpha=1,

						# coords
						nodeCoords = nodeCoords0,
						nodelayerCoords = nodeCoords0,

						layerColorDict={'Gene':'#DAA520', 'Protein':'#FF0000'},
						layerOrderDict={'SNP':4, 'Gene':3, 'miRNA':2, 'Protein':1},
						# node label style
						defaultNodeLabelStyle='oblique', # normal, italic or oblique
						layerPadding=0.25)
fig.savefig("Multilayer_Networks_cancer06_all.pdf")
