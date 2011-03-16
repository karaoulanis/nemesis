import os
import re
from time import time,ctime

def getFiles(root_path):
	""" Get files/paths from a given directory
	"""
	print "Searching in the src directory..."
	paths=[]
	names=[]
	for path, folders, files in os.walk(root_path):
		for name in files:
			filename,extension=os.path.splitext(name)
			if extension in ['.cpp','.h']:
				head,tail=os.path.split(path)
				paths.append(head+'/'+tail)
				names.append(name)
	return dict(zip(names, paths))
def findFunction(line,keywords):
	""" Find function declerations
	"""
	line=line.strip()
	if len(line)==0:
		return None
	words=line.replace('(',' ( ').replace(')',' ) ').split()
	if words[0] in keywords and line.find('(') and line.find(')'):
		ret=[line.replace(';','')]
		try:
			index=words.index('(')
			ret.append(words[index-1])
		except:
			ret=None
		return ret
if __name__=='__main__':
	""" Main function
	"""
	paths=getFiles('..'+'/src')
	names=paths.keys()
	names.sort()
	keywords=['int','double','boolean','const','void']
	dirlist=os.listdir('.')
	counter=0
	for name in names:
		filename,extension=os.path.splitext(name)
		if extension=='.h':
			#=======================================================
			# init
			#=======================================================
			class_name='%s'%filename
			test_name='%sTest'%filename
			cpp_name='%s_unittest.cpp'%filename
			if cpp_name in dirlist:
				continue
			keywords.append(class_name)
			counter+=1
			
			#=======================================================
			# open file
			#=======================================================
			path='%s/%s'%(paths[name],name)
			print path,"->",cpp_name
			ofile=open("raw/"+cpp_name,'w')
			ofile.write("// automatic generated unit test file	\n")
			ofile.write("// author  : Fotios E. Karaoulanis		\n")
			ofile.write('// date    : %s\n'%ctime(time()))
			ofile.write('// file    : %s\n'%path)
			ofile.write('\n')
		
			#=======================================================
			# write includes
			#=======================================================
			directory=paths[name].replace('../','')
			ofile.write('#include <%s/%s>\n'%(directory,name))
			ofile.write('#include <gtest/gtest.h>\n')
			ofile.write('\n')
			#~ incs=dfs(graph,name)
			#~ for inc in incs:
				#~ directory=paths[inc].replace('../','')
				#~ ofile.write('#include <%s/%s>\n'%(directory,inc))
			#~ ofile.write('\n')
			
			#=======================================================
			# read data from header file
			#=======================================================
			ifile=open(path,'r')
			lines=ifile.readlines()
			ifile.close()
			
			#=======================================================
			# write fixture
			#=======================================================
			ofile.write('class %s: public testing::Test {\n'%test_name)
			ofile.write('protected:\n')
			ofile.write('\t%s %s();\n'%(class_name,class_name.lower()))
			ofile.write('\tvirtual void SetUp()\n\t{\n\t}\n')
			ofile.write('\tvirtual void TearDown()\n\t{\n\t}\n')
			ofile.write('};\n\n')
			
			#=======================================================
			# write test for every (possible) function
			#=======================================================
			for line in lines:
				ret=findFunction(line,keywords)
				if ret!=None:
					(function,function_name)=ret
					ofile.write('// test: %s\n'%(function))
					ofile.write('TEST_F(%s,%s)\n{\n}\n\n'%(test_name,function_name))
					
			#=======================================================
			# close file
			#=======================================================
			ofile.close()
	raw_input('\n%i files left to be tested.\nPress <enter> to continue...'%counter)