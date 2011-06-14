import os
import tracker
import math

def xy(x,y,s,t):
	##############################################
	# Settings file
	##############################################
	# Create settings file
	try:
		setFile=open("plot.settings",'w')
	except IOError:
		print 'Cannot open settings file.' 
		sys.exit(0)
	# Select output type
	if(t=="PS"):
		setFile.write('set terminal postscript\n')
		setFile.write('set output \"%s.ps\"\n'%s)
	elif(t=="PNG"):
		setFile.write('set terminal png monochrome\n')
		setFile.write('set output \"%s.png\"\n'%s)
		setFile.write('set size 0.7,0.5\n')
		#setFile.write('set pointsize 0.7\n')
	else:
		print "Not supported file format!"
		sys.exit(0)
	# Select style
	setFile.write('set data style linespoints\n')
	# Finish settings file
	setFile.write('plot \"%s.dat"'%s)
	setFile.close()
	
	##############################################
	# Data file
	##############################################
	# Create data file
	try:
		datFile=open(s+'.dat','w')
	except IOError:
		print 'Cannot open data file.' 
		sys.exit(0)
	for i in range(len(x)):
		datFile.write('%f %f\n'%(x[i],y[i]))
	datFile.close()
	
	##############################################
	# Plot
	##############################################
	os.system("wgnupl32.exe plot.settings")
	#os.remove("plot.settings")

def time_data(id,k,s="plotdata",t="PS"):
	n=tracker.steps(id)
	x=[]
	y=[]
	for i in range(n):
		x.append(tracker.time(id,i))
		y.append(tracker.data(id,i)[k])
	xy(x,y,s,t)
def data_lambda(id,k,s="plotdata",t="PS"):
	n=tracker.steps(id)
	x=[]
	y=[]
	for i in range(n):
		x.append(tracker.data(id,i)[k])
		y.append(tracker.Lambda(id,i))
	xy(x,y,s,t)

def uy_lambda(id,s="plotdata",t="PS"):
	data_lambda(id,6,s,t)
def uz_lambda(id,s="plotdata",t="PS"):
	data_lambda(id,7,s,t)
def time_ux(id,s="plotdata",t="PS"):
	time_data(id,5,s,t)
def time_uy(id,s="plotdata",t="PS"):
	time_data(id,6,s,t)
def time_vx(id,s="plotdata",t="PS"):
	time_data(id,21,s,t)
def time_vy(id,s="plotdata",t="PS"):
	time_data(id,22,s,t)
def time_ax(id,s="plotdata",t="PS"):
	time_data(id,37,s,t)
def time_ay(id,s="plotdata",t="PS"):
	time_data(id,38,s,t)
