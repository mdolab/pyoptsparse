#!/usr/bin/env python
'''
pyOpt_history

Holds the Python Design Optimization History Class.

Copyright (c) 2008-2013 by pyOpt Developers
All rights reserved.
Revision: 1.0   $Date: 11/12/2009 21:00$


Developers:
-----------
- Dr. Ruben E. Perez (RP)
- Mr. Peter W. Jansen (PJ)

History
-------
	v. 1.0  - Initial Class Creation (PJ,RP 2009)
'''

__version__ = '$Revision: $'

'''
To Do:
	- shevling of optimizer
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys
import pdb
import array as ARRAY
import shelve

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================
#import extension

# =============================================================================
# Misc Definitions
# =============================================================================



# =============================================================================
# History Class
# =============================================================================
class History(object):
	
	'''
	Abstract Class for Optimizer History Object
	'''
	
	def __init__(self, filename, mode, optimizer=None, opt_prob=None, *args, **kwargs):
		
		'''
		Optimizer History Class Initialization
		
		**Arguments:**
		
		- filename  -> STR: Name for .bin and .cue file
		- mode      -> STR: Either read ('r') or write ('w') mode
		
		**Keyword arguments:**
		
		- optimizer -> INST: Opimizer class instance,  *Default* = None
		- opt_prob  -> STR: Optimization Problem Name, *Default* = None
		
		Documentation last updated:  April. 14, 2010 - Peter W. Jansen
		'''
		
		# 
		self.filename = filename
		self.mode = mode
		
		# 
		bin_name = filename + '.bin'
		cue_name = filename + '.cue'
		#opt_name = filename + '.opt'
		
		#
		if self.mode == 'w':
		
			# 
			if os.path.isfile(bin_name):
				os.remove(bin_name)
			#end
			if os.path.isfile(cue_name):
				os.remove(cue_name)
			#end
			#if os.path.isfile(opt_name):
			#	os.remove(opt_name)
			##end   
			
		else:
			
			if not os.path.isfile(bin_name):
				raise NameError('Error: filename %s.bin does not exist'%(filename))
			#end
			if not os.path.isfile(cue_name):
				raise NameError('Error: filename %s.cue does not exist'%(filename))
			#end
			#if not os.path.isfile(opt_name):
			#	raise NameError('Error: filename %s.opt does not exist'%(filename))
			##end
			
		#end
		
		# 
		self.bin_file = open(bin_name,mode+'b')
		self.cue_file = open(cue_name,mode)
		#self.opt_file = shelve.open(opt_name)
		
		#
		if self.mode == 'w':
			
			#
			if (optimizer == None):
				optname = 'None'
			else:
				optname = optimizer.name
			#end
			header = 'History for %s solving %s\n' %(optname,opt_prob)
			self.cue_file.write(header)
			#self.opt_file['optimizer'] = optimizer
			#self.opt_file.close()
			
		elif self.mode == 'r':
			
			#
			self.cues = {}
			self.icount = {}
#			if (self.opt_file['optimizer'].__dict__ != optimizer.__dict__):
#				print 'Warning Optimizer Instance and stored Optimizer do not match -- hot start aborted'
#				self.cue_file.close()
#				self.opt_file.close()
#				self.s_count = 0
#				return
#			#end
			lines = self.cue_file.readlines()			
			for line in lines[1:]:
				
				if len(line) < 3:
					break
				else:
					#read in positions
					tline = line.split()
					if self.cues.has_key(tline[2]):
						self.cues[tline[2]].append([int(tline[0]),int(tline[1])])
					else:
						self.cues[tline[2]] = [[int(tline[0]),int(tline[1])]]
						self.icount[tline[2]] = 0						
					#end
				#end
			#end 
			self.cue_file.close()
		#end
		
		# 
		self.s_count = 0
		
		
	def close(self):
		
		'''
		Close Optimizer History Files
		
		Documentation last updated:  December. 11, 2009 - Ruben E. Perez
		'''
		
		self.bin_file.close()
		if self.mode == 'w':
			self.cue_file.close()
		#end
		
		
	def read(self, index=[], ident=['obj']):
		
		'''
		Read Data from Optimizer History Files
		
		**Keyword arguments:**
		
		- index -> LIST,SCALAR: Index (list), [0,-1] for all, [] internal count, -1 for last, *Default* = []
		- ident -> STR: Indentifier, *Default* = 'obj'
		
		Documentation last updated:  April. 14, 2010 - Peter W. Jansen
		'''
		
		# index = [0,-1]
		bdata = {}
		hist_end = False
		
		for id in ident:
			bdata[id] = []
			if id in self.cues.keys():
				if isinstance(index,int):
					if (index == -1):
						index = len(self.cues[id])-1
					#end
					
					index = [index, index+1]
				elif isinstance(index,list):
					if (index == []):
						index = [self.icount[id], self.icount[id]+1]
						self.icount[id] += 1
					elif (index == [0,-1]):
						index = [0, len(self.cues[id])]
					#end
				else:
					raise ValueError('Index type not understood - must be either int or list')
				#end
			else:
				hist_end = True
				return (bdata,hist_end)
			#end
			for i in xrange(index[0],index[1]):
				
				#
				if (i >= len(self.cues[id])):
					hist_end = True
					return (bdata,hist_end)
				#end
				tvals = ARRAY.array('d')
				self.bin_file.seek(self.cues[id][i][0]*8,0)
				tvals.fromfile(self.bin_file,self.cues[id][i][1])
				bdata[id].append(numpy.array(tvals))
			#end
		#end
		
		return (bdata, hist_end)
		
		
	def write(self,bin_data,cue_data):
		
		'''
		Write Data to Optimizer History Files
		
		**Arguments:**
		
		- bin_data -> LIST/ARRAY: Data to be written to binary file
		- cue_data -> STR: Variable identifier for cue file
		
		Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
		'''        
		
		#
		bin_data = numpy.array(bin_data)
		tdata = ARRAY.array('d',bin_data.flatten())
		tdata.tofile(self.bin_file)
		self.bin_file.flush()
		
		# 
		self.cue_file.write('%d %d %s\n'%(self.s_count,len(bin_data.flatten()), cue_data))
		self.cue_file.flush()
		
		#
		self.s_count += len(bin_data.flatten())
		
		return
		
		
	def overwrite(self,bin_data,index):
		
		'''
		Overwrite Data on Optimizer History Files
		
		**Arguments:**
		
		- bin_data -> ARRAY: Data to overwrite old data
		- index -> INT: Starting index of old data
		
		Documentation last updated:  Feb. 03, 2011 - Peter W. Jansen
		'''        
		
		#
		bin_data = numpy.array(bin_data)
		tdata = ARRAY.array('d',bin_data.flatten())
		self.bin_file.seek(index*8,0)
		tdata.tofile(self.bin_file)
		self.bin_file.flush()
		self.bin_file.seek(0,2)
		
		return
	


#==============================================================================
# Optimizer History Test
#==============================================================================
if __name__ == '__main__':
	
	# Test Optimizer History
	print 'Testing Optimizer History...'
	hst = History()
	
