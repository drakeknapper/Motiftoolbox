#! /usr/bin/python

import posix
import time
import threading
import topParser as tp

class process():
	def __init__(self, pid, T=5, verbose=False):
		self.pid = pid
		self.T = T
		self.verbose = verbose
		self.info = tp.topParser()

		try: 
			self.info.parseProcs(pid=self.pid)

		except: 
			print "Process", pid, "not found."
			exit(-1)

		self.name = self.info.getInfo(self.pid, "COMMAND")

		try: posix.mkdir("procs")
		except: pass
		self._logName = "procs/"+self.name+str(self.pid)+".txt"
		self.log=open(self._logName,'w')
		self.log.close()
		self.start()

	def write(self, text):
		if verbose:
	
			try:
				self.log = open(self._logName,'a')
				self.log.write(str(time.time()-self.time0)+"\t"+text+'\n')
				self.log.close()

			except:
				error("write: An error has occured!")
	

		else:
			pass

	def start(self):
		self.procD = threading.Thread(target=self.mainLoop)
		self.deamon = True
		self.time0 = time.time()
		self.procD.start()

	def stop(self):
		self.deamon = False
		self.procD.join()
		

	def mainLoop(self):

		while self.deamon:
			
			try: 
				self.info.parseProcs(self.pid)

			except:
				print "process.mainLoop: Process",self.pid,"exited. Now exiting..."
				self.deamon = False
				break

			self.write(str(self.info.getInfo(self.pid, "CPU"))+"\t"+str(self.info.getInfo(self.pid, "MEM")))
			time.sleep(self.T)
