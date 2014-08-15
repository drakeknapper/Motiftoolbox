#! /usr/bin/python

import sys
import threading
import os
import posix
import time
import cPickle
import multiprocessing


start = time.time()


def sum(f, parameter, param_values, max_threads=None, kwargs={}):
	s = summer(max_threads=max_threads)
	return s.sum(function=f, parameter=parameter, param_values=param_values, kwargs=kwargs)


def distribute(f, parameter, param_values, kwargs={}, max_threads=None, nice=19, verbose=0):
	d = distributor(max_threads=max_threads, nice=nice, verbose=verbose)
	return d.distribute(function=f, parameter=parameter, param_values=param_values, kwargs=kwargs)


#===========================================================================

def error(text):
	errLog=open("ERROR_LOG.txt", 'w')
	errLog.write(text+'\n')
	errLog.close()

#===========================================================================

class master(object):

	def __init__(self, max_threads=None, nice=19, verbose=0):

		if max_threads == None:
			self.max_threads = multiprocessing.cpu_count()

		else:
			self.max_threads = max_threads

		self.nice = nice
		self.verbose = verbose

		name = sys.argv[0].split('.')[1][1:]

		if self.verbose:
			self.log_name = name+"_log.dat"
			self.log = open(self.log_name, 'w')
			self.log.close()

		self.ctrl_name = name+"_ctrl.dat"
		ctrl = open(self.ctrl_name, 'w')
		ctrl.write(str(self.max_threads)+','+str(self.nice))
		ctrl.close()

		self.writeLock = threading.Lock()

		if self.verbose:
			self.log_lock = threading.Lock()

		self.sigPause = threading.Event()

		self.threadBin = []
		self.runningThreads = 0


	def getCtrl(self):
		ctrl = open(self.ctrl_name, 'r')
		tmp = ctrl.read().split(',')
		self.max_threads = int(tmp[0])
		self.nice = int(tmp[1])
		ctrl.close()


	def write(self, text):
		self.log_lock.acquire()
		self.log = open(self.log_name, 'a')

		try:
			self.log.write(str(time.time()-start)+" "+text+'\n')
		
		except:
			error("write: An error has occured!")

		self.log.close()
		self.log_lock.release()


	def operate(self,**kwargs):
		self.runningThreads+=1
		r, w = os.pipe()
		pid = os.fork()
		if pid:	# PARENT
			os.close(w)

			tmpParam = kwargs[self.parameter]
			r = os.fdopen(r, 'rb')
			try: tmpResult = cPickle.load(r)
			except:
				try: posix.kill(pid, -1)				# Kill the child
				except: pass

				if self.verbose: self.write("An error occured in job "+str(pid))	# Tell the reader

				self.jobs.append(kwargs[self.parameter])		# Append the job to be done later
				self.runningThreads -= 1			# Tell your master that you've one it
				self.sigPause.set()
				return
			os.waitpid(pid, 0)

		else:	# CHILD
			posix.nice(self.nice)
			os.close(r)
			w = os.fdopen(w, 'wb')
			cPickle.dump(self.function(**kwargs), w, 2)
			w.close()
			sys.exit(0)


		self.writeLock.acquire()
		self.operate_return(tmpParam, tmpResult)
		self.runningThreads -= 1
		self.sigPause.set()
		self.writeLock.release()

#===========================================================================

class distributor(master):

	def __init__(self, max_threads=4, nice=19, verbose=0):
		master.__init__(self, max_threads, nice, verbose)


	def distribute(self, function, parameter, param_values, kwargs={}, raw=False): # e. g. mainloop
		if self.verbose: self.write("Master: Distributing the Paramter "+parameter+" among the net...")

		self.function = function
		self.parameter = parameter
		self.jobs = list(param_values)
		self.result = {}

		for j in self.jobs:
			committed = False

			while not committed:

				if self.runningThreads < self.max_threads:
					if self.verbose: self.write("Master: Distributing "+parameter+"="+str(j))
					kwargs[parameter] = j
					self.threadBin.append(threading.Thread(target=self.operate, kwargs=kwargs))
					self.threadBin[-1].start()
					committed=True
					time.sleep(0.01)

				else:
					if self.verbose: self.write("Master: Waiting for threads to finish...") 
					self.sigPause.wait()
					self.sigPause.clear()
					self.getCtrl()
		
		if self.verbose: self.write("====================FINISHING PROCESSES=====================")

		while self.runningThreads > 0:
			self.sigPause.wait()
			self.sigPause.clear()

		if raw:
			return self.result # Don't wait until all the results have been caught. ERROR!!!

		else:

			try:
				res=[]
				for i in range(len(self.result)):
					res.append(self.result[param_values[i]])

			except:
				if self.verbose: self.write("master: Not all kwargs found! Exiting...")

				for i in range(len(self.result)):
					print self.result.items()[i]
					print "ERROR: See file",self.log_name

				exit(-1)

			return res


	def operate_return(self, tmpParam, tmpResult):
		self.result[tmpParam] = tmpResult


#===========================================================================


class summer(master):

	def __init__(self, max_threads=4, verbose=0):
		master.__init__(self, max_threads, verbose)


	def sum(self, function, parameter, param_values, kwargs={}, raw=False): # e. g. mainloop
		if self.verbose: self.write("Master: Sum the output of function...")

		self.function = function
		self.parameter = parameter
		self.jobs = param_values
		self.result = None

		for j in self.jobs:
			committed=False
			while not committed:
				if self.runningThreads < self.max_threads:
					if self.verbose: self.write("Master: Distributing "+parameter+"="+str(j))
					kwargs[parameter] = j
					self.threadBin.append(threading.Thread(target=self.operate, kwargs=kwargs))
					self.threadBin[-1].start()
					committed=True
					time.sleep(0.01)
				else:
					if self.verbose: self.write("Master: Waiting for threads to finish...") 
					self.sigPause.wait()
					self.sigPause.clear()
		
		if self.verbose: self.write("====================FINISHING PROCESSES=====================")

		while self.runningThreads > 0:
			self.sigPause.wait()
			self.sigPause.clear()

		return self.result # Don't wait until all the results have been caught. ERROR!!!


	def operate_return(self, tmpParam, tmpResult):

		if self.result == None:
			self.result = tmpResult

		else:
			self.result += tmpResult
















