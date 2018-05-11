"""
Helper Functions Library

Preston Huft, 2018. Updated 11 May.

Contains functions I use often enough to really not want to rewrite them each 
time that I use them. Could also contain classes, I suppose.
"""

## LIBRARIES

from random import random as rn

## METHODS 
def nint(fnumber):
	""" Takes an float and rounds it naturally to an integer."""
	return int(0.5+fnumber)
	
def add_nested_lists(list1,list2):
	sumlist = []
	if (len(list1)==len(list2)):
		for i in range(0,len(list1)):
			sumlist.append([sum(x) for x in zip(list1[i],list2[i])])
	else:
		print("Lists aren't equal length.")
	return sumlist
	
def copy_nested_list(mylist):

		""" This assumes a depth of 1, i.e. a list of list of non-lists."""
		copy_list = []
		for l in mylist:
			copy_list.append(list(l))
		return copy_list
		
def rand_list(length,offset,scale):
	"""Returns a scalar list of len 'length' whose elements are random doubles
		centered around the double 'offset', between +/-'scale'/2.
	"""
	return [scale*(rn()-.5) + offset for i in range(0,length)]