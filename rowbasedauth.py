
"""
rowbasedauth.py
"""


import sys, os, string, time, re


BLACKLISTDIR = "/usr/local/Proasis2/Data/BLACKLIST"



def GetUsername(inpStr):
	"""
	Get username from cookie Proasis3User (or Proasis2User)
	"""

	inpStr = str(inpStr)
	c = dict(re.findall(r'\[(.+?)]=\[(.+?)]', inpStr))
	if "Proasis3User" in c.keys():
		return c["Proasis3User"]
	elif "Proasis2User" in c.keys():
		return c["Proasis2User"]
	else:
		return ""


def GetBlacklist(inpStr):
	"""
	for row based authentication
	"""

	retDict = {}

	# get curUser u from inpStr
	u = GetUsername(inpStr)

	#
	# blacklist file has filename BLACKLISTDIR/curUser.dat
	# - it is a python dictionary of type {regnoA:1, regnoB:1, ...}
	# - where hits with regnoA,regnoB cannot be viewed by curUser
	# these files created by separate cron job
	#
	curf = "%s/%s%s" % (BLACKLISTDIR, u, ".dat")

	if not os.path.isfile(curf):
		return retDict

	# get blacklist as python dictionary
	retDict = MyMethodToReadFile(curf)

	return retDict
