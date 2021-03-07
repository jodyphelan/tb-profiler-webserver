# import subprocess
# import sys
#
# def run_cmd(cmd):
# 	"""
# 	Wrapper to run a command using subprocess with 3 levels of verbosity and automatic exiting if command failed
# 	"""
# 	cmd = "set -u pipefail; " + cmd
# 	res = subprocess.call(cmd,shell=True)
# 	if res!=0:
# 		sys.stderr.write("Command Failed! Please Check!")
# 		exit(1)
#
#
# def cmd_out(cmd):
# 	cmd = "set -u pipefail; " + cmd
# 	try:
# 		res = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
# 		for l in res.stdout:
# 			yield l.decode().rstrip()
# 	except:
# 		sys.stderr.write("Command Failed! Please Check!")
# 		exit(1)
#
# def get_random_file(prefix = None,extension=None):
# 	randint = rand_generator.randint(1,999999)
# 	if prefix:
# 		if extension:
# 			return "%s.%s%s" % (prefix,randint,extension)
# 		else:
# 			return "%s.%s.txt" % (prefix,randint)
# 	else:
# 		if extension:
# 			return "%s.tmp%s" % (randint,extension)
# 		else:
# 			return "%s.tmp.txt" % (randint)
#
# def get_fastq_md5(filename):
# 	md5 = None
# 	for l in cmd_out("gunzip -c %s | awk 'NR%%4==2' | head -1000 | md5" % filename):
# 		row = l.rstrip().split()
# 		md5 = row[0]
# 	return md5
