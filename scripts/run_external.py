import subprocess

def run_external(command: list[str], shell=False, stdin=None) -> tuple[str, str]:
	"""run external command and return stout and stderr
	"""
	if stdin is None:
		result = subprocess.run(command, capture_output=True, text=True, shell=shell)
	else:
		result = subprocess.run(command, capture_output=True, text=True, input=stdin, shell=shell)

	return result.stdout, result.stderr
