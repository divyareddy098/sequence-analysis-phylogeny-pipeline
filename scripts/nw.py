# We'll use binary flags to encode direction: up, left, diagonal
D_FLAGS = {
	"up": 4,
	"left": 2,
	"diag": 1
}

def _initialize_matrix(num_rows: int, num_columns: int) -> list[list[int]]:
	""" Make a matrix of specified dimensions filled with 0s
	"""
	matrix = []

	for _ in range(num_rows):
		row = []
		for _ in range(num_columns):
			row.append(0)
		matrix.append(row)

	return matrix


def fill_matrix(
	seq_a: str,
	seq_b: str,
	match: int,
	mismatch: int,
	gap: int
) -> tuple[list[list[int]], list[list[int]]]:
	""" Fill a score and direction matrix using needleman wunsch algorithm approach
	"""
	# First make a matrix of the appropriate dimensions for score and direction
	n_rows = len(seq_a)+1
	n_cols = len(seq_b)+1
	score_mat = _initialize_matrix(n_rows, n_cols)
	direct_mat = _initialize_matrix(n_rows, n_cols)

	# fill first row and column
	for i in range(n_cols):
		score_mat[0][i] = gap*i
		direct_mat[0][i] = D_FLAGS["left"]

	for i in range(n_rows):
		score_mat[i][0] = gap*i
		direct_mat[i][0] = D_FLAGS["up"]

	# reset cell 0,0 in direct mat for tidyness as it's not used
	direct_mat[0][0] = 0

	# fill the matrices
	for i in range(1, n_rows):
		for j in range(1, n_cols):
			base_a = seq_a[i-1]
			base_b = seq_b[j-1]
			if base_a == base_b:
				match_mismatch = match
			else:
				match_mismatch = mismatch

			up = score_mat[i-1][j] + gap
			left = score_mat[i][j-1] + gap
			diag = score_mat[i-1][j-1] + match_mismatch

			best_score = max(up, left, diag)
			score_mat[i][j] = best_score

			if up == best_score:
				direct_mat[i][j] += D_FLAGS["up"]
			if left == best_score:
				direct_mat[i][j] += D_FLAGS["left"]
			if diag == best_score:
				direct_mat[i][j] += D_FLAGS["diag"]

	return score_mat, direct_mat


def traceback_single(
	seq_a: str,
	seq_b: str,
	score_mat: list[list[int]],
	direct_mat: list[list[int]]
) -> tuple[tuple[str, str], int]:
	i = len(seq_a)
	j = len(seq_b)

	aln_a = []
	aln_b = []
	score = score_mat[i][j]

	while i != 0 or j != 0:
		direction_bits = direct_mat[i][j]
		# use bitwise operations to check flags
		if D_FLAGS["diag"] & direction_bits:
			i -= 1
			j -= 1
			aln_a.append(seq_a[i])
			aln_b.append(seq_b[j])
			continue
		if D_FLAGS["up"] & direction_bits:
			i -= 1
			aln_a.append(seq_a[i])
			aln_b.append("-")
			continue
		if D_FLAGS["left"] & direction_bits:
			j -= 1
			aln_a.append("-")
			aln_b.append(seq_b[j])
			continue

	aln_a = "".join(aln_a[::-1])
	aln_b = "".join(aln_b[::-1])

	return (aln_a, aln_b), score


def needleman_wunsch(
	seq_a: str,
	seq_b: str,
	match: int,
	mismatch: int,
	gap: int,
) -> tuple[tuple[str, str], int]:
	
	score_mat, direct_mat = fill_matrix(seq_a, seq_b, match, mismatch, gap)
	alns, score = traceback_single(seq_a, seq_b, score_mat, direct_mat)
	
	return alns, score

def print_mats(matrix: list[list[str|int]], seqa: str, seqb: str, dir_mat: bool = False) -> None:
	if dir_mat:
		char_dict = {
			1: "\\", # diagonal
			2: "-", # left
			3: ">", # diag + left
			4: "'", # up
			5: "V", # diag + up
			6: "+", # up + left
			7: "*" # all
		}
	else:
		char_dict = {}
	# get max chars in cells
	max_chars = 0
	for row in matrix:
		for cell in row:
			if len(str(cell)) > max_chars:
				max_chars = len(str(cell))
	row_divider = "  +" + "+".join(["-"*max_chars for _ in range(len(matrix[0]))]) + "+"
	seq_delim = " "*max_chars
	top_line = seq_delim.join([" "]*2 + [i for i in seqb])
	to_print = [top_line]
	to_print.append(row_divider)
	for i, row in enumerate(matrix):
		if i == 0:
			seq = " "
		else:
			seq = seqa[i-1]
		row_contents = [f"{seq} |"]
		for j, cell in enumerate(row):
			contents_str = str(char_dict.get(cell, cell))
			cell_string = " " * (max_chars-len(contents_str)) + contents_str
			row_contents.append(f"{cell_string}|")
		to_print.append("".join(row_contents))
		to_print.append(row_divider)
	print("\n".join(to_print))


def print_alns(alignments):
	out_lines = [""]
	for n, a in enumerate(alignments):
		out_lines.append(f"Alignment: {n+1}")
		out_lines += _format_print(a)
	print("\n".join(out_lines))


def _format_print(alignments: tuple[str]):
	print_len = 30 # If you want to use this in a terminal application, use shutil.get_terminal_size().columns - 2 
	if len(alignments[0]) > print_len:
		prefixes = ["1 ", "  ", "2 "]
	else:
		prefixes = [""]*3
	aln_list = []
	for a,b in zip(*alignments):
		if a == b:
			aln_list.append("|")
		else:
			aln_list.append(" ")
	aln_string = "".join(aln_list)
	out_lines = []
	for i in range(0, len(aln_string), print_len):
		out_lines.append(prefixes[0] + alignments[0][i:i+print_len])
		out_lines.append(prefixes[1] + aln_string[i:i+print_len])
		out_lines.append(prefixes[2] + alignments[1][i:i+print_len])
		out_lines.append("")
	return out_lines