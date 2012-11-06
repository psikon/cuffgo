#include <R.h>
#include <Rdefines.h>    
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define LINELENGTH 350

SEXP parseCuffmerge(SEXP filename) {
	//FILE *cuffmerge = fopen(CHAR(STRING_ELT(filename, 0)), "r");
	FILE *cuffmerge = fopen(CHAR(STRING_ELT(filename, 0)), "r");

	char line[LINELENGTH];
	char endLine;
	int numLines = 0;
	char *token[27];
	int row = 0;
	int col;
	int *p_start;
	int *p_end;
	int *p_exon_number;
	int status = 0;

	// count lines for init
	if (cuffmerge) {
		while ((endLine = fgetc(cuffmerge)) != EOF) {
			if (endLine == '\n') {
				numLines++;
			}
		}
		if (numLines == 0) {
			Rprintf("ERROR reading file is empty : %s",
					CHAR(STRING_ELT(filename, 0)));
		}
		fclose(cuffmerge);
	} else {
		Rprintf("ERROR reading file: %s", CHAR(STRING_ELT(filename, 0)));
	}
	// init vectors for R
	SEXP chr, program, feature, start, end, score, strand, frame, gene_id,
			transcript_id, exon_number, gene_name, oId, contained_in,
			nearest_ref, class_code, tss_id, list, list_names;

	PROTECT(chr = NEW_CHARACTER(numLines));

	PROTECT(program = NEW_CHARACTER(numLines));

	PROTECT(feature = NEW_CHARACTER(numLines));

	PROTECT(start = NEW_INTEGER(numLines));
	p_start = INTEGER_POINTER(start);

	PROTECT(end = NEW_INTEGER(numLines));
	p_end = INTEGER_POINTER(end);

	PROTECT(score = NEW_CHARACTER(numLines));

	PROTECT(strand = NEW_CHARACTER(numLines));

	PROTECT(frame = NEW_CHARACTER(numLines));

	PROTECT(gene_id = NEW_CHARACTER(numLines));

	PROTECT(transcript_id = NEW_CHARACTER(numLines));

	PROTECT(exon_number = NEW_INTEGER(numLines));
	p_exon_number = INTEGER_POINTER(exon_number);

	PROTECT(gene_name = NEW_CHARACTER(numLines));

	PROTECT(oId = NEW_CHARACTER(numLines));

	PROTECT(contained_in = NEW_CHARACTER(numLines));

	PROTECT(nearest_ref = NEW_CHARACTER(numLines));

	PROTECT(class_code = NEW_CHARACTER(numLines));

	PROTECT(tss_id = NEW_CHARACTER(numLines));

	// names of the vectors
	char *names[17] = { "chr", "program", "feature", "start", "end", "score",
			"strand", "frame", "gene_id", "transcript_id", "exon_number",
			"gene_name", "oId", "contained_in", "nearest_ref", "class_code",
			"tss_id" };
	// open the file for reading only
	cuffmerge = fopen(CHAR(STRING_ELT(filename, 0)), "r");
	// read the line until the end of file
	while (fgets(line, sizeof line, cuffmerge)) {
		char status = 0; //0 - transcript, 1 - exon
		int pos = 0;
		col = 0;
		token[0] = strtok(line, "\t\"");
		while (token[col]) {

			// assigment of the values to separate vectors
			switch (col) {
			case 1: // case chr
				SET_STRING_ELT(chr, row, mkChar(token[0]));
				break;
			case 2: // case program
				SET_STRING_ELT(program, row, mkChar(token[1]));
				break;
			case 3: // case feature
				SET_STRING_ELT(feature, row, mkChar(token[2]));
				break;
			case 4: // case start
				p_start[row] = atof(token[3]);
				break;
			case 5: // case end
				p_end[row] = atoi(token[4]);
				break;
			case 6: // case score
				SET_STRING_ELT(score, row, mkChar(token[5]));
				break;
			case 7: // case strand
				SET_STRING_ELT(strand, row, mkChar(token[6]));
				break;
			case 8: // case frame
				SET_STRING_ELT(frame, row, mkChar(token[7]));
				break;
			case 9: // case gene_id annotation
				if (strstr(token[8], "gene_id") != 0) {
					SET_STRING_ELT(gene_id, row, mkChar(token[9]));
				} else {
					SET_STRING_ELT(gene_id, row, mkChar(""));
				}
				break;
			case 11: // case transcript_id annotation
				if (strstr(token[10], "transcript_id") != 0) {
					SET_STRING_ELT(transcript_id, row, mkChar(token[11]));
				} else {
					SET_STRING_ELT(transcript_id, row, mkChar(""));
				}
				break;
			case 13: // case exon number
				if (strstr(token[12], "exon_number") != 0) {
					p_exon_number[row] = atoi(token[13]);
				} else {
					p_exon_number[row] = 0;
				}
				break;
			case 15: //case gene name
				if (strstr(token[14], "gene_name") != 0) {
					SET_STRING_ELT(gene_name, row, mkChar(token[15]));
				} else {
					SET_STRING_ELT(gene_name, row, mkChar(""));
				}
				break;
			case 17: //case oID
				if (strstr(token[16], "oId") != 0) {
					SET_STRING_ELT(oId, row, mkChar(token[17]));
				} else {
					SET_STRING_ELT(oId, row, mkChar(""));
				}
				break;
			case 19: // case possible extra row contained_in
				if (strstr(token[18], "contained_in") != 0) {
					SET_STRING_ELT(contained_in, row, mkChar(token[19]));
					status = 1;
				} else {
					SET_STRING_ELT(contained_in, row, mkChar(""));
				}
				break;
			}
			if (status == 0) {
				switch (col) {
				case 20:
					if (strstr(token[18], "nearest_ref") != 0) {
						SET_STRING_ELT(nearest_ref, row, mkChar(token[19]));
					} else {
						SET_STRING_ELT(nearest_ref, row, mkChar(""));
					}
					break;
				case 21:
					if (strstr(token[20], "class_code") != 0) {
						SET_STRING_ELT(class_code, row, mkChar(token[21]));
					} else {
						SET_STRING_ELT(class_code, row, mkChar(""));
					}
					break;
				case 23:
					if (strstr(token[22], "tss_id") != 0) {
						SET_STRING_ELT(tss_id, row, mkChar(token[23]));
					} else {
						SET_STRING_ELT(tss_id, row, mkChar(""));
					}
					break;
				}
			}
			if (status == 1) {
				switch (col) {
				case 21:
					if (strstr(token[20], "nearest_ref") != 0) {
						SET_STRING_ELT(nearest_ref, row, mkChar(token[21]));
					} else {
						SET_STRING_ELT(nearest_ref, row, mkChar(""));
					}
					break;
				case 23:
					if (strstr(token[22], "class_code") != 0) {
						SET_STRING_ELT(class_code, row, mkChar(token[23]));
					} else {
						SET_STRING_ELT(class_code, row, mkChar(""));
					}
					break;
				case 25:
					if (strstr(token[24], "tss_id") != 0) {
						SET_STRING_ELT(tss_id, row, mkChar(token[25]));
					} else {
						SET_STRING_ELT(tss_id, row, mkChar(""));
					}
					break;
				}
			}

			//only for debugging
			//Rprintf("%d  %s \n", col, token[col]);
			col++; // for the switch case handling
			token[col] = strtok(NULL, "\t\""); // reset token
		}
		row++; //get next line
	}

// init the names of the vectors
	PROTECT(list_names = allocVector(STRSXP, 17));
// label the vectors with names (later used by r for reference)
	for (int i = 0; i < 17; i++) {
		SET_STRING_ELT(list_names, i, mkChar(names[i]));
	}
// build a list of the vectors
	PROTECT(list = allocVector(VECSXP, 17));
	SET_VECTOR_ELT(list, 0, chr);
	SET_VECTOR_ELT(list, 1, program);
	SET_VECTOR_ELT(list, 2, feature);
	SET_VECTOR_ELT(list, 3, start);
	SET_VECTOR_ELT(list, 4, end);
	SET_VECTOR_ELT(list, 5, score);
	SET_VECTOR_ELT(list, 6, strand);
	SET_VECTOR_ELT(list, 7, frame);
	SET_VECTOR_ELT(list, 8, gene_id);
	SET_VECTOR_ELT(list, 9, transcript_id);
	SET_VECTOR_ELT(list, 10, exon_number);
	SET_VECTOR_ELT(list, 11, gene_name);
	SET_VECTOR_ELT(list, 12, oId);
	SET_VECTOR_ELT(list, 13, contained_in);
	SET_VECTOR_ELT(list, 14, nearest_ref);
	SET_VECTOR_ELT(list, 15, class_code);
	SET_VECTOR_ELT(list, 16, tss_id);
	setAttrib(list, R_NamesSymbol, list_names);
	UNPROTECT(19);
// close the filestream
	fclose(cuffmerge);
	return list;
}

