#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define LINELENGTH 300

SEXP parseCufflinks(SEXP filename) {
	
	FILE *cufflinks = fopen(CHAR(STRING_ELT(filename, 0)), "r");

	char line[LINELENGTH];
	char endLine;
	int numLines = 0;
	char *token[25];
	int row = 0;
	int col;
	int *p_start;
	int *p_end;
	double *p_score;
	int *p_exon_number;
	double *p_fpkm;
	double *p_frac;
	double *p_conf_lo;
	double *p_conf_hi;
	double *p_cov;

	// count lines for init
	if (cufflinks) {
		while ((endLine = fgetc(cufflinks)) != EOF) {
			if (endLine == '\n') {
				numLines++;
			}
		}
		if (numLines == 0) {
			Rprintf("ERROR reading file is empty : %s",
					CHAR(STRING_ELT(filename, 0)));
		}
		fclose(cufflinks);
	} else {
		Rprintf("ERROR reading file: %s", CHAR(STRING_ELT(filename, 0)));
	}
	// init vectors for R
	SEXP chr, program, feature, start, end, score, strand, frame, gene_id,
			transcript_id, exon_number, fpkm, frac, conf_lo, conf_hi, cov, frs,
			list, list_names;

	PROTECT(chr = NEW_CHARACTER(numLines));

	PROTECT(program = NEW_CHARACTER(numLines));

	PROTECT(feature = NEW_CHARACTER(numLines));

	PROTECT(start = NEW_INTEGER(numLines));
	p_start = INTEGER_POINTER(start);

	PROTECT(end = NEW_INTEGER(numLines));
	p_end = INTEGER_POINTER(end);

	PROTECT(score = NEW_NUMERIC(numLines));
	p_score = NUMERIC_POINTER(score);

	PROTECT(strand = NEW_CHARACTER(numLines));

	PROTECT(frame = NEW_CHARACTER(numLines));

	PROTECT(gene_id = NEW_CHARACTER(numLines));

	PROTECT(transcript_id = NEW_CHARACTER(numLines));

	PROTECT(exon_number = NEW_INTEGER(numLines));
	p_exon_number = INTEGER_POINTER(exon_number);

	PROTECT(fpkm = NEW_NUMERIC(numLines));
	p_fpkm = NUMERIC_POINTER(fpkm);

	PROTECT(frac = NEW_NUMERIC(numLines));
	p_frac = NUMERIC_POINTER(frac);

	PROTECT(conf_lo = NEW_NUMERIC(numLines));
	p_conf_lo = NUMERIC_POINTER(conf_lo);

	PROTECT(conf_hi = NEW_NUMERIC(numLines));
	p_conf_hi = NUMERIC_POINTER(conf_hi);

	PROTECT(cov = NEW_NUMERIC(numLines));
	p_cov = NUMERIC_POINTER(cov);

	PROTECT(frs = NEW_CHARACTER(numLines));
	// names of the vectors
	char *names[17] = { "chr", "program", "feature", "start", "end", "score",
			"strand", "frame", "gene_id", "transcript_id", "exon_number",
			"fpkm", "frac", "conf_lo", "conf_hi", "cov", "frs" };
	// open the file for reading only
	cufflinks = fopen(CHAR(STRING_ELT(filename, 0)), "r");
	// read the line until the end of file
	while (fgets(line, sizeof line, cufflinks)) {
		char status = 0; //0 - transcript, 1 - exon
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
				if (strlen(token[2]) == 4) {
					status = 1; // make decision exon or transcript
				}
				break;
			case 4: // case start
				p_start[row] = atof(token[3]);
				break;
			case 5: // case end
				p_end[row] = atoi(token[4]);
				break;
			case 6: // case score
				p_score[row] = atof(token[5]);
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
			}
			if (status == 0) { // read line is a transcript
				p_exon_number[row] = 0; // transcript has no exon number
				switch (col) {
				case 13: //case fpkm
					if (strstr(token[12], "FPKM") != 0) {
						p_fpkm[row] = atof(token[13]);
					} else {
						p_fpkm[row] = 0;
					}
					break;
				case 15: //case frac
					if (strstr(token[14], "frac") != 0) {
						p_frac[row] = atof(token[15]);
					} else {
						p_frac[row] = 0;
					}
					break;
				case 17: //case conf_lo
					if (strstr(token[16], "conf_lo") != 0) {
						p_conf_lo[row] = atof(token[17]);
					} else {
						p_conf_lo[row] = 0;
					}
					break;
				case 19: //case conf_hi
					if (strstr(token[18], "conf_hi") != 0) {
						p_conf_hi[row] = atof(token[19]);
					} else {
						p_conf_hi[row] = 0;
					}
					break;
				case 21: //case coverage
					if (strstr(token[20], "cov") != 0) {
						p_cov[row] = atof(token[21]);
					} else {
						p_cov[row] = 0;
					}
					break;
				case 23: //case full read support
					if (strstr(token[22], "full_read_support") != 0) {
						SET_STRING_ELT(frs, row, mkChar(token[23]));
					} else {
						SET_STRING_ELT(frs, row, mkChar(""));
					}
					break;
				}
			}

			if (status == 1) {
				switch (col) {
				case 13: //case exon number
					if (strstr(token[12], "exon_number") != 0) {
						p_exon_number[row] = atof(token[13]);
					} else {
						p_exon_number[row] = 0;
					}
					break;
				case 15: //case fpkm
					if (strstr(token[14], "FPKM") != 0) {
						p_fpkm[row] = atof(token[15]);
					} else {
						p_fpkm[row] = 0;
					}
					break;
				case 17: //case frac
					if (strstr(token[16], "frac") != 0) {
						p_frac[row] = atof(token[17]);
					} else {
						p_frac[row] = 0;
					}
					break;
				case 19: //case conf_lo
					if (strstr(token[18], "conf_lo") != 0) {
						p_conf_lo[row] = atof(token[19]);
					} else {
						p_conf_lo[row] = 0;
					}
					break;
				case 21: // case conf_hi
					if (strstr(token[20], "conf_hi") != 0) {
						p_conf_hi[row] = atof(token[21]);
					} else {
						p_conf_hi[row] = 0;
					}
					break;
				case 23: // case coverage
					if (strstr(token[22], "cov") != 0) {
						p_cov[row] = atof(token[23]);
					} else {
						p_cov[row] = 0;
					}
					break;
				}
				SET_STRING_ELT(frs, row, mkChar("")); // exon has no full read support
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
	SET_VECTOR_ELT(list, 11, fpkm);
	SET_VECTOR_ELT(list, 12, frac);
	SET_VECTOR_ELT(list, 13, conf_lo);
	SET_VECTOR_ELT(list, 14, conf_hi);
	SET_VECTOR_ELT(list, 15, cov);
	SET_VECTOR_ELT(list, 16, frs);

	setAttrib(list, R_NamesSymbol, list_names);
	UNPROTECT(19);
	// close the filestream
	fclose(cufflinks);
	return list;
}

