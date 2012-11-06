#include <R.h>
#include <Rdefines.h>    
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define LINELENGTH 400

SEXP writeCufflinksGTF(SEXP filename, SEXP chr, SEXP program, SEXP feature,
		SEXP start, SEXP end, SEXP score, SEXP strand, SEXP frame, SEXP gene_id,
		SEXP transcript_id, SEXP exon_number, SEXP fpkm, SEXP frac,
		SEXP conf_lo, SEXP conf_hi, SEXP cov, SEXP frs, SEXP length) {

	FILE *gtf = fopen(CHAR(STRING_ELT(filename, 0)), "w");

	int *p_start, *p_end, *p_score, *p_exon_number, mylength;
	double *p_fpkm, *p_frac, *p_conf_lo, *p_conf_hi, *p_cov;

	PROTECT(start = AS_INTEGER(start));
	p_start = INTEGER_POINTER(start);

	PROTECT(end = AS_INTEGER(end));
	p_end = INTEGER_POINTER(end);

	PROTECT(score = AS_INTEGER(score));
	p_score = INTEGER_POINTER(score);

	PROTECT(exon_number = AS_INTEGER(exon_number));
	p_exon_number = INTEGER_POINTER(exon_number);

	PROTECT(fpkm = AS_NUMERIC(fpkm));
	p_fpkm = NUMERIC_POINTER(fpkm);

	PROTECT(frac = AS_NUMERIC(frac));
	p_frac = NUMERIC_POINTER(frac);

	PROTECT(conf_lo = AS_NUMERIC(conf_lo));
	p_conf_lo = NUMERIC_POINTER(conf_lo);

	PROTECT(conf_hi = AS_NUMERIC(conf_hi));
	p_conf_hi = NUMERIC_POINTER(conf_hi);

	PROTECT(cov = AS_NUMERIC(cov));
	p_cov = NUMERIC_POINTER(cov);

	PROTECT(length = AS_INTEGER(length));
	mylength = INTEGER_POINTER(length)[0];

	int i = 0;
	while (i < mylength) {

		if (strcmp(CHAR(STRING_ELT(feature, i)), "transcript") == 0) {
			fprintf(gtf,
					"%s\t%s\t%s\t%i\t%i\t%i\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; FPKM \"%f\"; "
							"frac \"%f\"; conf_lo \"%f\"; conf_hi \"%f\"; cov \"%f\"; full_read_support \"\" ;\n",
					CHAR(STRING_ELT(chr, i)), CHAR(STRING_ELT(program, i)),
					CHAR(STRING_ELT(feature, i)), p_start[i], p_end[i],
					p_score[i], CHAR(STRING_ELT(strand, i)),
					CHAR(STRING_ELT(frame, i)), CHAR(STRING_ELT(gene_id, i)),
					CHAR(STRING_ELT(transcript_id, i)), p_fpkm[i], p_frac[i],
					p_conf_lo[i], p_conf_hi[i], p_cov[i]);
		} else {
			fprintf(gtf,
					"%s\t%s\t%s\t%i\t%i\t%i\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%i\";"
							" FPKM \"%f\"; frac \"%f\"; conf_lo \"%f\"; conf_hi \"%f\"; cov \"%f\";\n",
					CHAR(STRING_ELT(chr, i)), CHAR(STRING_ELT(program, i)),
					CHAR(STRING_ELT(feature, i)), p_start[i], p_start[i],
					p_score[i], CHAR(STRING_ELT(strand, i)),
					CHAR(STRING_ELT(frame, i)), CHAR(STRING_ELT(gene_id, i)),
					CHAR(STRING_ELT(transcript_id, i)), p_exon_number[i],
					p_fpkm[i], p_frac[i], p_conf_lo[i], p_conf_hi[i], p_cov[i]);

		}
		i++;
	}

	fclose(gtf);

	UNPROTECT(10);
	return (R_NilValue);
}

SEXP writeCuffCompareGTF(SEXP filename, SEXP chr, SEXP program, SEXP feature,
		SEXP start, SEXP end, SEXP score, SEXP strand, SEXP frame, SEXP gene_id,
		SEXP transcript_id, SEXP exon_number, SEXP gene_name, SEXP oId,
		SEXP nearest_ref, SEXP class_code, SEXP tss_id, SEXP length) {

	FILE *gtf = fopen(CHAR(STRING_ELT(filename, 0)), "w");

	int *p_start, *p_end, *p_exon_number, mylength;

	PROTECT(start = AS_INTEGER(start));
	p_start = INTEGER_POINTER(start);

	PROTECT(end = AS_INTEGER(end));
	p_end = INTEGER_POINTER(end);

	PROTECT(exon_number = AS_INTEGER(exon_number));
	p_exon_number = INTEGER_POINTER(exon_number);

	PROTECT(length = AS_INTEGER(length));
	mylength = INTEGER_POINTER(length)[0];

	int i = 0;
	while (i < 50) {
		fprintf(gtf,
				"%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%i\"; "
						"gene_name \"%s\"; oId \"%s\"; nearest_ref \"%s\"; class_code \"%s\"; tss_id \"%s\" ;\n",
				CHAR(STRING_ELT(chr, i)), CHAR(STRING_ELT(program, i)),
				CHAR(STRING_ELT(feature, i)), p_start[i], p_end[i],
				CHAR(STRING_ELT(score, i)), CHAR(STRING_ELT(strand, i)),
				CHAR(STRING_ELT(frame, i)), CHAR(STRING_ELT(gene_id, i)),
				CHAR(STRING_ELT(transcript_id, i)), p_exon_number[i],
				CHAR(STRING_ELT(gene_name, i)), CHAR(STRING_ELT(oId, i)),
				CHAR(STRING_ELT(nearest_ref, i)),
				CHAR(STRING_ELT(class_code, i)), CHAR(STRING_ELT(tss_id, i)));
		i++;
	}
	fclose(gtf);

	UNPROTECT(4);
	return (R_NilValue);
}

SEXP writeCuffMergeGTF(SEXP filename, SEXP chr, SEXP program, SEXP feature,
		SEXP start, SEXP end, SEXP score, SEXP strand, SEXP frame, SEXP gene_id,
		SEXP transcript_id, SEXP exon_number, SEXP gene_name, SEXP oId,
		SEXP contained_in, SEXP nearest_ref, SEXP class_code, SEXP tss_id,
		SEXP length) {

	FILE *gtf = fopen(CHAR(STRING_ELT(filename, 0)), "w");

	int *p_start, *p_end, *p_exon_number, *p_length, mylength;

	PROTECT(start = AS_INTEGER(start));
	p_start = INTEGER_POINTER(start);

	PROTECT(end = AS_INTEGER(end));
	p_end = INTEGER_POINTER(end);

	PROTECT(exon_number = AS_INTEGER(exon_number));
	p_exon_number = INTEGER_POINTER(exon_number);

	PROTECT(length = AS_INTEGER(length));
	mylength = INTEGER_POINTER(length)[0];

	int i = 0;
	while (i < mylength) {
		if (strcmp(CHAR(STRING_ELT(contained_in, i)), "") != 0) {
			fprintf(gtf,
					"%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%i\"; "
							"gene_name \"%s\"; oId \"%s\"; contained_in \"%s\"; nearest_ref \"%s\"; class_code \"%s\"; tss_id \"%s\";\n",
					CHAR(STRING_ELT(chr, i)), CHAR(STRING_ELT(program, i)),
					CHAR(STRING_ELT(feature, i)), p_start[i], p_end[i],
					CHAR(STRING_ELT(score, i)), CHAR(STRING_ELT(strand, i)),
					CHAR(STRING_ELT(frame, i)), CHAR(STRING_ELT(gene_id, i)),
					CHAR(STRING_ELT(transcript_id, i)), p_exon_number[i],
					CHAR(STRING_ELT(gene_name, i)), CHAR(STRING_ELT(oId, i)),
					CHAR(STRING_ELT(contained_in, i)),
					CHAR(STRING_ELT(nearest_ref, i)),
					CHAR(STRING_ELT(class_code, i)),
					CHAR(STRING_ELT(tss_id, i)));
		} else {
			fprintf(gtf,
					"%s\t%s\t%s\t%i\t%i\t%s\t%s\t%s\tgene_id \"%s\"; transcript_id \"%s\"; exon_number \"%i\"; "
							"gene_name \"%s\"; oId \"%s\";  nearest_ref \"%s\"; class_code \"%s\"; tss_id \"%s\";\n",
					CHAR(STRING_ELT(chr, i)), CHAR(STRING_ELT(program, i)),
					CHAR(STRING_ELT(feature, i)), p_start[i], p_end[i],
					CHAR(STRING_ELT(score, i)), CHAR(STRING_ELT(strand, i)),
					CHAR(STRING_ELT(frame, i)), CHAR(STRING_ELT(gene_id, i)),
					CHAR(STRING_ELT(transcript_id, i)), p_exon_number[i],
					CHAR(STRING_ELT(gene_name, i)), CHAR(STRING_ELT(oId, i)),
					CHAR(STRING_ELT(nearest_ref, i)),
					CHAR(STRING_ELT(class_code, i)),
					CHAR(STRING_ELT(tss_id, i)));
		}
		i++;
	}
	fclose(gtf);
	UNPROTECT(4);
	return (R_NilValue);
}