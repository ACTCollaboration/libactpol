/*                           (C) 2002 C. Barth Netterfield */
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <assert.h>
#include <ctype.h>
#include <libgen.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <slimlib.h>
#include <zzip/zzip.h>

#include "actpol/getdata.h"

/***************************************************************/
/*                                                             */
/*    Structures for describing dirfile formats                */
/*                                                             */
/***************************************************************/

#define FIELD_LENGTH 50
#define MAX_FILENAME_LENGTH 250
#define MAX_LINE_LENGTH 250
#define MAX_LINCOM 3
#define MAX_IN_COLS 15
#define MAX_POLYNOM 6

struct RawEntryType {
  char field[FIELD_LENGTH+1];
  char file[MAX_FILENAME_LENGTH + FIELD_LENGTH + 2];
  char type;
  int size;
  int samples_per_frame;
};

struct FileHandle {
  ZZIP_FILE *fp;
  SLIMFILE *slim;
};

struct PolynomEntryType {
  char field[FIELD_LENGTH+1];
  char in_field[FIELD_LENGTH+1];
  int n_coeff;
  double a[MAX_POLYNOM];
};

struct LincomEntryType {
  char field[FIELD_LENGTH+1];
  int n_infields;
  char in_fields[MAX_LINCOM][FIELD_LENGTH+1];
  double m[MAX_LINCOM];
  double b[MAX_LINCOM];
};

struct LinterpEntryType {
  char field[FIELD_LENGTH+1];
  char raw_field[FIELD_LENGTH+1];
  char linterp_file[MAX_FILENAME_LENGTH];
  int n_interp;
  double *x;
  double *y;
};

struct MultiplyEntryType {
  char field[FIELD_LENGTH+1];
  char in_fields[2][FIELD_LENGTH+1];
};

struct MplexEntryType {
  char field[FIELD_LENGTH+1];
  char cnt_field[FIELD_LENGTH+1];
  char data_field[FIELD_LENGTH+1];
  int i;
  int max_i;
};

struct BitEntryType {
  char field[FIELD_LENGTH+1];
  char raw_field[FIELD_LENGTH+1];
  int bitnum;
  int numbits;
};

struct FormatType {
  char FileDirName[MAX_FILENAME_LENGTH];
  int frame_offset;
  struct RawEntryType first_field;
  struct RawEntryType *rawEntries;
  int n_raw;
  struct PolynomEntryType *polynomEntries;
  int n_polynom;
  struct LincomEntryType *lincomEntries;
  int n_lincom;
  struct LinterpEntryType *linterpEntries;
  int n_linterp;
  struct MultiplyEntryType *multiplyEntries;
  int n_multiply;
  struct MplexEntryType *mplexEntries;
  int n_mplex;
  struct BitEntryType *bitEntries;
  int n_bit;
};

const int  MAX_GETDATA_FILES_OPEN=128;

const char *GD_ERROR_CODES[15] = {"OK",
                            "Could not open Format file",
                            "Error in Format file",
                            "Could not open Data file",
                            "Field name too long",
                            "Field code not found in File Format",
                            "Unrecognized return type",
                            "Could not open field file",
                            "Could not open included Format file",
                            " ",
                            " ",
                            "Could not open Field File",
                            "Size mismatch in linear combination",
                            "Could not open interpolation file ",
                            "Too many levels of recursion"
};

static int DoField(int recurse_level,
    const struct FormatType *F, const char *field_code,
    int first_frame, int first_samp,
    int num_frames, int num_samp,
    char return_type, void *data_out,
    int *error_code);

/*
struct lines {
  int num_lines;
  char **line;
  char *_buf;
};
*/

static char *
read_file( ZZIP_FILE *fp )
{
    ZZIP_STAT stat;
    char *buf;
    size_t len;

    zzip_fstat( fp, &stat );
    len = stat.st_size;

    buf = (char *) malloc( len+1 );
    zzip_read( fp, buf, len );
    buf[len] = '\0';

    return buf;
}

/*
struct line_list {
  char *str;
  struct line_list *next;
}

static struct line_list *
read_lines_from_buffer( void *buf, size_t len )
{
  struct line_list *lines, *line, *prev;
  char *tok;
  char *s = (char *) buf;
  char *t;

  lines = NULL;

  while ( (tok = strsep(&s,"\n")) != NULL )
  {
    while ( isspace(*tok) ) tok++;
    t = strchr( tok, '#' );
    if ( t != NULL ) *t = '\0';

    if ( strlen(tok) == 0 ) continue;

    line = (struct line_list *) malloc( sizeof(struct line_list) );
    line->str = tok;
    line->next = NULL;

    if ( lines == NULL )
      lines = line;
    else
      prev->next = line;

    prev = line;
  }

  return lines;
}
*/

static char **
read_lines_from_buffer( char *buf, size_t *nlines )
{
  char **lines;
  size_t n, max_nlines;
  const size_t chunk = 64;
  char *s, *t, *tok;

  max_nlines = chunk;
  lines = (char **) malloc( max_nlines*sizeof(char *) );
  n = 0;

  s = buf;
  while ( (tok = strsep(&s,"\n")) != NULL )
  {
    // strip leading whitespace
    while ( isspace(*tok) ) tok++;

    // strip comments
    t = strchr( tok, '#' );
    if ( t != NULL ) *t = '\0';

    if ( strlen(tok) == 0 ) continue;

    if ( n == max_nlines )
    {
      max_nlines += chunk;
      lines = (char **) realloc( lines, max_nlines*sizeof(char *) ); 
    }

    lines[n++] = tok;
  }

  if ( n < max_nlines )
    lines = (char **) realloc( lines, n*sizeof(char *) );

  *nlines = n;
  return lines;
}

static bool
file_exists( const char *filename )
{
  ZZIP_FILE *fp = zzip_fopen( filename, "r" );
  if ( fp == NULL ) return false;
  zzip_fclose( fp );
  return true;
}

static size_t
file_size( const char *filename )
{
  printf( "file_size: %s\n", filename );
    ZZIP_FILE *fp = zzip_fopen( filename, "r" );
    assert( fp != NULL ); 
    ZZIP_STAT stat;
    zzip_fstat( fp, &stat );
    zzip_fclose( fp );
    return stat.st_size;
}

/***************************************************************************/
/*                                                                         */
/*    GetLine: read non-comment line from format file                      */
/*        The line is placed in *line.                                     */
/*        Returns 1 if successful, 0 if unsuccessful                       */
/*                                                                         */
/***************************************************************************/
static int GetLine(FILE *fp, char *line) {
  char *ret_val;
  int first_char;
  int i, len;

  do {
    ret_val = fgets(line, MAX_LINE_LENGTH, fp);
    first_char = 0;
    while (line[first_char] == ' ' || line[first_char] == '\t') ++first_char;
    line += first_char;
  } while (ret_val && (line[0] == '#' || line[0] == 0 || line[1] == 0));


  if (ret_val) {
    /* truncate comments from end of lines */
    len = strlen(line);
    for (i = 0; i < len; i++) {
      if (line[i]=='#')
        line[i] = '\0';
    }

    return(1); /* a line was read */
  }
  return(0);  /* there were no valid lines */
}


/***************************************************************************/
/*                                                                         */
/*   FreeF: free any entries that have been allocated in F                 */
/*                                                                         */
/***************************************************************************/
static void FreeF(struct FormatType *F) {
  if (F->n_raw > 0) free(F->rawEntries);
  if (F->n_polynom > 0) free(F->polynomEntries);
  if (F->n_lincom > 0) free(F->lincomEntries);
  if (F->n_multiply > 0) free(F->multiplyEntries);
  if (F->n_linterp >0) free(F->linterpEntries);
  if (F->n_mplex > 0) free(F->mplexEntries);
  if (F->n_bit > 0) free(F->bitEntries);
}

/***************************************************************************/
/*                                                                         */
/*   ParseRaw: parse a RAW data type in the formats file                   */
/*                                                                         */
/***************************************************************************/
static void ParseRaw(char in_cols[MAX_IN_COLS][MAX_LINE_LENGTH],
    struct RawEntryType *R, const char* subdir, int *error_code) {
  strcpy(R->field, in_cols[0]); /* field */
  if (strcmp(subdir, ".") == 0)
  {
    snprintf(R->file, MAX_FILENAME_LENGTH + FIELD_LENGTH + 2, "%s",
      in_cols[0]); /* filename */
  }
  else
  {
    snprintf(R->file, MAX_FILENAME_LENGTH + FIELD_LENGTH + 2, "%s/%s", subdir,
      in_cols[0]); /* path and filename */
  }
  //R->fp = -1; /* file not opened yet */
  //R->slim = NULL; /* No SLIMFILE either */
  switch (in_cols[2][0]) {
    case 'c':
      R->size = 1;
      break;
    case 's': case 'u':
      R->size = 2;
      break;
    case 'S': case 'U': case 'f': case 'i':
      R->size = 4;
      break;
    case 'd':
      R->size = 8;
      break;
    default:
      printf("bad raw %c\n", in_cols[2][0]);
      *error_code = GD_E_FORMAT;
      return;
  }
  R->type = in_cols[2][0];
  R->samples_per_frame = atoi(in_cols[3]);
  if (R->samples_per_frame<=0) {
    *error_code = GD_E_FORMAT;
    return;
  }
}

/***************************************************************************/
/*                                                                         */
/*  ParsePolynom: parse a POLYNOM data type in the formats file              */
/*                                                                         */
/***************************************************************************/
static void ParsePolynom(char in_cols[MAX_IN_COLS][MAX_LINE_LENGTH],
    int n_cols, struct PolynomEntryType *P, int *error_code) {
  int i;
  P->n_coeff = n_cols - 3;
  if ((P->n_coeff < 2) || (P->n_coeff > MAX_POLYNOM)) {
    *error_code = GD_E_FORMAT;
    return;
  }
  strcpy(P->field, in_cols[0]); /* field */
  strncpy(P->in_field, in_cols[2], FIELD_LENGTH);
  for (i=0; i<MAX_POLYNOM; i++) {
    P->a[i] = (i<P->n_coeff) ? atof(in_cols[i+3]) : 0.;
  }
}

/***************************************************************************/
/*                                                                         */
/*  ParseLincom: parse a LINCOM data type in the formats file              */
/*                                                                         */
/***************************************************************************/
static void ParseLincom(char in_cols[MAX_IN_COLS][MAX_LINE_LENGTH],
    int n_cols, struct LincomEntryType *L, int *error_code) {
  int i;
  strcpy(L->field, in_cols[0]); /* field */
  L->n_infields = atoi(in_cols[2]);
  if ((L->n_infields<1) || (L->n_infields>MAX_LINCOM) ||
      (n_cols < L->n_infields*3 + 3)) {
    *error_code = GD_E_FORMAT;
    return;
  }
  for (i=0; i<L->n_infields; i++) {
    strncpy(L->in_fields[i], in_cols[i*3+3], FIELD_LENGTH);
    L->m[i] = atof(in_cols[i*3+4]);
    L->b[i] = atof(in_cols[i*3+5]);
  }
}
/***************************************************************************/
/*                                                                         */
/*  ParseLinterp: parse a LINTERP data type in the formats file            */
/*                                                                         */
/***************************************************************************/
static void ParseLinterp(char in_cols[MAX_IN_COLS][MAX_LINE_LENGTH],
                         const char* linterp_prefix,
                         struct LinterpEntryType *L) {
  char full_name[MAX_LINE_LENGTH];
  char *tokptr, *last_tokptr; 

  strcpy(L->field, in_cols[0]); /* field */
  strncpy(L->raw_field, in_cols[2], FIELD_LENGTH);
  L->n_interp = -1; /* linterp file not read yet */

  if(linterp_prefix != NULL) {
    strcpy(full_name, in_cols[3]);
    strcpy(L->linterp_file, linterp_prefix);

    // tokenize on / until file is reached 
    tokptr = strtok (full_name,"//");
    while ( tokptr != NULL) {
      last_tokptr = tokptr;
      tokptr=strtok(NULL, "//");
    } 
    strcat(L->linterp_file,last_tokptr);

  } else {
    // copy the name directly over
    strcpy(L->linterp_file, in_cols[3]);
  }
}

/***************************************************************************/
/*                                                                         */
/*   ParseMultiply: parse MULTIPLY data type entry in formats file         */
/*                                                                         */
/***************************************************************************/
static void ParseMultiply(char in_cols[MAX_IN_COLS][MAX_LINE_LENGTH],
    int n_cols, struct MultiplyEntryType *M,
    int *error_code) {

  if (n_cols<4) {
    *error_code = GD_E_FORMAT;
    return;
  }

  strcpy(M->field, in_cols[0]); /* field */

  strncpy(M->in_fields[0], in_cols[2], FIELD_LENGTH);
  strncpy(M->in_fields[1], in_cols[3], FIELD_LENGTH);
}

/***************************************************************************/
/*                                                                         */
/*   ParseMplex: parse MPLEX data type entry in formats file               */
/*                                                                         */
/***************************************************************************/
static void ParseMplex(char in_cols[MAX_IN_COLS][MAX_LINE_LENGTH],
    int n_cols, struct MplexEntryType *M,
    int *error_code) {

  if (n_cols<6) {
    *error_code = GD_E_FORMAT;
    return;
  }

  strcpy(M->field, in_cols[0]); /* field */
  strncpy(M->cnt_field, in_cols[2], FIELD_LENGTH);
  strncpy(M->data_field, in_cols[3], FIELD_LENGTH);
  M->i = atoi(in_cols[4]);
  M->max_i = atoi(in_cols[5]);
  if ((M->max_i<1) || (M->max_i < M->i)) {
    *error_code = GD_E_FORMAT;
    return;
  }
}

/***************************************************************************/
/*                                                                         */
/*   ParseBit: parse BIT data type entry in formats file                   */
/*                                                                         */
/***************************************************************************/
static void ParseBit(char in_cols[MAX_IN_COLS][MAX_LINE_LENGTH],
    int n_cols,
    struct BitEntryType *B,
    int *error_code) {
  int error = 0;

  strcpy(B->field, in_cols[0]); /* field */
  strncpy(B->raw_field, in_cols[2], FIELD_LENGTH); /* field */

  B->bitnum=atoi(in_cols[3]);
  if (n_cols>4) {
    B->numbits=atoi(in_cols[4]);
  } else {
    B->numbits=1;
  }

  if (B->numbits<1) error = 1;
  if (B->bitnum<0) error = 1;
  if (B->bitnum + B->numbits - 1 > 31) error = 1;

  if (error) {
    *error_code=GD_E_FORMAT;
    return;
  }
}

/***************************************************************************/
/*                                                                         */
/*   Compare functions for sorting the lists (using stdlib qsort)          */
/*                                                                         */
/***************************************************************************/
static int RawCmp(const void *A, const void *B) {
  return (strcmp(((struct RawEntryType *)A)->field,
        ((struct RawEntryType *)B)->field));
}

static int PolynomCmp(const void *A, const void *B) {
  return (strcmp(((struct PolynomEntryType *)A)->field,
        ((struct PolynomEntryType *)B)->field));
}

static int LincomCmp(const void *A, const void *B) {
  return (strcmp(((struct LincomEntryType *)A)->field,
        ((struct LincomEntryType *)B)->field));
}

static int LinterpCmp(const void *A, const void *B) {
  return (strcmp(((struct LinterpEntryType *)A)->field,
        ((struct LinterpEntryType *)B)->field));
}

static int MultiplyCmp(const void *A, const void *B) {
  return (strcmp(((struct MultiplyEntryType *)A)->field,
        ((struct MultiplyEntryType *)B)->field));
}

static int MplexCmp(const void *A, const void *B) {
  return (strcmp(((struct MplexEntryType *)A)->field,
        ((struct MplexEntryType *)B)->field));
}

static int BitCmp(const void *A, const void *B) {
  return (strcmp(((struct BitEntryType *)A)->field,
        ((struct BitEntryType *)B)->field));
}

/***************************************************************************/
/*                                                                         */
/*  ParseFormatFile: Perform the actual parsing of the format file.  This  */
/*     function is called from GetFormat once for the main format file and */
/*     once for each included file.                                        */
/*                                                                         */
/***************************************************************************/
static int ParseFormatFile(ZZIP_FILE* fp, struct FormatType *F, const char* filedir,
    const char* subdir, const char* linterp_prefix, char*** IncludeList, int *i_include)
{
  //char instring[MAX_LINE_LENGTH];
  char in_cols[MAX_IN_COLS][MAX_LINE_LENGTH];
  int n_cols, error_code = GD_E_OK;

  size_t i, nlines;
  char *buf = read_file( fp );
  char **lines = read_lines_from_buffer( buf, &nlines );

  /***** start parsing ****/
  for ( i = 0; i < nlines; i++ )
  {
    //printf( "%s\n", lines[i] );
    /* ok, brute force parse...  slow and ugly but convenient... */
    n_cols = sscanf(lines[i], "%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s",
        in_cols[0], in_cols[1], in_cols[2], in_cols[3],
        in_cols[4], in_cols[5], in_cols[6], in_cols[7],
        in_cols[8], in_cols[9], in_cols[10], in_cols[11],
        in_cols[12], in_cols[13], in_cols[14]);

    //printf( "%s: %d\n", lines[i], n_cols );
    if (n_cols<2) {
      error_code = GD_E_FORMAT;
    } else if (strlen(in_cols[0])>FIELD_LENGTH) {
      error_code = GD_E_FIELD;
    } else if (strcmp(in_cols[1], "RAW")==0) {
      F->n_raw++;
      F->rawEntries =
        realloc(F->rawEntries, F->n_raw*sizeof(struct RawEntryType));
      ParseRaw(in_cols, F->rawEntries+F->n_raw - 1, subdir, &error_code);
    } else if (strcmp(in_cols[1], "POLYNOM")==0) {
      F->n_polynom++;
      F->polynomEntries =
        realloc(F->polynomEntries,
            F->n_polynom*sizeof(struct PolynomEntryType));
      ParsePolynom(in_cols, n_cols, F->polynomEntries+F->n_polynom - 1,
          &error_code);
    } else if (strcmp(in_cols[1], "LINCOM")==0) {
      F->n_lincom++;
      F->lincomEntries =
        realloc(F->lincomEntries,
            F->n_lincom*sizeof(struct LincomEntryType));
      ParseLincom(in_cols, n_cols, F->lincomEntries+F->n_lincom - 1,
          &error_code);
    } else if (strcmp(in_cols[1], "LINTERP")==0) {
      F->n_linterp++;
      F->linterpEntries =
        realloc(F->linterpEntries,
            F->n_linterp*sizeof(struct LinterpEntryType));
      ParseLinterp(in_cols, linterp_prefix, F->linterpEntries+F->n_linterp - 1);
    } else if (strcmp(in_cols[1], "MULTIPLY")==0) {
      F->n_multiply++;
      F->multiplyEntries =
        realloc(F->multiplyEntries,
            F->n_multiply*sizeof(struct MultiplyEntryType));
      ParseMultiply(in_cols, n_cols, F->multiplyEntries+F->n_multiply - 1,
          &error_code);
    } else if (strcmp(in_cols[1], "MPLEX")==0) {
      F->n_mplex++;
      F->mplexEntries =
        realloc(F->mplexEntries,
            F->n_mplex*sizeof(struct MplexEntryType));
      ParseMplex(in_cols, n_cols, F->mplexEntries+F->n_mplex - 1,
          &error_code);
    } else if (strcmp(in_cols[1], "BIT")==0) {
      F->n_bit++;
      F->bitEntries =
        realloc(F->bitEntries,
            F->n_bit*sizeof(struct BitEntryType));
      ParseBit(in_cols, n_cols, F->bitEntries+F->n_bit - 1,
          &error_code);
    } else if (strcmp(in_cols[0], "FRAMEOFFSET")==0) {
      F->frame_offset = atoi(in_cols[1]);
    } else if (strcmp(in_cols[0], "INCLUDE")==0) {
      int i, found = 0;
      char format_file[MAX_FILENAME_LENGTH + 6];
      char new_subdir[MAX_FILENAME_LENGTH + 1];
      ZZIP_FILE* new_fp = NULL;

      /* Run through the include list to see if we've already included this
       * file */
      for (i = 0; i < *i_include; ++i)
        if (strcmp(in_cols[1], (*IncludeList)[i]) == 0) {
          found = 1;
          break;
        }

      /* If we found the file, we won't reopen it.  Conitnue parsing. */
      if (found)
        continue;

      /* Otherwise, try to open the file */
      snprintf(format_file, MAX_FILENAME_LENGTH + 6, "%s/%s/%s", filedir,
          subdir, in_cols[1]);
      new_fp = zzip_fopen(format_file, "r");

      /* If opening the file failed, set the error code and abort parsing. */
      if (new_fp == NULL) {
	printf("%s\n", format_file);
        error_code = GD_E_OPEN_INCLUDE;
        break;
      }

      /* If we got here, we managed to open the inlcuded file; parse it */
      *IncludeList = realloc(*IncludeList, ++(*i_include) * sizeof(char*));
      (*IncludeList)[*i_include - 1] = strdup(in_cols[1]);

      /* extract the subdirectory name - dirname both returns a volatile string
       * and modifies its argument, ergo strcpy */
      strcpy(format_file, in_cols[1]);
      if (strcmp(subdir, ".") == 0)
        strcpy(new_subdir, dirname(format_file));
      else
        snprintf(new_subdir, MAX_FILENAME_LENGTH, "%s/%s", subdir,
            dirname(format_file));

      error_code = ParseFormatFile(new_fp, F, filedir, new_subdir, linterp_prefix,
          IncludeList, i_include);
      zzip_fclose(new_fp);
    } else {
      printf("bad format a: %s %s\n", in_cols[0], in_cols[1]);
      error_code = GD_E_FORMAT;
    }

    /* break out of loop (so we can return) if we've encountered an error */
    if (error_code != GD_E_OK)
      break;
  }

  free( lines );
  free( buf );

  return error_code;
}

/***************************************************************************/
/*                                                                         */
/*  GetDataClose: close all open file handles.                             */
/*                                                                         */
/***************************************************************************/

void GetDataClose(struct FormatType *F) {
  FreeF(F);
}

/***************************************************************************/
/*                                                                         */
/*   GetFormat: Read format file and fill structure.  The format           */
/*      is cached.                                                         */
/*                                                                         */
/***************************************************************************/
struct FormatType *GetFormat(const char *filedir, const char *linterp_prefix, int *error_code) {
  int i_format, i;
  //struct stat statbuf;
  ZZIP_FILE *fp;
  char format_file[MAX_FILENAME_LENGTH+6];
  struct FormatType *F;
  char raw_data_filename[MAX_FILENAME_LENGTH+FIELD_LENGTH+2];
  char **IncludeList = NULL;
  int i_include;

  /***** open the format file (if there is one) ******/
  snprintf(format_file, MAX_FILENAME_LENGTH+6, "%s/format", filedir);
  fp = zzip_fopen(format_file, "r");
  if (fp == NULL) {
    printf("%s\n", format_file);
    *error_code = GD_E_OPEN_FORMAT;
    return (NULL);
  }

  F = (struct FormatType *) malloc( sizeof(struct FormatType) );

  strcpy(F->FileDirName, filedir);
  F->n_raw = F->n_polynom = F->n_lincom = F->n_multiply = F->n_linterp = F->n_mplex = F->n_bit
    = 0;
  F->frame_offset = 0;
  F->rawEntries = NULL;
  F->polynomEntries = NULL;
  F->lincomEntries = NULL;
  F->multiplyEntries = NULL;
  F->linterpEntries = NULL;
  F->mplexEntries = NULL;
  F->bitEntries = NULL;

  /* Parse the file.  This will take care of any necessary inclusions */
  i_include = 1;
  IncludeList = malloc(sizeof(char*));
  IncludeList[0] = strdup("format");
  *error_code = ParseFormatFile(fp, F, filedir, ".", linterp_prefix, &IncludeList, &i_include);
  zzip_fclose(fp);

  /* Clean up IncludeList.  We don't need it anymore */
  for (i = 0; i < i_include; ++i)
    free(IncludeList[i]);
  free(IncludeList);

  if (*error_code!=GD_E_OK) {
    FreeF(F);
    return(NULL);
  }

  /** Now sort the lists */
  if (F->n_raw > 1) {
    for (i=0; i<F->n_raw; i++) {
      snprintf(raw_data_filename, MAX_FILENAME_LENGTH+FIELD_LENGTH+2, 
          "%s/%s", filedir, F->rawEntries[i].file);
      //if (stat(raw_data_filename, &statbuf) >=0) {
      if ( file_exists(raw_data_filename) ) {
        F->first_field = F->rawEntries[i];
        break;
      }
      snprintf(raw_data_filename, MAX_FILENAME_LENGTH+FIELD_LENGTH+2, 
          "%s/%s.slm", filedir, F->rawEntries[i].file);
      //if (stat(raw_data_filename, &statbuf) >=0) {
      if ( file_exists(raw_data_filename) ) {
        F->first_field = F->rawEntries[i];
        break;
      }
    }

    qsort(F->rawEntries, F->n_raw, sizeof(struct RawEntryType),
        RawCmp);
  }

  if (F->n_polynom > 1) {
    qsort(F->polynomEntries, F->n_polynom, sizeof(struct PolynomEntryType),
        PolynomCmp);
  }
  if (F->n_lincom > 1) {
    qsort(F->lincomEntries, F->n_lincom, sizeof(struct LincomEntryType),
        LincomCmp);
  }
  if (F->n_linterp > 1) {
    qsort(F->linterpEntries, F->n_linterp, sizeof(struct LinterpEntryType),
        LinterpCmp);
  }
  if (F->n_multiply > 1) {
    qsort(F->multiplyEntries, F->n_multiply, sizeof(struct MultiplyEntryType),
        MultiplyCmp);
  }
  if (F->n_mplex > 1) {
    qsort(F->mplexEntries, F->n_mplex, sizeof(struct MplexEntryType),
        MplexCmp);
  }
  if (F->n_bit > 1) {
    qsort(F->bitEntries, F->n_bit, sizeof(struct BitEntryType),
        BitCmp);
  }
  return(F);
}

/***************************************************************************/
/*                                                                         */
/*     File File Frame numbers into dataout                                */
/*                                                                         */
/***************************************************************************/
static void FillFileFrame(void *dataout, char rtype, int s0, int n) {
  int i;

  switch(rtype) {
    case 'c':
      for (i=0; i<n; i++) {
        ((char*)dataout)[i] = (char)i+s0;
      }
      break;
    case 'i': /* for compatibility with creaddata. (deprecated) */
    case 'S':
      for (i=0; i<n; i++) {
        ((int*)dataout)[i] = (int)i+s0;
      }
      break;
    case 's':
      for (i=0; i<n; i++) {
        ((short*)dataout)[i] = (short)i+s0;
      }
      break;
    case 'U':
      for (i=0; i<n; i++) {
        ((unsigned int *)dataout)[i] = (unsigned int)i+s0;
      }
      break;
    case 'u':
      for (i=0; i<n; i++) {
        ((unsigned short *)dataout)[i] = (unsigned short)i+s0;
      }
      break;
    case 'f':
      for (i=0; i<n; i++) {
        ((float*)dataout)[i] = (float)i+s0;
      }
      break;
    case 'd':
      for (i=0; i<n; i++) {
        ((double*)dataout)[i] = (double)i+s0;
      }
      break;
  }
}
/***************************************************************************/
/*                                                                         */
/*    ConvertType: copy data to output buffer while converting type        */
/*           Returns error code                                            */
/*                                                                         */
/***************************************************************************/
static int ConvertType(unsigned char *data_in, char in_type,
                       void *data_out, char out_type, int n) {
  int i;

  if (out_type=='n') { /* null return type: don't return data */
    return(0);
  }

  switch (in_type) {
    case 'c':
      switch (out_type) {
        case 'c':
          for (i=0;i<n;i++) ((unsigned char*) data_out)[i]=data_in[i];
          break;
        case 's':
          for (i=0;i<n;i++) ((short*)data_out)[i]=data_in[i];
          break;
        case 'u':
          for (i=0;i<n;i++) ((unsigned short*)data_out)[i]=data_in[i];
          break;
        case 'i': case 'S':
          for (i=0;i<n;i++) ((int*)data_out)[i]=data_in[i];
          break;
        case 'U':
          for (i=0;i<n;i++)  ((unsigned int*)data_out)[i]=data_in[i];
          break;
        case 'f':
          for (i=0;i<n;i++) ((float*)data_out)[i]=data_in[i];
          break;
        case 'd':
          for (i=0;i<n;i++) ((double*)data_out)[i]=data_in[i];
          break;
        default:
          return (GD_E_BAD_RETURN_TYPE);
      }
      break;
    case 's':
      switch (out_type) {
        case 'c':
          for (i=0;i<n;i++) ((unsigned char*) data_out)[i]=((short *)data_in)[i];
          break;
        case 's':
          for (i=0;i<n;i++) ((short*)data_out)[i]=((short *)data_in)[i];
          break;
        case 'u':
          for (i=0;i<n;i++) ((unsigned short*)data_out)[i]=((short *)data_in)[i];
          break;
        case 'S': case 'i':
          for (i=0;i<n;i++) ((int*)data_out)[i]=((short *)data_in)[i];
          break;
        case 'U':
          for (i=0;i<n;i++) ((unsigned int*)data_out)[i]=((short *)data_in)[i];
          break;
        case 'f':
          for (i=0;i<n;i++) ((float*)data_out)[i]=((short *)data_in)[i];
          break;
        case 'd':
          for (i=0;i<n;i++) ((double*)data_out)[i]=((short *)data_in)[i];
          break;
        default:
          return (GD_E_BAD_RETURN_TYPE);
      }
      break;
    case 'u':
      switch (out_type) {
        case 'c':
          for (i=0;i<n;i++) ((unsigned char*) data_out)[i]=
            ((unsigned short *)data_in)[i];
          break;
        case 's':
          for (i=0;i<n;i++) ((short*)data_out)[i]=
            ((unsigned short *)data_in)[i];
          break;
        case 'u':
          for (i=0;i<n;i++) ((unsigned short*)data_out)[i]=
            ((unsigned short *)data_in)[i];
          break;
        case 'i': case 'S':
          for (i=0;i<n;i++) ((int*)data_out)[i]=
            ((unsigned short *)data_in)[i];
          break;
        case 'U':
          for (i=0;i<n;i++)  ((unsigned int*)data_out)[i]=
            ((unsigned short *)data_in)[i];
          break;
        case 'f':
          for (i=0;i<n;i++) ((float*)data_out)[i]=
            ((unsigned short *)data_in)[i];
          break;
        case 'd':
          for (i=0;i<n;i++) ((double*)data_out)[i]=
            ((unsigned short *)data_in)[i];
          break;
        default:
          return (GD_E_BAD_RETURN_TYPE);
      }
      break;
    case 'i':
    case 'S':
      switch (out_type) {
        case 'c':
          for (i=0;i<n;i++) ((unsigned char*) data_out)[i]=((int *)data_in)[i];
          break;
        case 's':
          for (i=0;i<n;i++) ((short*)data_out)[i]=((int *)data_in)[i];
          break;
        case 'u':
          for (i=0;i<n;i++) ((unsigned short*)data_out)[i]=((int *)data_in)[i];
          break;
        case 'i': case 'S':
          for (i=0;i<n;i++) ((int*)data_out)[i]=((int *)data_in)[i];
          break;
        case 'U':
          for (i=0;i<n;i++)  ((unsigned int*)data_out)[i]=((int *)data_in)[i];
          break;
        case 'f':
          for (i=0;i<n;i++) ((float*)data_out)[i]=((int *)data_in)[i];
          break;
        case 'd':
          for (i=0;i<n;i++) ((double*)data_out)[i]=((int *)data_in)[i];
          break;
        default:
          return (GD_E_BAD_RETURN_TYPE);
      }
      break;
    case 'U':
      switch (out_type) {
        case 'c':
          for (i=0;i<n;i++) ((unsigned char*) data_out)[i]=
            ((unsigned *)data_in)[i];
          break;
        case 's':
          for (i=0;i<n;i++) ((short*)data_out)[i]=
            ((unsigned *)data_in)[i];
          break;
        case 'u':
          for (i=0;i<n;i++) ((unsigned short*)data_out)[i]=
            ((unsigned *)data_in)[i];
          break;
        case 'i': case 'S':
          for (i=0;i<n;i++) ((int*)data_out)[i]=
            ((unsigned *)data_in)[i];
          break;
        case 'U':
          for (i=0;i<n;i++)  ((unsigned int*)data_out)[i]=
            ((unsigned *)data_in)[i];
          break;
        case 'f':
          for (i=0;i<n;i++) ((float*)data_out)[i]=
            ((unsigned *)data_in)[i];
          break;
        case 'd':
          for (i=0;i<n;i++) ((double*)data_out)[i]=
            ((unsigned *)data_in)[i];
          break;
        default:
          return (GD_E_BAD_RETURN_TYPE);
      }
      break;
    case 'f':
      switch (out_type) {
        case 'c':
          for (i=0;i<n;i++) ((unsigned char*) data_out)[i]=((float *)data_in)[i];
          break;
        case 's':
          for (i=0;i<n;i++) ((short*)data_out)[i]=((float *)data_in)[i];
          break;
        case 'u':
          for (i=0;i<n;i++) ((unsigned short*)data_out)[i]=((float *)data_in)[i];
          break;
        case 'i': case 'S':
          for (i=0;i<n;i++) ((int*)data_out)[i]=((float *)data_in)[i];
          break;
        case 'U':
          for (i=0;i<n;i++)  ((unsigned int*)data_out)[i]=((float *)data_in)[i];
          break;
        case 'f':
          for (i=0;i<n;i++) ((float*)data_out)[i]=((float *)data_in)[i];
          break;
        case 'd':
          for (i=0;i<n;i++) ((double*)data_out)[i]=((float *)data_in)[i];
          break;
        default:
          return (GD_E_BAD_RETURN_TYPE);
      }
      break;
    case 'd':
      switch (out_type) {
        case 'c':
          for (i=0;i<n;i++) ((unsigned char*) data_out)[i]=((double *)data_in)[i];
          break;
        case 's':
          for (i=0;i<n;i++) ((short*)data_out)[i]=((double *)data_in)[i];
          break;
        case 'u':
          for (i=0;i<n;i++) ((unsigned short*)data_out)[i]=((double *)data_in)[i];
          break;
        case 'i': case 'S':
          for (i=0;i<n;i++) ((int*)data_out)[i]=((double *)data_in)[i];
          break;
        case 'U':
          for (i=0;i<n;i++)  ((unsigned int*)data_out)[i]=((double *)data_in)[i];
          break;
        case 'f':
          for (i=0;i<n;i++) ((float*)data_out)[i]=((double *)data_in)[i];
          break;
        case 'd':
          for (i=0;i<n;i++) ((double*)data_out)[i]=((double *)data_in)[i];
          break;
        default:
          return (GD_E_BAD_RETURN_TYPE);
          break;
      }
      break;
    default:
      printf("internal getdata bug: unknown type shouldn't make it here!\n");
      return (GD_E_BAD_RETURN_TYPE);
  }

  return(GD_E_OK);
}


/***************************************************************************/
/*                                                                         */
/*      FillZero: fill data buffer with zero of the appropriate type       */
/*        used if s0<0 - fill up to 0, or up to ns+s0, whichever is less   */
/*                                                                         */
/***************************************************************************/
static int FillZero(char *databuffer, char type, int s0, int ns) {
  int nz;

  if (s0>=0) return 0;

  if (s0+ns>0) nz = -s0;
  else nz = ns;

  switch (type) {
    case 'c':
      memset(databuffer, 0, nz);
      break;
    case 's':
    case 'u':
      memset(databuffer, 0, nz*sizeof(short));
      break;
    case 'i':
    case 'S':
    case 'U':
      memset(databuffer, 0, nz*sizeof(int));
      break;
    case 'f':
      memset(databuffer, 0, nz*sizeof(float));
      break;
    case 'd':
      memset(databuffer, 0, nz*sizeof(double));
      break;
  }

  return (nz);

}
/***************************************************************************/
/*                                                                         */
/*   Get samples per frame for field, given FormatType *F                  */
/*                                                                         */
/***************************************************************************/
static int GetSPF(int recurse_level, const char *field_code, const struct FormatType *F, int *error_code) {
  struct RawEntryType tR;
  struct RawEntryType *R;
  struct PolynomEntryType tP;
  struct PolynomEntryType *P;
  struct LincomEntryType tL;
  struct LincomEntryType *L;
  struct MultiplyEntryType tM;
  struct MultiplyEntryType *M;
  struct BitEntryType tB;
  struct BitEntryType *B;
  struct LinterpEntryType tI;
  struct LinterpEntryType *I;
  int spf;

  if (!F) { /* don't crash */
    return(0);
  }

  if (recurse_level > 10) {
    *error_code = GD_E_RECURSE_LEVEL;
    return(0);
  }

  if ((strcmp(field_code,"FILEFRAM")==0) ||
      (strcmp(field_code,"INDEX")==0)) {
    return(1);
  }

  /***************************************/
  /** Check to see if it is a raw entry **/
  /* binary search for the field */
  /* make a RawEntry we can compare to */
  strncpy(tR.field, field_code, FIELD_LENGTH);
  /** use the stdlib binary search */
  R = bsearch(&tR, F->rawEntries, F->n_raw,
      sizeof(struct RawEntryType), RawCmp);
  if (R!=NULL) {
    spf = R->samples_per_frame;
    return(spf);
  }

  /***************************************/
  /** Check to see if it is a polynom entry **/
  /* binary search for the field */
  /* make a RawEntry we can compare to */
  strncpy(tP.field, field_code, FIELD_LENGTH);
  /** use the stdlib binary search */
  P = bsearch(&tP, F->polynomEntries, F->n_polynom,
      sizeof(struct PolynomEntryType), PolynomCmp);
  if (P!=NULL) {
    spf = GetSPF(recurse_level+1, P->in_field, F, error_code);
    return(spf);
  }

  /***************************************/
  /** Check to see if it is a lincom entry **/
  /* binary search for the field */
  /* make a RawEntry we can compare to */
  strncpy(tL.field, field_code, FIELD_LENGTH);
  /** use the stdlib binary search */
  L = bsearch(&tL, F->lincomEntries, F->n_lincom,
      sizeof(struct LincomEntryType), LincomCmp);
  if (L!=NULL) {
    spf = GetSPF(recurse_level+1, L->in_fields[0], F, error_code);
    return(spf);
  }

  /***************************************/
  /** Check to see if it is a multiply entry **/
  /* binary search for the field */
  /* make a RawEntry we can compare to */
  strncpy(tM.field, field_code, FIELD_LENGTH);
  /** use the stdlib binary search */
  M = bsearch(&tM, F->multiplyEntries, F->n_multiply,
      sizeof(struct MultiplyEntryType), MultiplyCmp);
  if (M != NULL) {
    int spf2;
    spf = GetSPF(recurse_level+1,M->in_fields[0], F, error_code);
    spf2 = GetSPF(recurse_level+1,M->in_fields[1], F, error_code);
    if (spf2 > spf)
      spf = spf2;
    return(spf);
  }

  /***************************************/
  /** Check to see if it is a bit entry **/
  /* binary search for the field */
  /* make a BitEntry we can compare to */
  strncpy(tB.field, field_code, FIELD_LENGTH);
  /** use the stdlib binary search */
  B = bsearch(&tB, F->bitEntries, F->n_bit,
      sizeof(struct BitEntryType), BitCmp);
  if (B!=NULL) {
    spf = GetSPF(recurse_level+1,B->raw_field, F, error_code);
    return(spf);
  }

  /***************************************/
  /** Check to see if it is a linterp entry **/
  /* binary search for the field */
  /* make a LinterpEntry we can compare to */
  strncpy(tI.field, field_code, FIELD_LENGTH);
  /** use the stdlib binary search */
  I = bsearch(&tI, F->linterpEntries, F->n_linterp,
      sizeof(struct LinterpEntryType), LinterpCmp);
  if (I!=NULL) {
    spf = GetSPF(recurse_level+1, I->raw_field, F, error_code);
    return(spf);
  }

  printf("field_code = %s\n", field_code);
  *error_code = GD_E_BAD_CODE;
  return(0);
}

static inline void open_raw(struct FileHandle *R, const char *FileDirName,
                           const char *ChannelName) {

  char datafilename[2 * MAX_FILENAME_LENGTH + FIELD_LENGTH + 2];

  R->slim = NULL;

  /* Try to open raw file */
  snprintf(datafilename, 2 * MAX_FILENAME_LENGTH + FIELD_LENGTH + 2, 
           "%s/%s", FileDirName, ChannelName);
  R->fp = zzip_open(datafilename, O_RDONLY);
  if (R->fp)
    return;

  /* Try to open slimmed raw file */
  snprintf(datafilename, 2 * MAX_FILENAME_LENGTH + FIELD_LENGTH + 2, 
           "%s/%s.slm", FileDirName, ChannelName);
  R->slim = slimopen(datafilename, "r");
}

static inline off_t seek_wrap(struct FileHandle *R, off_t offset, int whence) {
  if (R->fp)
    return zzip_seek(R->fp, offset, whence);
  if (R->slim)
    return slimseek(R->slim, offset, whence);
  return (off_t)-1;
}

static inline ssize_t read_wrap(struct FileHandle *R, void *buf, size_t count) {
  if (R->fp)
    return zzip_read(R->fp, buf, count);
  if (R->slim)
    return slimread(buf, 1, count, R->slim); // syntax matches fread, not read.
  return (ssize_t)-1;
}

/***************************************************************************/
/*                                                                         */
/*   Look to see if the field code belongs to a raw.  If so, parse it.     */
/*                                                                         */
/***************************************************************************/
static int DoIfRaw(const struct FormatType *F, const char *field_code,
    int first_frame, int first_samp,
    int num_frames, int num_samp,
    char return_type, void *data_out,
    int *error_code, int *n_read) {

  struct RawEntryType tR;
  struct RawEntryType *R;
  int s0, ns, bytes_read;
  unsigned char *databuffer;
  struct FileHandle FH;

  /******* binary search for the field *******/
  /* make a RawEntry we can compare to */
  strncpy(tR.field, field_code, FIELD_LENGTH);
  /** use the stdlib binary search */
  R = bsearch(&tR, F->rawEntries, F->n_raw,
      sizeof(struct RawEntryType), RawCmp);
  if (R==NULL) return(0);

  /** if we got here, we found the field! **/
  s0 = first_samp + first_frame*R->samples_per_frame;
  ns = num_samp + num_frames*R->samples_per_frame;

  /** open the file */
  open_raw(&FH, F->FileDirName, R->file);
  if (FH.fp<0 && FH.slim == NULL) {
    *n_read = 0;
    *error_code = GD_E_OPEN_RAWFIELD;
    return(1);
  }

  databuffer = (unsigned char *)malloc(ns*R->size);

  *n_read = 0;
  if (s0 < 0) {
    *n_read = FillZero((char *)databuffer, R->type, s0, ns);
    ns -= *n_read;
    s0 = 0;
  }

  if (ns>0) {
    seek_wrap(&FH, s0*R->size, SEEK_SET);
    bytes_read = read_wrap(&FH, databuffer + *n_read*R->size, ns*R->size);
    *n_read += bytes_read/R->size;
  }

  *error_code =
    ConvertType(databuffer, R->type, data_out, return_type, *n_read);

  free(databuffer);

  if (FH.fp >= 0)
    zzip_close(FH.fp);
  if (FH.slim != NULL)
    slimclose(FH.slim);

  return(1);
}


/***************************************************************************/
/*                                                                         */
/*            AllocTmpbuff: allocate a buffer of the right type and size   */
/*                                                                         */
/***************************************************************************/
static void *AllocTmpbuff(char type, int n) {
  assert(n > 0);
  void *buff=NULL;
  switch(type) {
    case 'n':
      buff = NULL;
      break;
    case 'c':
      buff = malloc(n*sizeof(char));
      break;
    case 'i':
    case 'S':
    case 'U':
      buff = malloc(n*sizeof(int));
      break;
    case 's':
    case 'u':
      buff = malloc(n*sizeof(short));
      break;
    case 'f':
      buff = malloc(n*sizeof(float));
      break;
    case 'd':
      buff = malloc(n*sizeof(double));
      break;
    default:
      printf("Unexpected bad type error in AllocTmpbuff (%c)\n",type);
      abort();
      break;
  }
  if ((type != 'n') && (buff==NULL)) {
    printf("Memory Allocation error in AllocTmpbuff\n");
  }
  return(buff);
}

/***************************************************************************/
/*                                                                         */
/*   ScaleData: out = m*in+b                                               */
/*                                                                         */
/***************************************************************************/
static void ScaleData(void *data, char type, int npts, double m, double b) {
  char *data_c;
  short *data_s;
  unsigned short *data_u;
  unsigned *data_U;
  int *data_i;
  float *data_f;
  double *data_d;

  int i;

  switch(type) {
    case 'n':
      break;
    case 'c':
      data_c = (char *)data;
      for (i=0; i<npts; i++) {
        data_c[i] =(char)((double)data_c[i] * m + b);
      }
      break;
    case 's':
      data_s = (short *)data;
      for (i=0; i<npts; i++) {
        data_s[i] =  (short)((double)data_s[i] * m + b);
      }
      break;
    case 'u':
      data_u = (unsigned short *)data;
      for (i=0; i<npts; i++) {
        data_u[i] = (unsigned short)((double)data_u[i] * m + b);
      }
      break;
    case 'S': case 'i':
      data_i = (int *)data;
      for (i=0; i<npts; i++) {
        data_i[i] = (int)((double)data_i[i] * m + b);
      }
      break;
    case 'U':
      data_U = (unsigned*)data;
      for (i=0; i<npts; i++) {
        data_U[i] = (unsigned)((double)data_U[i] * m + b);
      }
      break;
    case 'f':
      data_f = (float *)data;
      for (i=0; i<npts; i++) {
        data_f[i] = (float)((double)data_f[i] * m + b);
      }
      break;
    case 'd':
      data_d = (double *)data;
      for (i=0; i<npts; i++) {
        data_d[i] = data_d[i] * m + b;
      }
      break;
    default:
      printf("Another impossible error\n");
      abort();
      break;
  }
}

/***************************************************************************/
/*                                                                         */
/*   PolynomData: out = a[0] + a[1]*in + ...                               */
/*                                                                         */
/***************************************************************************/

/* Macros to reduce tangly code */
#define POLYNOM5(t,npts) \
  for (i = 0; i < npts; i++) ((t*)data)[i] = (t)( \
      ((t*)data)[i] * ((t*)data)[i] * ((t*)data)[i] \
      * ((t*)data)[i] * ((t*)data)[i] * a[5] \
      + ((t*)data)[i] * ((t*)data)[i] * ((t*)data)[i] * ((t*)data)[i] * a[4] \
      + ((t*)data)[i] * ((t*)data)[i] * ((t*)data)[i] * a[3] \
      + ((t*)data)[i] * ((t*)data)[i] * a[2] \
      + ((t*)data)[i] * a[1] + a[0] \
      )

#define POLYNOM4(t,npts) \
  for (i = 0; i < npts; i++) ((t*)data)[i] = (t)( \
      ((t*)data)[i] * ((t*)data)[i] * ((t*)data)[i] * ((t*)data)[i] * a[4] \
      + ((t*)data)[i] * ((t*)data)[i] * ((t*)data)[i] * a[3] \
      + ((t*)data)[i] * ((t*)data)[i] * a[2] \
      + ((t*)data)[i] * a[1] + a[0] \
      )

#define POLYNOM3(t,npts) \
  for (i = 0; i < npts; i++) ((t*)data)[i] = (t)( \
      ((t*)data)[i] * ((t*)data)[i] * ((t*)data)[i] * a[3] \
      + ((t*)data)[i] * ((t*)data)[i] * a[2] \
      + ((t*)data)[i] * a[1] + a[0] \
      )

#define POLYNOM2(t,npts) \
  for (i = 0; i < npts; i++) ((t*)data)[i] = (t)( \
      ((t*)data)[i] * ((t*)data)[i] * a[2] \
      + ((t*)data)[i] * a[1] + a[0] \
      )

#define POLYNOM(t) \
  switch (n) { \
    case 2: POLYNOM2(t,npts); break; \
    case 3: POLYNOM3(t,npts); break; \
    case 4: POLYNOM4(t,npts); break; \
    case 5: POLYNOM5(t,npts); break; \
  }

static void PolynomData(void *data, char type, int npts, int n, const double *a) {
  int i;

  if (n == 1) {
    /* no need to duplicate this case */
    ScaleData(data, type, npts, a[1], a[0]);
  } else {
    switch (type) {
      case 'n':                              break;
      case 'c':     POLYNOM(          char); break;
      case 's':     POLYNOM(         short); break;
      case 'u':     POLYNOM(unsigned short); break;
      case 'S':     POLYNOM(           int); break;
      case 'U':     POLYNOM(  unsigned int); break;
      case 'f':     POLYNOM(         float); break;
      case 'd':     POLYNOM(        double); break;
      default:
        printf("Another impossible error\n");
        abort();
        break;
    }
  }
}

/***************************************************************************/
/*                                                                         */
/*            AddData: add B to A.  B is unchanged                         */
/*                                                                         */
/***************************************************************************/
static void AddData(void *A, int spfA, void *B, int spfB, char type, int n) {
  int i;

  switch(type) {
    case 'n': /* null read */
      break;
    case 'c':
      for (i=0; i<n; i++) {
        ((char*)A)[i] += ((char*)B)[i * spfB / spfA];
      }
      break;
    case 'S': case 'i':
      for (i=0; i<n; i++) {
        ((int*)A)[i] += ((int*)B)[i * spfB / spfA];
      }
      break;
    case 's':
      for (i=0; i<n; i++) {
        ((short*)A)[i] += ((short*)B)[i * spfB / spfA];
      }
      break;
    case 'u':
      for (i=0; i<n; i++) {
        ((unsigned short*)A)[i] += ((unsigned short*)B)[i * spfB / spfA];
      }
      break;
    case 'U':
      for (i=0; i<n; i++) {
        ((unsigned*)A)[i] += ((unsigned*)B)[i * spfB / spfA];
      }
      break;
    case 'f':
      for (i=0; i<n; i++) {
        ((float*)A)[i] += ((float*)B)[i * spfB / spfA];
      }
      break;
    case 'd':
      for (i=0; i<n; i++) {
        ((double*)A)[i] += ((double*)B)[i * spfB / spfA];
      }
      break;
    default:
      printf("Unexpected bad type error in AddData\n");
      abort();
      break;
  }
}

/***************************************************************************/
/*                                                                         */
/*            MultiplyData: multiply B by A.  B is unchanged               */
/*                                                                         */
/***************************************************************************/
static void MultiplyData(void *A, int spfA, void *B, int spfB, char type, int n)
{
  int i;

  switch(type) {
    case 'n': /* null read */
      break;
    case 'c':
      for (i=0; i<n; i++) {
        ((char*)A)[i] *= ((char*)B)[i * spfB / spfA];
      }
      break;
    case 'S': case 'i':
      for (i=0; i<n; i++) {
        ((int*)A)[i] *= ((int*)B)[i * spfB / spfA];
      }
      break;
    case 's':
      for (i=0; i<n; i++) {
        ((short*)A)[i] *= ((short*)B)[i * spfB / spfA];
      }
      break;
    case 'u':
      for (i=0; i<n; i++) {
        ((unsigned short*)A)[i] *= ((unsigned short*)B)[i * spfB / spfA];
      }
      break;
    case 'U':
      for (i=0; i<n; i++) {
        ((unsigned*)A)[i] *= ((unsigned*)B)[i * spfB / spfA];
      }
      break;
    case 'f':
      for (i=0; i<n; i++) {
        ((float*)A)[i] *= ((float*)B)[i * spfB / spfA];
      }
      break;
    case 'd':
      for (i=0; i<n; i++) {
        ((double*)A)[i] *= ((double*)B)[i * spfB / spfA];
      }
      break;
    default:
      printf("Unexpected bad type error in MultiplyData\n");
      abort();
      break;
  }
}

/***************************************************************************/
/*                                                                         */
/*   Look to see if the field code belongs to a polynom.  If so, parse it.  */
/*                                                                         */
/***************************************************************************/
static int DoIfPolynom(int recurse_level,
    const struct FormatType *F, const char *field_code,
    int first_frame, int first_samp,
    int num_frames, int num_samp,
    char return_type, void *data_out,
    int *error_code, int *n_read) {
  struct PolynomEntryType tP;
  struct PolynomEntryType *P;
  void *tmpbuf;
  int i;
  int spf1, spf2, num_samp2, first_samp2;
  int n_read2;

  /******* binary search for the field *******/
  /* make a PolynomEntry we can compare to */
  strncpy(tP.field, field_code, FIELD_LENGTH);
  /** use the stdlib binary search */
  P = bsearch(&tP, F->polynomEntries, F->n_polynom,
      sizeof(struct PolynomEntryType), PolynomCmp);
  if (P==NULL) return(0);

  /*****************************************/
  /** if we got here, we found the field! **/
  /** read into dataout and scale the first element **/
  spf1 = GetSPF(recurse_level+1,P->in_field, F, error_code);
  if (*error_code != GD_E_OK) return(1);

  /* read and scale the first field and record the number of samples
   * returned */
  *n_read = DoField(recurse_level+1, F, P->in_field,
      first_frame, first_samp,
      num_frames, num_samp,
      return_type, data_out,
      error_code);

  if (*error_code != GD_E_OK)
    return(1);

  /* Nothing to polynomise */
  if (*n_read == 0)
    return 1;

  PolynomData(data_out, return_type, *n_read, P->n_coeff, P->a);

  return(1);
}

/***************************************************************************/
/*                                                                         */
/*   Look to see if the field code belongs to a lincom.  If so, parse it.  */
/*                                                                         */
/***************************************************************************/
static int DoIfLincom(int recurse_level,
    const struct FormatType *F, const char *field_code,
    int first_frame, int first_samp,
    int num_frames, int num_samp,
    char return_type, void *data_out,
    int *error_code, int *n_read) {
  struct LincomEntryType tL;
  struct LincomEntryType *L;
  void *tmpbuf;
  int i;
  int spf1, spf2, num_samp2, first_samp2;
  int n_read2;

  /******* binary search for the field *******/
  /* make a LincomEntry we can compare to */
  strncpy(tL.field, field_code, FIELD_LENGTH);
  /** use the stdlib binary search */
  L = bsearch(&tL, F->lincomEntries, F->n_lincom,
      sizeof(struct LincomEntryType), LincomCmp);
  if (L==NULL) return(0);

  /*****************************************/
  /** if we got here, we found the field! **/
  /** read into dataout and scale the first element **/
  spf1 = GetSPF(recurse_level+1,L->in_fields[0], F, error_code);
  if (*error_code != GD_E_OK) return(1);

  /* read and scale the first field and record the number of samples
   * returned */
  *n_read = DoField(recurse_level+1, F, L->in_fields[0],
      first_frame, first_samp,
      num_frames, num_samp,
      return_type, data_out,
      error_code);

  if (*error_code != GD_E_OK)
    return(1);

  /* Nothing to lincomise */
  if (*n_read == 0)
    return 1;

  ScaleData(data_out, return_type, *n_read, L->m[0], L->b[0]);

  if (L->n_infields > 1) {
    for (i=1; i<L->n_infields; i++) {

      /* find the samples per frame of the next field */
      spf2 = GetSPF(recurse_level+1, L->in_fields[i], F, error_code);
      if (*error_code != GD_E_OK) return(1);

      /* calculate the first sample and number of samples to read of the
       * next field */
      num_samp2 = (int)ceil((double)*n_read * spf2 / spf1);
      first_samp2 = (first_frame * spf2 + first_samp * spf2 / spf1);

      /* Allocate a temporary buffer for the next field */
      tmpbuf = AllocTmpbuff(return_type, num_samp2);
      if (!tmpbuf && return_type != 'n') {
        return(0);
      }

      /* read the next field */
      n_read2 = DoField(recurse_level+1, F, L->in_fields[i],
          0, first_samp2,
          0, num_samp2,
          return_type, tmpbuf,
          error_code);
      if (*error_code != GD_E_OK) {
        free(tmpbuf);
        return(1);
      }

      ScaleData(tmpbuf, return_type, n_read2, L->m[i], L->b[i]);

      if (n_read2 > 0 && n_read2 * spf1 != *n_read * spf2) {
        *n_read = n_read2 * spf1 / spf2;
      }

      AddData(data_out, spf1, tmpbuf, spf2, return_type, *n_read);

      free(tmpbuf);
    }
  }

  return(1);
}

/***************************************************************************/
/*                                                                         */
/*  Look to see if the field code belongs to a multiply.  If so, parse it. */
/*                                                                         */
/***************************************************************************/
static int DoIfMultiply(int recurse_level,
    const struct FormatType *F, const char *field_code,
    int first_frame, int first_samp,
    int num_frames, int num_samp,
    char return_type, void *data_out,
    int *error_code, int *n_read) {
  struct MultiplyEntryType tM;
  struct MultiplyEntryType *M;
  void *tmpbuf;
  int spf1, spf2, num_samp2, first_samp2;
  int n_read2;

  /******* binary search for the field *******/
  /* make a MultiplyEntry we can compare to */
  strncpy(tM.field, field_code, FIELD_LENGTH);
  /** use the stdlib binary search */
  M = bsearch(&tM, F->multiplyEntries, F->n_multiply,
      sizeof(struct MultiplyEntryType), MultiplyCmp);
  if (M==NULL) return(0);

  /*****************************************/
  /** if we got here, we found the field! **/
  /** read into dataout and scale the first element **/

  /* find the samples per frame of the first field */
  spf1 = GetSPF(recurse_level+1, M->in_fields[0], F, error_code);
  if (*error_code != GD_E_OK) return(1);

  /* read the first field and record the number of samples
   * returned */
  *n_read = DoField(recurse_level+1, F, M->in_fields[0],
      first_frame, first_samp,
      num_frames, num_samp,
      return_type, data_out,
      error_code);

  if (*error_code != GD_E_OK)
    return 1;

  /* Nothing to multiply */
  if (*n_read == 0)
    return 1;

  /* find the samples per frame of the second field */
  spf2 = GetSPF(recurse_level+1, M->in_fields[1], F, error_code);
  if (*error_code != GD_E_OK) return(1);

  /* calculate the first sample and number of samples to read of the
   * second field */
  num_samp2 = (int)ceil((double)*n_read * spf2 / spf1);
  first_samp2 = (first_frame * spf2 + first_samp * spf2 / spf1);

  /* Allocate a temporary buffer for the second field */
  tmpbuf = AllocTmpbuff(return_type, num_samp2);
  if (!tmpbuf && return_type != 'n') {
    return(0);
  }

  /* read the second field */
  n_read2 = DoField(recurse_level+1, F, M->in_fields[1],
      0, first_samp2,
      0, num_samp2,
      return_type, tmpbuf,
      error_code);
  if (*error_code != GD_E_OK) {
    free(tmpbuf);
    return(1);
  }

  if (n_read2 > 0 && n_read2 * spf1 < *n_read * spf2) {
    *n_read = n_read2 * spf1 / spf2;
  }
  MultiplyData(data_out, spf1, tmpbuf, spf2, return_type, *n_read);
  free(tmpbuf);

  return(1);
}

/***************************************************************************/
/*                                                                         */
/*   Look to see if the field code belongs to a bitfield.  If so, parse it.*/
/*                                                                         */
/***************************************************************************/
static int DoIfBit(int recurse_level,
    const struct FormatType *F, const char *field_code,
    int first_frame, int first_samp,
    int num_frames, int num_samp,
    char return_type, void *data_out,
    int *error_code, int *n_read) {
  struct BitEntryType tB;
  struct BitEntryType *B;
  unsigned *tmpbuf;
  int i;
  int spf;
  int ns;
  unsigned mask;

  /******* binary search for the field *******/
  /* make a BitEntry we can compare to */
  strncpy(tB.field, field_code, FIELD_LENGTH);
  /** use the stdlib binary search */
  B = bsearch(&tB, F->bitEntries, F->n_bit,
      sizeof(struct BitEntryType), BitCmp);
  if (B==NULL) return(0);

  /*****************************************/
  /** if we got here, we found the field! **/
  spf = GetSPF(recurse_level+1, B->raw_field, F, error_code);
  if (*error_code!=GD_E_OK) {
    *n_read = 0;
    return(1);
  }

  ns = num_samp + num_frames*spf;
  tmpbuf = (unsigned *)malloc(ns*sizeof(unsigned));

  *n_read = DoField(recurse_level+1, F, B->raw_field,
      first_frame, first_samp,
      num_frames, num_samp,
      'U', tmpbuf,
      error_code);
  if (*error_code!=GD_E_OK) {
    free(tmpbuf);
    return(1);
  }

  if (B->numbits==32) mask = 0xffffffff;
  else mask = (unsigned)(pow(2,B->numbits)-0.9999);

  for (i=0; i<*n_read; i++) {
    tmpbuf[i] = (tmpbuf[i]>>B->bitnum) & mask;
  }

  *error_code = ConvertType((unsigned char *)tmpbuf, 'U',
      data_out, return_type, *n_read);
  free(tmpbuf);

  return(1);
}

/***************************************************************************/
/*                                                                         */
/*   ReadLinterpFile: Read in the linterp data for this field              */
/*                                                                         */
/***************************************************************************/
void MakeDummyLinterp(struct LinterpEntryType *E); /* prototype => no warning... */
void MakeDummyLinterp(struct LinterpEntryType *E) {
  E->n_interp = 2;
  E->x = (double *)malloc(2*sizeof(double));
  E->y = (double *)malloc(2*sizeof(double));
  E->x[0] = 0;
  E->y[0] = 0;
  E->x[1] = 1;
  E->y[1] = 1;
}

static int ReadLinterpFile(struct LinterpEntryType *E) {
  ZZIP_FILE *fp;
  int i;
  char line[255];

  fp = zzip_fopen(E->linterp_file, "r");
  if (fp==NULL) {
    MakeDummyLinterp(E);
    return (GD_E_OPEN_LINFILE);
  }

  size_t nlines;
  char *buf = read_file( fp );
  char **lines = read_lines_from_buffer( buf, &nlines );
  zzip_fclose(fp);

  /* first read the file to see how big it is */
  /*
  i=0;
  while (GetLine(fp, line)) {
    i++;
  }
  */
  if ( nlines < 2 ) {
    MakeDummyLinterp(E);
    free( lines );
    free( buf );
    return (GD_E_OPEN_LINFILE);
  }
  E->n_interp = nlines;
  E->x = (double *)malloc(nlines*sizeof(double));
  E->y = (double *)malloc(nlines*sizeof(double));

  for (i=0; i<E->n_interp; i++)
    sscanf(lines[i], "%lg %lg",&(E->x[i]), &(E->y[i]));

  free( lines );
  free( buf );

  return (GD_E_OK);
}

/***************************************************************************/
/*                                                                         */
/*   GetIndex: just linearly search - we are probably right to start with  */
/*                                                                         */
/***************************************************************************/
static int GetIndex(double x, double lx[], int idx, int n) {
  /* increment until we are bigger */
  while ((idx < n-2) && (x > lx[idx])) {
    idx++;
  }
  /* decrement until we are smaller */
  while ((idx > 0) && (x < lx[idx])) {
    idx--;
  }

  return(idx);
}

/***************************************************************************/
/*                                                                         */
/*   LinterpData: calibrate data using lookup table lx and ly              */
/*                                                                         */
/***************************************************************************/
static void LinterpData(void *data, char type, int npts,
    double *lx, double *ly, int n_ln) {
  int i, idx=0;
  double x;

  for (i=0; i<npts; i++) {
    switch(type) {
      case 'n':
        return;
        break;
      case 'c':
        x = ((char *)data)[i];
        idx = GetIndex(x, lx, idx, n_ln);
        ((char *)data)[i] = (char)(ly[idx] + (ly[idx+1]-ly[idx])/
                                   (lx[idx+1]-lx[idx]) * (x-lx[idx]));
        break;
      case 's':
        x = ((short *)data)[i];
        idx = GetIndex(x, lx, idx, n_ln);
        ((short *)data)[i] = (short)(ly[idx] + (ly[idx+1]-ly[idx])/
                                     (lx[idx+1]-lx[idx]) * (x-lx[idx]));
        break;
      case 'u':
        x = ((unsigned short *)data)[i];
        idx = GetIndex(x, lx, idx,n_ln);
        ((unsigned short *)data)[i] =
          (unsigned short)(ly[idx] + (ly[idx+1]-ly[idx])/
                           (lx[idx+1]-lx[idx]) * (x-lx[idx]));
        break;
      case 'S': case 'i':
        x = ((int *)data)[i];
        idx = GetIndex(x, lx, idx, n_ln);
        ((int *)data)[i] = (int)(ly[idx] + (ly[idx+1]-ly[idx])/
                                 (lx[idx+1]-lx[idx]) * (x-lx[idx]));
        break;
      case 'U':
        x = ((unsigned *)data)[i];
        idx = GetIndex(x, lx, idx, n_ln);
        ((unsigned *)data)[i] =
          (unsigned)(ly[idx] + (ly[idx+1]-ly[idx])/
                     (lx[idx+1]-lx[idx]) * (x-lx[idx]));
        break;
      case 'f':
        x = ((float *)data)[i];
        idx = GetIndex(x, lx, idx, n_ln);
        ((float *)data)[i] = (float)(ly[idx] + (ly[idx+1]-ly[idx])/
                                     (lx[idx+1]-lx[idx]) * (x-lx[idx]));
        break;
      case 'd':
        x = ((double *)data)[i];
        idx = GetIndex(x, lx, idx, n_ln);
        ((double *)data)[i] = (double)(ly[idx] + (ly[idx+1]-ly[idx])/
                                       (lx[idx+1]-lx[idx]) * (x-lx[idx]));
        break;
      default:
        printf("Another impossible error\n");
        abort();
        break;
    }
  }
}

/***************************************************************************/
/*                                                                         */
/*   Look to see if the field code belongs to a bitfield.  If so, parse it.*/
/*                                                                         */
/***************************************************************************/
static int DoIfLinterp(int recurse_level,
    const struct FormatType *F, const char *field_code,
    int first_frame, int first_samp,
    int num_frames, int num_samp,
    char return_type, void *data_out,
    int *error_code, int *n_read) {
  struct LinterpEntryType tI;
  struct LinterpEntryType *I;

  /******* binary search for the field *******/
  /* make a LinterpEntry we can compare to */
  strncpy(tI.field, field_code, FIELD_LENGTH);
  /** use the stdlib binary search */
  I = bsearch(&tI, F->linterpEntries, F->n_linterp,
      sizeof(struct LinterpEntryType), LinterpCmp);
  if (I==NULL) return(0);

  /*****************************************/
  /** if we got here, we found the field! **/
  if (I->n_interp<0) {
    *error_code = ReadLinterpFile(I);
    if (*error_code != GD_E_OK) {
      *n_read = 0;
      return(1);
    }
  }
  *n_read = DoField(recurse_level+1, F, I->raw_field,
      first_frame, first_samp,
      num_frames, num_samp,
      return_type, data_out,
      error_code);
  if (*error_code!=GD_E_OK) return(1);
  LinterpData(data_out, return_type, *n_read, I->x, I->y, I->n_interp);
  return(1);

}

/***************************************************************************/
/*                                                                         */
/*  DoField: Doing one field once F has been identified                    */
/*                                                                         */
/***************************************************************************/
static int DoField(int recurse_level,
    const struct FormatType *F, const char *field_code,
    int first_frame, int first_samp,
    int num_frames, int num_samp,
    char return_type, void *data_out,
    int *error_code) {
  int n_read;

  if (recurse_level>10) {
    *error_code = GD_E_RECURSE_LEVEL;
    return(0);
  }


  /********************************************/
  /* if Asking for "FILEFRAM" or "INDEX", just return it */
  if ((strcmp(field_code,"FILEFRAM")==0) ||
      (strcmp(field_code,"INDEX")==0)) {
    n_read = num_frames + num_samp;
    if (data_out!=NULL) {
      FillFileFrame(data_out, return_type, first_frame+first_samp+
          F->frame_offset, n_read);
    }
    *error_code=GD_E_OK;
    return(n_read);
  }

  if (DoIfRaw(F, field_code,
        first_frame, first_samp,
        num_frames, num_samp,
        return_type, data_out,
        error_code, &n_read)) {
    return(n_read);
  } else if (DoIfPolynom(recurse_level,
        F, field_code,
        first_frame, first_samp,
        num_frames, num_samp,
        return_type, data_out,
        error_code, &n_read)) {
    return(n_read);
  } else if (DoIfLincom(recurse_level,
        F, field_code,
        first_frame, first_samp,
        num_frames, num_samp,
        return_type, data_out,
        error_code, &n_read)) {
    return(n_read);
  } else if (DoIfBit(recurse_level,
        F, field_code,
        first_frame, first_samp,
        num_frames, num_samp,
        return_type, data_out,
        error_code, &n_read)) {
    return(n_read);
  } else if (DoIfLinterp(recurse_level,
        F, field_code,
        first_frame, first_samp,
        num_frames, num_samp,
        return_type, data_out,
        error_code, &n_read)) {
    return(n_read);
  } else if (DoIfMultiply(recurse_level,
        F, field_code,
        first_frame, first_samp,
        num_frames, num_samp,
        return_type, data_out,
        error_code, &n_read)) {
    return(n_read);
  } else {
    *error_code = GD_E_BAD_CODE;
    return(0);
  }
}

/***************************************************************************/
/*                                                                         */
/*  GetData: read BLAST format files.                                      */
/*    filename_in: the name of the file directory (raw files are in here)  */
/*    field_code: the name of the field you want to read                   */
/*    first_frame, first_samp: the first sample read is                    */
/*              first_samp + samples_per_frame*first_frame                 */
/*    num_frames, num_samps: the number of samples read is                 */
/*              num_samps + samples_per_frame*num_frames                   */
/*    return_type: data type of *data_out.  's': 16 bit signed             */
/*              'u' 16bit unsiged. 'S' 32bit signed 'U' 32bit unsigned     */
/*              'c' 8 bit signed                                           */
/*    void *data_out: array to put the data                                */
/*    *error_code: error code is returned here. If error_code==null,       */
/*               GetData prints the error message and exits                */
/*                                                                         */
/*    return value: returns number of samples actually read into data_out  */
/*                                                                         */
/***************************************************************************/
int GetData(const struct FormatType *F, const char *field_code,
    int first_frame, int first_samp,
    int num_frames, int num_samp,
    char return_type, void *data_out,
    int *error_code) {

  int n_read=0;

  *error_code = GD_E_OK;

  first_frame -= F->frame_offset;

  n_read = DoField(0, F, field_code,
      first_frame, first_samp,
      num_frames, num_samp,
      return_type, data_out,
      error_code);

  return(n_read);
}

/***************************************************************************/
/*                                                                         */
/*    Get the number of frames available                                   */
/*                                                                         */
/***************************************************************************/
int GetNFrames(const struct FormatType *F, int *error_code, const char *in_field) {
  char raw_data_filename[2 * MAX_FILENAME_LENGTH + FIELD_LENGTH + 2];
  //struct stat statbuf;
  long st_size;
  int nf;

  *error_code = GD_E_OK;

  if (!F || F->n_raw==0) {
    *error_code = GD_E_FORMAT;
    return(0);
  }

  /* load the first valid raw field, either as a regular or a slim file */
  snprintf(raw_data_filename, 2 * MAX_FILENAME_LENGTH + FIELD_LENGTH + 2, 
      "%s/%s", F->FileDirName, F->first_field.file);
  //if (stat(raw_data_filename, &statbuf) < 0) {
  if ( !file_exists(raw_data_filename) ) {
    snprintf(raw_data_filename, 2 * MAX_FILENAME_LENGTH + FIELD_LENGTH + 2, 
             "%s/%s.slm", F->FileDirName, F->first_field.file);
    st_size = slimrawsize(raw_data_filename);
    if (st_size < 0)
      return(0);
  } else
  {
    st_size = file_size( raw_data_filename );
  }

  nf = st_size/
    (F->first_field.size*F->first_field.samples_per_frame);
  nf += F->frame_offset;

  return(nf);
}

/***************************************************************************/
/*                                                                         */
/*    Get the number of samples for each frame for the given field         */
/*                                                                         */
/***************************************************************************/
int GetSamplesPerFrame(const struct FormatType *F, const char *field_name, int *error_code) {

  *error_code = GD_E_OK;

  if (!F || F->n_raw==0) {
    *error_code = GD_E_FORMAT;
    return(0);
  }

  return GetSPF(0, field_name, F, error_code);
}


/* vim: ts=2 sw=2 et
*/
