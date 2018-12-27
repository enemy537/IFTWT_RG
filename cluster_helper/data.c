/* The C clustering library.
 * Copyright (C) 2002 Michiel Jan Laurens de Hoon.
 *
 * This library was written at the Laboratory of DNA Information Analysis,
 * Human Genome Center, Institute of Medical Science, University of Tokyo,
 * 4-6-1 Shirokanedai, Minato-ku, Tokyo 108-8639, Japan.
 * Contact: michiel.dehoon 'AT' riken.jp
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 * 
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 * 
 */

/* This file contains C code needed for Cluster 3.0, particularly file reading
 * and data handling. It is platform-independent; platform-dependent code is
 * located in windows/gui.c (Microsoft Windows), in mac/Controller.m (Mac OS X),
 * and in x11/gui.c (X11 using Motif).
 * 
 * Michiel de Hoon, (michiel.dehoon 'AT' riken.jp).
 * University of Tokyo, Human Genome Center.
 * 2003.01.10.
*/

/*============================================================================*/
/* Header files                                                               */
/*============================================================================*/

/* Standard C header files */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Local header files */
#include "data.h"
#include "cluster.h" /* The C clustering library */


/*============================================================================*/
/* Data declaration                                                           */
/*============================================================================*/

static int _rows = 0;
static int _columns = 0;
static double* _geneweight = NULL;
static double* _arrayweight = NULL;
static double* _geneorder = NULL;  /* Saves gene order in the data file */
static double* _arrayorder = NULL; /* Saves array order in the data file */
static int* _geneindex = NULL;   /* Set by clustering methods for file output */
static int* _arrayindex = NULL;  /* Set by clustering methods for file output */
static char* _uniqID = NULL;     /* Stores UNIQID identifier in the data file */
static char** _geneuniqID = NULL;
static char** _genename = NULL;
static char** _arrayname = NULL;
static double** _data = NULL;
static int** _mask = NULL;

/*============================================================================*/
/* Utility routines                                                           */
/*============================================================================*/

static char* GetLine(FILE* inputfile)
/* The function GetLine reads one line from the inputfile, and returns it as a
 * null-terminated string. If inputfile is at EOF, an empty string is returned.
 * Empty lines are skipped.
 * If this function fails due to memory allocation error, it returns NULL.
 * The calling routine should free the char* returned by GetLine.
 */
{ int c;
  int n = 0;
  int size = 1023;
  char* temp;
  char* line = malloc((size+1)*sizeof(char));
  if (!line) return NULL;
  while (n==0)
  { while ((c = getc(inputfile))!=EOF && c!='\r' && c!='\n')
    { if (n == size)
      { size *= 2;
        temp = realloc(line,(size+1)*sizeof(char));
        if (!temp)
        { free(line);
          return NULL;
        }
        line = temp;
      }
      line[n] = (char)c;
      n++;
    }
    if (c=='\r')
    { c = getc(inputfile);
      if (c!='\n' && c!=EOF) ungetc(c,inputfile);
    }
    if (c==EOF) break;
  }
  line[n] = '\0';
  temp = realloc(line,(n+1)*sizeof(char));
  if (!temp)
  /* This should not happen, as temp is smaller than line.
   * But let's check to make sure anyway. */
  { free(line);
    return NULL;
  }
  return temp;
}

static char* tokenize(char* s)
{
  char* p = s;
  while (1)
  {
    if (*p=='\0') return NULL;
    if (*p=='\t')
    {
      *p = '\0';
      return p+1;
    }
    p++;
  }
}

static char* MakeID(const char* name, int i)
{ int n;
  char* ID;
  int ndigits = 1;
  int remainder = i;
  while (remainder/=10) ndigits++; /* Count how many digits there are in i */
  n = strlen(name) + ndigits + 2;
  /* One more for the X, and one more for the \0 termination character */
  ID = malloc(n*sizeof(char));
  if (ID) sprintf(ID, "%s%dX",name,i);
  return ID;
}

static int SetClusterIndex(char which, int k, int* clusterid)
{ int i;
  int cluster;
  int counter = 0;
  int* index = NULL;
  if (which=='g')
  { index = malloc(_rows*sizeof(int));
    if (!index) return 0;
    for (i=0; i<_rows; i++) index[i] = i;
    sort(_rows, _geneorder, index);
    for (cluster = 0; cluster < k; cluster++)
    { for (i = 0; i < _rows; i++)
      { const int j = index[i];
        if (clusterid[j]==cluster)
        { _geneindex[counter] = j;
          counter++;
        }
      }
    }
  }
  if (which=='a')
  { index = malloc(_columns*sizeof(int));
    if (!index) return 0;
    for (i=0; i<_columns; i++) index[i] = i;
    sort(_columns, _arrayorder, index);
    for (cluster = 0; cluster < k; cluster++)
    { for (i = 0; i < _columns; i++)
      { const int j = index[i];
        if (clusterid[j]==cluster)
        { _arrayindex[counter] = j;
          counter++;
        }
      }
    }
  }
  free(index);
  return 1;
}

static int
PerformGeneSOM(FILE* file, int XDim, int YDim, int iterations, double tau,
               char metric)
{ int i = 0;
  int j = 0;
  int k;
  int ok;

  int (*Group)[2] = malloc(_rows*sizeof(int[2]));
  double*** Nodes = malloc(XDim*sizeof(double**));
  int* clusterid = malloc(_rows*sizeof(int));
  int* index = malloc(_columns*sizeof(int));
  if (Nodes)
  { for (i = 0; i < XDim; i++)
    { Nodes[i] = malloc(YDim*sizeof(double*));
      j = 0;
      if (!Nodes[i]) break;
      for ( ; j < YDim; j++)
      { Nodes[i][j] = malloc(_columns*sizeof(double));
        if (!Nodes[i][j]) break;
      }
      if (j < YDim) break;
    }
  }
  if (!Group || !clusterid || !index || !Nodes || i < XDim || j < YDim)
  { if (Group) free(Group);
    if (clusterid) free(clusterid);
    if (index) free(index);
    if (Nodes)
    { if (i < XDim)
      { while (j--) free(Nodes[i][j]);
        free(Nodes[i]);
      }
      while (i--)
      { for (j = 0; j < YDim; j++) free(Nodes[i][j]);
        free(Nodes[i]);
      }
      free(Nodes);
    }
    return 0;
  }

  somcluster(_rows, _columns, _data, _mask, _arrayweight, 0,
    XDim, YDim, tau, iterations, metric, Nodes, Group);

  for (i=0; i<_rows; i++) clusterid[i] = Group[i][0] * YDim + Group[i][1];
  free(Group);

  for (k=0; k<_columns; k++) index[k] = k;
  sort(_columns, _arrayorder, index);
  fputs("NODE", file);
  for (i=0; i<_columns; i++) fprintf(file, "\t%s", _arrayname[index[i]]);
  putc('\n', file);
  for (i=0; i<XDim; i++)
  { for (j=0; j<YDim; j++)
    { fprintf(file, "NODE(%d,%d)", i, j);
      for (k=0; k<_columns; k++) fprintf(file, "\t%f", Nodes[i][j][index[k]]);
      putc('\n', file);
    }
  }
  free(index);

  for (i=0;i<XDim;i++)
  { for (j=0; j<YDim; j++) free(Nodes[i][j]);
    free(Nodes[i]);
  }
  free(Nodes);

  ok = SetClusterIndex('g', XDim * YDim, clusterid);
  free(clusterid);
  return ok;
}

static int
PerformArraySOM(FILE* file, int XDim, int YDim, int iterations, double tau,
                char metric)
{ int i = 0;
  int j = 0;
  int k;
  int ok;
  int (*Group)[2] = malloc(_columns*sizeof(int[2]));
  double*** Nodes = malloc(XDim*sizeof(double**));
  int* clusterid = malloc(_columns*sizeof(int));
  if (Nodes)
  { for (i = 0; i < XDim; i++)
    { Nodes[i] = malloc(YDim*sizeof(double*));
      j = 0;
      if (!Nodes[i]) break;
      for ( ; j < YDim; j++)
      { Nodes[i][j] = malloc(_rows*sizeof(double));
        if (!Nodes[i][j]) break;
      }
      if (j < YDim) break;
    }
  }
  if (!Group || !clusterid || !Nodes || i < XDim || j < YDim)
  { if (Group) free(Group);
    if (Nodes)
    { if (i < XDim)
      { while (j--) free(Nodes[i][j]);
        free(Nodes[i]);
      }
      while (i--)
      { for (j = 0; j < YDim; j++) free(Nodes[i][j]);
        free(Nodes[i]);
      }
      free(Nodes);
    }
    free(clusterid);
    return 0;
  }

  somcluster(_rows, _columns, _data, _mask, _geneweight, 1,
    XDim, YDim, tau, iterations, metric, Nodes, Group);

  for (i=0; i<_columns; i++)
    clusterid[i] = Group[i][0] * YDim + Group[i][1];
  free(Group);

  fprintf(file, "%s\t", _uniqID);
  for (i=0; i<XDim; i++)
    for (j=0; j<YDim; j++) fprintf(file, "\tNODE(%d,%d)", i, j);
  putc('\n', file);

  for (k=0;k<_rows;k++)
  { int index = _geneindex[k];
    fprintf(file, "%s\t", _geneuniqID[index]);
    if (_genename[index]) fputs(_genename[index], file);
    else fputs(_geneuniqID[index], file);
    for (i=0; i<XDim; i++)
      for (j=0; j<YDim; j++) fprintf(file, "\t%f", Nodes[i][j][index]);
    putc('\n', file);
  }

  for (i=0;i<XDim;i++)
  { for (j=0; j<YDim; j++) free(Nodes[i][j]);
    free(Nodes[i]);
  }
  free(Nodes);
  ok = SetClusterIndex('a', XDim * YDim, clusterid);
  free(clusterid);
  return ok;
}

/*============================================================================*/
/* Data handling routines                                                     */
/*============================================================================*/

void Free(void)
{ int row, column;
  if (_data)
  { for (row = 0; row < _rows; row++)
      if (_data[row]) free(_data[row]);
    free(_data);
  }
  if (_mask)
  { for (row = 0; row < _rows; row++)
      if (_mask[row]) free(_mask[row]);
    free(_mask);
  }
  if (_geneuniqID)
  { for (row = 0; row < _rows; row++)
      if (_geneuniqID[row]) free(_geneuniqID[row]);
    free(_geneuniqID);
  }
  if (_genename)
  { for (row = 0; row < _rows; row++)
      if (_genename[row]) free(_genename[row]);
    free(_genename);
  }
  if (_arrayname)
  { for (column = 0; column < _columns; column++)
      if (_arrayname[column]) free(_arrayname[column]);
    free(_arrayname);
  }
  if (_geneorder) free(_geneorder);
  if (_arrayorder) free(_arrayorder);
  if (_geneindex) free(_geneindex);
  if (_arrayindex) free(_arrayindex);
  if (_geneweight) free(_geneweight);
  if (_arrayweight) free(_arrayweight);
  if (_uniqID) free(_uniqID);
  _genename = NULL;
  _geneuniqID = NULL;
  _geneweight = NULL;
  _geneorder = NULL;
  _geneindex = NULL;
  _arrayname = NULL;
  _arrayweight = NULL;
  _arrayorder = NULL;
  _arrayindex = NULL;
  _data = NULL;
  _mask = NULL;
  _uniqID = NULL;
  _rows = 0;
  _columns = 0;
}

int GetRows(void)
{ return _rows;
}

int GetColumns(void)
{ return _columns;
}

char* Load(FILE* file)
/* Load in data from tab-delimited text file.
 * If an error occurs, an error message is returned.
 * If no error occurs, the string "ok" is returned.
 * In case of insufficient memory, NULL is returned.
 * All messages are allocated with malloc, and should be
 * freed by the calling routine. */
{ int row, column;           /* Counters for data matrix */
  int fileRow, fileColumn;   /* Counters for rows and columns in the file */
  int n;
  int nFileColumns;
  char* line;
  char* s;

  int geneNameColumn = -1;
  int geneWeightColumn = -1;
  int geneOrderColumn = -1;
  int arrayWeightRow = -1;
  int arrayOrderRow = -1;

  /* Deallocate previously allocated space */
  Free();

  /* Parse header line (first line) to find out what the columns are */
  line = GetLine(file);
  if (!line) return NULL;
  if(line[0]=='\0')
  { const char text[] = "Attempting to read empty file";
    const int m = strlen(text) + 1;
    char* error = malloc(m*sizeof(char));
    free(line);
    if (!error) return NULL;
    strcpy(error, text);
    return error;
  }
  s = tokenize(line); /* Skip the first column UNIQID */
  fileColumn = 1;
  _columns = 0;
  while (s)
  { char* token = s;
    s = tokenize(s);
    if (!strcmp(token,"NAME")) geneNameColumn = fileColumn;
    else if (!strcmp(token,"GWEIGHT")) geneWeightColumn = fileColumn;
    else if (!strcmp(token,"GORDER")) geneOrderColumn = fileColumn;
    else _columns++;
    fileColumn++;
  }
  free(line);
  nFileColumns = fileColumn;
  if (nFileColumns < 2) 
  { const char text[] = "Less than two columns found in data file";
    const int m = strlen(text) + 1;
    char* error = malloc(m*sizeof(char));
    if (!error) return NULL;
    strcpy(error, text);
    _columns = 0;
    return error;
  }

  /* Check if the other rows in the file have the same number of columns */
  fileRow = 1;
  while (1)
  { line = GetLine(file);
    if (!line) return NULL;
    if (line[0]=='\0')
    { free(line);
      break; /* Reached the end of the file */
    }
    else
    /* Parse the first column to find out what the rows contain */
    { fileColumn = 1; /* One more columns than tabs */
      for (s=line; (*s)!='\0'; s++) if(*s=='\t') fileColumn++;
      s = tokenize(line);
      if(!strcmp(line,"EWEIGHT")) arrayWeightRow=fileRow;
      else if (!strcmp(line,"EORDER")) arrayOrderRow=fileRow;
      else if (line[0]=='\0') s = NULL; /* no gene name found */ 
      else _rows++;
      free(line);
      fileRow++;
      if (s==NULL)
      { char* temp;
        char* text;
        n = 1024;
        text = malloc(n*sizeof(char));
        if (!text) return NULL;
        sprintf(text, "Gene name is missing in line %d", fileRow);
	n = strlen(text) + 1;
        temp = realloc(text,n*sizeof(char));
        if (!temp)
        { free(text);
          return NULL;
        }
        _rows = 0;
        _columns = 0;
	return temp;
      }
      if (fileColumn < nFileColumns)
      { char* temp;
        char* text;
        int n = 1024;
        text = malloc(n*sizeof(char));
        if (!text) return NULL;
        sprintf(text,
                "Error reading line %d: only %d columns found (%d expected)",
                fileRow, fileColumn, nFileColumns);
	n = strlen(text) + 1;
        temp = realloc(text,n*sizeof(char));
        if (!temp)
        { free(text);
          return NULL;
        }
        _rows = 0;
        _columns = 0;
	return temp;
      }
      if (fileColumn > nFileColumns)
      { char* temp;
        char* text;
        n = 1024;
        text = malloc(n*sizeof(char));
        if (!text) return NULL;
        sprintf(text,
                "Error reading line %d: %d columns given (%d expected)",
                fileRow, fileColumn, nFileColumns);
	n = strlen(text) + 1;
        temp = realloc(text,n*sizeof(char));
        if (!temp)
        { free(text);
          return NULL;
        }
        _rows = 0;
        _columns = 0;
	return text;
      }
    }
  }

  /* Read the first line into a string */
  rewind(file);
  line = GetLine(file);
  if (!line) return NULL;
  
  /* Save which word the user used instead of UniqID */
  s = tokenize(line);
  n = strlen(line);
  _uniqID = malloc((n+1)*sizeof(char));
  if (!_uniqID)
  { free(line);
    goto exit;
  }
  strcpy(_uniqID, line);

  /* Allocate space for array names (experiment names) and save them */
  _arrayname = calloc(_columns, sizeof(char*));
  if (!_arrayname)
  { free(line);
    goto exit;
  }
  column = 0;
  fileColumn = 1;
  while (column < _columns)
  { char* token = s;
    s = tokenize(s);
    n = strlen(token);
    if (fileColumn!=0 &&
        fileColumn!=geneNameColumn &&
        fileColumn!=geneWeightColumn &&
        fileColumn!=geneOrderColumn)
    { _arrayname[column] = malloc((n+1)*sizeof(char));
      if (!_arrayname[column])
      { free(line);
        goto exit;
      }
      strcpy(_arrayname[column], token);
      column++;
    }
    fileColumn++;
  }
  free(line);

  /* Allocate space for array weights */
  _arrayweight = malloc(_columns*sizeof(double));
  if (!_arrayweight) goto exit;
  _arrayorder = malloc(_columns*sizeof(double));
  if (!_arrayorder) goto exit;
  _arrayindex = malloc(_columns*sizeof(int));
  if (!_arrayindex) goto exit;
  for (column = 0; column < _columns; column++)
  { _arrayweight[column] = 1.;
    _arrayorder[column] = column;
  }

  /* Allocate space for data */
  _data = calloc(_rows, sizeof(double*));
  if (!_data) goto exit;
  _mask = calloc(_rows, sizeof(int*));
  if (!_mask) goto exit;
  for (row = 0; row < _rows; row++)
  { _data[row] = malloc(_columns*sizeof(double));
    _mask[row] = malloc(_columns*sizeof(int));
    if (!_data[row] || !_mask[row]) goto exit;
  }

  /* Allocate space for gene quantities */
  _geneweight = malloc(_rows*sizeof(double));
  _geneorder = malloc(_rows*sizeof(double));
  _geneindex = malloc(_rows*sizeof(int));
  _geneuniqID = calloc(_rows, sizeof(char*));
  _genename = calloc(_rows, sizeof(char*));
  if (!_geneweight || !_geneorder || !_geneindex || !_geneuniqID || !_genename)
    goto exit;

  /* Unless a GWEIGHT column exists, 
   * fill the gene weights with the default value */
  if (geneWeightColumn == -1)
    for (row = 0; row < _rows; row++) _geneweight[row] = 1.;
  /* Unless a GORDER column exist, set the gene order to the default value */
  if (geneOrderColumn == -1)
    for (row = 0; row < _rows; row++) _geneorder[row] = row;

  /* Read in gene data */
  row = 0;
  fileRow = 1;
  while (1)
  { line = GetLine(file); /* Reached end of file */
    if (!line) goto exit;
    if (line[0]=='\0')
    { free(line);
      break;
    }
    if (fileRow==arrayWeightRow)
    { column = 0;
      fileColumn = 1;
      /* Skipping UNIQID column */
      s = tokenize(line);
      while (column < _columns)
      { char* error = NULL;
        char* token = s;
        s = tokenize(s);
        if (fileColumn!=geneNameColumn &&
            fileColumn!=geneWeightColumn &&
            fileColumn!=geneOrderColumn)
        {
          _arrayweight[column] = 0; /* Default value */
          if(token[0]!='\0')
	  { const double number = strtod(token, &error);
	    if (!(*error)) _arrayweight[column] = number;
          }
          column++;
        }
        fileColumn++;
      }
    }
    else if (fileRow==arrayOrderRow)
    { column = 0;
      fileColumn = 1;
      /* Skipping UNIQID column */
      s = tokenize(line);
      while (column < _columns)
      { char* error = NULL;
        char* token = s;
        s = tokenize(s);
        if (fileColumn!=geneNameColumn &&
            fileColumn!=geneWeightColumn &&
            fileColumn!=geneOrderColumn)
        {
          _arrayorder[column] = 0; /* Default value */
          if(token[0]!='\0')
          { const double number = strtod(token, &error);
            if (!(*error)) _arrayorder[column] = number;
          }
          column++;
        }
        fileColumn++;
      }
    }
    else
    { column = 0;
      fileColumn = 0;
      s = line;
      while (s)
      { char* token = s;
        s = tokenize(s);
        if (fileColumn==0)
        { n = strlen(token) + 1;
          _geneuniqID[row] = malloc(n*sizeof(char));
          if (!_geneuniqID[row])
          { free(line);
            goto exit;
          }
          strcpy(_geneuniqID[row], token);
        }
        else if (fileColumn==geneNameColumn)
        { n = strlen(token) + 1;
          _genename[row] = malloc(n*sizeof(char));
          if (!_genename[row])
          { free(line);
            goto exit;
          }
          strcpy(_genename[row],token);
        }
        else if (fileColumn==geneWeightColumn)
        { char* error = NULL;
          double number = strtod(token, &error);
          if (!(*error)) _geneweight[row] = number;
          else _geneweight[row] = 0.;
        }
        else if (fileColumn==geneOrderColumn)
        { char* error = NULL;
          double number = strtod(token, &error);
          if (!(*error)) _geneorder[row] = number;
          else _geneorder[row] = 0.;
        }
        else
        { char* error = NULL;
          _data[row][column] = 0;
          _mask[row][column] = 0;
          if (token[0]!='\0') /* Otherwise it is a missing value */
          { double number = strtod(token, &error);
            if (!(*error))
            { _data[row][column] = number;
              _mask[row][column] = 1;
            }
          }
          column++;
        }
        fileColumn++;
      }
      row++;
    }
    fileRow++;
    free(line);
  }
  sort(_rows, _geneorder, _geneindex);
  sort(_columns, _arrayorder, _arrayindex);
  return "ok";
exit:
  Free();
  return NULL;
}

int Save(FILE* outputfile, int geneID, int arrayID)
{ int row, column;
  if (geneID) fputs("GID\t", outputfile);
  fputs(_uniqID, outputfile);
  fputs("\tNAME\tGWEIGHT", outputfile);
  /* Now add headers for data columns */
  for (column = 0; column < _columns; column++)
  { putc('\t', outputfile);
    fputs(_arrayname[_arrayindex[column]], outputfile);
  }
  putc('\n', outputfile);

  if (arrayID)
  { fputs("AID", outputfile);
    if (geneID) putc('\t',outputfile);
    fputs("\t\t", outputfile);
    for (column = 0; column < _columns; column++)
    { char* ID = MakeID("ARRY",_arrayindex[column]);
      if (!ID) return 0;
      putc('\t', outputfile);
      fputs(ID, outputfile);
      free(ID);
    }
    putc('\n', outputfile);
  }

  fputs("EWEIGHT", outputfile);
  if (geneID) putc('\t', outputfile);
  fputs("\t\t", outputfile);
  for (column = 0; column < _columns; column++)
    fprintf(outputfile, "\t%f", _arrayweight[_arrayindex[column]]);
  putc('\n', outputfile);

  for (row = 0; row < _rows; row++)
  { int index = _geneindex[row];
    if (geneID)
    { char* ID = MakeID("GENE",index);
      if (!ID) return 0;
      fputs(ID, outputfile);
      free(ID);
      putc('\t', outputfile);
    }

    fputs(_geneuniqID[index], outputfile);
    putc('\t', outputfile);
    if (_genename[index]) fputs(_genename[index], outputfile);
    else fputs(_geneuniqID[index], outputfile);
    fprintf(outputfile, "\t%f", _geneweight[index]);

    for (column = 0; column < _columns; column++)
    { int columnindex = _arrayindex[column];
      putc('\t', outputfile);
      if (_mask[index][columnindex])
        fprintf(outputfile, "%f", _data[index][columnindex]);
    }
    putc('\n', outputfile);
  }
  return 1;
}

int SelectSubset(int n, const int use[])
{ int row;
  double** data = malloc(n*sizeof(double*));
  int** mask = malloc(n*sizeof(int*));
  char** geneuniqID = malloc(n*sizeof(char*));
  char** genename = malloc(n*sizeof(char*));
  double* geneorder = malloc(n*sizeof(double));
  double* geneweight = malloc(n*sizeof(double));
  
  if (!data || !mask || !geneuniqID || !genename || !geneorder || !geneweight)
  { if (data) free(data);
    if (mask) free(mask);
    if (geneuniqID) free(geneuniqID);
    if (genename) free(genename);
    if (geneorder) free(geneorder);
    if (geneweight) free(geneweight);
    return 0;
  }

  n = 0;
  for (row = 0; row < _rows; row++)
  { if (use[row])
    { data[n] = _data[row];
      mask[n] = _mask[row];
      geneuniqID[n] = _geneuniqID[row];
      genename[n] = _genename[row];
      geneorder[n] = _geneorder[row];
      geneweight[n] = _geneweight[row];
      n++;
    }
    else
    { free(_data[row]);
      free(_mask[row]);
      free(_geneuniqID[row]);
      if (_genename[row]) free(_genename[row]);
    }
  }
  free(_data);
  free(_mask);
  free(_geneuniqID);
  free(_genename);
  free(_geneorder);
  free(_geneweight);
  _rows = n;
  _data = data;
  _mask = mask;
  _geneuniqID = geneuniqID;
  _genename = genename;
  _geneorder = geneorder;
  _geneweight = geneweight;
  sort(_rows, _geneorder, _geneindex);
  return 1;
}

void LogTransform(void)
{ int row, column;
  for (row = 0; row < _rows; row++)
  { /* Log transformation */
    for (column = 0; column < _columns; column++)
    { if (_mask[row][column] && _data[row][column] > 0)
        _data[row][column] = log(_data[row][column])/log(2.);
      else _mask[row][column]=0;
    }
  }
  return;
}

int AdjustGenes(int MeanCenter, int MedianCenter, int Normalize)
{ int row, column;
  for (row = 0; row < _rows; row++)
  { /* Center genes */
    if (MeanCenter || MedianCenter)
    { int counter = 0;
      double* temp = malloc(_columns*sizeof(double));
      if (!temp) return 0;
      for (column = 0; column < _columns; column++)
      { if (_mask[row][column])
        { temp[counter] = _data[row][column];
          counter++;
        }
      }
      if (counter > 0)
      { if (MeanCenter)
        { double rowmean = mean(counter, temp);
          for (column = 0; column < _columns; column++) 
            if (_mask[row][column]) _data[row][column] -= rowmean;
        }
        else if (MedianCenter)
        { double rowmedian = median(counter, temp);
          for (column = 0; column < _columns; column++)
            if (_mask[row][column]) _data[row][column] -= rowmedian;
        }
      }
      free(temp);
    }
    /* Normalize genes */
    if (Normalize)
    { double ssqu = 0;
      for (column = 0; column < _columns; column++)
      { if (_mask[row][column])
        { double term = _data[row][column];
          ssqu += term*term;
        }
      }
      if (ssqu > 0) /* Avoid dividing by zero */
      { double std = sqrt(ssqu);
        for (column = 0; column < _columns; column++)
          if (_mask[row][column]) _data[row][column] /= std;
      }
    }
  }
  return 1;
}

int AdjustArrays(int MeanCenter, int MedianCenter, int Normalize)
{ int row, column;
  /* Center Arrays */
  if (MeanCenter || MedianCenter)
  { double* temp = malloc(_rows*sizeof(double));
    if (!temp) return 0;
    for (column = 0; column < _columns; column++)
    { int counter = 0;
      for (row = 0; row < _rows; row++)
      { if (_mask[row][column])
        { temp[counter] = _data[row][column];
          counter++;
        }
      }
      if (counter > 0)
      { if (MeanCenter)
        { double columnmean = mean(counter,temp);
          for (row = 0; row < _rows; row++)
            if (_mask[row][column])
              _data[row][column] -= columnmean;
        }
        else if (MedianCenter)
        { double columnmedian = median(counter,temp);
          for (row = 0; row < _rows; row++)
            if (_mask[row][column])
              _data[row][column] -= columnmedian;
        }
      }
    }
    free(temp);
  }

  /* Normalize arrays */
  if (Normalize)
  { for (column = 0; column < _columns; column++)
    { double ssqu = 0;
      for (row = 0; row < _rows; row++)
        if (_mask[row][column])
        { double term = _data[row][column];
          ssqu += term * term;
        }
      if (ssqu > 0) /* Avoid dividing by zero */
      { double std = sqrt(ssqu);
        for (row = 0; row < _rows; row++)
          if (_mask[row][column]) _data[row][column] /= std;
      }
    }
  }
  return 1;
}

int PerformSOM(FILE* GeneFile, int GeneXDim, int GeneYDim, int GeneIters,
  double GeneTau, char GeneMetric, FILE* ArrayFile, int ArrayXDim,
  int ArrayYDim, int ArrayIters, double ArrayTau, char ArrayMetric)
{ int ok = 1;
  if (GeneIters>0)
  { ok = PerformGeneSOM(GeneFile,
                        GeneXDim,
                        GeneYDim,
                        GeneIters,
                        GeneTau,
                        GeneMetric);
    if (!ok) return 0;
  }
  else sort(_rows, _geneorder, _geneindex);
  if (ArrayIters>0)
  { ok = PerformArraySOM(ArrayFile,
                         ArrayXDim,
                         ArrayYDim,
                         ArrayIters,
                         ArrayTau,
                         ArrayMetric);
    if (!ok) return 0;
  }
  else sort(_columns, _arrayorder, _arrayindex);
  return ok;
}

int FilterRow(int Row, int bStd, int bPercent, int bAbsVal, int bMaxMin,
  double absVal, double percent, double std, int numberAbs, double maxmin)
{ int Count = 0;
  int CountAbs = 0;
  double Sum = 0;
  double Sum2 = 0;
  double Min = 10000000;
  double Max = -10000000;
  /* Compute some row stats */
  int Column;
  for (Column = 0; Column < _columns; Column++)
  { if (_mask[Row][Column])
    { double value = _data[Row][Column];
      Sum += value;
      Sum2 += value*value;
      Count++;
      Min = min(value,Min);
      Max = max(value,Max);
      if (fabs(value) >= absVal) CountAbs++;
    }
  }
  /* Filter based on percent values present;
   * remove rows with too many missing values.
   */
  if (bPercent)
  { int number = (int) ceil(percent*_columns/100);
    if (Count < number) return 0;
  }
  /* Remove rows with low SD */
  if (bStd)
  { if (Count > 1)
    { double Ave = Sum / (double) Count;
      double Var = (Sum2 - 2 * Ave * Sum + Count * Ave * Ave)/ (Count-1);
      if (sqrt(Var) < std) return 0;
    }
    else return 0;
  }
  /* Remove rows with too few extreme values */
  if (bAbsVal && CountAbs < numberAbs) return 0;
  /* Remove rows with too small Max-Min */
  if (bMaxMin && Max - Min < maxmin) return 0;
  return 1;
}

const char*
CalculateWeights(double GeneCutoff, double GeneExponent, char GeneDist,
  double ArrayCutoff, double ArrayExponent, char ArrayDist)
{ double* geneweight = NULL;
  double* arrayweight = NULL;
  if (GeneCutoff && GeneExponent && GeneDist)
  { geneweight = calculate_weights(_rows, _columns, _data, _mask,
                                   _arrayweight, 0, GeneDist,
                                   GeneCutoff, GeneExponent);
    if (!geneweight)
      return "Insufficient memory to calculate the row weights";
  }
  if (ArrayCutoff && ArrayExponent && ArrayDist)
  { arrayweight = calculate_weights(_rows, _columns, _data, _mask,
                                    _geneweight, 1, ArrayDist,
                                    ArrayCutoff, ArrayExponent);
    if (!arrayweight)
    { if (geneweight) free(geneweight);
      return "Insufficient memory to calculate the column weights";
    }
  }
  if (geneweight)
  { free(_geneweight);
    _geneweight = geneweight;
  }
  if (arrayweight)
  { free(_arrayweight);
    _arrayweight = arrayweight;
  }
  return NULL;
}


int HierarchicalCluster(FILE* file, char metric, int transpose, char method)
{ int i;
  int ok = 0;
  const int nNodes = (transpose ? _columns : _rows) - 1;
  int* index = (transpose==0) ? _geneindex : _arrayindex;
  const double* order = (transpose==0) ? _geneorder : _arrayorder;
  double* weight = (transpose==0) ? _arrayweight : _geneweight;
  const char* keyword = (transpose==0) ? "GENE" : "ARRY";
 
  char** nodeID = calloc(nNodes, sizeof(char*));
  /* Perform hierarchical clustering. */
  Node* tree = treecluster(_rows, _columns, _data, _mask, weight, transpose,
                           metric, method, NULL);
  if (!tree || !nodeID)
  { if (tree) free(tree);
    if (nodeID) free(nodeID);
    return 0;
  }

  if (metric=='e' || metric=='b')
  /* Scale all distances such that they are between 0 and 1 */
  { double scale = 0.0;
    for (i = 0; i < nNodes; i++)
      if (tree[i].distance > scale) scale = tree[i].distance;
    if (scale) for (i = 0; i < nNodes; i++) tree[i].distance /= scale;
  }

  /* Now we join nodes */
  for (i = 0; i < nNodes; i++)
  { int min1 = tree[i].left;
    int min2 = tree[i].right;
    /* min1 and min2 are the elements that are to be joined */
    char* ID1;
    char* ID2;
    nodeID[i] = MakeID("NODE",i+1);
    if (!nodeID[i]) break;
    if (min1 < 0)
    { int index1 = -min1-1;
      ID1 = nodeID[index1];
      tree[i].distance = max(tree[i].distance, tree[index1].distance);
    }
    else
      ID1 = MakeID(keyword, min1);
    if (min2 < 0)
    { int index2 = -min2-1;
      ID2 = nodeID[index2];
      tree[i].distance = max(tree[i].distance, tree[index2].distance);
    }
    else
      ID2 = MakeID(keyword, min2);
 
    if (ID1 && ID2)
    { fprintf(file, "%s\t%s\t%s\t", nodeID[i], ID1, ID2);
      fprintf(file, "%f\n", 1.0-tree[i].distance);
    }
    if (ID1 && min1>=0) free(ID1);
    if (ID2 && min2>=0) free(ID2);
    if (!ID1 || !ID2) break;
  }

  /* Now set up order based on the tree structure */
  if (i==nNodes) /* Otherwise we encountered the break */
    ok = sorttree(nNodes, tree, order, index);
  for (i = 0; i < nNodes; i++) if (nodeID[i]) free(nodeID[i]);
  free(nodeID);
  free(tree);

  return ok;
}

int GeneKCluster(int k, int nTrials, char method, char dist, int* NodeMap)
{ int ifound = 0;
  double error;
  int ok;
  kcluster(k, _rows, _columns, _data, _mask,
    _arrayweight, 0, nTrials, method, dist, NodeMap, &error, &ifound);
  ok = SetClusterIndex('g', k, NodeMap);
  if (ok) return ifound;
  return -1;
}

int ArrayKCluster(int k, int nTrials, char method, char dist, int* NodeMap)
{ int ifound = 0;
  double error;
  int ok;
  kcluster(k, _rows, _columns, _data, _mask,
    _geneweight, 1, nTrials, method, dist, NodeMap, &error, &ifound);
  ok = SetClusterIndex('a', k, NodeMap);
  if (ok) return ifound;
  return -1;
}

int SaveGeneKCluster(FILE* file, int k, const int* NodeMap)
{ int i, cluster;
  int* geneindex = malloc(_rows*sizeof(int));
  if (!geneindex) return 0;
  fprintf(file, "%s\tGROUP\n", _uniqID);
  for (i=0; i<_rows; i++) geneindex[i] = i;
  sort(_rows,_geneorder,geneindex);
  for (cluster = 0; cluster < k; cluster++)
  { for (i = 0; i < _rows; i++)
    { const int j = geneindex[i];
      if (NodeMap[j]==cluster)
        fprintf(file, "%s\t%d\n", _geneuniqID[j], NodeMap[j]);
    }
  }
  free(geneindex);
  return 1;
}

int SaveArrayKCluster(FILE* file, int k, const int* NodeMap)
{ int i, cluster;
  int* arrayindex = malloc(_columns*sizeof(int));
  if (!arrayindex) return 0;
  fputs("ARRAY\tGROUP\n", file);
  for (i=0; i<_columns; i++) arrayindex[i] = i;
  sort(_columns,_arrayorder,arrayindex);
  for (cluster = 0; cluster < k; cluster++)
  { for (i = 0; i < _columns; i++)
    { const int j = arrayindex[i];
      if (NodeMap[j]==cluster)
        fprintf(file, "%s\t%d\n", _arrayname[j], NodeMap[j]);
    }
  }
  free(arrayindex);
  return 1;
}

const char* PerformGenePCA(FILE* coordinatefile, FILE* pcfile)
{
  int i = 0;
  int j = 0;
  const int nmin = min(_rows,_columns);
  double** u = malloc(_rows*sizeof(double*));
  double** v = malloc(nmin*sizeof(double*));
  double* w = malloc(nmin*sizeof(double));
  double* m = malloc(_columns*sizeof(double));
  if (u)
  { for (i = 0; i < _rows; i++)
    { u[i] = malloc(_columns*sizeof(double));
      if (!u[i]) break;
    }
  }
  if (v)
  { for (j = 0; j < nmin; j++)
    { v[j] = malloc(nmin*sizeof(double));
      if (!v[j]) break;
    }
  }
  if (!u || !v || !w || !m || i < _rows || j < nmin)
  { if (u)
    { while (i--) free(u[i]);
      free(u);
    }
    if (v)
    { while (j--) free(v[j]);
      free(v);
    }
    if (w) free(w);
    if (m) free(m);
    return "Insufficient Memory for PCA calculation";
  }
  for (j = 0; j < _columns; j++)
  { double value;
    m[j] = 0.0;
    for (i = 0; i < _rows; i++)
    { value = _data[i][j];
      u[i][j] = value;
      m[j] += value;
    }
    m[j] /= _rows;
    for (i = 0; i < _rows; i++) u[i][j] -= m[j];
  }
  pca(_rows, _columns, u, v, w);
  fprintf(coordinatefile, "%s\tNAME\tGWEIGHT", _uniqID);
  for (j=0; j < nmin; j++)
    fprintf(coordinatefile, "\t%f", w[j]);
  putc('\n', coordinatefile);
  fprintf(pcfile, "EIGVALUE");
  for (j=0; j < _columns; j++)
    fprintf(pcfile, "\t%s", _arrayname[j]);
  putc('\n', pcfile);
  fprintf(pcfile, "MEAN");
  for (j=0; j < _columns; j++)
    fprintf(pcfile, "\t%f", m[j]);
  putc('\n', pcfile);
  if (_rows>_columns)
  { for (i=0; i<_rows; i++)
    { fprintf(coordinatefile, "%s\t",_geneuniqID[i]);
      if (_genename[i]) fputs(_genename[i], coordinatefile);
      else fputs(_geneuniqID[i], coordinatefile);
      fprintf(coordinatefile, "\t%f", _geneweight[i]);
      for (j=0; j<_columns; j++)
        fprintf(coordinatefile, "\t%f", u[i][j]);
      putc('\n', coordinatefile);
    }
    for (i = 0; i < nmin; i++)
    { fprintf(pcfile, "%f", w[i]);
      for (j=0; j < _columns; j++)
        fprintf(pcfile, "\t%f", v[i][j]);
      putc('\n', pcfile);
    }
  }
  else
  { for (i=0; i<_rows; i++)
    { fprintf(coordinatefile, "%s\t",_geneuniqID[i]);
      if (_genename[i]) fputs(_genename[i], coordinatefile);
      else fputs(_geneuniqID[i], coordinatefile);
      fprintf(coordinatefile, "\t%f", _geneweight[i]);
      for (j=0; j<nmin; j++)
        fprintf(coordinatefile, "\t%f", v[i][j]);
      putc('\n', coordinatefile);
    }
    for (i = 0; i < _rows; i++)
    { fprintf(pcfile, "%f", w[i]);
      for (j=0; j < _columns; j++)
        fprintf(pcfile, "\t%f", u[i][j]);
      putc('\n', pcfile);
    }
  }
  for (i = 0; i < _rows; i++) free(u[i]);
  for (i = 0; i < nmin; i++) free(v[i]);
  free(u);
  free(v);
  free(w); 
  free(m);
  return NULL;
}

const char* PerformArrayPCA(FILE* coordinatefile, FILE* pcfile)
{
  int i = 0;
  int j = 0;
  const int nmin = min(_rows,_columns);
  double** u = malloc(_columns*sizeof(double*));
  double** v = malloc(nmin*sizeof(double*));
  double* w = malloc(nmin*sizeof(double));
  double* m = malloc(_rows*sizeof(double));
  if (u)
  { for (i = 0; i < _columns; i++)
    { u[i] = malloc(_rows*sizeof(double));
      if (!u[i]) break;
    }
  }
  if (v)
  { for (j = 0; j < nmin; j++)
    { v[j] = malloc(nmin*sizeof(double));
      if (!v[j]) break;
    }
  }
  if (!u || !v || !w || !m || i < _columns || j < nmin)
  { if (u)
    { while (i--) free(u[i]);
      free(u);
    }
    if (v)
    { while (j--) free(v[j]);
      free(v);
    }
    if (w) free(w);
    if (m) free(m);
    return "Insufficient Memory for PCA calculation";
  }
  for (j = 0; j < _rows; j++)
  { double value;
    m[j] = 0.0;
    for (i = 0; i < _columns; i++)
    { value = _data[j][i];
      u[i][j] = value;
      m[j] += value;
    }
    m[j] /= _columns;
    for (i = 0; i < _columns; i++) u[i][j] -= m[j];
  }
  pca(_columns, _rows, u, v, w);
  fprintf(coordinatefile, "EIGVALUE");
  for (j=0; j < _columns; j++)
    fprintf(coordinatefile, "\t%s", _arrayname[j]);
  putc('\n', coordinatefile);
  fprintf(coordinatefile, "EWEIGHT");
  for (j=0; j < _columns; j++)
    fprintf(coordinatefile, "\t%f", _arrayweight[j]);
  putc('\n', coordinatefile);
  fprintf(pcfile, "%s\tNAME\tMEAN", _uniqID);
  for (j=0; j < nmin; j++)
    fprintf(pcfile, "\t%f", w[j]);
  putc('\n', pcfile);
  if (_rows>_columns)
  { for (i = 0; i < nmin; i++)
    { fprintf(coordinatefile, "%f", w[i]);
      for (j=0; j<_columns; j++)
        fprintf(coordinatefile, "\t%f", v[j][i]);
      putc('\n', coordinatefile);
    }
    for (i = 0; i < _rows; i++)
    { fprintf(pcfile, "%s\t",_geneuniqID[i]);
      if (_genename[i]) fputs(_genename[i], pcfile);
      else fputs(_geneuniqID[i], pcfile);
      fprintf(pcfile, "\t%f", m[i]);
      for (j=0; j<_columns; j++)
        fprintf(pcfile, "\t%f", u[j][i]);
      putc('\n', pcfile);
    }
  }
  else /* _rows < _columns */
  { for (i=0; i<_rows; i++)
    { fprintf(coordinatefile, "%f", w[i]);
      for (j=0; j<_columns; j++)
        fprintf(coordinatefile, "\t%f", u[j][i]);
      putc('\n', coordinatefile);
    }
    for (i = 0; i < _rows; i++)
    { fprintf(pcfile, "%s\t",_geneuniqID[i]);
      if (_genename[i]) fputs(_genename[i], pcfile);
      else fputs(_geneuniqID[i], pcfile);
      fprintf(pcfile, "\t%f", m[i]);
      for (j=0; j<nmin; j++)
        fprintf(pcfile, "\t%f", v[j][i]);
      putc('\n', pcfile);
    }
  }
  for (i = 0; i < _columns; i++) free(u[i]);
  for (i = 0; i < nmin; i++) free(v[i]);
  free(u);
  free(v);
  free(w); 
  free(m);
  return NULL;
}
